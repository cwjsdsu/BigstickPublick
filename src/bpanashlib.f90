!  ADDED 7.11.3 started April, 2025 by CWJ@SDSU
! library of routines useful for constructing a basis in PANASh
! mostly through organizing wave vectors
! 
! Two new options: 'jp' to project out states of good J
!              and 'ro' to orthonormalize sets of vectors read in
!
module wfn_organize
	use bmpi_mod
	use nodeinfo
	implicit none
	
	integer :: nkeep_tmp ! target nkeep, though not all may be kept
	
contains
	
!.... BOSS ROUTINE FOR OPTION 'jp' TO PROJECT OUT STATES OF GOOD J
!.... NOTE THAT J^2 has already been set up in main program
!
!	
	subroutine Jproject_boss
		use io
		use lanczos_info
		use wfn_mod
		use bvectorlib_mod
		use localvectors
		use mod_reorthog
		use para_util_mod
	    use apply_ham
	    use apply_ham_omp
		use lanczos_util
		implicit none
		
		integer :: nvectors_in ! # of vectors read in
		integer :: op_rank ! rank of operator, i.e., 0, 1, 2...
		     ! note: may have to generalize to half-integers, in the case of adding/removing a single nucleon
		real(4) :: tol_def = 0.05
		real(4) :: tol
		integer :: ikeep
		integer :: istatestart, istatestop
		integer :: nkeep_tmp,dim_tmp,nkeep_new
		integer :: i, ivec,k
		real(kind=4) :: ee,xxjj,xxt,xxtt
		real(kind=8) :: dnorm0,da,db
		logical,allocatable :: keep(:)
		real, allocatable :: xjlist(:),elist(:)
	    real(kind = egv_prec),allocatable ::  work3(:)  ! work(np)
	    real(kind=egv_prec), allocatable ::  jlvec(:,:),xje(:)
		real(kind=4) :: xjj
		real(kind=4) :: beta_prev  ! used for convergence test
		
		integer :: ierr
		
		integer :: info
!...... open file..... get # of vectors available
		
		if(nproc==1)then
			storelanczosincore1 = .true.
		else
			storelanczosincoreMPI = .true.
		end if
	    call distribute_lanczos_pieces
		
	    call wfn_read_nkeep(oldwfnfile, ikeep) ! does BCAST
		
!        do i = 1,nkeep
!              ! new interface - we say which vec to read, it checks
!              ! KSM:  This will be very slow, only need to read head part of each vector
!              call wfn_readeigenvec(oldwfnfile,frag1, fcomm1_index, vec1,i,xe,xj,xt2)  
!              if(iproc==0)print*,i,xe,xj,xt2
!        enddo


!...... how many vectors to read in?

        if(iproc==0)then
			print*,' '
			print*,' There are ', ikeep,' initial wavefunctions '
			print*,' '
            print*,' Enter start, stop for initial states '	   	
            print*,' (Enter 0,0  to read all )'
            read*,istatestart,istatestop
			if(istatestart==0)istatestart = 1
			if(istatestop==0)istatestop = ikeep
			nvectors_in = istatestop-istatestart+1
			print*,' Reading in ',nvectors_in,' vectors, starting with # ',istatestart		
				
!...... what is rank of operator, i.e., how many for each vector			
			print*,' What is the rank (i.e, angular momentum) of operator? '
			read*,op_rank
			print*,' Operator rank = ',op_rank
			
!..... determine tolerance			
			tol = tol_def ! here one could read in the tolerance
			print*,' Tolerance for keeping a projected vector is ',tol
			
		end if
!..... BCAST		
#ifdef _MPI
    call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(istatestart,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(istatestop,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(nvectors_in,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(tol,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(op_rank,1,0,MPI_COMM_WORLD,ierr)
#endif        

!...... set up arrays.....

        nkeep_tmp = nvectors_in * (2*op_rank+1)  ! this is the max # of vectors we can generate; might be fewer
		niter = nkeep_tmp
		
	    call setup_localvectors
		
		allocate(alpha(2*op_rank+1),beta(2*op_rank+1))
		allocate(elist(nvectors_in*(2*op_rank+1)))
		allocate(xjlist(nvectors_in*(2*op_rank+1)))
		allocate(keep(nvectors_in*(2*op_rank+1)))
		allocate(work3(3*(2*op_rank+1)))
		allocate(jlvec(2*op_rank+1,2*op_rank+1))
		allocate(xje(2*op_rank+1))
		
		elist = 0.0
		xjlist = 0.0
		keep = .false.
		
!.....	loop over initial vectors....
		call br_set_histpos(0)
		call br_reset_histbuf()
		ivec = 1
        do i = istatestart,istatestop
   	    ! frag1 says that the slicing of the vector follows vec1
          ! new interface, we say which vec to read and it checks
           call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,ee,xxjj,xxtt)
!           xt =(-0.5 + sqrt(xtt+0.25))
           xxt = 0.0
		   
!           if(iproc==0)write(filenumber,*)e,xj,xt
           if(iproc==0)print*,i,ee,xxjj,xxt
!.............. NORMALIZE VECTOR, just in case it's not normalized		   
           call br_grab_vec1()! push vec1 -> br_reg		
		   call br_normalize(dnorm0) !normalizes br_reg
!		   call br_restore_vec1()       
          call br_set_histpos(ivec-1)

           call br_add2hist(ivec)	
		   keep(ivec) = .true.
		   elist(ivec) = ee
		   xjlist(ivec) = xxjj 
		   ivec = ivec + 2*op_rank+1 ! need to skip	   
        enddo		

	    call   procOP_clock(0,'set','all')

!..... apply J^2 in Lanczos....
        ivec = 1
        do i = 1,nvectors_in
!.......... read vector from history			

			alpha = 0.d0
			beta = 0.d0
			dim_tmp = 1
			call br_set_histpos(ivec)
			call br_retrieve_hist(ivec)
			call br_restore_vec1()

		    do k = 1,2*op_rank+1

	            call initialize_final('n')
				
!		        if (useHZSomp .and. num_threads_global>1) then
!		              call applyHbundled_omp('n')    ! THE CENTRAL ENGINE: APPLY MAT-VEC MULTIPLY
!		        else
		               call applyHbundled_g('n')    ! THE CENTRAL ENGINE: APPLY MAT-VEC MULTIPLY
!		        end if  
				
		        call br_grab_vec2()   ! push vec2 -> br_reg
		        call br_remove_prev_overlap(da)
		        call br_orthogonalize_restricted(ivec,1)! 1 says ignore last history entry (prev)

		        call br_normalize(db)
		        alpha(k) = da !real(da,kind(0.0e0))
		        beta(k) = db !real(db, kind(0.0e0))
!				print*,k, alpha(k),beta(k)
				! check size of db
!				print*,abs(db),tol

                if(k==1)then
					beta_prev = 1.0
				else
					beta_prev = beta(k-1)
				end if
	
				if(abs(db) < tol*beta_prev  .or. k==2*op_rank+1)then
					dim_tmp = k
					exit
				end if
				dim_tmp = k
				call br_restore_vec1()
				call br_add2hist(ivec+k)
			
			  
		    end do !k

			! NOW SOLVE 
			jlvec = 0.0
			do k = 1,dim_tmp
				jlvec(k,k)=alpha(k)
				if(k < dim_tmp)then
					jlvec(k,k+1)=beta(k)
					jlvec(k+1,k)=beta(k)
				end if
				
			end do
!			print*,' ALPHAS ',alpha(1:dim_tmp)
!			print*,' BETAS ',beta(1:dim_tmp)
            call DSYEV( 'V','U', dim_tmp, jlvec, 2*op_rank+1, xje, WORK3, 3*(2*op_rank+1), INFO )
			if(iproc==0)print*,' J(J+1) vals ',xje(1:dim_tmp)
			
			! TRANFORM THE SELECTED VECTORS AND EXTRACT J-VALUES, and write to history
			call br_transform_basis_restricted(2*op_rank+1,jlvec,ivec,ivec+dim_tmp-1)
			do k = 1,dim_tmp
				xjj = real( -0.5 + sqrt(xje(k) + 0.25), kind=4)

				xjlist(ivec+k-1)=xjj
				keep(ivec +k-1) = .true.
				elist(ivec+k-1) = elist(ivec)
				
			end do
			
 		    ivec = ivec + 2*op_rank+1 ! need to skip	   
			
		end do

!   a test as to whether vectors are properly placed		
!		do ivec = 1,niter
!			call br_set_histpos(ivec)			
!			call br_retrieve_hist(ivec)
!			call br_normalize(db)
!			print*,' Checking normalization of temp vectors ',db
!		end do

!...... write to file

!........ count up how many to keep

        nkeep_new = 0
		do ivec = 1,niter
			if(keep(ivec))nkeep_new = nkeep_new +1
			
		end do
		if(iproc==0)print*,' Keeping ',nkeep_new, 'vectors (out of ',nkeep_tmp,')'

!....... NOW WRITE NKEEP TO NEW FILE....
        call wfn_write_nkeep(nkeep_new) ! write number of vectors to wfn file
		
!......  WRITE OUT PROJECTED STATES......

        ivec = 0		

		do k = 1,niter
			call br_set_histpos(k)
			
			call br_retrieve_hist(k)
			if(.not.keep(k))cycle
			ivec = ivec + 1
			call br_restore_vec1()
			call wfn_writeeigenvec(wfnfile,frag1, vec1,ivec,elist(k),xjlist(k),0.0)			

		end do
        call wfn_close_file(wfnfile)
		
		
		return
	end subroutine Jproject_boss

!=================================================================================	
!...... BOSS ROUTINE FOR OPTION 'ro' TO READ IN AND ORTHONORMALIZE WAVE VECTORS....
	
	subroutine read_n_ortho_boss(orthonorm)
		use io
		use lanczos_info
		use wfn_mod
		use bvectorlib_mod
		use localvectors
		use mod_reorthog
		
		implicit none
		logical :: orthonorm
		integer :: nkeep_tmp ! # of vectors read in
		integer :: nkeep_final
		real(4) :: tol_def = 0.001
		real(4) :: tol
		logical :: full,first_file
		integer :: ikeep
		integer :: istatestart, istatestop
		integer :: nvectors_in ! # of vectors read in
		integer :: sum_vectors
		integer :: i,ivec
		real(4) :: e,xj,xtt,xt
		real(8) :: dnorm0
		logical,allocatable :: keep(:)
		real, allocatable :: xjlist(:),elist(:)
		integer :: ierr

!...... how many vectors to read in?

        if(iproc==0)then
			print*,' '
			print*,' How many vectors to read in? '
			read*,nkeep_tmp
			print*,' Reading in ',nkeep_tmp,' vectors '		
			
!..... determine tolerance			
			tol = tol_def ! here one could read in the tolerance
			print*,' Tolerance for keeping a vector is ',tol
			
		end if
!..... BCAST		
#ifdef _MPI
    call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(nkeep_tmp,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(tol,1,0,MPI_COMM_WORLD,ierr)
#endif 		
!..... reserve rooom for nkeep vectors...
		
		if(nproc==1)then
			storelanczosincore1 = .true.
		else
			storelanczosincoreMPI = .true.
		end if
	    call overlaptribution
		
        niter = nkeep_tmp
		allocate(keep(nkeep_tmp),xjlist(nkeep_tmp),elist(nkeep_tmp))

        call setup_localvectors

!...... read in files until nvectors_in saved
        full = .false.
		first_file = .true.
		ikeep = 0
		sum_vectors = 0
		call br_set_histpos(0)
		call br_reset_histbuf()
		ivec = 0
        do while (.not.full)
			if(first_file)then
				first_file = .false.
			    call wfn_read_nkeep(oldwfnfile, ikeep) ! does BCAST
				if(iproc==0)print*,' The first file has ',ikeep,' wave vectors '
			else
				if(iproc==0)print*,' - - -  NEXT FILE - - - '
			    call wfn_ropen_file(oldwfnfile)
			    call read_wfn_header(oldwfnfile,.true.)		
			    call wfn_read_nkeep(oldwfnfile, ikeep) ! does BCAST
				call checkbasisdim
				if(iproc==0)print*,' This file has ',ikeep,' wave vectors ' 
			end if
			if(iproc==0)then
                print*,' Enter start, stop for initial states '	   	
                print*,' (Enter 0,0  to read all )'
                read*,istatestart,istatestop
			end if
!--- BCAST---------------			
#ifdef _MPI
    call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(ikeep,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(istatestart,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(istatestop,1,0,MPI_COMM_WORLD,ierr)
#endif 
			
			if(istatestart==0)istatestart = 1
			if(istatestop==0)istatestop = ikeep
			nvectors_in = istatestop-istatestart+1
			
			if(sum_vectors+nvectors_in > nkeep_tmp)then
				nvectors_in = nkeep_tmp-sum_vectors +1
				istatestop = istatestart +nvectors_in -1
				
				if(iproc==0)print*,' Final file, taking only ',istatestart,' through ',istatestop
			end if
				
!............. LOOP OVER WAVE VECTORS AND STORE IN HISTORY...............				
			
            do i = istatestart,istatestop
	            call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,e,xj,xtt)
!           xt =(-0.5 + sqrt(xtt+0.25))
	            xt = 0.0
                call br_normalize(dnorm0)
		   
!           if(iproc==0)write(filenumber,*)e,xj,xt
	            if(iproc==0)print*,i,e,xj,xt,dnorm0
!.............. NORMALIZE VECTOR, just in case it's not normalized
  	            ivec = ivec+1
				if(ivec > nkeep_tmp)then
					print*,' problem, error in reading too many vectors ',ivec
					stop 891
				end if
	            call br_grab_vec1()! push vec1 -> br_reg		   
	            call br_add2hist(ivec)					
				xjlist(ivec)=xj
				elist(ivec)=e
			end do ! i
			
            call wfn_close_file(oldwfnfile)
			oldwfnfile=14 ! A KLUGE, because normally this is set  = 0; check if it works with MPI
			
			sum_vectors = sum_vectors+nvectors_in
			if(ivec /= sum_vectors)then
				print*,' miscount of vectors ',ivec,sum_vectors
				stop 899
			end if
			if(sum_vectors== nkeep_tmp)then
				if(iproc==0)print*,' All vectors read in '
				if(iproc==0 .and. orthonorm)print*,' Next: orthonormalize! '
				full = .true.
			end if
		end do

!......Gram-Schmidt orthogonalization
   if(orthonorm)then
        nkeep_final = 1
		call br_set_histpos(0)
		keep(1)=.true.

        do i = 2,nkeep_tmp
	        call br_orthogonalize_in_place(i,1,i-1) ! 0 says include last history entry (prev)
            call br_normalize_in_place(i,dnorm0)			
			if(iproc==0)print*,i,'norm ',dnorm0
			if(dnorm0 > tol)then
				nkeep_final = nkeep_final + 1
				keep(i)=.true.
			else
				keep(i)=.false.
			end if			
		end do
		if(iproc==0)print*,' Keeping ',nkeep_final,' states '
	else
		nkeep_final = nkeep_tmp
		keep=.true.
	end if
!...... write to disk.....
        call wfn_write_nkeep(nkeep_final) ! write number of vectors to wfn file

        ivec = 0
		
        do i = 1,nkeep_tmp
			if(keep(i))then
				ivec = ivec + 1
				call br_set_histpos(ivec)
				call br_retrieve_hist(ivec)
				call br_restore_vec1()
				call wfn_writeeigenvec(wfnfile,frag1, vec1,ivec,elist(ivec),xjlist(ivec),0.0)			
			end if		
		end do
        call wfn_close_file(wfnfile)
		
!...... close out
        return
	end subroutine read_n_ortho_boss	
	
	
end module wfn_organize