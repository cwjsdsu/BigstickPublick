
!=============================================================
!   subroutine expectator_p
!       master subroutine to control calculation of expectation values
!
!  NOTE: There will have to be modifications when sections of the lanczos vector are broken into pieces
!
!  CALLS
!   setup_localvectors
!   read_wfn_nkeep
!   dnormvec_p
!   wfn_readeigenvec
!    dnormvec_p
!  br_load_vec2_from_vec1
!   applyobsbundled
subroutine expectator_p

  use precisions
  use basis
  use io
  use localvectors
  use obs
  use fragments
  use nodeinfo
  use flags3body
  use mod_reorthog
  use wfn_mod
  use bmpi_mod
  use bvectorlib_mod
  use lanczos_util
  use apply_obs_mod
  implicit none

  integer :: ikeep
  integer :: i,j,n
  real :: xe,xj,xt
  real (kind = lanc_prec) :: exph
  real (kind = 8) :: dnorm
  logical :: zeroflag
  integer :: ierr

  !.....NOTE: LIMITED TO 2-BODY OBSERVABLES....CAN BE FIXED BUT MUST ADD ROUTINES TO applybsbundled 

  if(iproc==0 .and. threebody)print*,' Expectation values not yet set up for 3-body forces'
!........... TEMP...... WILL NEED TO FIX LATER....
 
  call setup_localvectors
  call wfn_read_nkeep(oldwfnfile, nkeep) ! does BCAST
  if(iproc==0)then
       print*,nkeep,' states '
       if(writeout)write(resultfile,'(''nkeep ='',1x,i10)')nkeep
       print*,' '
       print*,' STATE       E          J            <H >        (norm)'
       if(writeout)write(resultfile,*)'  STATE      E            J           T^2           <H >      (norm)'
  end if
  twoobsflag = .false.  ! may change
  do ikeep = 1,nkeep
     i = ikeep
     ! new interface, we say which vector to read - it checks
     call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,xe,xj,xt) ! KSM: updated
     call dnormvec_p('n','i',dnorm,zeroflag)
     call br_load_vec2_from_vec1()
     xj2 = 0.0e0_lanc_prec
     call applyobsbundled(1)
     if(iproc==0)then
         write(6,101)i,xe,xj,xj2, dnorm
!         if(writeout)write(resultfile,101)i,e,xj,xj2,dnorm
         if(writeout)write(resultfile,101)i,xe,xj,xt,xj2,dnorm
101  format(i5,2x,2(1x,f11.4),2(1x,f13.6),1x,f10.5)
     end if
  end do ! ikeep
  
  return
end subroutine expectator_p

!==================================================================
! 
! added 7.11.1 July 2024
! for option '(h)'
!
! compute matrix elements of a scalar Hamiltonian-like 1+2-body operator
! in a basis of previously computed states
!
subroutine hmatrixelements_boss
    use localvectors
    use bvectorlib_mod
	
    use localblocks
    use mod_reorthog
	use para_util_mod
	use lanczos_info
	use io
	use nodeinfo
	use lanczos_util
    use wfn_mod
	use apply_obs_mod
    use apply_ham
    use apply_ham_omp
    use apply_ham_block
    use pocc_mod
	
	implicit none
	integer :: i,f
	real, allocatable :: hme(:,:)
	real(8) :: dsclrprod = 0.d0
	
    call wfn_read_nkeep(oldwfnfile, nkeep)  ! does BCAST
	if(iproc==0)then
		print*,' There are ',nkeep,' states '
	end if
	niter = nkeep+1
    call setup_localvectors

	call distribute_lanczos_pieces
	
	allocate(etx(nkeep),j2tx(nkeep),t2tx(nkeep),partx(nkeep))
	call read_pivot_block(.true.)
	if(iproc==0)write(resultfile,'(" ")')

	call pocc_write_table_header() !  E   Ex ...
	if(iproc==0)then
	   do i = 1,nkeep

          write(resultfile,111)i,etx(i),etx(i)-etx(1),j2tx(i),t2tx(i),partx(i)
	   
        end do
	end if

11 format(i5,3x,2f10.5,2x,2f8.3)
111 format(i5,3x,2f10.5,2x,2f8.3,3x,i2)
  
	
	allocate(hme(nkeep,nkeep))
	hme = 0.0
	
!	call br_set_histpos(0)
	call procOP_clock(iproc,'set','all')
	
	do i = 1,nkeep
		call br_retrieve_hist(i)
		call br_restore_vec1()
        call initialize_final('n')
        call applyhbundled_g('n')			
		do f = 1,i			
			

			call br_retrieve_hist(f)
					
	        call br_restore_vec1()		
				
			dsclrprod = 0.d0
	 	    call doverlapvec(dsclrprod,.true.)  ! modified in 7.7.0 so as to work with fragments
			hme(i,f)= dsclrprod
			hme(f,i)= dsclrprod
!			br_histpos = br_histpos+1
			
		end do
!		br_histpos = br_histpos+1
	end do
	
	if(iproc==0)print*,' all done '
	
	if(iproc==0)then
		write(resultfile,'(" ")')
		do i = 1,nkeep
			do f = 1,i
				if(abs(j2tx(i)-j2tx(f)) > 0.0001)cycle
				if(abs(hme(i,f))> 1e-6)write(resultfile,*)i,f,hme(i,f)
				
				
			end do
		end do
		write(resultfile,*)-1,-1,0.0
		
	end if
!	print*,hme
		
	return
	
end subroutine hmatrixelements_boss



!==================================================================
