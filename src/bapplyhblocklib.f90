!===================================================================
!
!  file BAPPLYHLIB.f90
!
!  routines to control application of the Hamiltonian
!  for parallel calculations with (optional) fragmentation of the basis and vectors
!
!  started 9/2011 by CWJ
!
!  NOTE: As of version 7.0.0, "obsolete" libraries of apply H have been discarded
!
!===================================================================
!
!  added in 7.5.1; allows for contiguous pointers
!  this attribute is from fortran 2003;
!  if your compiler does not include fortran 2003, 
!  comment out the statements with "contiguous" and 
!  uncomment the lines without
!
!  added in 7.8.3: skip 'annexed' opbundles; removed in 7.8.4 
!


module apply_ham_block
	use localblocks
contains

!
! subroutine applyHbundled
!
! note: default is going from vecin to vecout but this can be reversed depending on hchar
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)
!   vchar = 'n' (normal), 'r' (reverse)
!      (fragments) of lanczos vectors stored in module localvectors
!        in vec1 and vec2; if normal  H vec1 = vec2
!                          if reverse H vec2 = vec1

subroutine applyHbundled_block
  use nodeinfo
  use flagger
  use flags3body
  use precisions
  use opbundles
  use fragments
  use interaction
  use basis
  use localvectors
  use localblocks
  use bmpi_mod
  use lanczos_info
  use coupledmatrixelements,only:call_spe
  implicit none

  integer iprocs, procstart,procstop
  real(kind=8) :: sum
  integer :: tid
  integer(kind=basis_prec) :: i, j
  integer :: ierr

!........ OPTION TO SIMULATE MPI ON A SINGLE "CORE".....

  if(distributeMPI .and. nproc == 1)then
      procstart = 0
      procstop  = nprocs -1
  else
      procstart = iproc
      procstop  = iproc
  end if

  call clocker('hmu','sta')

  if(noisy0) print *, "Starting applyHbundled_g"
  
!  write(99,*)iproc,' starting ',blockvec_1
  do iprocs = procstart,procstop
     call proc_clock(iprocs,'sta')

     if(.not.threebody)then
!........... PP .................

     if(noisy0) print *, "Starting applyHbundled_g PP"
     call clocker('ppo','sta')
     call procOP_clock(iprocs,'sta','PPO')
     call applyhPPbundled_block('f',opbundlestart(iprocs), opbundleend(iprocs))
     call procOP_clock(iprocs,'end','PPO')

     call clocker('ppo','end')
     call clocker('ppb','sta')
     call procOP_clock(iprocs,'sta','PPB')
     call applyhPPbundled_block('b',opbundlestart(iprocs), opbundleend(iprocs))
     call procOP_clock(iprocs,'end','PPB')

     call clocker('ppb','end')
!............NN ..................
     call clocker('nno','sta')
     call procOP_clock(iprocs,'sta','NNO')

     if(noisy0) print *, "Starting applyHbundled_g NN"
     call applyhNNbundled_block('f',opbundlestart(iprocs), opbundleend(iprocs))
     call applyhNNbundled_block('b',opbundlestart(iprocs), opbundleend(iprocs))
     call clocker('nno','end')
     call procOP_clock(iprocs,'end','NNO')


!........... PN....................................
     call clocker('one','sta')

     if(noisy0) print *, "Starting applyHbundled_g PN"

	     call clocker('pno','sta')
	     call procOP_clock(iprocs,'sta','PNO')
        call applyhPNbundled_block('f',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNbundled_block('h',opbundlestart(iprocs), opbundleend(iprocs))
		call clocker('pno','end')
        call procOP_clock(iprocs,'end','PNO')
     call procOP_clock(iprocs,'sta','PNB')

        call clocker('pnb','sta')
        call applyhPNbundled_block('b',opbundlestart(iprocs), opbundleend(iprocs))
		call clocker('pnb','end')
        call procOP_clock(iprocs,'end','PNB')
		
     call clocker('one','end')

     else
!............PPP.......................................
    call clocker('ppp','sta')
     call procOP_clock(iprocs,'sta','PPP')

    call applyhPPPbundled_block('f',opbundlestart(iprocs), opbundleend(iprocs))
    call applyhPPPbundled_block('b',opbundlestart(iprocs), opbundleend(iprocs))
    call clocker('ppp','end')
     call procOP_clock(iprocs,'end','PPP')

!............PPN.......................................
!    call clocker('ppn','sta')
 !    call procOP_clock(iprocs,'sta','PPN')

        call clocker('pXf','sta')
	    call procOP_clock(iprocs,'sta','PXF')

        call applyhPPNbundled_block('f',opbundlestart(iprocs), opbundleend(iprocs))
        call clocker('pXf','end')
        call procOP_clock(iprocs,'end','PXF')
	    call procOP_clock(iprocs,'sta','PXB')
        call clocker('pXb','sta')
        call applyhPPNbundled_block('b',opbundlestart(iprocs), opbundleend(iprocs))
        call clocker('pXb','end')

     call procOP_clock(iprocs,'end','PXB')

!    call clocker('ppn','end')
!............PNN.......................................
!    call clocker('pnn','sta')
call procOP_clock(iprocs,'sta','PYF')

     call clocker('pYf','sta')
        call applyhPNNbundled_block('f',opbundlestart(iprocs), opbundleend(iprocs))
        call clocker('pYf','end')
        call clocker('pYb','sta')
        call procOP_clock(iprocs,'end','PYF')
	    call procOP_clock(iprocs,'sta','PYB')
		
        call applyhPNNbundled_block('b',opbundlestart(iprocs), opbundleend(iprocs))
        call clocker('pYb','end')

!		call clocker('pnn','end')
     call procOP_clock(iprocs,'end','PYB')

!...........NNN......................................
    call clocker('nnn','sta')
     call procOP_clock(iprocs,'sta','NNN')

    call applyhNNNbundled_block('f',opbundlestart(iprocs), opbundleend(iprocs))
    call applyhNNNbundled_block('b',opbundlestart(iprocs), opbundleend(iprocs))
    call clocker('nnn','end')
     call procOP_clock(iprocs,'end','NNN')


    end if  ! if threebody


    call proc_clock(iprocs,'end')

!........... NOW REDUCE....................
   call clocker('hmu','end')


   if(noisy0) print *, "applyhbundled_g: Doing Reduce"
         ! 'n' corresponds with vec2 here because we are looking at the output vector
          ! call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
          ! Do reduce only onto root node.  We will be sending this data to 
          ! the slices from the isfragroot nodes (rank=0 in fcomm1, fcomm2, hcomm),
          ! so other nodes don't need it.
		  
		  call clocker('wait','sta')
		  
!		  call BMPI_Barrier(icomm,ierr)   ! REMOVE WHEN DONE
		  call clocker('wai','end')
		  
    call clocker('blr','sta')
		  
!	write(99,*)iproc,' VECTOR ',blockvec_2
#ifdef _MPI
    call BMPI_REDUCE(blockvec_2, size(blockvec_2), MPI_SUM, 0, fcomm2, ierr) ! in place
#endif
    call clocker('blr','end')
	
  end do  ! iprocs

  if(noisy0) print *, "applyhbundled_g: Returning"
  return

end subroutine applyHbundled_block
!==================================================

!  subroutine applyhPPbundled
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)
!   vchar = 'n' (normal), 'r' (reverse)
!      (fragments) of lanczos vectors stored in module localvectors
!        in vec1 and vec2; if normal  H vec1 = vec2
!                          if reverse H vec2 = vec1
!
! cleaned up in 7.7.9: removed some loops that have become disused
!
!===================================================================
subroutine applyhPPbundled_block (hchar,startbundle,endbundle )

  use nodeinfo
  use localvectors
  use system_parameters
!  use sectors
  use jumpNbody
  use precisions
  use interaction
!  use lanczos_info
  use opbundles
  use fragments
  use basis
  use lanczos_info
  use flagger
  use bmpi_mod
  use butil_mod
  use contigpointervectors, only :  p2b_1sd,p2b_2sd
  implicit none

  logical :: pinfo
  integer :: ibundle
  character(1) :: hchar
  integer :: startbundle,endbundle
  integer :: ivec

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=4)   xme, prod
  integer(kind=8) :: statei, statef,nsd,psdi,psdf
  integer(kind=8) :: locationi,locationf
  integer(kind=basis_prec) :: statefoff, stateistart, stateistop
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(4) :: mythread


!..............................................................
!
!  SET UP POINTERS


! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
        if( hchar == 'f')then
           p2b_1sd => p2b_isd
           p2b_2sd => p2b_fsd
        else
           p2b_1sd => p2b_fsd
           p2b_2sd => p2b_isd
        endif

  do ibundle = startbundle,endbundle

     if(opbundle(ibundle)%optype /= 'PP')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
	 
	 if(diagonalsectorsonly .and. opbundle(ibundle)%isector /= opbundle(ibundle)%fsector)cycle
	 
	 call bundle_clock(ibundle,'sta')

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
  csdstart = opbundle(ibundle)%nxstart
  csdend   = opbundle(ibundle)%nxend
  xjmpstart = opbundle(ibundle)%pxstart
  xjmpend   = opbundle(ibundle)%pxend
  cstride   = opbundle(ibundle)%cstride

  ncstates = (csdend +cstride -csdstart)/cstride

!--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
!          this makes for simple OpenMP threading
!$omp parallel private(xjmp, Xoplabel, xme, num_threads, mythread,nsd, pinfo)    &
!$omp          private(istart, iend, chunk, csd, csd_index, ivec,statef,statei)  &
!$omp          private(statefoff, stateistart, stateistop, prod)  &
!$omp          firstprivate(cstride, ncstates, xjmpstart,xjmpend)  &
!$omp          shared(p2b_op, p2b_1sd, p2b_2sd, p2b_phase, hmatpp,block_1,block_2)
  num_threads =  omp_get_num_threads()
  mythread = omp_get_thread_num()
  ! thread local vec2, reduce at end


! KSM:  chunks are guarenteed not to overlap on statef, so we
! KSM:  don't have to worry about about collisions between threads.
! KSM:  each thread gets a different range of neutron SDs.
  pinfo = ncstates > 10
  chunk = (ncstates + num_threads - 1)/num_threads
  istart = mythread*chunk + 1
  iend = bmin((mythread + 1)*chunk,ncstates)
  csd_index = csdstart + (istart - 1)*cstride - cstride
  if(istart <= iend)then

!......... THERE ARE TWO VERSIONS.....
!          1st way: store in jumps index to PP matrix element;
!          this has faster setup
!          2nd way: store PP matrix elements directly;
!          slower set up, but on MPI nodes reduced memory load
	
  do csd = istart, iend
     csd_index = csd_index + cstride
     nsd = nstart(csd_index)

!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
          Xoplabel = p2b_op(xjmp)
          xme = hmatpp(Xoplabel)
!--------- GET PHASE.........................................
          xme = xme*p2b_phase(xjmp)
!---------- GET INITIAL, FINAL SDs and place in basis..............

          statei = p2b_1sd(xjmp)+ nsd !csd_index
          statef = p2b_2sd(xjmp)+nsd !csd_index
		  locationi = (statei-1)*dimblock
		  locationf = (statef-1)*dimblock
		  
		  do ivec = 1,dimblock
             blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
		  end do
      end do  ! xjmp
   end do  ! csd
end if
!$omp end parallel
  call bundle_clock(ibundle,'end')

  end do ! ibundle
  return
end subroutine applyhPPbundled_block

!===================================================================
!  subroutine applyhNNbundled
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., NN between two (sub) sectors, etc)
!   vchar = 'n' (normal), 'r' (reverse)
!      (fragments) of lanczos vectors stored in module localvectors
!        in vec1 and vec2; if normal  H vec1 = vec2
!                          if reverse H vec2 = vec1
!   
!
!===================================================================
subroutine applyhNNbundled_block (hchar,startbundle,endbundle )
   use localvectors
   use nodeinfo
   use system_parameters
   use jumpNbody
   use precisions
   use interaction
   use opbundles
   use fragments
   use basis
   use lanczos_info
   use flagger
   use bmpi_mod
   use butil_mod
   use contigpointervectors, only :  n2b_1sd,n2b_2sd
   implicit none

   ! arguments
   character(1),intent(in) :: hchar
   integer,intent(in) :: startbundle,endbundle

! --------- NOTE: basestart, basestop stored in module fragments

!------------------------------------------------------------

   integer :: ibundle
   integer(kind=8) csdstart, csdend,csd, csd_index,cstride,ncstates, pstride
   integer(kind=8) xjmp,xjmpstart,xjmpend
   integer(kind=8):: Xoplabel
   real(kind=4)   xme
   integer(kind=8) :: statei, statef,psd,nsdi,nsdf,statef_start,statef_end
   integer(kind=8) :: locationi,locationf
   
   integer(kind=8) :: istart, iend, chunk

!-------- OpenMP functions ---------------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: num_threads
   integer :: mythread
   integer :: ivec

!..............................................................
!..............................................................
!
!  SET UP POINTERS

         if( hchar == 'f' )then
            n2b_1sd => n2b_isd
            n2b_2sd => n2b_fsd
         else
            n2b_1sd => n2b_fsd
            n2b_2sd => n2b_isd
         endif


      do ibundle = startbundle,endbundle
         if(opbundle(ibundle)%optype /= 'NN')cycle
         if(opbundle(ibundle)%hchar /= hchar )cycle
!		 if(opbundle(ibundle)%annexed)cycle
		 
         if(diagonalsectorsonly .and. opbundle(ibundle)%isector /= opbundle(ibundle)%fsector)cycle
	     call bundle_clock(ibundle,'sta')
		 
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
         csdstart = opbundle(ibundle)%pxstart
         csdend   = opbundle(ibundle)%pxend
         cstride  = opbundle(ibundle)%cstride   !
         xjmpstart = opbundle(ibundle)%nxstart
         xjmpend   = opbundle(ibundle)%nxend
         ncstates = (csdend +cstride -csdstart)/cstride
!--------- OUTER LOOP OVER CONJUGATE PROTON SDs---------
!          this makes for simple OpenMP threading
!       NOTE CSTRIDE OVER PROTON SDs 

! firstprivate gives each thread its own copy, but initializes it
!    better than shared for read-only vars
! private gives each thread its own copy
!$omp parallel private(xjmp, Xoplabel, xme, num_threads, mythread,psd)         &
!$omp          private(istart, iend, chunk, csd, csd_index, ivec,statef,statei)  &
!$omp          firstprivate(cstride, ncstates, xjmpstart, xjmpend)  &
!$omp          shared(n2b_op, n2b_1sd, n2b_2sd, n2b_phase, hmatnn,block_1,block_2)
         num_threads =  omp_get_num_threads()
         mythread = omp_get_thread_num()
         ! thread local vec2, reduce at end

         chunk = (ncstates + num_threads - 1)/num_threads
         istart = mythread*chunk + 1
         iend = bmin((mythread + 1)*chunk,ncstates)
         csd_index = csdstart + (istart - 1)*cstride - cstride
         if(istart <= iend)then

!..... THIS FOLLOWING IS TO TRY TO FIND AN OPTIMAL ORDERING OF LOOPS
            pstride = pstridecut +1 
            if(iend-istart > 0)pstride   = pstart(csdstart+cstride)-pstart(csdstart)


                  do csd = istart, iend
                     csd_index = csd_index + cstride
                     psd = pstart(csd_index)
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
                     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
                        Xoplabel = n2b_op(xjmp)
!		if(Xoplabel==0)print*,' ZERO LABEL ',iproc,ibundle,xjmp
                        xme = hmatnn(Xoplabel)
!--------- GET PHASE.........................................
                        xme = xme*n2b_phase(xjmp)
!---------- GET INITIAL, FINAL SDs and place in basis..............
                        statei = n2b_1sd(xjmp)+psd ! csd_index
                        statef = n2b_2sd(xjmp)+psd  !csd_index
			  		  locationi = (statei-1)*dimblock
			  		  locationf = (statef-1)*dimblock
			  		    do ivec = 1,dimblock
! 			               block_2(ivec,statef) = block_2(ivec,statef) + xme*block_1(ivec,statei)
                           blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
			  		    end do
                     end do  ! xjmp
                  end do  ! csd

			  end if

!$omp end parallel
call bundle_clock(ibundle,'end')

   end do ! ibundle
   return
end subroutine applyhNNbundled_block

!=================================================================
!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPNbundled_block (hchar,startbundle,endbundle )
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use interaction
  use opbundles
  use fragments
  use lanczos_info
  use flagger
  use bmpi_mod
  implicit none

  integer :: ibundle,startbundle,endbundle
  character(1) :: hchar

!------------------------------------------------------------
  integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer :: a,b,c,d
  integer(kind=8) :: coplabel,doplabel
  integer :: phasep,phasen
  real(kind=4) ::   xme
  integer(kind=basis_prec) :: statei, statef,locationi,locationf
  integer(kind=4) num_threads

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  integer :: ivec

  if(applyXXonly)return    ! don't do PN
!..............................................................
!..............................................................
!
!  SET UP POINTERS


  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
	 
	 if(diagonalsectorsonly .and. opbundle(ibundle)%isector /= opbundle(ibundle)%fsector)cycle
     call bundle_clock(ibundle,'sta')

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads

! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( hchar /= 'b' )then

!$omp parallel do private( mythread,startp_thread,npjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)   &
!$omp          private(coplabel,doplabel,xme,ivec,statei,statef)                   &
!$omp          firstprivate(njmpstart, njmpend)       &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)              &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cpnpair,dpnpair,block_1,block_2)
     do mythread = 0,numpthreads -1


     startp_thread = opbundle(ibundle)%startp_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     npjmps_thread = opbundle(ibundle)%startp_thread(mythread+1) - startp_thread

! KSM:  start/stop set up so that each proton final state appears on only one thread
! KSM:  prevents collison over update of voutp(statef) below
!---------   Forward direction ------------------------------
     do pjmp = startp_thread + 1, startp_thread + npjmps_thread
        psdi = p1b_isd(pjmp)       ! initial proton slater determinant
        psdf = p1b_fsd(pjmp)       ! final proton SD
        phasep = p1b_phase(pjmp)   ! phase of proton jumps
        a = p1b_cop(pjmp)     ! KSM: Proton 1-body creation label
        c = p1b_dop(pjmp)     ! KSM: Proton 1-body destruction label
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
        do njmp = njmpstart,njmpend
!----------- FIND MATRIX ELEMTN --------------------------------------------
           b = n1b_cop(njmp)  ! KSM: Neutron 1-body creation label
           d = n1b_dop(njmp)  ! KSM: Neutron 1-body destruction label
           phasen = n1b_phase(njmp)
           coplabel = cpnpair(b,a)
           doplabel = dpnpair(d,c)
           xme = hmatpn(coplabel + doplabel)   ! get matrix element
           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
           nsdi = n1b_isd(njmp)
           nsdf = n1b_fsd(njmp)
           statei = nsdi + psdi                ! initial state in combined basis
           statef = nsdf + psdf                ! final state in combined basis
		   locationi = (statei-1)*dimblock
		   locationf = (statef-1)*dimblock
		   do ivec = 1,dimblock
               blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
		   end do
!		  do ivec = 1,dimblock
!             block_2(ivec,statef) = block_2(ivec,statef) + xme*block_1(ivec,statei)
! 		  end do
        end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do
else
!---- Backward direction using hermiticity ------------------- 

!$omp parallel do private(mythread,startn_thread,nnjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)   &
!$omp          private(coplabel,doplabel,xme,ivec,statei,statef)                   &
!$omp          firstprivate(pjmpstart, pjmpend)                     &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)    &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)    &
!$omp          shared(cpnpair,dpnpair,block_1,block_2)
  do mythread = 0, numnthreads-1
     ! thread local vec2, reduce at end


     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread

!...... OPTION TO SWITCH ORDER OF LOOPS WHEN NO OpenMP.......

     if(numnthreads > 1 .or. disableNoOMPloopswitch)then

     do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
        nsdi = n1b_isd(njmp)
        nsdf = n1b_fsd(njmp)
        b  = n1b_cop(njmp)
        d  = n1b_dop(njmp)
        phasen = n1b_phase(njmp)
        do pjmp = pjmpstart,pjmpend
           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
           psdf = p1b_fsd(pjmp)       ! final proton SD
           phasep = p1b_phase(pjmp)   ! phase of proton jumps
           a  = p1b_cop(pjmp) 
           c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
!----------- FIND MATRIX ELEMENT -------------------------------------------
           coplabel = cpnpair(b,a)
           doplabel = dpnpair(d,c)		   
           xme = hmatpn(coplabel + doplabel)     ! get matrix element
           xme = xme*phasep*phasen               ! multiply matrix element by jump phases
           statei = nsdi + psdi                  ! initial state in combined basis
           statef = nsdf + psdf                  ! final state in combined basis
		   locationi = (statei-1)*dimblock
		   locationf = (statef-1)*dimblock
		   do ivec = 1,dimblock
               blockvec_2(ivec+locationi) = blockvec_2(ivec+locationi) + xme*blockvec_1(ivec+locationf)
		   end do
		   
! 		  do ivec = 1,dimblock
!              block_2(ivec,statei) = block_2(ivec,statei) + xme*block_1(ivec,statef)
! 		  end do
        end do  ! pjmp
     end do  ! njmp

     else


	     startp_thread = opbundle(ibundle)%startp_thread(mythread)     !  starting position for proton 1-body jumps for this thread
	     npjmps_thread = opbundle(ibundle)%startp_thread(mythread+1) - startp_thread
	     do pjmp = startp_thread + 1, startp_thread + npjmps_thread
	        psdi = p1b_isd(pjmp)       ! initial proton slater determinant
	        psdf = p1b_fsd(pjmp)       ! final proton SD
	        phasep = p1b_phase(pjmp)   ! phase of proton jumps
	        a = p1b_cop(pjmp)     ! KSM: Proton 1-body creation label
	        c = p1b_dop(pjmp)     ! KSM: Proton 1-body destruction label
	!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
	        do njmp = njmpstart,njmpend
	!----------- FIND MATRIX ELEMTN --------------------------------------------
	           b = n1b_cop(njmp)  ! KSM: Neutron 1-body creation label
	           d = n1b_dop(njmp)  ! KSM: Neutron 1-body destruction label
	           phasen = n1b_phase(njmp)
	           coplabel = cpnpair(b,a)
	           doplabel = dpnpair(d,c)
	           xme = hmatpn(coplabel + doplabel)   ! get matrix element
	           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
	           nsdi = n1b_isd(njmp)
	           nsdf = n1b_fsd(njmp)
	           statei = nsdi + psdi                ! initial state in combined basis
	           statef = nsdf + psdf                ! final state in combined basis
			   locationi = (statei-1)*dimblock
			   locationf = (statef-1)*dimblock
			   do ivec = 1,dimblock
	               blockvec_2(ivec+locationi) = blockvec_2(ivec+locationi) + xme*blockvec_1(ivec+locationf)
			   end do
!	 		  do ivec = 1,dimblock
!	              block_2(ivec,statei) = block_2(ivec,statei) + xme*block_1(ivec,statef)
!	 		  end do
	        end do  ! njmp
	     end do  ! pjmp        




     end if
  end do
!$omp end parallel do


  end if
  call bundle_clock(ibundle,'end')

  end do  ! ibundle

  return
end subroutine applyhPNbundled_block

!===============================================================
!
!  apply 3-body jumps for PPP
!
subroutine applyhPPPbundled_block (hchar,startbundle,endbundle )


  use nodeinfo
  use localvectors
  use system_parameters
  use precisions
  use jump3body
  use interactions3body
  use opbundles
  use fragments
  use basis
  use flagger
  use bmpi_mod
  use butil_mod
  use contigpointervectors, only :  p3b_1sd,p3b_2sd
  implicit none

  integer :: ibundle,startbundle,endbundle
  character(1) :: hchar

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=4)   xme
  integer(kind=8) :: statei, statef,nsd,psdi,psdf
  integer(kind=8) :: locationi,locationf
  
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(4) :: mythread
  integer(4) :: ivec

!..............................................................
!
!  SET UP POINTERS

!
        if( hchar == 'f')then
           p3b_1sd => p3b_isd
           p3b_2sd => p3b_fsd
        else
           p3b_1sd => p3b_fsd
           p3b_2sd => p3b_isd
        endif

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PPP')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle

     call bundle_clock(ibundle,'sta')

!...... EXTRACT INFORMATION FROM OPBUNDLE ........

  csdstart = opbundle(ibundle)%nxstart
  csdend   = opbundle(ibundle)%nxend
  xjmpstart = opbundle(ibundle)%pxstart
  xjmpend   = opbundle(ibundle)%pxend
  cstride   = opbundle(ibundle)%cstride

  ncstates = (csdend +cstride -csdstart)/cstride

!--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
!          this makes for simple OpenMP threading
!$omp parallel private(xjmp, Xoplabel, xme, num_threads, mythread)      &
!$omp          private(istart, iend, chunk, csd, csd_index, ivec,statef,statei,nsd)  &
!$omp          shared( cstride,  ncstates,xjmpstart,xjmpend)  &
!$omp          shared(p3b_op, p3b_1sd, p3b_2sd, p3b_phase, hmatppp,block_1,block_2)
  num_threads =  omp_get_num_threads()
  mythread = omp_get_thread_num()
  ! thread local vec2, reduce at end


  chunk = (ncstates + num_threads - 1)/num_threads
  istart = mythread*chunk + 1
  iend = bmin((mythread + 1)*chunk,ncstates)
  csd_index = csdstart + (istart - 1)*cstride - cstride
  if(istart <= iend)then

  if(num_threads > 1 .or. disableNoOMPloopswitch)then
  do csd = istart, iend
     csd_index = csd_index + cstride
     nsd = nstart(csd_index)
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
          Xoplabel = p3b_op(xjmp)
          xme = hmatppp(Xoplabel)
!--------- GET PHASE.........................................
          xme = xme*p3b_phase(xjmp)
!---------- GET INITIAL, FINAL SDs and place in basis..............

          statei = p3b_1sd(xjmp)+nsd !csd_index
          statef = p3b_2sd(xjmp)+nsd !csd_index
	      locationi = (statei-1)*dimblock
	      locationf = (statef-1)*dimblock
	      do ivec = 1,dimblock
              blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
	      end do
! 		  do ivec = 1,dimblock
!              block_2(ivec,statef) = block_2(ivec,statef) + xme*block_1(ivec,statei)
! 		  end do
      end do  ! xjmp
   end do  ! csd

   else  ! ONLY 1 THREAD - ADDED VERSION 7.1.6 JULY 2012................

!--------- LOOP OVER 3-BODY JUMPS IN THIS SECTOR JUMPS.............
     csd_index = nstart(csd_index+1)-1
     cstride = 1
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
          Xoplabel = p3b_op(xjmp)
          xme = hmatppp(Xoplabel)
!--------- GET PHASE.........................................
          xme = xme*p3b_phase(xjmp)
          psdi = p3b_1sd(xjmp)
          psdf = p3b_2sd(xjmp)
          nsd = csd_index
       
!---------- GET INITIAL, FINAL SDs and place in basis..............
          do csd = istart, iend
              nsd = nsd + cstride
              statei = psdi+ nsd 
              statef = psdf+nsd 
   	          locationi = (statei-1)*dimblock
   		      locationf = (statef-1)*dimblock
   		      do ivec = 1,dimblock
                  blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
   		      end do
!	 		  do ivec = 1,dimblock
!	              block_2(ivec,statef) = block_2(ivec,statef) + xme*block_1(ivec,statei)
!	 		  end do
          end do
      end do  ! xjmp

   end if
   end if
!$omp end parallel
call bundle_clock(ibundle,'end')

   end do  ! ibundle
   return
   end subroutine applyhPPPbundled_block
!============================================================
!
!  apply 3-body jumps for NNN
!
!

subroutine applyhNNNbundled_block (hchar,startbundle,endbundle )
   use nodeinfo
   use localvectors
   use system_parameters
   use precisions
   use jump3body
   use interactions3body
   use opbundles
   use fragments
   use basis
   use flagger
   use bmpi_mod
   use butil_mod
   use contigpointervectors, only :  n3b_1sd,n3b_2sd
   implicit none

   integer :: ibundle,startbundle,endbundle
   character(1) :: hchar

!------------------------------------------------------------

   integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index,pstride
   integer(kind=8) xjmp,xjmpstart,xjmpend
   integer(kind=8):: Xoplabel
   real(kind=4)   xme
   integer(kind=8) :: statei, statef,psd,nsdi,nsdf,statef_start,statef_end
   integer(kind=8) :: locationi,locationf
   
!-------- OpenMP functions ---------------------------------
   integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
   integer(kind=4) :: num_threads
   integer(kind=8) :: istart, iend, chunk
   integer(4) :: mythread
   integer(4) :: ivec

!..............................................................
!
!  SET UP POINTERS

         if( hchar == 'f')then
            n3b_1sd => n3b_isd
            n3b_2sd => n3b_fsd
         else
            n3b_1sd => n3b_fsd
            n3b_2sd => n3b_isd
         endif


      do ibundle = startbundle,endbundle
         if(opbundle(ibundle)%optype /= 'NNN')cycle
         if(opbundle(ibundle)%hchar /= hchar )cycle
	     call bundle_clock(ibundle,'sta')
		 
!		 if(opbundle(ibundle)%annexed)cycle

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!        csdstart = pstart(opbundle(ibundle)%pxstart)
!        csdend   = pstart(opbundle(ibundle)%pxend)
         csdstart = opbundle(ibundle)%pxstart
         csdend   = opbundle(ibundle)%pxend
         xjmpstart = opbundle(ibundle)%nxstart
         xjmpend   = opbundle(ibundle)%nxend
         cstride   = opbundle(ibundle)%cstride

         ncstates = (csdend +cstride -csdstart)/cstride

!--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
!          this makes for simple OpenMP threading
!$omp parallel private(xjmp, Xoplabel, xme, num_threads, mythread)      &
!$omp          private(istart, iend, chunk, csd, csd_index, ivec,statef,statei)  &
!$omp          private(statef_start,statef_end,psd) & 
!$omp          shared(cstride, ncstates,xjmpstart,xjmpend,pstride)  &
!$omp          shared(n3b_op, n3b_1sd, n3b_2sd, n3b_phase, hmatnnn,block_1,block_2)
         num_threads =  omp_get_num_threads()
         mythread = omp_get_thread_num()
         ! thread local vec2, reduce at end


         chunk = (ncstates + num_threads - 1)/num_threads
         istart = mythread*chunk + 1
         iend = bmin((mythread + 1)*chunk,ncstates)
         csd_index = csdstart + (istart - 1)*cstride - cstride
         if(istart <= iend)then

!..... THIS FOLLOWING IS TO TRY TO FIND AN OPTIMAL ORDERING OF LOOPS
            pstride = pstridecut +1 
            if(iend-istart > 0)pstride   = pstart(csdstart+cstride)-pstart(csdstart)
            if(num_threads > 1 .or. pstride > pstridecut .or. disableNoOMPloopswitch)then

            do csd = istart, iend
               csd_index = csd_index + cstride
               psd = pstart(csd_index)
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
               do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
                  Xoplabel = n3b_op(xjmp)
                  xme = hmatnnn(Xoplabel)
!--------- GET PHASE.........................................
                  xme = xme*n3b_phase(xjmp)
!---------- GET INITIAL, FINAL SDs and place in basis..............
                  statei = n3b_1sd(xjmp)+psd !csd_index
                  statef = n3b_2sd(xjmp)+psd !csd_index
	   		   locationi = (statei-1)*dimblock
	   		   locationf = (statef-1)*dimblock
	   		   do ivec = 1,dimblock
	                  blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
	   		   end do
!		 		  do ivec = 1,dimblock
!		              block_2(ivec,statef) = block_2(ivec,statef) + xme*block_1(ivec,statei)
!		 		  end do
               end do  ! xjmp
            end do  ! csd

         else  ! only 1 thread  ADDED V7.1.6 July 2013
!--------- LOOP OVER 3-BODY JUMPS IN THIS SECTOR JUMPS.............
            csd_index = pstart(csdstart)     
            do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
               Xoplabel = n3b_op(xjmp)
               xme = hmatnnn(Xoplabel)
!--------- GET PHASE.........................................
               xme = xme*n3b_phase(xjmp)
               statei = n3b_1sd(xjmp) + csd_index
               statef_start = n3b_2sd(xjmp) +csd_index
!---------- GET INITIAL, FINAL SDs and place in basis..............
               statef_end = statef_start+pstride*(iend-istart)
               do statef = statef_start,statef_end,pstride
				   locationi = (statei-1)*dimblock
				   locationf = (statef-1)*dimblock
				   do ivec = 1,dimblock
		               blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
				   end do
!		  		  do ivec = 1,dimblock
!		               block_2(ivec,statef) = block_2(ivec,statef) + xme*block_1(ivec,statei)
!		  		  end do
                  statei = statei+pstride
               end do  !
            end do  ! xjmp
         end if
      end if
!$omp end parallel
call bundle_clock(ibundle,'end')

   end do !ibundle
   return
end subroutine applyhNNNbundled_block

!==========================================================
!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPNNbundled_block (hchar,startbundle,endbundle )
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use jump3body
  use interactions3body
  use opbundles
  use fragments
  use flagger
  use bmpi_mod
  use contigpointervectors, only :  n2b_1sd,n2b_2sd
  
  implicit none

  integer :: ibundle,startbundle,endbundle
   character(1) :: hchar

!------------------------------------------------------------
  integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer a,bb,c,dd
  integer(kind=8) :: coplabel,doplabel
  integer :: phasep,phasen
  real(kind=4)   xme
  integer(kind=8) :: statei, statef
  integer(kind=8) :: locationi,locationf
  
  integer(kind=4) num_threads
  integer(4)  :: ivec

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  logical    :: launched(0:3)


!  print*,iproc,' proc ',hchar,startbundle,endbundle

!..............................................................

!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if(  hchar /= 'b')then

  do ibundle = startbundle,endbundle
	  
!	  print*,iproc,ibundle,opbundle(ibundle)%optype,opbundle(ibundle)%hchar
     if(opbundle(ibundle)%optype /= 'PNN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
	 
     call bundle_clock(ibundle,'sta')
	 
!	 if(opbundle(ibundle)%annexed)cycle
	 
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

!$omp parallel do private( mythread,startp_thread,npjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,ivec,statei,statef)                   &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)              &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(block_1,block_2)   &
!$omp          shared(n2b_isd,n2b_fsd,n2b_phase,n2b_cop,n2b_dop)              &
!$omp          shared(cpnntriplet,dpnntriplet)
   do mythread = 0,numpthreads -1
     ! thread local vec2, reduce at end
 
     startp_thread = opbundle(ibundle)%startp_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     npjmps_thread = opbundle(ibundle)%startp_thread(mythread+1) - startp_thread

!---------   Forward direction ------------------------------
    do pjmp = startp_thread + 1, startp_thread + npjmps_thread
        psdi = p1b_isd(pjmp)       ! initial proton slater determinant
        psdf = p1b_fsd(pjmp)       ! final proton SD
        phasep = p1b_phase(pjmp)   ! phase of proton jumps
        a = p1b_cop(pjmp) 
        c = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
        do njmp = njmpstart,njmpend
!----------- FIND MATRIX ELEMTN --------------------------------------------

           phasen = n2b_phase(njmp)
           bb = n2b_cop(njmp)
           dd = n2b_dop(njmp)

           coplabel = cpnntriplet(bb,a)
           doplabel = dpnntriplet(dd,c)

           xme = hmatpnn(coplabel + doplabel)   ! get matrix element
           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
           nsdi = n2b_isd(njmp)
           nsdf = n2b_fsd(njmp)
           statei = nsdi + psdi                ! initial state in combined basis
           statef = nsdf + psdf                ! final state in combined basis
		   locationi = (statei-1)*dimblock
		   locationf = (statef-1)*dimblock
!		   if(ibundle==7 .and.hchar=='b')print*,' bundle7a ',locationi,locationf,xme,blockvec_1(1+locationi)
		   
!if(blockvec_1(1+locationi)/=0.0)print*,ibundle,locationi,locationf,xme*blockvec_1(1+locationi)
		   
		   do ivec = 1,dimblock
               blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
		   end do
!  		  do ivec = 1,dimblock
 !              block_2(ivec,statef) = block_2(ivec,statef) + xme*block_1(ivec,statei)
  !		  end do
        end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do
call bundle_clock(ibundle,'end')

  end do ! ibundle
else

  do ibundle = startbundle,endbundle
!	  print*,iproc,ibundle,opbundle(ibundle)%optype,opbundle(ibundle)%hchar
	  
     if(opbundle(ibundle)%optype /= 'PNN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
	 call bundle_clock(ibundle,'sta')
	 
!	 if(opbundle(ibundle)%annexed)cycle
	 
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

!---- Backward direction using hermiticity ------------------- 

!$omp parallel do private( mythread,startn_thread,nnjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,ivec,statei,statef)                   &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)    &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(n2b_isd,n2b_fsd,n2b_phase,n2b_cop,n2b_dop)    &
!$omp          shared(cpnntriplet,dpnntriplet)         &
!$omp          shared(block_1,block_2)
    do mythread = 0, numnthreads-1
     ! thread local vec2, reduce at end


     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread
!...... OPTION TO SWITCH ORDER OF LOOPS WHEN NO OpenMP.......

!     if(numnthreads > 1 .or. disableNoOMPloopswitch)then
     do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
        nsdi = n2b_isd(njmp)
        nsdf = n2b_fsd(njmp)
        bb  = n2b_cop(njmp)
        dd  = n2b_dop(njmp)
        phasen = n2b_phase(njmp)
        do pjmp = pjmpstart,pjmpend
           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
           psdf = p1b_fsd(pjmp)       ! final proton SD
           phasep = p1b_phase(pjmp)   ! phase of proton jumps
           a  = p1b_cop(pjmp) 
           c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
!----------- FIND MATRIX ELEMENT -------------------------------------------
           coplabel = cpnntriplet(bb,a)
           doplabel = dpnntriplet(dd,c)
           xme = hmatpnn(coplabel + doplabel)     ! get matrix element
           xme = xme*phasep*phasen               ! multiply matrix element by jump phases
           statei = nsdi + psdi                  ! initial state in combined basis
           statef = nsdf + psdf                  ! final state in combined basis
		   locationi = (statei-1)*dimblock
		   locationf = (statef-1)*dimblock
!		   if(ibundle==7.and.hchar=='b')print*,' bundle7b ',locationf,locationi,xme,blockvec_1(1+locationf)
		   
!		   if(blockvec_1(1+locationf)/=0.0)print*,ibundle,locationf,locationi,xme*blockvec_1(1+locationf)
		   
		   do ivec = 1,dimblock
               blockvec_2(ivec+locationi) = blockvec_2(ivec+locationi) + xme*blockvec_1(ivec+locationf)
		   end do
!  		  do ivec = 1,dimblock
 !              block_2(ivec,statei) = block_2(ivec,statei) + xme*block_1(ivec,statef)
 ! 		  end do
        end do  ! pjmp
     end do  ! njmp
!     else    ! REVERSE ORDER OF LOOPS
!........***WARNING*** SOME BUG WITH THE FOLLOWING.....

!     do pjmp = pjmpstart,pjmpend
!           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
!           psdf = p1b_fsd(pjmp)       ! final proton SD
!           phasep = p1b_phase(pjmp)   ! phase of proton jumps
!           a  = p1b_cop(pjmp) 
!           c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
!           do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
!              nsdi = n2b_isd(njmp)
!              nsdf = n2b_fsd(njmp)
!              bb  = n2b_cop(njmp)
!              dd  = n2b_dop(njmp)
!              phasen = n2b_phase(njmp)

!----------- FIND MATRIX ELEMENT -------------------------------------------
!              coplabel = cpnntriplet(bb,a)
!              doplabel = dpnntriplet(dd,c)
!              xme = hmatpnn(coplabel + doplabel)     ! get matrix element
!              xme = xme*phasep*phasen               ! multiply matrix element by jump phases
!              statei = nsdi + psdi                  ! initial state in combined basis
!              statef = nsdf + psdf                  ! final state in combined basis
!   		   locationi = (statei-1)*dimblock
!   		   locationf = (statef-1)*dimblock
		   
!   		   do ivec = 1,dimblock
!                  blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
!   		   end do
!	 		  do ivec = 1,dimblock
!	              block_2(ivec,statei) = block_2(ivec,statei) + xme*block_1(ivec,statef)
!	 		  end do
!           end do
!      end do  ! pjmp
!      end if
  end do
!$omp end parallel do
call bundle_clock(ibundle,'end')

  end do ! ibundle
  end if
  return
end subroutine applyhPNNbundled_block

!==========================================================

!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPPNbundled_block (hchar,startbundle,endbundle )
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use jump3body
  use interactions3body
  use opbundles
  use fragments
  use flagger
  use bmpi_mod
  use contigpointervectors, only :  p2b_1sd,p2b_2sd
 
  implicit none

  integer :: ibundle,startbundle,endbundle
   character(1) :: hchar,vchar

!------------------------------------------------------------
  integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer a,bb,c,dd
  integer(kind=8) :: coplabel,doplabel
  integer :: phasep,phasen
  real(kind=4)   xme
  integer(kind=8) :: statei, statef
  integer(kind=8) :: locationi,locationf
  
  integer(kind=4) num_threads
  integer(4) :: ivec
  
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  logical    :: launched(0:3)

!..............................................................
!

! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( hchar /= 'b' )then

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PPN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
call bundle_clock(ibundle,'sta')
	 
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

!$omp parallel do private( mythread,startp_thread,npjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,ivec,statei,statef)                   &
!$omp          shared(p2b_isd,p2b_fsd,p2b_phase,p2b_cop,p2b_dop)              &
!$omp          shared(ibundle,opbundle)                                       &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cppntriplet,dppntriplet)                   &
!$omp          shared(block_1,block_2)
!     mythread = omp_get_thread_num()


    do mythread = 0,numpthreads -1
     ! thread local vec2, reduce at end


     startp_thread = opbundle(ibundle)%startp_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     npjmps_thread = opbundle(ibundle)%startp_thread(mythread+1) - startp_thread
!---------   Forward direction ------------------------------
!      do pjmp = pjmpstart,pjmpend
     do pjmp = startp_thread + 1, startp_thread + npjmps_thread
        psdi = p2b_isd(pjmp)       ! initial proton slater determinant
        psdf = p2b_fsd(pjmp)       ! final proton SD
        phasep = p2b_phase(pjmp)   ! phase of proton jumps
        bb = p2b_cop(pjmp) 
        dd = p2b_dop(pjmp)
!           print*,' (F) ',bb,dd,pjmp,startp_thread,npjmps_thread,pjmpstart,pjmpend

!--------- LOOP OVER NEUTRON JUMPS --------------------------------
        do njmp = njmpstart,njmpend
!----------- FIND MATRIX ELEMENT ---------------------------

           phasen = n1b_phase(njmp)
           a = n1b_cop(njmp)
           c = n1b_dop(njmp)
           coplabel = cppntriplet(bb,a)
           doplabel = dppntriplet(dd,c)

           xme = hmatppn(coplabel + doplabel)   ! get matrix element
           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
           nsdi = n1b_isd(njmp)
           nsdf = n1b_fsd(njmp)
           statei = nsdi + psdi                ! initial state in combined basis
           statef = nsdf + psdf                ! final state in combined basis
		   locationi = (statei-1)*dimblock
		   locationf = (statef-1)*dimblock
		   do ivec = 1,dimblock
               blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
		   end do
!  		  do ivec = 1,dimblock
 !              block_2(ivec,statef) = block_2(ivec,statef) + xme*block_1(ivec,statei)
  !		  end do
       end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do
call bundle_clock(ibundle,'end')

  end do ! ibundle
else

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PPN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
	 call bundle_clock(ibundle,'sta')
	 
!	 if(opbundle(ibundle)%annexed)cycle
	 
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

!---- Backward direction using hermiticity ------------------- 

!$omp parallel do private( mythread,startn_thread,nnjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,ivec,statei,statef)                   &
!$omp          shared(p2b_isd,p2b_fsd,p2b_phase,p2b_cop,p2b_dop)              &
!$omp          shared(ibundle,opbundle)                                       &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cppntriplet,dppntriplet)                   &
!$omp          shared(block_1,block_2)
    do mythread = 0, numnthreads-1
     ! thread local vec2, reduce at end

     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread

!...... OPTION TO SWITCH ORDER OF LOOPS WHEN NO OpenMP.......

!     if(numnthreads > 1 .or. disableNoOMPloopswitch)then

     do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
        nsdi = n1b_isd(njmp)
        nsdf = n1b_fsd(njmp)
        a  = n1b_cop(njmp)
        c  = n1b_dop(njmp)
        phasen = n1b_phase(njmp)
        do pjmp = pjmpstart,pjmpend
           psdi = p2b_isd(pjmp)       ! initial proton slater determinant
           psdf = p2b_fsd(pjmp)       ! final proton SD
           phasep = p2b_phase(pjmp)   ! phase of proton jumps
           bb  = p2b_cop(pjmp) 
           dd  = p2b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS ---------------------------
!----------- FIND MATRIX ELEMENT ------------------------------

           coplabel = cppntriplet(bb,a)
           doplabel = dppntriplet(dd,c)
           xme = hmatppn(coplabel + doplabel)     ! get matrix element
           xme = xme*phasep*phasen               ! multiply matrix element by jump phases
           statei = nsdi + psdi                  ! initial state in combined basis
           statef = nsdf + psdf                  ! final state in combined basis
		   locationi = (statei-1)*dimblock
		   locationf = (statef-1)*dimblock
		   do ivec = 1,dimblock
               blockvec_2(ivec+locationi) = blockvec_2(ivec+locationi) + xme*blockvec_1(ivec+locationf)
		   end do
!  		  do ivec = 1,dimblock
!               block_2(ivec,statei) = block_2(ivec,statei) + xme*block_1(ivec,statef)
!  		  end do
        end do  ! pjmp
     end do  ! njmp

	 !........***WARNING*** SOME BUG WITH THE FOLLOWING.....
	 
!     else
!     do pjmp = pjmpstart,pjmpend
!           psdi = p2b_isd(pjmp)       ! initial proton slater determinant
!           psdf = p2b_fsd(pjmp)       ! final proton SD
!           phasep = p2b_phase(pjmp)   ! phase of proton jumps
!           bb  = p2b_cop(pjmp) 
!           dd = p2b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
!           do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
!              nsdi = n1b_isd(njmp)
!              nsdf = n1b_fsd(njmp)
!              a  = n1b_cop(njmp)
!              c  = n1b_dop(njmp)
!              phasen = n1b_phase(njmp)

!----------- FIND MATRIX ELEMENT -------------------------------------------
!              coplabel = cppntriplet(bb,a)
!              doplabel = dppntriplet(dd,c)
!              xme = hmatppn(coplabel + doplabel)     ! get matrix element
!              xme = xme*phasep*phasen               ! multiply matrix element by jump phases
!              statei = nsdi + psdi                  ! initial state in combined basis
!              statef = nsdf + psdf                  ! final state in combined basis
!   		   locationi = (statei-1)*dimblock
!   		   locationf = (statef-1)*dimblock
!   		   do ivec = 1,dimblock
!                  blockvec_2(ivec+locationf) = blockvec_2(ivec+locationf) + xme*blockvec_1(ivec+locationi)
!   		   end do
!	 		  do ivec = 1,dimblock
!	              block_2(ivec,statei) = block_2(ivec,statei) + xme*block_1(ivec,statef)
!	 		  end do
!           end do
!      end do  ! pjmp
!     end if
  end do
!$omp end parallel do
call bundle_clock(ibundle,'end')

  end do ! ibundle

  end if
  return
end subroutine applyhPPNbundled_block
!==========================================================



end module apply_ham_block
