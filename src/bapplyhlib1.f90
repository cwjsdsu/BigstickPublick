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
module contigpointervectors
	use precisions
	implicit none
! --------- NOTE: basestart, basestop stored in module fragments
	real (kind = lanc_prec), pointer, contiguous :: vecin(: ) 
	real (kind = lanc_prec), pointer, contiguous :: vecout(: ) 
    integer(kind=basis_prec), pointer, contiguous :: p2b_1sd(:), p2b_2sd(:)
    integer(kind=basis_prec), pointer, contiguous :: n2b_1sd(:), n2b_2sd(:)
    integer(kind=basis_prec), pointer, contiguous :: p3b_1sd(:), p3b_2sd(:)
    integer(kind=basis_prec), pointer, contiguous :: n3b_1sd(:), n3b_2sd(:)
!  real (kind = lanc_prec), pointer :: vecin(: ) 
!  real (kind = lanc_prec), pointer :: vecout(: ) 
!  integer(kind=basis_prec), pointer :: p2b_1sd(:), p2b_2sd(:)
!  integer(kind = basis_prec), pointer :: n2b_1sd(:), n2b_2sd(:)
!  integer(kind=basis_prec), pointer :: p3b_1sd(:), p3b_2sd(:)
!  integer(kind = basis_prec), pointer :: n3b_1sd(:), n3b_2sd(:)
	
end module contigpointervectors

module apply_ham
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

subroutine applyHbundled_g(vchar)
  use nodeinfo
  use flagger
  use flags3body
  use precisions
  use opbundles
  use fragments
  use interaction
  use basis
  use localvectors
  use bmpi_mod
  use lanczos_info
  use coupledmatrixelements,only:call_spe
  implicit none

  integer iprocs, procstart,procstop
  character(1) :: vchar  
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


  if(noisy0) print *, "Starting applyHbundled_g"
  do iprocs = procstart,procstop
     call proc_clock(iprocs,'sta')
	
!............SPE......................................
     if(call_spe)then
        call clocker('spe','sta')
        call procOP_clock(iprocs,'sta','SPE')
        if(noisy0) print *, "Starting applyHbundled_g SPE"
        call applySPEbundled_g(vchar,opbundlestart(iprocs), opbundleend(iprocs))
        call clocker('spe','end')
        call procOP_clock(iprocs,'end','SPE')
     endif
    if(applypntrace)call applyPNavgdiag(vchar,opbundlestart(iprocs), opbundleend(iprocs))

     if(.not.threebody)then
!........... PP .................

     if(noisy0) print *, "Starting applyHbundled_g PP"
     call clocker('ppo','sta')
     call procOP_clock(iprocs,'sta','PPO')
     call applyhPPbundled_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
     call procOP_clock(iprocs,'end','PPO')

     call clocker('ppo','end')
     call clocker('ppb','sta')
     call procOP_clock(iprocs,'sta','PPB')
     call applyhPPbundled_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))
     call procOP_clock(iprocs,'end','PPB')

     call clocker('ppb','end')
!............NN ..................
     call clocker('nno','sta')
     call procOP_clock(iprocs,'sta','NNO')

     if(noisy0) print *, "Starting applyHbundled_g NN"
     call applyhNNbundled_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
     call applyhNNbundled_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))
     call clocker('nno','end')
     call procOP_clock(iprocs,'end','NNO')


!........... PN....................................
     call clocker('one','sta')

     if(noisy0) print *, "Starting applyHbundled_g PN"
     if(useTR)then
        call applyhPNbundledTR_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNbundledTR_g(vchar,'h',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNbundledTR_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))
     else
	     call clocker('pno','sta')
	     call procOP_clock(iprocs,'sta','PNO')
        call applyhPNbundled_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNbundled_g(vchar,'h',opbundlestart(iprocs), opbundleend(iprocs))
		call clocker('pno','end')
        call procOP_clock(iprocs,'end','PNO')
     call procOP_clock(iprocs,'sta','PNB')

        call clocker('pnb','sta')
        call applyhPNbundled_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))
		call clocker('pnb','end')
        call procOP_clock(iprocs,'end','PNB')
		
     end if
     call clocker('one','end')

     else
!............PPP.......................................
    call clocker('ppp','sta')
     call procOP_clock(iprocs,'sta','PPP')

    call applyhPPPbundled_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
    call applyhPPPbundled_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))
    call clocker('ppp','end')
     call procOP_clock(iprocs,'end','PPP')

!............PPN.......................................
    call clocker('ppn','sta')
     call procOP_clock(iprocs,'sta','PPN')

    if(useTR)then
        call applyhPPNbundledTR_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPPNbundledTR_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))
    else
        call applyhPPNbundled_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPPNbundled_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))

    endif
     call procOP_clock(iprocs,'end','PPN')

    call clocker('ppn','end')
!............PNN.......................................
    call clocker('pnn','sta')
     call procOP_clock(iprocs,'sta','PNN')

    if(useTR)then
        call applyhPNNbundledTR_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNNbundledTR_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))

    else
        call applyhPNNbundled_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNNbundled_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))

    endif
    call clocker('pnn','end')
     call procOP_clock(iprocs,'end','PNN')

!...........NNN......................................
    call clocker('nnn','sta')
     call procOP_clock(iprocs,'sta','NNN')

    call applyhNNNbundled_g(vchar,'f',opbundlestart(iprocs), opbundleend(iprocs))
    call applyhNNNbundled_g(vchar,'b',opbundlestart(iprocs), opbundleend(iprocs))
    call clocker('nnn','end')
     call procOP_clock(iprocs,'end','NNN')


    end if  ! if threebody

    ! If using thread local output storage we have to reduce to vec2
    ! unconverted bundle types are still dumping in vec2 so we just add to it
    ! from vec2thread
    if(useVec2Thread) then
! the schedule here is static with a chunksize.  This make sense
! because each chunk will have predictable runtime
!$omp parallel do                      &
!$omp    private(i, j, tid, sum)          &
!$omp    shared(vec2threadflat, v2s, v2e)  &
!$omp    schedule(static, 1024)
      do i = v2s, v2e
        sum = vec2(i)
        !! do tid = 0, ompNumThreads-1
        !!    sum = sum + vec2thread(i, tid)
        !!    vec2thread(i, tid) = 0.0
        !! end do
        
        ! sum the contributions to vec2(i)
        do j = (i - v2s), vec2threadend, vec2threadchunk
            sum = sum + vec2threadflat(j)
            vec2threadflat(j) = 0.0;
        end do
        vec2(i) = real(sum,kind=lanc_prec)
      end do
!$omp end parallel do
    end if

    call proc_clock(iprocs,'end')

!........... NOW REDUCE....................
!            WHEN LANCZOS VECTOR BROKEN, 
!            THIS WILL GET MORE COMPLICATED
   if(noisy0) print *, "applyhbundled_g: Doing Reduce"
   if(useNewReorthog) then
       if(vchar == 'n') then
         ! 'n' corresponds with vec2 here because we are looking at the output vector
          ! call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
          ! Do reduce only onto root node.  We will be sending this data to 
          ! the slices from the isfragroot nodes (rank=0 in fcomm1, fcomm2, hcomm),
          ! so other nodes don't need it.
          call BMPI_REDUCE(vec2, size(vec2), MPI_SUM, 0, fcomm2, ierr) ! in place
       else
          ! call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, fcomm1, ierr) ! in place
          call BMPI_REDUCE(vec1, size(vec1), MPI_SUM, 0, fcomm1, ierr) ! in place
       end if
   else
       if(nproc > 1 .and. vchar == 'n')then
          call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
       endif
       if(nproc > 1 .and. vchar == 'r')then
          call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, fcomm1, ierr) ! in place
       end if
   end if

  end do  ! iprocs
  if(noisy0) print *, "applyhbundled_g: Returning"
  return

end subroutine applyHbundled_g
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
subroutine applyhPPbundled_orig (vchar,hchar,startbundle,endbundle )

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
  use contigpointervectors, only : vecin,vecout, p2b_1sd,p2b_2sd
  implicit none

  logical :: pinfo
  integer :: ibundle
  character(1) :: hchar,vchar
  integer :: startbundle,endbundle

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=4)   xme, prod
  integer(kind=8) :: statei, statef,nsd,psdi,psdf
  integer(kind=basis_prec) :: statefoff, stateistart, stateistop
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(kind=8) :: vs
  integer(4) :: mythread
  real(kind=lanc_prec), pointer :: voutp(:)


!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
  select case(vchar)
     case ('n')

        vecin  => vec1
        vecout => vec2
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


     case ('r')

        vecin  => vec2
        vecout => vec1
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
        if( hchar == 'b')then  !reversed from above
           p2b_1sd => p2b_isd
           p2b_2sd => p2b_fsd
        else
           p2b_1sd => p2b_fsd
           p2b_2sd => p2b_isd
        endif
     case default
        print *, "bad vchar=", vchar
        stop
  end select

!  do ibundle = endbundle,startbundle,-1
  do ibundle = startbundle,endbundle

     if(opbundle(ibundle)%optype /= 'PP')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
	 
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
!$omp parallel private(vs,xjmp, Xoplabel, xme, num_threads, mythread,nsd, pinfo)    &
!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!$omp          private(statefoff, stateistart, stateistop, prod)  &
!$omp          private(voutp) &
!$omp          firstprivate(cstride, ncstates, xjmpstart,xjmpend)  &
!$omp          shared(vecin, vecout)  &
!$omp          shared(vec2threadchunkm1) &
!$omp          shared(p2b_op, p2b_1sd, p2b_2sd, p2b_phase, hmatpp)
  num_threads =  omp_get_num_threads()
  mythread = omp_get_thread_num()
  ! thread local vec2, reduce at end
  if(useVec2Thread) then
    !! voutp(v2s:v2e) => vec2thread(:, mythread)
    vs = mythread * vec2threadchunk;
    voutp(v2s:v2e) => vec2threadflat(vs : vs + vec2threadchunkm1)
  else
    voutp(v2s:v2e) => vecout
  end if

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
          voutp(statef) = voutp(statef) + xme*vecin(statei)
      end do  ! xjmp
   end do  ! csd
end if
!$omp end parallel
!--------------OR DO HERMITIAN/BACKWARDS APPLICATION----------
  call bundle_clock(ibundle,'end')

  end do ! ibundle
  return
end subroutine applyhPPbundled_orig

subroutine applyhPPbundled_thread (vchar,hchar,startbundle,endbundle )
  use nodeinfo
  use localvectors
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
  use contigpointervectors, only : vecin,vecout, p2b_1sd,p2b_2sd
  implicit none

  integer :: ibundle
  character(1) :: hchar,vchar
  integer :: startbundle,endbundle

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=4)   xme, prod
  integer(kind=8) :: statei, statef,nsd,psdi,psdf
  integer(kind=basis_prec) :: statefoff, stateistart, stateistop
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(kind=8) :: vs
  integer(4) :: mythread
  real(kind=lanc_prec), pointer :: voutp(:)


!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
  select case(vchar)
     case ('n')

        vecin  => vec1
        vecout => vec2
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


     case ('r')

        vecin  => vec2
        vecout => vec1
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
        if( hchar == 'b')then  !reversed from above
           p2b_1sd => p2b_isd
           p2b_2sd => p2b_fsd
        else
           p2b_1sd => p2b_fsd
           p2b_2sd => p2b_isd
        endif
     case default
        print *, "bad vchar=", vchar
        stop
  end select
! We use dynamic scheduling here because the iterations have
! extemely variable runtime.   Sometimes we just cycle.
! chunk size is one bundle for maximum balance. 
!$omp parallel do private(vs, xjmp, Xoplabel, xme, nsd)    &
!$omp          private(mythread, num_threads) &
!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!$omp          private(psdi, psdf)  &
!$omp          private(statefoff, stateistart, stateistop, prod)  &
!$omp          private(voutp) &
!$omp          private(csdstart, csdend, xjmpstart, xjmpend, cstride) &
!$omp          private(ncstates) &
!$omp          shared(opbundle)  &
!$omp          shared(vecin, vecout)  &
!$omp          shared(vec2threadchunkm1) &
!$omp          shared(p2b_op, p2b_1sd, p2b_2sd, p2b_phase, hmatpp) &
!$omp     schedule(dynamic,1)
   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'PP')cycle
      if(opbundle(ibundle)%hchar /= hchar )cycle
! 	  if(opbundle(ibundle)%annexed)cycle
	  
	  call bundle_clock(ibundle,'sta')
	  
      num_threads =  omp_get_num_threads()
      mythread = omp_get_thread_num()

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
      csdstart = opbundle(ibundle)%nxstart
      csdend   = opbundle(ibundle)%nxend
      xjmpstart = opbundle(ibundle)%pxstart
      xjmpend   = opbundle(ibundle)%pxend
      cstride   = opbundle(ibundle)%cstride

      ncstates = (csdend +cstride -csdstart)/cstride

!--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
      ! thread local vec2, reduce at end
      ! voutp(v2s:v2e) => vec2thread(:, mythread)
      vs = mythread * vec2threadchunk
      voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)

      chunk = ncstates ! instead of / num_threads
      istart = 1
      iend = ncstates
      csd_index = csdstart + (istart - 1)*cstride - cstride
      if(istart > iend) cycle

!......... THERE ARE TWO VERSIONS.....
!          1st way: store in jumps index to PP matrix element;
!          this has faster setup
!          2nd way: store PP matrix elements directly;
!          slower set up, but on MPI nodes reduced memory load

      if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO PP MATRIX ELEMENTS
         ! KSM - Test both ways!!
         if(.false.) then
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
                  voutp(statef) = voutp(statef) + xme*vecin(statei)
               end do  ! xjmp
            end do  ! csd
         else
            ! Trial speedup
            ! result so for is mysteriously slower
            do xjmp = xjmpstart,xjmpend
!------------------ FETCH MATRIX ELEMENT...............................
               Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
               xme = hmatpp(Xoplabel)
!--------- GET PHASE.........................................
               xme = xme*p2b_phase(xjmp)
               psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
               psdf = p2b_2sd(xjmp)   ! KSM: final P SD

!============================================
               csd_index = csdstart + (istart - 1)*cstride - cstride
               do csd = istart, iend
                  csd_index = csd_index + cstride
                  nsd = nstart(csd_index)
                  statei = psdi + nsd 
                  statef = psdf + nsd
                  voutp(statef) = voutp(statef) + xme*vecin(statei)
               end do ! csd

!============================================
! KSM speedup.  Everything out of the loop that can go
!!               statefoff = psdf - psdi;
!!               stateistart = csd_index + psdi + cstride;
!!               stateistop = stateistart + (iend - istart)*cstride
!!               ! really iterating over chunk of neutrons.   Input
!!               ! and output have the same state offset due to neutrons
!!               ! because they are unaffected by PP
!!               ! There is no way to generate overlap
!!               do statei=stateistart, stateistop, cstride
!!                  statef = statei + statefoff;
!!                  prod = xme * vecin(statei)
!!                  voutp(statef) = voutp(statef) + prod
!!               end do
            end do  ! xjmp
         end if
      else    ! USE STORED PP MATRIX ELEMENTS DIRECTLY........
         do csd = istart, iend
            csd_index = csd_index + cstride
            nsd = nstart(csd_index)

!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
            do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
               xme = p2b_me(xjmp)
!                 Xoplabel = p2b_op(xjmp)
!                 xme = hmatpp(Xoplabel)
!--------- GET PHASE.........................................
               xme = xme*p2b_phase(xjmp)
!---------- GET INITIAL, FINAL SDs and place in basis..............
               statei = p2b_1sd(xjmp)+ nsd !csd_index
               statef = p2b_2sd(xjmp)+nsd !csd_index
               voutp(statef) = voutp(statef) + xme*vecin(statei)
            end do  ! xjmp
         end do  ! csd
      end if
!--------------OR DO HERMITIAN/BACKWARDS APPLICATION----------
     call bundle_clock(ibundle,'end')

   end do ! ibundle
!$omp end parallel do
   return
end subroutine applyhPPbundled_thread

subroutine applyhPPbundled_g (vchar,hchar,startbundle,endbundle )
   use localvectors
   implicit none
   character(1), intent(in) :: hchar,vchar
   integer, intent(in) :: startbundle,endbundle

   if(useVec2Thread) then
      call applyhPPbundled_thread(vchar, hchar, startbundle, endbundle)
   else
      call applyhPPbundled_orig(vchar, hchar, startbundle, endbundle)
   end if
end subroutine applyhPPbundled_g

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
subroutine applyhNNbundled_orig (vchar,hchar,startbundle,endbundle )
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
   use contigpointervectors, only : vecin,vecout, n2b_1sd,n2b_2sd
   implicit none

   ! arguments
   character(1),intent(in) :: hchar,vchar
   integer,intent(in) :: startbundle,endbundle

! --------- NOTE: basestart, basestop stored in module fragments

!------------------------------------------------------------

   integer :: ibundle
   integer(kind=8) csdstart, csdend,csd, csd_index,cstride,ncstates, pstride
   integer(kind=8) xjmp,xjmpstart,xjmpend
   integer(kind=8):: Xoplabel
   real(kind=4)   xme
   integer(kind=8) :: statei, statef,psd,nsdi,nsdf,statef_start,statef_end
   integer(kind=8) :: istart, iend, chunk
   integer(kind=8) :: vs

!-------- OpenMP functions ---------------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: num_threads
   integer :: mythread
   real(kind=lanc_prec), pointer :: voutp(:)

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
   select case(vchar)
      case ('n')
         vecin  => vec1
         vecout => vec2
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
         if( hchar == 'f' )then
            n2b_1sd => n2b_isd
            n2b_2sd => n2b_fsd
         else
            n2b_1sd => n2b_fsd
            n2b_2sd => n2b_isd
         endif

      case ('r')
         vecin  => vec2
         vecout => vec1
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
         if( hchar == 'b')then  !reversed from above
            n2b_1sd => n2b_isd
            n2b_2sd => n2b_fsd
         else
            n2b_1sd => n2b_fsd
            n2b_2sd => n2b_isd
         endif
      case default
         print *, "bad vchar=", vchar
         stop
      end select

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
!$omp parallel private(vs, xjmp, Xoplabel, xme, num_threads, mythread,psd)         &
!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!$omp          private(voutp) &
!$omp          firstprivate(cstride, ncstates, xjmpstart, xjmpend)  &
!$omp          shared(vecin, vecout)  &
!$omp          shared(vec2threadchunkm1) &
!$omp          shared(n2b_op, n2b_1sd, n2b_2sd, n2b_phase, hmatnn)
         num_threads =  omp_get_num_threads()
         mythread = omp_get_thread_num()
         ! thread local vec2, reduce at end
         if(useVec2Thread) then
            !! voutp(v2s:v2e) => vec2thread(:, mythread)
            vs = mythread * vec2threadchunk
            voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
         else
            voutp(v2s:v2e) => vecout
         end if
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
                        voutp(statef) = voutp(statef) +  xme*vecin(statei)
                     end do  ! xjmp
                  end do  ! csd

			  end if

!$omp end parallel
call bundle_clock(ibundle,'end')

   end do ! ibundle
   return
end subroutine applyhNNbundled_orig

! This version has each thread write into a different output
! buffer that is reduced with other threads at the end.
subroutine applyhNNbundled_thread (vchar,hchar,startbundle,endbundle )
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
   use contigpointervectors, only : vecin,vecout, n2b_1sd,n2b_2sd
   implicit none

   ! arguments
   character(1),intent(in) :: hchar,vchar
   integer,intent(in) :: startbundle,endbundle

! --------- NOTE: basestart, basestop stored in module fragments

!------------------------------------------------------------

   integer :: ibundle
   integer(kind=8) csdstart, csdend,csd, csd_index,cstride,ncstates, pstride
   integer(kind=8) xjmp,xjmpstart,xjmpend
   integer(kind=8):: Xoplabel
   real(kind=4)   xme
   integer(kind=8) :: statei, statef,psd,nsdi,nsdf,statef_start,statef_end
   integer(kind=8) :: istart, iend
   integer(kind=8) :: vs

!-------- OpenMP functions ---------------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: mythread
   real(kind=lanc_prec), pointer :: voutp(:)

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
   select case(vchar)
   case ('n')
      vecin  => vec1
      vecout => vec2
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
      if( hchar == 'f' )then
         n2b_1sd => n2b_isd
         n2b_2sd => n2b_fsd
      else
         n2b_1sd => n2b_fsd
         n2b_2sd => n2b_isd
      endif

   case ('r')
      vecin  => vec2
      vecout => vec1
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
      if( hchar == 'b')then  !reversed from above
         n2b_1sd => n2b_isd
         n2b_2sd => n2b_fsd
      else
         n2b_1sd => n2b_fsd
         n2b_2sd => n2b_isd
      endif
   case default
      print *, "bad vchar=", vchar
      stop
   end select

! firstprivate gives each thread its own copy, but initializes it
! this is   better than shared for read-only vars

! schedule is dynamic because runtime per step is quite variable.
!$omp parallel do                                         &
!$omp         private(vs, ibundle, mythread, voutp)       &
!$omp         private(csdstart, csdend, cstride)          &
!$omp         private(xjmpstart, xjmpend, ncstates)       &
!$omp         private(istart, iend, csd_index)            &
!$omp         private(pstride)                            &
!$omp         private(xjmp, Xoplabel, xme)                &
!$omp         private(statei, statef)                     &
!$omp         private(statef_start, statef_end)           &
!$omp         firstprivate(startbundle, endbundle)        &
!$omp         shared(vecin)                               &
!$omp         shared(vec2threadchunkm1)                   &
!$omp         shared(opbundle,pstart,hmatnn)              &
!$omp         shared(n2b_op, n2b_phase,n2b_1sd,n2b_2sd)    &
!$omp         shared(storeXXmesjumps)                     &
!$omp   schedule(dynamic, 1)
   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'NN')cycle
      if(opbundle(ibundle)%hchar /= hchar )cycle
! 	  if(opbundle(ibundle)%annexed)cycle
	  
	  call bundle_clock(ibundle,'sta')
	  
!....... Thread specific setup
      mythread = omp_get_thread_num()
      ! thread local vec2, reduce at end
      if(useVec2Thread) then
         ! voutp(v2s:v2e) => vec2thread(:, mythread)
         vs = mythread * vec2threadchunk
         voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
      else
         voutp(v2s:v2e) => vecout
      end if
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

      istart = 1
      iend = ncstates
      csd_index = csdstart + (istart - 1)*cstride - cstride
      if(istart > ncstates) cycle
!..... THIS FOLLOWING IS TO TRY TO FIND AN OPTIMAL ORDERING OF LOOPS
      pstride = pstridecut +1 
      if(iend-istart > 0) pstride = pstart(csdstart+cstride)-pstart(csdstart)

!............ THERE ARE TWO VERSIONS.....
!             1st way: store in jumps index to NN matrix element;
!             this has faster setup
!             2nd way: store NN matrix elements directly;
!             slower set up, but on MPI nodes reduced memory load
      if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO NN MATRIX ELEMENTS
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
         csd_index = pstart(csdstart)

         do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
            Xoplabel = n2b_op(xjmp)
!		            if(Xoplabel==0)print*,' ZERO LABEL ',iproc,ibundle,xjmp
            xme = hmatnn(Xoplabel)
!--------- GET PHASE.........................................
            xme = xme*n2b_phase(xjmp)
            statei = n2b_1sd(xjmp) + csd_index
            statef_start = n2b_2sd(xjmp) +csd_index
!---------- GET INITIAL, FINAL SDs and place in basis..............
            statef_end = statef_start+pstride*(iend-istart)
            do statef = statef_start,statef_end,pstride
               voutp(statef) = voutp(statef) +  xme*vecin(statei)
               statei = statei+pstride
            end do  !
         end do  ! xjmp
      else
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
         csd_index = pstart(csdstart)
         do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
            xme = n2b_me(xjmp)
!--------- GET PHASE.........................................
            xme = xme*n2b_phase(xjmp)
            statei = n2b_1sd(xjmp) + csd_index
            statef_start = n2b_2sd(xjmp) +csd_index
!---------- GET INITIAL, FINAL SDs and place in basis..............
            statef_end = statef_start+pstride*(iend-istart)
            do statef = statef_start,statef_end,pstride
               voutp(statef) = voutp(statef) +  xme*vecin(statei)
               statei = statei+pstride
            end do  !
         end do  ! xjmp
      end if
	  call bundle_clock(ibundle,'end')
   end do ! ibundle
!$omp end parallel do
   return
end subroutine applyhNNbundled_thread
subroutine applyhNNbundled_g (vchar,hchar,startbundle,endbundle )
   use localvectors
   implicit none
   ! arguments
   character(1),intent(in) :: hchar,vchar
   integer,intent(in) :: startbundle,endbundle

   if(useVec2Thread) then
      call applyhNNbundled_thread(vchar,hchar,startbundle,endbundle )
   else
      call applyhNNbundled_orig(vchar,hchar,startbundle,endbundle )
   end if
end subroutine applyhNNbundled_g
!=================================================================
!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPNbundled_orig (vchar,hchar,startbundle,endbundle )
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
  use contigpointervectors, only : vecin,vecout
  implicit none

  integer :: ibundle,startbundle,endbundle
  character(1) :: hchar,vchar

!------------------------------------------------------------
  integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer :: a,b,c,d
  integer(kind=8) :: coplabel,doplabel
  integer :: phasep,phasen
  real(kind=4) ::   xme
  integer(kind=basis_prec) :: statei, statef
  integer(kind=4) num_threads

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  integer(kind=basis_prec) :: vs

  real(kind=lanc_prec), pointer :: voutp(:)


  if(applyXXonly)return    ! don't do PN
!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
  select case(vchar)
     case ('n')

        vecin  => vec1
        vecout => vec2


     case ('r')

        vecin  => vec2
        vecout => vec1

     case default
        print *, "bad vchar=", vchar
        stop
  end select


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
  if( (vchar == 'n' .and. hchar /= 'b') .or. & 
       (vchar == 'r' .and. hchar == 'b') )then

!$omp parallel do private(vs, mythread,startp_thread,npjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)   &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          private(voutp) &
!$omp          shared(vec2threadchunkm1) &
!$omp          firstprivate(njmpstart, njmpend)       &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)              &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cpnpair,dpnpair,vecin,vecout)
     do mythread = 0,numpthreads -1
        ! thread local vec2, reduce at end
        if(useVec2Thread) then
          ! voutp(v2s:v2e) => vec2thread(:, mythread)
          vs = mythread * vec2threadchunk
          voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
        else
          voutp(v2s:v2e) => vecout
        end if

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
           voutp(statef) = voutp(statef) + xme*vecin(statei)
        end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do
else
!---- Backward direction using hermiticity ------------------- 

!$omp parallel do private(vs, mythread,startn_thread,nnjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)   &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          private(voutp) &
!$omp          firstprivate(pjmpstart, pjmpend)                     &
!$omp          shared(vec2threadchunkm1)                            &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)    &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)    &
!$omp          shared(cpnpair,dpnpair,vecin,vecout)
  do mythread = 0, numnthreads-1
     ! thread local vec2, reduce at end
     if(useVec2Thread) then
       ! voutp(v2s:v2e) => vec2thread(:, mythread)
       vs = mythread * vec2threadchunk
       voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
     else
       voutp(v2s:v2e) => vecout
     end if

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
           voutp(statei) = voutp(statei) + xme*vecin(statef)
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
	           voutp(statei) = voutp(statei) + xme*vecin(statef)
	        end do  ! njmp
	     end do  ! pjmp        

!     do pjmp = pjmpstart,pjmpend
!           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
!           psdf = p1b_fsd(pjmp)       ! final proton SD
!           phasep = p1b_phase(pjmp)   ! phase of proton jumps
!           a  = p1b_cop(pjmp) 
!           c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
!           do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
!              nsdi = n1b_isd(njmp)
!              nsdf = n1b_fsd(njmp)
 !             b  = n1b_cop(njmp)
 !             d  = n1b_dop(njmp)
 !             phasen = n1b_phase(njmp)

!----------- FIND MATRIX ELEMENT -------------------------------------------
!              coplabel = cpnpair(b,a)
!              doplabel = dpnpair(d,c)
!              xme = hmatpn(coplabel + doplabel)     ! get matrix element
!              xme = xme*phasep*phasen               ! multiply matrix element by jump phases
!              statei = nsdi + psdi                  ! initial state in combined basis
!              statef = nsdf + psdf                  ! final state in combined basis
!              voutp(statei) = voutp(statei) + xme*vecin(statef)
!           end do
!      end do  ! pjmp


     end if
  end do
!$omp end parallel do


  end if
  call bundle_clock(ibundle,'end')

  end do  ! ibundle

  return
end subroutine applyhPNbundled_orig

subroutine applyhPNbundled_thread (vchar,hchar,startbundle,endbundle )
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
   use contigpointervectors, only : vecin,vecout
   implicit none

!-- Arguments -----------------------------------------------
   integer, intent(in) :: startbundle,endbundle
   character(1),intent(in) :: hchar,vchar

!------------------------------------------------------------
   integer :: ibundle
   integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

   integer(kind=8) pjmp,pjmpstart,pjmpend
   integer(kind=8) njmp,njmpstart,njmpend
   integer :: a,b,c,d
   integer(kind=8) :: coplabel,doplabel
   integer :: phasep,phasen
   real(kind=4) ::   xme
   integer(kind=basis_prec) :: statei, statef
   integer(kind=4) num_threads
   integer(kind=basis_prec) :: vs

!-------- OpenMP functions/vars -----------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: mythread,numpthreads,numnthreads
   integer(8) :: startp_thread, npjmps_thread
   integer(8) :: startn_thread, nnjmps_thread

   real(kind=lanc_prec), pointer :: voutp(:)


!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
      select case(vchar)
      case ('n')
         vecin  => vec1
         vecout => vec2

      case ('r')
         vecin  => vec2
         vecout => vec1

      case default
         print *, "bad vchar=", vchar
         stop
      end select

! schedule is dynamic because runtime per step is variable
!$omp parallel do                                               &
!$omp     private(ibundle)                                      &
!$omp     private(vs, mythread, voutp)                          &
!$omp     private(pjmpstart, pjmpend, njmpstart, njmpend)       &
!$omp     private(numpthreads, numnthreads)                     &
!$omp     private(startp_thread, npjmps_thread)                 &
!$omp     private(pjmp, njmp, psdi, psdf)                       &
!$omp     private(a, b, c, d)                                   &
!$omp     private(phasep, phasen, coplabel, doplabel)           &
!$omp     private(xme, nsdi, nsdf, statei, statef)              &
!$omp     shared(useVec2Thread)                                 &
!$omp     shared(vec2threadchunkm1)                             &
!$omp     shared(vchar, hchar)                                  &
!$omp     schedule(dynamic,1)
      do ibundle = startbundle,endbundle
         if(opbundle(ibundle)%optype /= 'PN')cycle
         if(opbundle(ibundle)%hchar /= hchar )cycle
!		 if(opbundle(ibundle)%annexed)cycle
		 
	     call bundle_clock(ibundle,'sta')
		 
         mythread = omp_get_thread_num()
         ! thread local vec2, reduce at end of applyHbundled_g
         if(useVec2Thread) then
            !! voutp(v2s:v2e) => vec2thread(:, mythread)
            vs = mythread * vec2threadchunk
            voutp(v2s:v2e)=>vec2threadflat(vs: vs+vec2threadchunkm1)
         else
            voutp(v2s:v2e) => vecout
         end if

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
         pjmpstart = opbundle(ibundle)%pxstart
         pjmpend   = opbundle(ibundle)%pxend
         njmpstart = opbundle(ibundle)%nxstart
         njmpend   = opbundle(ibundle)%nxend
         numpthreads = opbundle(ibundle)%numpthreads
         numnthreads = opbundle(ibundle)%numnthreads

!  print*,numberthreads,' threads '
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
         if( (vchar == 'n' .and. hchar /= 'b') .or. & 
             (vchar == 'r' .and. hchar == 'b') )then

            !  starting position for proton 1-body jumps for this thread
            ! bit of a hack here to keep this routine similar to _orig version
            ! we take the complete thread range and ignore the divisions
            startp_thread = opbundle(ibundle)%startp_thread(0)     
            npjmps_thread = opbundle(ibundle)%startp_thread(numpthreads) - startp_thread

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
                  voutp(statef) = voutp(statef) + xme*vecin(statei)
               end do  ! njmp
            end do  ! pjmp        
         else
!---- Backward direction using hermiticity ------------------- 

         !  starting position for proton 1-body jumps for this thread
         ! bit of a hack here to keep this routine similar to _orig version
         ! we take the complete thread range and ignore the divisions
         startn_thread = opbundle(ibundle)%startn_thread(0)     
         nnjmps_thread = opbundle(ibundle)%startn_thread(numnthreads) - startn_thread

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
!		               if(coplabel+doplabel==0)then
!			               print*,iproc,'(b)',njmp,coplabel,doplabel
!			               print*,a,b,c,d,ibundle
!		               end if		   
                  xme = hmatpn(coplabel + doplabel)     ! get matrix element
                  xme = xme*phasep*phasen               ! multiply matrix element by jump phases
                  statei = nsdi + psdi                  ! initial state in combined basis
                  statef = nsdf + psdf                  ! final state in combined basis
                  voutp(statei) = voutp(statei) + xme*vecin(statef)
               end do  ! pjmp
            end do  ! njmp

         else

            do pjmp = pjmpstart,pjmpend
               psdi = p1b_isd(pjmp)       ! initial proton slater determinant
               psdf = p1b_fsd(pjmp)       ! final proton SD
               phasep = p1b_phase(pjmp)   ! phase of proton jumps
               a  = p1b_cop(pjmp) 
               c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
               do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
                  nsdi = n1b_isd(njmp)
                  nsdf = n1b_fsd(njmp)
                  b  = n1b_cop(njmp)
                  d  = n1b_dop(njmp)
                  phasen = n1b_phase(njmp)

!----------- FIND MATRIX ELEMENT -------------------------------------------
                  coplabel = cpnpair(b,a)
                  doplabel = dpnpair(d,c)
!   		            if(coplabel+doplabel==0)then
!   			            print*,iproc,'(c)',njmp,coplabel,doplabel
!   		            end if
                  xme = hmatpn(coplabel + doplabel)     ! get matrix element
                  xme = xme*phasep*phasen               ! multiply matrix element by jump phases
                  statei = nsdi + psdi                  ! initial state in combined basis
                  statef = nsdf + psdf                  ! final state in combined basis
                  voutp(statei) = voutp(statei) + xme*vecin(statef)
               end do
            end do  ! pjmp
         end if
      end if
	  call bundle_clock(ibundle,'end')
	  
   end do  ! ibundle
!$omp end parallel do

   return
end subroutine applyhPNbundled_thread

subroutine applyhPNbundled_g (vchar,hchar,startbundle,endbundle )
   use localvectors
   implicit none
   integer, intent(in) :: startbundle,endbundle
   character(1),intent(in) :: hchar,vchar

   if(useVec2Thread) then
      call applyhPNbundled_thread (vchar,hchar,startbundle,endbundle )
   else
      call applyhPNbundled_orig (vchar,hchar,startbundle,endbundle )
   end if
   return
end subroutine applyhPNbundled_g

!==========================================================
!=================================================================
!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPNbundledTR_g (vchar,hchar,startbundle,endbundle )
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use interaction
  use opbundles
  use fragments
  use bmpi_mod
  use contigpointervectors, only : vecin,vecout
  implicit none

  integer :: ibundle,startbundle,endbundle
   character(1) :: hchar,vchar

!------------------------------------------------------------
  integer(8) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer a,b,c,d
  integer(kind=8) :: coplabel,doplabel
  integer :: phasep,phasen
  real(kind=4)   xme
  integer(kind=8) :: statei, statef
  integer(kind=4) num_threads
  integer(kind=basis_prec) :: vs

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  logical    :: launched(0:3)

  real(kind=lanc_prec), pointer :: voutp(:)

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
  select case(vchar)
     case ('n')

        vecin  => vec1
        vecout => vec2


     case ('r')

        vecin  => vec2
        vecout => vec1

     case default
        print *, "bad vchar=", vchar
        stop
  end select

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
	 
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads

!  print*,numberthreads,' threads '
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( (vchar == 'n' .and. opbundle(ibundle)%hchar /= 'b') .or. & 
       (vchar == 'r' .and. opbundle(ibundle)%hchar == 'b') )then

!$omp parallel do private(vs, mythread,startp_thread,npjmps_thread)                  &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)          &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          private(voutp) &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)              &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(vec2threadchunkm1)                                      &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cpnpair,dpnpair,vecin,vecout)
!     mythread = omp_get_thread_num()
     do mythread = 0,numpthreads -1
        ! thread local vec2, reduce at end
        if(useVec2Thread) then
          !! voutp(v2s:v2e) => vec2thread(:, mythread)
          vs = mythread * vec2threadchunk
          voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
        else
          voutp(v2s:v2e) => vecout
        end if

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
              b = n1b_cop(njmp)
              d = n1b_dop(njmp)
              phasen = n1b_phase(njmp) *phasepnpair(b,a)*phasepnpair(d,c)
              coplabel = cpnpair(b,a)
              doplabel = dpnpair(d,c)
              xme = hmatpn(coplabel + doplabel)   ! get matrix element
              xme = xme*phasep*phasen             ! multiply matrix element by jump phases
   !           write(87,*)a,b,c,d,xme, phasepnpair(b,a),phasepnpair(d,c)
              nsdi = n1b_isd(njmp)
              nsdf = n1b_fsd(njmp)
              statei = nsdi + psdi                ! initial state in combined basis
              statef = nsdf + psdf                ! final state in combined basis
              voutp(statef) = voutp(statef) + xme*vecin(statei)
           end do  ! njmp
        end do  ! pjmp        
     end do
!$omp end parallel do
else
!---- Backward direction using hermiticity ------------------- 

!$omp parallel do private(vs, mythread,startn_thread,nnjmps_thread)           &
!$omp          private(voutp) &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)   &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)    &
!$omp          shared(ibundle,opbundle)                             &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)    &
!$omp          shared(vec2threadchunkm1)                            &
!$omp          shared(cpnpair,dpnpair,vecin,vecout)
!     mythread = omp_get_thread_num()
        do mythread = 0, numnthreads-1
        ! thread local vec2, reduce at end
        if(useVec2Thread) then
          !! voutp(v2s:v2e) => vec2thread(:, mythread)
          vs = mythread * vec2threadchunk
          voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
        else
          voutp(v2s:v2e) => vecout
        end if
        startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
        nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread
        do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
           nsdi = n1b_isd(njmp)
           nsdf = n1b_fsd(njmp)
           b  = n1b_cop(njmp)
           d  = n1b_dop(njmp)
           phasen = n1b_phase(njmp)
           do pjmp = pjmpstart,pjmpend
              psdi = p1b_isd(pjmp)       ! initial proton slater determinant
              psdf = p1b_fsd(pjmp)       ! final proton SD
              a  = p1b_cop(pjmp) 
              c  = p1b_dop(pjmp)
              phasep = p1b_phase(pjmp)*phasepnpair(b,a)*phasepnpair(d,c)   ! phase of proton jumps

   !--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
   !----------- FIND MATRIX ELEMENT -------------------------------------------
              coplabel = cpnpair(b,a)
              doplabel = dpnpair(d,c)
              xme = hmatpn(coplabel + doplabel)     ! get matrix element
              xme = xme*phasep*phasen               ! multiply matrix element by jump phases
              statei = nsdi + psdi                  ! initial state in combined basis
              statef = nsdf + psdf                  ! final state in combined basis
              voutp(statei) = voutp(statei) + xme*vecin(statef)
           end do  ! pjmp
        end do  ! njmp
     end do
!$omp end parallel do

  end if

  end do ! ibundle
  return
end subroutine applyhPNbundledTR_g

!==========================================================

!==================================================================
!
!  subroutine applyspes
!
!  applies single-particle energies -- purely diagonal
!
!====================================================================
subroutine applySPEbundled_orig (vchar,startbundle,endbundle )
   use basis
   use sectors
   use diagh
   use precisions
   use lanczos_info
   use localvectors
   use nodeinfo
   use system_parameters
   use opbundles
   use fragments
   use contigpointervectors, only : vecin,vecout
   implicit none

   integer :: ibundle,startbundle,endbundle
   character(1) :: vchar
!------------------------------------------------------------

   integer(kind=8) nsdstart,nsdend,psdstart,psdend
   integer(kind=8) xjmp,xjmpstart,xjmpend

   integer(kind=8) :: ibasis
   integer(kind=4) num_threads
!-------- OpenMP functions ---------------------------------
   integer(kind=4) :: omp_get_thread_num, omp_get_num_threads

   real(4)    ::  pspe,nspe
   integer(kind=8) :: ip,in

   if(vchar == 'n')then
      vecin  => vec1
      vecout => vec2
   else
      vecin  => vec2
      vecout => vec1
   end if


   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'SPE')cycle
! 	  if(opbundle(ibundle)%annexed)cycle
	  
	  call bundle_clock(ibundle,'sta')

      psdstart = opbundle(ibundle)%pxstart
      psdend   = opbundle(ibundle)%pxend
      nsdstart = opbundle(ibundle)%nxstart
      nsdend   = opbundle(ibundle)%nxend

! note: first private says that each such private copy is initialized with 
! the value before the parallel pragma
!................ LOOP OVER PROTON SDs in that sector..........
!$omp parallel do private(ip,in,pspe,nspe,ibasis) & 
!$omp  firstprivate(nsdstart, nsdend) &
!$omp  shared(psdstart,psdend,vecin,vecout) & 
!$omp  shared(pstart,nstart,pspe_h,nspe_h)
      do ip = psdstart,psdend
           pspe = pspe_h(ip)   ! the proton contribution to the single-particle energies
           ibasis = pstart(ip) + nstart(nsdstart)
!............... LOOP OVER NEUTRON SDS in sector jsc........................
           do in = nsdstart,nsdend
              nspe = nspe_h(in)  ! neutron contributions to s.p.e.
              vecout(ibasis) = vecout(ibasis) + vecin(ibasis)*( pspe + nspe )  ! add spes
              ibasis = ibasis+1   ! neutron SDs are contiguous
           end do  !in
      end do  !ip
!$omp end parallel do
call bundle_clock(ibundle,'sta')

   end do ! ibundle
   return
end subroutine applySPEbundled_orig

! Alternate version that threads across bundles, not slices
subroutine applySPEbundled_thread(vchar,startbundle,endbundle )
   use basis
   use sectors
   use diagh
   use precisions
   use lanczos_info
   use localvectors
   use nodeinfo
   use system_parameters
   use opbundles
   use fragments
   use contigpointervectors, only : vecin,vecout
   implicit none

   integer :: ibundle,startbundle,endbundle
   character(1) :: vchar
!------------------------------------------------------------

   integer(kind=8) nsdstart,nsdend,psdstart,psdend
   integer(kind=8) xjmp,xjmpstart,xjmpend

   integer(kind=8) :: ibasis
   real(kind=lanc_prec), pointer :: voutp(:)
!-------- OpenMP functions ---------------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: num_threads
   integer :: mythread
   integer(kind=basis_prec) :: vs

   real(4)    ::  pspe,nspe
   integer(kind=8) :: ip,in

   if(vchar == 'n')then
      vecin  => vec1
      vecout => vec2
   else
      vecin  => vec2
      vecout => vec1
   end if

! note: first private says that each such private copy is initialized with 
! the value before the parallel pragma

!$omp parallel do private(ibundle,ip,in,pspe,nspe,ibasis) & 
!$omp  private(psdstart, psdend, nsdstart, nsdend)        &
!$omp  private(voutp)                                     &
!$omp  private(vs, mythread)                              &
!$omp  shared(vecin,vecout)                               & 
!$omp  shared(vec2threadchunkm1)                          &
!$omp  shared(pstart,nstart,pspe_h,nspe_h)                &
!$omp     schedule(dynamic,1)
   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'SPE')cycle
	  call bundle_clock(ibundle,'sta')
! 	  if(opbundle(ibundle)%annexed)cycle
	  
      mythread = omp_get_thread_num()
      !! voutp(v2s:v2e) => vec2thread(:, mythread)
      vs = mythread * vec2threadchunk
      voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)

      psdstart = opbundle(ibundle)%pxstart
      psdend   = opbundle(ibundle)%pxend
      nsdstart = opbundle(ibundle)%nxstart
      nsdend   = opbundle(ibundle)%nxend

!.... LOOP OVER PROTON SDs in that sector..........
      do ip = psdstart,psdend
           pspe = pspe_h(ip)   ! the proton contribution to the single-particle energies
           ibasis = pstart(ip) + nstart(nsdstart)
!......... LOOP OVER NEUTRON SDS in sector jsc........................
           do in = nsdstart,nsdend
              nspe = nspe_h(in)  ! neutron contributions to s.p.e.
              voutp(ibasis) = voutp(ibasis) + vecin(ibasis)*( pspe + nspe )  ! add spes
              ibasis = ibasis+1   ! neutron SDs are contiguous
           end do  !in
      end do  !ip
	  call bundle_clock(ibundle,'end')
   end do ! ibundle
!$omp end parallel do
   return
end subroutine applySPEbundled_thread

subroutine applySPEbundled_g(vchar,startbundle,endbundle )
   use localvectors
   implicit none
   character(1), intent(in) :: vchar
   integer, intent(in) :: startbundle,endbundle

   if(useVec2Thread) then
      call applySPEbundled_thread(vchar, startbundle, endbundle)
   else
      call applySPEbundled_orig(vchar, startbundle, endbundle)
   end if
end subroutine applySPEbundled_g

!===============================================================
!
!  apply 3-body jumps for PPP
!
subroutine applyhPPPbundled_g (vchar,hchar,startbundle,endbundle )


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
  use contigpointervectors, only : vecin,vecout, p3b_1sd,p3b_2sd
  implicit none

  integer :: ibundle,startbundle,endbundle
  character(1) :: hchar,vchar

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=4)   xme
  integer(kind=8) :: statei, statef,nsd,psdi,psdf
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(4) :: mythread
  integer(kind=basis_prec) :: vs

  real(kind=lanc_prec), pointer :: voutp(:)

!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
  select case(vchar)
     case ('n')

        vecin  => vec1
        vecout => vec2
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
        if( hchar == 'f')then
           p3b_1sd => p3b_isd
           p3b_2sd => p3b_fsd
        else
           p3b_1sd => p3b_fsd
           p3b_2sd => p3b_isd
        endif

     case ('r')

        vecin  => vec2
        vecout => vec1
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
        if( hchar == 'b')then  !reversed from above
           p3b_1sd => p3b_isd
           p3b_2sd => p3b_fsd
        else
           p3b_1sd => p3b_fsd
           p3b_2sd => p3b_isd
        endif
     case default
        print *, "bad vchar=", vchar
        stop
  end select

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PPP')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle

!...... EXTRACT INFORMATION FROM OPBUNDLE ........

  csdstart = opbundle(ibundle)%nxstart
  csdend   = opbundle(ibundle)%nxend
  xjmpstart = opbundle(ibundle)%pxstart
  xjmpend   = opbundle(ibundle)%pxend
  cstride   = opbundle(ibundle)%cstride

  ncstates = (csdend +cstride -csdstart)/cstride

!--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
!          this makes for simple OpenMP threading
!$omp parallel private(xjmp, Xoplabel, xme, num_threads, mythread, vs)      &
!$omp          private(voutp) &
!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei,nsd)  &
!$omp          shared(vecin, cstride, vecout,  ncstates,xjmpstart,xjmpend)  &
!$omp          shared(vec2threadchunkm1)                                    &
!$omp          shared(p3b_op, p3b_1sd, p3b_2sd, p3b_phase, hmatppp)
  num_threads =  omp_get_num_threads()
  mythread = omp_get_thread_num()
  ! thread local vec2, reduce at end
  if(useVec2Thread) then
    !! voutp(v2s:v2e) => vec2thread(:, mythread)
    vs = mythread * vec2threadchunk
    voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
  else
    voutp(v2s:v2e) => vecout
  end if

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
          voutp(statef) = voutp(statef) + xme*vecin(statei)
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
              voutp(statef) = voutp(statef) + xme*vecin(statei)
          end do
      end do  ! xjmp

   end if
   end if
!$omp end parallel

   end do  ! ibundle
   return
   end subroutine applyhPPPbundled_g
!============================================================
!
!  apply 3-body jumps for NNN
!
!

subroutine applyhNNNbundled_g (vchar,hchar,startbundle,endbundle )
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
   use contigpointervectors, only : vecin,vecout, n3b_1sd,n3b_2sd
   implicit none

   integer :: ibundle,startbundle,endbundle
   character(1) :: hchar,vchar

!------------------------------------------------------------

   integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index,pstride
   integer(kind=8) xjmp,xjmpstart,xjmpend
   integer(kind=8):: Xoplabel
   real(kind=4)   xme
   integer(kind=8) :: statei, statef,psd,nsdi,nsdf,statef_start,statef_end
!-------- OpenMP functions ---------------------------------
   integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
   integer(kind=4) :: num_threads
   integer(kind=8) :: istart, iend, chunk
   integer(4) :: mythread
   integer(kind=basis_prec) :: vs

   real(kind=lanc_prec), pointer :: voutp(:)

!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
   select case(vchar)
      case ('n')
         vecin  => vec1
         vecout => vec2
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
         if( hchar == 'f')then
            n3b_1sd => n3b_isd
            n3b_2sd => n3b_fsd
         else
            n3b_1sd => n3b_fsd
            n3b_2sd => n3b_isd
         endif


      case ('r')
         vecin  => vec2
         vecout => vec1
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
         if( hchar == 'b')then  !reversed from above
            n3b_1sd => n3b_isd
            n3b_2sd => n3b_fsd
         else
            n3b_1sd => n3b_fsd
            n3b_2sd => n3b_isd
         endif
      case default
         print *, "bad vchar=", vchar
         stop
      end select
      do ibundle = startbundle,endbundle
         if(opbundle(ibundle)%optype /= 'NNN')cycle
         if(opbundle(ibundle)%hchar /= hchar )cycle
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
!$omp parallel private(xjmp, Xoplabel, xme, num_threads, mythread,vs)      &
!$omp          private(voutp) &
!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!$omp          private(statef_start,statef_end,psd) & 
!$omp          shared(vecin, cstride, vecout,  ncstates,xjmpstart,xjmpend,pstride)  &
!$omp          shared(vec2threadchunkm1)            &
!$omp          shared(n3b_op, n3b_1sd, n3b_2sd, n3b_phase, hmatnnn)
         num_threads =  omp_get_num_threads()
         mythread = omp_get_thread_num()
         ! thread local vec2, reduce at end
         if(useVec2Thread) then
            !! voutp(v2s:v2e) => vec2thread(:, mythread)
            vs = mythread * vec2threadchunk
            voutp(v2s:v2e) => vec2threadflat(vs: vs+vec2threadchunkm1)
         else
            voutp(v2s:v2e) => vecout
         end if

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
                  voutp(statef) = voutp(statef) + xme*vecin(statei)
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
                  voutp(statef) = voutp(statef) +  xme*vecin(statei)
                  statei = statei+pstride
               end do  !
            end do  ! xjmp
         end if
      end if
!$omp end parallel
   end do !ibundle
   return
end subroutine applyhNNNbundled_g
!============================================================
!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPNNbundledTR_g (vchar,hchar,startbundle,endbundle )
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
  use contigpointervectors, only : vecin,vecout, n2b_1sd,n2b_2sd
  implicit none

  integer :: ibundle,startbundle,endbundle
   character(1) :: hchar,vchar
!------------------------------------------------------------
  integer(kind=8) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer a,bb,c,dd
  integer(kind=8) :: coplabel,doplabel
  integer :: phasep,phasen
  real(kind=4)   xme
  integer(kind=8) :: statei, statef
  integer(kind=4) num_threads
  
  real(kind=lanc_prec), pointer :: voutp(:)

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(kind=basis_prec) :: vs
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  logical    :: launched(0:3)

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
  select case(vchar)
     case ('n')

        vecin  => vec1
        vecout => vec2

     case ('r')

        vecin  => vec2
        vecout => vec1

     case default
        print *, "bad vchar=", vchar
        stop

  end select

! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( (vchar == 'n' .and. hchar /= 'b') .or. & 
       (vchar == 'r' .and. hchar == 'b') )then

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PNN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
	 
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

!$omp parallel do private(vs, mythread,startp_thread,npjmps_thread)           &
!$omp          private(voutp) &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)              &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(vec2threadchunkm1)   &
!$omp          shared(n2b_isd,n2b_fsd,n2b_phase,n2b_cop,n2b_dop)              &
!$omp          shared(cpnntriplet,dpnntriplet,vecin,vecout)
!     mythread = omp_get_thread_num()

  do mythread = 0,numpthreads -1
     ! thread local vec2, reduce at end
     if(useVec2Thread) then
       !! voutp(v2s:v2e) => vec2thread(:, mythread)
       vs = mythread * vec2threadchunk
       voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
     else
       voutp(v2s:v2e) => vecout
     end if

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

           bb = n2b_cop(njmp)
           dd = n2b_dop(njmp)
           phasen = n2b_phase(njmp)*phasepnn(bb,a)*phasepnn(dd,c)

           coplabel = cpnntriplet(bb,a)
           doplabel = dpnntriplet(dd,c)

           xme = hmatpnn(coplabel + doplabel)   ! get matrix element
           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
           nsdi = n2b_isd(njmp)
           nsdf = n2b_fsd(njmp)
           statei = nsdi + psdi                ! initial state in combined basis
           statef = nsdf + psdf                ! final state in combined basis


           voutp(statef) = voutp(statef) + xme*vecin(statei)

        end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do

  end do ! ibundle

else

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PNN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
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

!$omp parallel do private(vs, mythread,startn_thread,nnjmps_thread)           &
!$omp          private(voutp) &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)    &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(vec2threadchunkm1)   &
!$omp          shared(n2b_isd,n2b_fsd,n2b_phase,n2b_cop,n2b_dop)    &
!$omp          shared(cpnntriplet,dpnntriplet,vecin,vecout)
   do mythread = 0, numnthreads-1
     ! thread local vec2, reduce at end
     if(useVec2Thread) then
       !! voutp(v2s:v2e) => vec2thread(:, mythread)
       vs = mythread * vec2threadchunk
       voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
     else
       voutp(v2s:v2e) => vecout
     end if

     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread

!...... OPTION TO SWITCH ORDER OF LOOPS WHEN NO OpenMP.......

     if(numnthreads > 1 .or. disableNoOMPloopswitch)then

     do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
        nsdi = n2b_isd(njmp)
        nsdf = n2b_fsd(njmp)
        bb  = n2b_cop(njmp)
        dd  = n2b_dop(njmp)
        phasen = n2b_phase(njmp)
        do pjmp = pjmpstart,pjmpend
           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
           psdf = p1b_fsd(pjmp)       ! final proton SD
           a  = p1b_cop(pjmp) 
           c  = p1b_dop(pjmp)
           phasep = p1b_phase(pjmp)*phasepnn(bb,a)*phasepnn(dd,c)   ! phase of proton jumps

!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
!----------- FIND MATRIX ELEMENT -------------------------------------------
           coplabel = cpnntriplet(bb,a)
           doplabel = dpnntriplet(dd,c)
           xme = hmatpnn(coplabel + doplabel)     ! get matrix element
           xme = xme*phasep*phasen               ! multiply matrix element by jump phases
           statei = nsdi + psdi                  ! initial state in combined basis
           statef = nsdf + psdf                  ! final state in combined basis

           voutp(statei) = voutp(statei) + xme*vecin(statef)

        end do  ! pjmp
     end do  ! njmp
     else    ! REVERSE ORDER OF LOOPS
     do pjmp = pjmpstart,pjmpend
           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
           psdf = p1b_fsd(pjmp)       ! final proton SD
           phasep = p1b_phase(pjmp)   ! phase of proton jumps
           a  = p1b_cop(pjmp) 
           c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
           do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
              nsdi = n2b_isd(njmp)
              nsdf = n2b_fsd(njmp)
              bb  = n2b_cop(njmp)
              dd  = n2b_dop(njmp)
              phasen = n2b_phase(njmp)

!----------- FIND MATRIX ELEMENT -------------------------------------------
              coplabel = cpnntriplet(bb,a)
              doplabel = dpnntriplet(dd,c)
              xme = hmatpnn(coplabel + doplabel)     ! get matrix element
              xme = xme*phasep*phasen               ! multiply matrix element by jump phases
              statei = nsdi + psdi                  ! initial state in combined basis
              statef = nsdf + psdf                  ! final state in combined basis
              voutp(statei) = voutp(statei) + xme*vecin(statef)
           end do
      end do  ! pjmp
      end if

  end do
!$omp end parallel do
  end do ! ibundle
  end if
  return
end subroutine applyhPNNbundledTR_g

!==========================================================
!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPNNbundled_g (vchar,hchar,startbundle,endbundle )
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
  use contigpointervectors, only : vecin,vecout, n2b_1sd,n2b_2sd
  
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
  integer(kind=4) num_threads

  real(kind=lanc_prec), pointer :: voutp(:)

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(kind=basis_prec) :: vs
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  logical    :: launched(0:3)

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
  select case(vchar)
     case ('n')

        vecin  => vec1
        vecout => vec2


     case ('r')

        vecin  => vec2
        vecout => vec1
     case default
        print *, "bad vchar=", vchar
        stop

  end select

! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( (vchar == 'n' .and. hchar /= 'b') .or. & 
       (vchar == 'r' .and. hchar == 'b') )then

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PNN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
	 
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

!$omp parallel do private(vs, mythread,startp_thread,npjmps_thread)           &
!$omp          private(voutp) &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)              &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(vec2threadchunkm1)   &
!$omp          shared(n2b_isd,n2b_fsd,n2b_phase,n2b_cop,n2b_dop)              &
!$omp          shared(cpnntriplet,dpnntriplet,vecin,vecout)
   do mythread = 0,numpthreads -1
     ! thread local vec2, reduce at end
     if(useVec2Thread) then
       ! voutp(v2s:v2e) => vec2thread(:, mythread)
       vs = mythread * vec2threadchunk
       voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
     else
       voutp(v2s:v2e) => vecout
     end if
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
           voutp(statef) = voutp(statef) + xme*vecin(statei)

        end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do

  end do ! ibundle
else

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PNN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
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

!$omp parallel do private(vs, mythread,startn_thread,nnjmps_thread)           &
!$omp          private(voutp) &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)    &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(n2b_isd,n2b_fsd,n2b_phase,n2b_cop,n2b_dop)    &
!$omp          shared(cpnntriplet,dpnntriplet,vecin,vecout)         &
!$omp          shared(vec2threadchunkm1)
    do mythread = 0, numnthreads-1
     ! thread local vec2, reduce at end
     if(useVec2Thread) then
       ! voutp(v2s:v2e) => vec2thread(:, mythread)
       vs = mythread * vec2threadchunk
       voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
     else
       voutp(v2s:v2e) => vecout
     end if

     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread
!...... OPTION TO SWITCH ORDER OF LOOPS WHEN NO OpenMP.......

     if(numnthreads > 1 .or. disableNoOMPloopswitch)then
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
           voutp(statei) = voutp(statei) + xme*vecin(statef)

        end do  ! pjmp
     end do  ! njmp
     else    ! REVERSE ORDER OF LOOPS
     do pjmp = pjmpstart,pjmpend
           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
           psdf = p1b_fsd(pjmp)       ! final proton SD
           phasep = p1b_phase(pjmp)   ! phase of proton jumps
           a  = p1b_cop(pjmp) 
           c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
           do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
              nsdi = n2b_isd(njmp)
              nsdf = n2b_fsd(njmp)
              bb  = n2b_cop(njmp)
              dd  = n2b_dop(njmp)
              phasen = n2b_phase(njmp)

!----------- FIND MATRIX ELEMENT -------------------------------------------
              coplabel = cpnntriplet(bb,a)
              doplabel = dpnntriplet(dd,c)
              xme = hmatpnn(coplabel + doplabel)     ! get matrix element
              xme = xme*phasep*phasen               ! multiply matrix element by jump phases
              statei = nsdi + psdi                  ! initial state in combined basis
              statef = nsdf + psdf                  ! final state in combined basis
              voutp(statei) = voutp(statei) + xme*vecin(statef)
           end do
      end do  ! pjmp
      end if
  end do
!$omp end parallel do
  end do ! ibundle
  end if
  return
end subroutine applyhPNNbundled_g

!==========================================================

!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPPNbundled_g (vchar,hchar,startbundle,endbundle )
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
  use contigpointervectors, only : vecin,vecout, p2b_1sd,p2b_2sd
 
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
  integer(kind=4) num_threads
  
  real(kind=lanc_prec), pointer :: voutp(:)

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(kind=basis_prec) :: vs
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  logical    :: launched(0:3)

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
  select case(vchar)
     case ('n')

        vecin  => vec1
        vecout => vec2


     case ('r')

        vecin  => vec2
        vecout => vec1

     case default
        print *, "bad vchar=", vchar
        stop
  end select

! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( (vchar == 'n' .and. hchar /= 'b') .or. & 
       (vchar == 'r' .and. hchar == 'b') )then

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PPN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
	 
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

!$omp parallel do private(vs, mythread,startp_thread,npjmps_thread)           &
!$omp          private(voutp)                                                 &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(p2b_isd,p2b_fsd,p2b_phase,p2b_cop,p2b_dop)              &
!$omp          shared(ibundle,opbundle)                                       &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cppntriplet,dppntriplet,vecin,vecout)                   &
!$omp          shared(vec2threadchunkm1)
!     mythread = omp_get_thread_num()


    do mythread = 0,numpthreads -1
     ! thread local vec2, reduce at end
     if(useVec2Thread) then
       !! voutp(v2s:v2e) => vec2thread(:, mythread)
       vs = mythread * vec2threadchunk
       voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
     else
       voutp(v2s:v2e) => vecout
     end if

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
           voutp(statef) = voutp(statef) + xme*vecin(statei)
       end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do
  end do ! ibundle
else

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PPN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
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

!$omp parallel do private(vs, mythread,startn_thread,nnjmps_thread)           &
!$omp          private(voutp)                                                 &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(p2b_isd,p2b_fsd,p2b_phase,p2b_cop,p2b_dop)              &
!$omp          shared(ibundle,opbundle)                                       &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cppntriplet,dppntriplet,vecin,vecout)                   &
!$omp          shared(vec2threadchunkm1)
    do mythread = 0, numnthreads-1
     ! thread local vec2, reduce at end
     if(useVec2Thread) then
       !! voutp(v2s:v2e) => vec2thread(:, mythread)
       vs = mythread * vec2threadchunk
       voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
     else
       voutp(v2s:v2e) => vecout
     end if

     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread

!...... OPTION TO SWITCH ORDER OF LOOPS WHEN NO OpenMP.......

     if(numnthreads > 1 .or. disableNoOMPloopswitch)then

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
           voutp(statei) = voutp(statei) + xme*vecin(statef)
        end do  ! pjmp
     end do  ! njmp

     else
     do pjmp = pjmpstart,pjmpend
           psdi = p2b_isd(pjmp)       ! initial proton slater determinant
           psdf = p2b_fsd(pjmp)       ! final proton SD
           phasep = p2b_phase(pjmp)   ! phase of proton jumps
           bb  = p2b_cop(pjmp) 
           dd = p2b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
           do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
              nsdi = n1b_isd(njmp)
              nsdf = n1b_fsd(njmp)
              a  = n1b_cop(njmp)
              c  = n1b_dop(njmp)
              phasen = n1b_phase(njmp)

!----------- FIND MATRIX ELEMENT -------------------------------------------
              coplabel = cppntriplet(bb,a)
              doplabel = dppntriplet(dd,c)
              xme = hmatppn(coplabel + doplabel)     ! get matrix element
              xme = xme*phasep*phasen               ! multiply matrix element by jump phases
              statei = nsdi + psdi                  ! initial state in combined basis
              statef = nsdf + psdf                  ! final state in combined basis
              voutp(statei) = voutp(statei) + xme*vecin(statef)
           end do
      end do  ! pjmp
     end if
  end do
!$omp end parallel do
  end do ! ibundle

  end if
  return
end subroutine applyhPPNbundled_g
!==========================================================

!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPPNbundledTR_g (vchar,hchar,startbundle,endbundle )
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use jump3body
  use interactions3body
  use opbundles
  use fragments
  use bmpi_mod
  use contigpointervectors, only : vecin,vecout, p2b_1sd,p2b_2sd
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
  integer(kind=4) num_threads

  real(kind=lanc_prec), pointer :: voutp(:)

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(kind=basis_prec) :: vs
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  logical    :: launched(0:3)

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
  select case(vchar)
     case ('n')

        vecin  => vec1
        vecout => vec2


     case ('r')

        vecin  => vec2
        vecout => vec1

     case default
        print *, "bad vchar=", vchar
        stop
  end select

! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( (vchar == 'n' .and. hchar /= 'b') .or. & 
       (vchar == 'r' .and. hchar == 'b') )then

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PPN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
	 
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

!$omp parallel do private(vs, mythread,startp_thread,npjmps_thread)           &
!$omp          private(voutp)                                                 &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(p2b_isd,p2b_fsd,p2b_phase,p2b_cop,p2b_dop)              &
!$omp          shared(ibundle,opbundle)                                       &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cppntriplet,dppntriplet,vecin,vecout)                   &
!$omp          shared(vec2threadchunkm1)
!     mythread = omp_get_thread_num()


    do mythread = 0,numpthreads -1
     ! thread local vec2, reduce at end
     if(useVec2Thread) then
       !! voutp(v2s:v2e) => vec2thread(:, mythread)
       vs = mythread * vec2threadchunk
       voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
     else
       voutp(v2s:v2e) => vecout
     end if

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
           a = n1b_cop(njmp)
           c = n1b_dop(njmp)
           phasen = n1b_phase(njmp)*phaseppn(bb,a)*phaseppn(dd,c)

           coplabel = cppntriplet(bb,a)
           doplabel = dppntriplet(dd,c)

           xme = hmatppn(coplabel + doplabel)   ! get matrix element
           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
           nsdi = n1b_isd(njmp)
           nsdf = n1b_fsd(njmp)
           statei = nsdi + psdi                ! initial state in combined basis
           statef = nsdf + psdf                ! final state in combined basis
           voutp(statef) = voutp(statef) + xme*vecin(statei)

       end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do
  end do !ibundle

else

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PPN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
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

!$omp parallel do private(vs, mythread,startn_thread,nnjmps_thread)           &
!$omp          private(voutp)                                                 &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,bb,c,dd) &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(p2b_isd,p2b_fsd,p2b_phase,p2b_cop,p2b_dop)              &
!$omp          shared(ibundle,opbundle)                                       &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cppntriplet,dppntriplet,vecin,vecout)                   &
!$omp          shared(vec2threadchunkm1)
!     mythread = omp_get_thread_num()
    do mythread = 0, numnthreads-1
     ! thread local vec2, reduce at end
     if(useVec2Thread) then
       !! voutp(v2s:v2e) => vec2thread(:, mythread)
       vs = mythread * vec2threadchunk
       voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
     else
       voutp(v2s:v2e) => vecout
     end if
     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread
!     do njmp = njmpstart,njmpend
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
           phasep = p2b_phase(pjmp)*phaseppn(bb,a)*phaseppn(dd,c)   ! phase of proton jumps

!--------- LOOP OVER NEUTRON JUMPS ---------------------------
!----------- FIND MATRIX ELEMENT ------------------------------

           coplabel = cppntriplet(bb,a)
           doplabel = dppntriplet(dd,c)
           xme = hmatppn(coplabel + doplabel)     ! get matrix element
           xme = xme*phasep*phasen               ! multiply matrix element by jump phases
           statei = nsdi + psdi                  ! initial state in combined basis
           statef = nsdf + psdf                  ! final state in combined basis
           voutp(statei) = voutp(statei) + xme*vecin(statef)


        end do  ! pjmp
     end do  ! njmp
  end do
!$omp end parallel do
  end do ! ibundle

  end if
  return
end subroutine applyhPPNbundledTR_g
!--------------------------------------------------------
!
! ADDED IN 7.6.3
! instead of applying full PN piece, just apply diagonal part
! to get correct application, only apply where SPEs are applied 
! note SPEs are diagonal in proton space and in neutron space, so no conflict
!
subroutine applyPNavgdiag(vchar,startbundle,endbundle)
    use basis
    use sectors
    use precisions
    use lanczos_info
    use localvectors
    use nodeinfo
    use system_parameters
    use opbundles
    use fragments
    use contigpointervectors, only : vecin,vecout
	use tracy
    implicit none

    integer :: ibundle,startbundle,endbundle
    character(1) :: vchar
  !------------------------------------------------------------

    integer(kind=4) num_threads
  !-------- OpenMP functions ---------------------------------
    integer(kind=4) :: omp_get_thread_num, omp_get_num_threads

    real(4)    ::  pspe,nspe
    integer(kind=8) :: istate, ibasis,fbasis
	integer :: is    ! which sector
	real(kind=4) :: xme

    if(vchar == 'n')then
       vecin  => vec1
       vecout => vec2
    else
       vecin  => vec2
       vecout => vec1
    end if


    do ibundle = startbundle,endbundle
       if(opbundle(ibundle)%optype /= 'SPE')cycle
!  	   if(opbundle(ibundle)%annexed)cycle
	   
       is=opbundle(ibundle)%isector
	   if(is < 0 .or. is > nsectors(1))then
		   print*,' some problem with sectors in pn diagonal trace '
		   print*,is
		   stop
	   end if
	   ibasis = xsd(1)%sector(is)%basisstart
	   fbasis = xsd(1)%sector(is)%basisend
	   xme = sectortrace(is)
	   do istate = ibasis,fbasis
		   vecout(istate)=vecout(istate)+xme*vecin(istate)
	   end do

    end do ! ibundle
    return
end subroutine applyPNavgdiag

end module apply_ham
