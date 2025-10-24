!===================================================================
!
!  file BAPPLYHLIBOMP_block.f90
!  
!  shan: using destination distribution as hint for work load distribution
!  07/05/2016
!
!  subroutines are adapted from routines found in bapplyhlibOMP.f90
!  to work with block Lanczos using openmp.
!  -modified by RMZ SDSU 3/2022 
!
!  routines to control application of the Hamiltonian
!  for parallel calculations with (optional) fragmentation of the basis and vectors
!
!  started 9/2011 by CWJ
!
!  NOTE: As of version 7.0.0, "obsolete" libraries of apply H have been discarded
!
!===================================================================

module applyhamblockomp
	use apply_ham_omp ! use same data structs for single vector case
	implicit none
	
contains

!
! subroutine applyHbundled_block_omp
!- Shan's omp algorithm applied to block Lanczos
! - Non-optiomized versions of these routines can be found in bapplyhlibOMP.f90
! - Modified by RMZ 3/2022 SDSU

! note: default is going from vecin to vecout but this can be reversed depending on hchar
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)
!   vchar = 'n' (normal), 'r' (reverse)
!      (fragments) of lanczos vectors stored in module localvectors
!        in vec1 and vec2; if normal  H vec1 = vec2
!                          if reverse H vec2 = vec1

!subroutine applyHbundled_omp(vchar)
subroutine applyHbundled_block_omp
  use nodeinfo
  use flagger
  use flags3body
  use precisions
  use opbundles
  use fragments
!  use shampoo
  use interaction
  use basis
  use localvectors
  use bmpi_mod
  use lanczos_info
  use contigpointervectors
  use sectors
  use diagh
  use jumpNbody
  use timing_parallel
  use timing
  use coupledmatrixelements, only:call_spe
  use system_parameters
!  use sectors
  use localblocks

  implicit none

  character(1) :: vchar, hchar 
  real(kind=8) :: sum
  integer :: tid
  integer(kind=basis_prec) :: i, j
  integer :: ierr

  integer :: ibundle
  integer :: mythread
  integer :: fs

  integer(kind=basis_prec) :: bleft, bright
  integer(kind=basis_prec) :: psdstart, psdend, nsdstart, nsdend, csdstart, csdend
  integer(kind=basis_prec) :: ip, in, statef, statei

  integer(kind=basis_prec) :: cstride, xjmpstart, xjmpend, xjmp
  integer(kind=basis_prec) :: istart, iend, csd, psd
  integer(kind=basis_prec) :: Xoplabel
  
  integer(kind=basis_prec) :: ncstates, csd_index, nsd, psdi, psdf, nsdf, nsdi

  integer(kind=basis_prec) :: pjmpstart, pjmpend, njmpstart, njmpend, pjmp, njmp
  integer :: phasep, phasen, a, b, c, d
  integer(kind=basis_prec) :: coplabel, doplabel
        
  real(kind=4) :: xme
  real(4) :: pspe, nspe
  
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer :: ivec
  integer(kind=8) :: locationi,locationf
  integer, save :: numiter = 0
  double precision :: timer1, t12

  double precision :: timers(8),get_Wtime

  double precision :: thrtimers(0:256), elapsed


  vchar = 'n'

  if (vchar /= 'n') then
     print *, "Vchar NN can only be n not ", vchar
     stop 4
  end if
 
  call proc_clock(iproc,'sta')

  vecin  => blockvec_1
  vecout => blockvec_2

  timers = 0.0
  thrtimers = 0.0
  
  numiter = numiter + 1
  timer1 = get_Wtime()
  if (verbose_ompham) print *, "Entering OMP      ", timer1, " iter ", numiter

  mythread = 0
  call clocker('hmu','sta')
  
!print*, "made it here a1"
!$omp parallel &
!$omp private(ibundle, mythread, fs, bleft, bright)  &
!$omp private(psdstart, psdend, nsdstart, nsdend, csdstart, csdend)  &
!$omp private(ip, ivec, pspe, statef, in, nspe, statei) & !statei
!$omp private(cstride, xjmpstart, xjmpend, xjmp) &
!$omp private(istart, iend, csd, psd) &
!$omp private(Xoplabel, xme) &
!$omp private(hchar) & 
!$omp private(locationi,locationf)&
!$omp private(ncstates, csd_index, nsd, psdi, psdf, nsdf, nsdi) &
!$omp private(pjmpstart, pjmpend, njmpstart, njmpend, pjmp, njmp) &
!$omp private(phasep, phasen, a, b, c, d) &
!$omp private(coplabel, doplabel) &
!$omp private(p2b_1sd, p2b_2sd, n2b_1sd, n2b_2sd) &
!$omp private(i, t12, elapsed) &
!$omp reduction(+ : timers) 
!!$omp private(locationi,locationf) & 
  
  mythread = omp_get_thread_num()
  do ibundle = opbundlestart(iproc), opbundleend(iproc)
!	 if(opbundle(ibundle)%annexed)cycle
   
     t12 = get_Wtime()
     fs = opbundle(ibundle)%fsector
     bleft = xsd(1)%sector(fs)%basisstart
     bright = xsd(1)%sector(fs)%basisend

     select case(opbundle(ibundle)%optype)
     case ('SPE')

        psdstart = opbundle(ibundle)%pxstart
        psdend   = opbundle(ibundle)%pxend
        nsdstart = opbundle(ibundle)%nxstart
        nsdend   = opbundle(ibundle)%nxend

        do ip = threadStart(ibundle, mythread), threadStop(ibundle, mythread)
           pspe = pspe_h(ip)   ! the proton contribution to the single-particle energies
           statei = pstart(ip) 
           statef = pstart(ip) + nstart(nsdstart)
           
          
           do in = nsdstart,nsdend
              nspe = nspe_h(in)  ! neutron contributions to s.p.e.
              ! if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) then
              !   print *, "Wrong exe SPE ", iproc, ibundle, statef, in, mythread, mybasisstart(mythread), mybasisstop(mythread)
              !   stop 11
              ! end if
             
              locationi = (statei-1)*dimblock
              locationf = (statef-1)*dimblock

            do ivec = 1,dimblock
                vecout(ivec+locationf) =vecout(ivec+locationf) + (pspe + nspe)*vecin(ivec+locationf)
   
            end do
              !vecout(statef) = vecout(statef) + vecin(statef)* (pspe + nspe)  ! add spes
              statef = statef+1   ! neutron SDs are contiguous
           end do  !in
        end do  !ip
        elapsed = get_Wtime() - t12
        timers(1) = timers(1) + elapsed
        thrtimers(mythread) = thrtimers(mythread) + elapsed
        
     case ('NN')
        hchar = opbundle(ibundle)%hchar
        if( hchar == 'b' )then
           n2b_1sd => n2b_fsd
           n2b_2sd => n2b_isd
        else
           n2b_1sd => n2b_isd
           n2b_2sd => n2b_fsd
        endif

        psdstart = opbundle(ibundle)%pxstart
        psdend   = opbundle(ibundle)%pxend
        cstride  = opbundle(ibundle)%cstride   !
        xjmpstart = opbundle(ibundle)%nxstart
        xjmpend   = opbundle(ibundle)%nxend
        
        istart = threadStart(ibundle, mythread)
        iend   = threadStop(ibundle, mythread)
        if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO NN MATRIX ELEMENTS
           do csd = istart, iend, cstride
              psd = pstart(csd)
              do xjmp = xjmpstart,xjmpend
                 Xoplabel = n2b_op(xjmp)
                 xme = hmatnn(Xoplabel)
                 xme = xme*n2b_phase(xjmp)
                 statei = n2b_1sd(xjmp) + psd
                 statef = n2b_2sd(xjmp) + psd
		           locationi = (statei-1)*dimblock
		           locationf = (statef-1)*dimblock
		  
                 do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
                 end do
                 !vecout(statef) = vecout(statef) +  xme*vecin(statei)
              end do  !
           end do  ! xjmp
           
        else

           !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
           do csd = istart, iend
              psd = pstart(csd)
              do xjmp = xjmpstart,xjmpend
                 !--------- FETCH MATRIX ELEMENT...............................
                 xme = n2b_me(xjmp)
                 !--------- GET PHASE.........................................
                 xme = xme*n2b_phase(xjmp)
                 statei = n2b_1sd(xjmp) + psd
                 statef = n2b_2sd(xjmp) + psd
!                 if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) then
!                    print *, "Wrong exe NN1 ", iproc, ibundle, statef, n2b_2sd(xjmp), mythread, mybasisstart(mythread), mybasisstop(mythread)
!                    stop 11
!                 end if
                 locationi = (statei-1)*dimblock
		           locationf = (statef-1)*dimblock
		  
                 do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
                 end do
                 !vecout(statef) = vecout(statef) +  xme*vecin(statei)
              end do  !
           end do  ! xjmp
        end if

        elapsed = get_Wtime() - t12
        timers(2) = timers(2) + elapsed
        thrtimers(mythread) = thrtimers(mythread) + elapsed

     case ('PP')
        hchar = opbundle(ibundle)%hchar
        if( hchar == 'b')then
           p2b_1sd => p2b_fsd
           p2b_2sd => p2b_isd
        else
           p2b_1sd => p2b_isd
           p2b_2sd => p2b_fsd
        endif
        csdstart = opbundle(ibundle)%nxstart
        csdend   = opbundle(ibundle)%nxend
        cstride   = opbundle(ibundle)%cstride
        xjmpstart = opbundle(ibundle)%pxstart
        xjmpend   = opbundle(ibundle)%pxend
        istart = threadStart(ibundle, mythread)
        iend   = threadStop(ibundle, mythread)

        if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO PP MATRIX ELEMENTS

           do xjmp = istart, iend
              !------------------ FETCH MATRIX ELEMENT...............................
              Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
              xme = hmatpp(Xoplabel)
              !--------- GET PHASE.........................................
              xme = xme*p2b_phase(xjmp)
              psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
              psdf = p2b_2sd(xjmp)   ! KSM: final P SD

              do csd = csdstart, csdend, cstride
                 nsd = nstart(csd)
                 statef = psdf + nsd
                 statei = psdi + nsd
                 locationi = (statei-1)*dimblock
		           locationf = (statef-1)*dimblock
		  
                 do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
                 end do
                 !vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do ! csd
           end do  ! xjmp

           if (ompNumThreads <= 1) then
              elapsed = get_Wtime() - t12
              thrtimers(mythread) = thrtimers(mythread) + elapsed

              if(hchar == 'b')then
                  timers(4) = timers(4) + elapsed
              else
                  timers(3) = timers(3) + elapsed
              end if
              cycle     ! may not necessary
           end if

           ! do irreglar: left
           istart = xjmpstart
           do i = mythread -1, 0, -1
              if (threadStart(ibundle, i) > 0) then
                 istart = threadStop(ibundle, i) + 1
                 exit
              end if
           end do
           if (threadStart(ibundle, mythread) > 0) then
              iend = threadStart(ibundle, mythread) -1
           else
              iend = xjmpend
              do i = mythread + 1, ompNumThreads-1
                 if (threadStart(ibundle, i) > 0) then
                    iend = threadStart(ibundle, i) -1
                    exit
                 end if
              end do
           end if

           do xjmp = istart, iend
              
              psdf = p2b_2sd(xjmp)   ! KSM: final P SD

              if (psdf < mybasisstart(mythread) -1) cycle
              if (psdf >= mybasisstop(mythread)) cycle 

!              i = 0
!              do while (psdf >= mybasisstop(i))
!                 i = i + 1
!              end do
!              if (i /= mythread) cycle

              Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
              xme = hmatpp(Xoplabel)
              !--------- GET PHASE.........................................
              xme = xme*p2b_phase(xjmp)
              psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
              
              do csd = csdstart, csdend, cstride
                 nsd = nstart(csd)
                 statef = psdf + nsd
                 statei = psdi + nsd
                 locationi = (statei-1)*dimblock
		           locationf = (statef-1)*dimblock
		  
                 do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
                 end do
                 !vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do ! csd
           end do  ! xjmp

           if (threadStart(ibundle, mythread) > 0 ) then
              istart = threadStop(ibundle, mythread) +1
              iend = xjmpend
              do i = mythread + 1, ompNumThreads -1
                 if (threadStart(ibundle, i) > 0) then
                    iend = threadStart(ibundle, i) -1 
                 end if
              end do
           else
              istart = 0
              iend   = -1
           end if
           
           do xjmp = istart, iend
              
              psdf = p2b_2sd(xjmp)   ! KSM: final P SD

              if (psdf < mybasisstart(mythread) -1) cycle
              if (psdf >= mybasisstop(mythread) ) cycle
              
!              i = 0
!              do while (psdf >= mybasisstop(i))
!                 i = i + 1
!              end do
!              if (i /= mythread) cycle

              Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
              xme = hmatpp(Xoplabel)
              !--------- GET PHASE.........................................
              xme = xme*p2b_phase(xjmp)
              psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
              
              do csd = csdstart, csdend, cstride
                 nsd = nstart(csd)
                 statef = psdf + nsd
                 statei = psdi + nsd
                 locationi = (statei-1)*dimblock
		           locationf = (statef-1)*dimblock
		  
                 do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
                 end do
                 !vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do ! csd
           end do  ! xjmp

           
        else    ! USE STORED PP MATRIX ELEMENTS DIRECTLY........ 
           do xjmp = xjmpstart,xjmpend
              psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
              psdf = p2b_2sd(xjmp)   ! KSM: final P SD
              xme = p2b_me(xjmp)
              xme = xme*p2b_phase(xjmp)
                               
i = 0
do while (psdf >= mybasisstop(i))
  i = i + 1
end do
if (i /= mythread) cycle

              do csd = csdstart, csdend, cstride
                 nsd = nstart(csd)
                 statef = psdf+nsd !csd_index
!!!                 if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
                 
                 statei = psdi+ nsd !csd_index
                 locationi = (statei-1)*dimblock
		           locationf = (statef-1)*dimblock
		  
                 do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
                 end do
                 !vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do  ! xjmp
           end do  ! csd
        end if

        elapsed = get_Wtime() - t12
        thrtimers(mythread) = thrtimers(mythread) + elapsed

        if( hchar == 'b')then
           timers(4) = timers(4) + elapsed
        else
           timers(3) = timers(3) + elapsed
        end if
     case ('PN')

        pjmpstart = opbundle(ibundle)%pxstart
        pjmpend   = opbundle(ibundle)%pxend
        njmpstart = opbundle(ibundle)%nxstart
        njmpend   = opbundle(ibundle)%nxend

        hchar = opbundle(ibundle)%hchar
        istart = threadStart(ibundle, mythread)
        iend   = threadStop(ibundle, mythread)
        if (hchar /= 'b') then
           do pjmp = istart, iend
              psdi = p1b_isd(pjmp)       ! initial proton slater determinant
              psdf = p1b_fsd(pjmp)       ! final proton SD
              phasep = p1b_phase(pjmp)   ! phase of proton jumps
              a = p1b_cop(pjmp)     ! KSM: Proton 1-body creation label
              c = p1b_dop(pjmp)     ! KSM: Proton 1-body destruction label
              !--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
              do njmp = njmpstart,njmpend
                 nsdf = n1b_fsd(njmp)
                 statef = nsdf + psdf                ! final state in combined basis
!                 if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle

                 !----------- FIND MATRIX ELEMTN --------------------------------------------
                 b = n1b_cop(njmp)  ! KSM: Neutron 1-body creation label
                 d = n1b_dop(njmp)  ! KSM: Neutron 1-body destruction label
                 phasen = n1b_phase(njmp)
                 coplabel = cpnpair(b,a)
                 doplabel = dpnpair(d,c)
                 xme = hmatpn(coplabel + doplabel)   ! get matrix element
                 xme = xme*phasep*phasen             ! multiply matrix element by jump phases
                 nsdi = n1b_isd(njmp)
                 statei = nsdi + psdi                ! initial state in combined basis
                 locationi = (statei-1)*dimblock
		           locationf = (statef-1)*dimblock
		  
                 do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
                 end do
                 !vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do  ! njmp
           end do  ! pjmp

          elapsed = get_Wtime() - t12
          thrtimers(mythread) = thrtimers(mythread) + elapsed
          timers(5) = timers(5) + elapsed
           
      else
         !---- Backward direction using hermiticity -------------------
         ! do regular first

         istart = threadStart(ibundle, mythread)
         iend   = threadStop(ibundle, mythread)
         
         do pjmp = istart, iend
            psdf = p1b_fsd(pjmp)       ! final proton SD
            psdi = p1b_isd(pjmp)       ! initial proton slater determinant
            phasep = p1b_phase(pjmp)   ! phase of proton jumps
            a  = p1b_cop(pjmp)
            c  = p1b_dop(pjmp)

            do njmp = njmpstart,njmpend
               nsdi = n1b_isd(njmp)
               nsdf = n1b_fsd(njmp)
               b  = n1b_cop(njmp)
               d  = n1b_dop(njmp)
               phasen = n1b_phase(njmp)

               statef = nsdf + psdf                  ! final state in combined basis
               statei = nsdi + psdi                  ! initial state in combined basis
!!!               if (statei < mybasisstart(mythread) .or. statei > mybasisstop(mythread)) cycle
               
               coplabel = cpnpair(b,a)
               doplabel = dpnpair(d,c)
               xme = hmatpn(coplabel + doplabel)     ! get matrix element
               xme = xme*phasep*phasen               ! multiply matrix element by jump phases
               locationi = (statei-1)*dimblock
               locationf = (statef-1)*dimblock
      
               do ivec = 1,dimblock
                  vecout(ivec+locationi) =vecout(ivec+locationi) + xme*vecin(ivec+locationf)
               end do
               !vecout(statei) = vecout(statei) + xme*vecin(statef)
            end do  ! pjmp
         end do  ! njmp

         if (ompNumThreads <= 1) then 
           elapsed = get_Wtime() - t12
           thrtimers(mythread) = thrtimers(mythread) + elapsed
           timers(6) = timers(6) + elapsed
           cycle ! may not necessary
         end if

         ! do irreglar: left
         istart = pjmpstart
         do i = mythread -1, 0, -1
            if (threadStart(ibundle, i) > 0) then
               istart = threadStop(ibundle, i) + 1
               exit
            end if
         end do
         if (threadStart(ibundle, mythread) > 0) then
            iend = threadStart(ibundle, mythread) -1
         else
            iend = pjmpend
            do i = mythread + 1, ompNumThreads-1
               if (threadStart(ibundle, i) > 0) then
                  iend = threadStart(ibundle, i) -1
                  exit
               end if
            end do
         end if
         
!          write(*,117) ibundle, iproc, ' thread ', mythread, '  left ', istart, iend
117 format("Bundle ", i6, i6, a8, i4, a7, 2i12)

         do pjmp = istart, iend

            psdi = p1b_isd(pjmp)       ! initial proton slater determinant
            if (psdi < mybasisstart(mythread) -1) cycle
            if (psdi >= mybasisstop(mythread)) cycle

! i = 0
! do while (psdi >= mybasisstop(i))
!   i = i + 1
! end do
! if (i /= mythread) cycle

            psdf = p1b_fsd(pjmp)       ! final proton SD
            phasep = p1b_phase(pjmp)   ! phase of proton jumps
            a  = p1b_cop(pjmp)
            c  = p1b_dop(pjmp)

            do njmp = njmpstart,njmpend
               nsdi = n1b_isd(njmp)
               nsdf = n1b_fsd(njmp)
               b  = n1b_cop(njmp)
               d  = n1b_dop(njmp)
               phasen = n1b_phase(njmp)

               statef = nsdf + psdf                  ! final state in combined basis
               statei = nsdi + psdi                  ! initial state in combined basis
!!!               if (statei < mybasisstart(mythread) .or. statei > mybasisstop(mythread)) cycle
               
               coplabel = cpnpair(b,a)
               doplabel = dpnpair(d,c)
               xme = hmatpn(coplabel + doplabel)     ! get matrix element
               xme = xme*phasep*phasen               ! multiply matrix element by jump phases
               locationi = (statei-1)*dimblock
               locationf = (statef-1)*dimblock
      
               do ivec = 1,dimblock
                  vecout(ivec+locationi) =vecout(ivec+locationi) + xme*vecin(ivec+locationf)
               end do
               !vecout(statei) = vecout(statei) + xme*vecin(statef)
            end do  ! pjmp
         end do  ! njmp

         ! do irregular: right

         if (threadStart(ibundle, mythread) > 0 ) then
            istart = threadStop(ibundle, mythread) +1
            iend = pjmpend
            do i = mythread + 1, ompNumThreads -1
               if (threadStart(ibundle, i) > 0) then
                  iend = threadStart(ibundle, i) -1 
               end if
            end do
         else
            istart = 0
            iend   = -1
         end if

!          write(*,118) ibundle, iproc, ' thread ', mythread, ' right ', istart, iend
118 format("Bundle ", i6, i6, a8, i4, a7, 2i12)

         do pjmp = istart, iend

            psdi = p1b_isd(pjmp)       ! initial proton slater determinant
            if (psdi < mybasisstart(mythread) -1) cycle
            if (psdi >= mybasisstop(mythread)) cycle

! i = 0
! do while (psdi >= mybasisstop(i))
!   i = i + 1
! end do
! if (i /= mythread) cycle

            
            psdf = p1b_fsd(pjmp)       ! final proton SD
            phasep = p1b_phase(pjmp)   ! phase of proton jumps
            a  = p1b_cop(pjmp)
            c  = p1b_dop(pjmp)

            do njmp = njmpstart,njmpend
               nsdi = n1b_isd(njmp)
               nsdf = n1b_fsd(njmp)
               b  = n1b_cop(njmp)
               d  = n1b_dop(njmp)
               phasen = n1b_phase(njmp)

               statef = nsdf + psdf                  ! final state in combined basis
               statei = nsdi + psdi                  ! initial state in combined basis
!!!               if (statei < mybasisstart(mythread) .or. statei > mybasisstop(mythread)) cycle
               
               coplabel = cpnpair(b,a)
               doplabel = dpnpair(d,c)
               xme = hmatpn(coplabel + doplabel)     ! get matrix element
               xme = xme*phasep*phasen               ! multiply matrix element by jump phases
               locationi = (statei-1)*dimblock
               locationf = (statef-1)*dimblock
      
               do ivec = 1,dimblock
                  vecout(ivec+locationi) =vecout(ivec+locationi) + xme*vecin(ivec+locationf)
               end do
               !vecout(statei) = vecout(statei) + xme*vecin(statef)
            end do  ! pjmp
         end do  ! njmp

        elapsed = get_Wtime() - t12
        thrtimers(mythread) = thrtimers(mythread) + elapsed
        timers(6) = timers(6) + elapsed
         
      end if
      
     case default
        print *, "wrong opbundle type exe ", opbundle(ibundle)%optype
        stop 30
        
     end select

  end do

!$omp end parallel
  
  time_procSPE(iproc) = time_procSPE(iproc) + timers(1)
  time_procNN(iproc)  = time_procNN(iproc)  + timers(2)
  time_procPP(iproc)  = time_procPP(iproc)  + timers(3)
  time_procPPb(iproc) = time_procPPb(iproc) + timers(4)
  time_procPN(iproc)  = time_procPN(iproc)  + timers(5)
  time_procPNb(iproc)  = time_procPNb(iproc) + timers(6)
  
  time_spe = time_spe + timers(1)
  time_nn  = time_nn  + timers(2)
  time_pp  = time_pp  + timers(3) + timers(4)
  time1body = time1body + timers(5) + timers(6)
  
  if (verbose_ompham) write(*, 230) "Finishing OMP     ", iproc, get_Wtime() - timer1, " iter ", numiter
 230 format(a19, i6, f12.4, a6, i6)
  if (verbose_ompham) print *, iproc, "Thread timing dist ", thrtimers(0:ompNumThreads)

if (.false.) then
  !............SPE......................................
  if(call_spe)then
     call clocker('spe','sta')
     call procOP_clock(iproc,'sta','SPE')
     if(noisy0) print *, "Starting applyHbundled_block_g SPE"
     call applySPEbundled_block_omp(vchar,opbundlestart(iproc), opbundleend(iproc))
     if(applypntrace)call applyPNavgdiag(vchar,opbundlestart(iproc), opbundleend(iproc))
     call clocker('spe','end')
     call procOP_clock(iproc,'end','SPE')
  end if

  if(.not.threebody)then
     !........... PP .................
     call clocker('ppo','sta')
     call procOP_clock(iproc,'sta','PPO')
     if(noisy0) print *, "Starting applyHbundled_block_g PP"
     call applyhPPbundled_block_omp(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
     call procOP_clock(iproc,'end','PPO')
     
     call procOP_clock(iproc,'sta','PPB')
     call applyhPPbundled_block_omp(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     call procOP_clock(iproc,'end','PPB')

     call clocker('ppo','end')
     !............NN ..................
     call clocker('nno','sta')
     call procOP_clock(iproc,'sta','NNO')

     if(noisy0) print *, "Starting applyHbundled_g NN"
     call applyhNNbundled_block_omp(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
     call applyhNNbundled_block_omp(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     call clocker('nno','end')
     call procOP_clock(iproc,'end','NNO')

     !........... PN....................................
     call clocker('one','sta')
     
     if(noisy0) print *, "Starting applyHbundled_g PN"
     if(useTR)then
        call procOP_clock(iproc,'sta','PNO')
        call applyhPNbundledTR_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPNbundledTR_g(vchar,'h',opbundlestart(iproc), opbundleend(iproc))
        call procOP_clock(iproc,'end','PNO')
        
        call procOP_clock(iproc,'sta','PNB')
        call applyhPNbundledTR_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        call procOP_clock(iproc,'end','PNB')
     else
        call procOP_clock(iproc,'sta','PNO')
        call applyhPNbundled_block_omp(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPNbundled_block_omp(vchar,'h',opbundlestart(iproc), opbundleend(iproc))
        call procOP_clock(iproc,'end','PNO')

        call procOP_clock(iproc,'sta','PNB') 
        call applyhPNbundled_block_omp(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        call procOP_clock(iproc,'end','PNB')
     end if
     call clocker('one','end')

  else
   ! three-body is not optimized for openMP
     !............PPP.......................................
     call clocker('ppp','sta')
     call procOP_clock(iproc,'sta','PPP')
     
     call applyhPPPbundled_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
     call applyhPPPbundled_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     call clocker('ppp','end')
     call procOP_clock(iproc,'end','PPP')

     !............PPN.......................................
     call clocker('ppn','sta')
     call procOP_clock(iproc,'sta','PPN')
     
     if(useTR)then
        call applyhPPNbundledTR_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPPNbundledTR_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     else
        call applyhPPNbundled_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPPNbundled_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        
     endif
     call procOP_clock(iproc,'end','PPN')
     
     call clocker('ppn','end')
     !............PNN.......................................
     call clocker('pnn','sta')
     call procOP_clock(iproc,'sta','PNN')
     
     if(useTR)then
        call applyhPNNbundledTR_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPNNbundledTR_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        
     else
        call applyhPNNbundled_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPNNbundled_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        
     endif
     call clocker('pnn','end')
     call procOP_clock(iproc,'end','PNN')
     
     !...........NNN......................................
     call clocker('nnn','sta')
     call procOP_clock(iproc,'sta','NNN')
     
     call applyhNNNbundled_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
     call applyhNNNbundled_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     call clocker('nnn','end')
     call procOP_clock(iproc,'end','NNN')
     
     
  end if  ! if threebody

end if

  ! If using thread local output storage we have to reduce to vec2
  ! unconverted bundle types are still dumping in vec2 so we just add to it
  ! from vec2thread
!  if(useVec2Thread) then
     ! the schedule here is static with a chunksize.  This make sense
     ! because each chunk will have predictable runtime
!!$omp parallel do                      &
!!$omp    private(i, j, tid, sum)          &
!!$omp    shared(vec2threadflat, v2s, v2e)  &
!!$omp    schedule(static, 1024)
!     do i = v2s, v2e
!        sum = vec2(i)
        !! do tid = 0, ompNumThreads-1
        !!    sum = sum + vec2thread(i, tid)
        !!    vec2thread(i, tid) = 0.0
        !! end do
        
        ! sum the contributions to vec2(i)
!        do j = (i - v2s), vec2threadend, vec2threadchunk
!           sum = sum + vec2threadflat(j)
!           vec2threadflat(j) = 0.0;
!        end do
!        vec2(i) = real(sum,kind=lanc_prec)
!     end do
!     !$omp end parallel do
!  end if
  
  call proc_clock(iproc,'end')
  call clocker('hmu','end')
  

  !........... NOW REDUCE....................
  !            WHEN LANCZOS VECTOR BROKEN, 
  !            THIS WILL GET MORE COMPLICATED
  if(noisy0) print *, "applyhbundled_g: Doing Reduce"
#ifdef _MPI  
  if(useNewReorthog) then
     if(vchar == 'n') then
        ! 'n' corresponds with vec2 here because we are looking at the output vector
        ! call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
        ! Do reduce only onto root node.  We will be sending this data to 
        ! the slices from the isfragroot nodes (rank=0 in fcomm1, fcomm2, hcomm),
        ! so other nodes don't need it.
        call clocker('blr', 'sta')
        call BMPI_REDUCE(blockvec_2, size(blockvec_2), MPI_SUM, 0, fcomm2, ierr) ! in place
        call clocker('blr', 'end')
     else
        ! call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, fcomm1, ierr) ! in place
        call BMPI_REDUCE(blockvec_1, size(blockvec_1), MPI_SUM, 0, fcomm1, ierr) ! in place
     end if
  else
     if(nproc > 1 .and. vchar == 'n')then
        call BMPI_ALLREDUCE(blockvec_2, size(blockvec_2), MPI_SUM, fcomm2, ierr) ! in place
     endif
     if(nproc > 1 .and. vchar == 'r')then
        call BMPI_ALLREDUCE(blockvec_1, size(blockvec_1), MPI_SUM, fcomm1, ierr) ! in place
     end if
  end if
#endif
  
  if (numiter == 12.and.verbose_ompham) write(*, 220)  "Finishing Reduce   ", iproc, get_Wtime() - timer1, " iter ", numiter
220 format(a19, i6, f12.4, a6, i6)

  if(noisy0) print *, "applyhbundled_g: Returning"
  return

end subroutine applyHbundled_block_omp

subroutine applyhPPbundled_block_omp (vchar,hchar,startbundle,endbundle )
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
  use contigpointervectors
  use localblocks
  use sectors

  implicit none

  integer :: ibundle
  character(1) :: hchar,vchar
  integer :: startbundle,endbundle

  integer :: fs
  integer(kind=8) :: bleft, bright

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
  logical, save :: yes = .true.

  integer :: ivec
  integer(kind=8) :: locationi,locationf

  if (vchar /= 'n') then
     print *, "Vchar can only be n not ", vchar
     stop 4
  end if
  
  vecin  => blockvec_1
  vecout => blockvec_2

  if( hchar == 'f')then
     p2b_1sd => p2b_isd
     p2b_2sd => p2b_fsd
  else
     p2b_1sd => p2b_fsd
     p2b_2sd => p2b_isd
  endif

!$omp parallel private(ibundle, fs, bleft, bright)  &
!$omp          private(xjmp, Xoplabel, xme, nsd)    &
!$omp          private(mythread) &
!$omp          private(istart,ivec, iend, csd, csd_index, statef, statei)  &
!$omp          private(psdi, psdf)  &
!$omp          private(locationi,locationf)&
!$omp          private(csdstart, csdend, xjmpstart, xjmpend, cstride) &
!$omp          private(ncstates) 
  do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'PP')cycle
      if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
      
      mythread = omp_get_thread_num()

      fs = opbundle(ibundle)%fsector
      bleft = xsd(1)%sector(fs)%basisstart
      bright = xsd(1)%sector(fs)%basisend
!!      if (bleft > mybasisstop(mythread) .or. bright < mybasisstart(mythread)) cycle
      
      !...... EXTRACT INFORMATION FROM OPBUNDLE ........
      csdstart = opbundle(ibundle)%nxstart
      csdend   = opbundle(ibundle)%nxend
      xjmpstart = opbundle(ibundle)%pxstart
      xjmpend   = opbundle(ibundle)%pxend
      cstride   = opbundle(ibundle)%cstride
      
      ncstates = (csdend +cstride -csdstart)/cstride
      istart = 1
      iend = ncstates
      csd_index = csdstart + (istart - 1)*cstride - cstride

      if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO PP MATRIX ELEMENTS
         ! KSM - Test both ways!!
         if(.false.) then
            do csd = istart, iend
               csd_index = csd_index + cstride
               nsd = nstart(csd_index)

               !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
               do xjmp = xjmpstart,xjmpend
                  statef = p2b_2sd(xjmp)+ nsd !csd_index
                  if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
                  !--------- FETCH MATRIX ELEMENT...............................
                  Xoplabel = p2b_op(xjmp)
                  xme = hmatpp(Xoplabel)
                  !--------- GET PHASE.........................................
                  xme = xme*p2b_phase(xjmp)
                  !---------- GET INITIAL, FINAL SDs and place in basis..............
                  statei = p2b_1sd(xjmp)+ nsd !csd_index
                  locationi = (statei-1)*dimblock
                  locationf = (statef-1)*dimblock
              
                  do ivec = 1,dimblock
                       vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
                  end do
                  !vecout(statef) = vecout(statef) + xme*vecin(statei)
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
                  statef = psdf + nsd
                  if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
                  
                  statei = psdi + nsd 
                  locationi = (statei-1)*dimblock
                  locationf = (statef-1)*dimblock
              
                  do ivec = 1,dimblock
                       vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
                  end do
                  !vecout(statef) = vecout(statef) + xme*vecin(statei)
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
               statef = p2b_2sd(xjmp)+nsd !csd_index
               if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle

               !--------- FETCH MATRIX ELEMENT...............................
               xme = p2b_me(xjmp)
               !                 Xoplabel = p2b_op(xjmp)
               !                 xme = hmatpp(Xoplabel)
               !--------- GET PHASE.........................................
               xme = xme*p2b_phase(xjmp)
               !---------- GET INITIAL, FINAL SDs and place in basis..............
               statei = p2b_1sd(xjmp)+ nsd !csd_index
               locationi = (statei-1)*dimblock
               locationf = (statef-1)*dimblock
              
               do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
               end do
               !vecout(statef) = vecout(statef) + xme*vecin(statei)
            end do  ! xjmp
         end do  ! csd
      end if
      !--------------OR DO HERMITIAN/BACKWARDS APPLICATION----------
   end do ! ibundle
   !$omp end parallel
   return
 end subroutine applyhPPbundled_block_omp


! This version has each thread write into a different output
! buffer that is reduced with other threads at the end.
subroutine applyhNNbundled_block_omp(vchar,hchar,startbundle,endbundle )
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
   use contigpointervectors
   use localblocks
   use sectors

   implicit none

   ! arguments
   character(1),intent(in) :: hchar,vchar
   integer,intent(in) :: startbundle,endbundle
   integer :: ivec
   integer(kind=8) :: locationi,locationf
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

   integer :: fs
   integer(kind=8) :: bleft, bright

   if (vchar /= 'n') then
      print *, "Vchar NN can only be n not ", vchar
      stop 4
   end if

   vecin  => blockvec_1
   vecout => blockvec_2

   if( hchar == 'f' )then
      n2b_1sd => n2b_isd
      n2b_2sd => n2b_fsd
   else
      n2b_1sd => n2b_fsd
      n2b_2sd => n2b_isd
   endif

! firstprivate gives each thread its own copy, but initializes it
! this is   better than shared for read-only vars

!$omp parallel private(ibundle, fs, bleft, bright)        &
!$omp         private(mythread)       &
!$omp         private(csdstart, csdend, cstride)          &
!$omp         private(xjmpstart, xjmpend, ncstates)       &
!$omp         private(istart, iend, ivec, csd, csd_index, psd)            &
!$omp         private(pstride)                            &
!$omp         private(xjmp, Xoplabel, xme)                &
!$omp         private(statei, statef)                     &
!$omp          private(locationi,locationf)&
!$omp         private(statef_start, statef_end)           &
!$omp         firstprivate(startbundle, endbundle)        

   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'NN')cycle
      if(opbundle(ibundle)%hchar /= hchar )cycle
!	  if(opbundle(ibundle)%annexed)cycle
	  
      !....... Thread specific setup
      mythread = omp_get_thread_num()

      fs = opbundle(ibundle)%fsector
      bleft = xsd(1)%sector(fs)%basisstart
      bright = xsd(1)%sector(fs)%basisend
!!!      if (bleft > mybasisstop(mythread) .or. bright < mybasisstart(mythread)) cycle

      !...... EXTRACT INFORMATION FROM OPBUNDLE ........
      csdstart = opbundle(ibundle)%pxstart
      csdend   = opbundle(ibundle)%pxend
      cstride  = opbundle(ibundle)%cstride   !
      xjmpstart = opbundle(ibundle)%nxstart
      xjmpend   = opbundle(ibundle)%nxend

      istart = threadStart(ibundle, mythread)
      iend   = threadStop(ibundle, mythread)
      if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO NN MATRIX ELEMENTS
         !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
         do csd = istart, iend, cstride 
            psd = pstart(csd)
            do xjmp = xjmpstart,xjmpend
              Xoplabel = n2b_op(xjmp)
              xme = hmatnn(Xoplabel)
              xme = xme*n2b_phase(xjmp)
              statei = n2b_1sd(xjmp) + psd
              statef = n2b_2sd(xjmp) + psd
              !              if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
!              if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) then
!                 print *, "Wrong exe NN ", iproc, ibundle, statef, istart, iend, cstride, csd, n2b_2sd(xjmp), mythread, mybasisstart(mythread), mybasisstop(mythread)
!                 print *, "Wrong exe NN2 ", pstart(116627), pstart(116840)
!                 stop 11
!              end if
              
              !vecout(statef) = vecout(statef) +  xme*vecin(statei)
              locationi = (statei-1)*dimblock
              locationf = (statef-1)*dimblock
      
               do ivec = 1,dimblock
                    vecout(ivec+locationf) = 	vecout(ivec+locationf) + xme*vecin(ivec+locationi)
               end do
            end do  !
         end do  ! xjmp
      else
         !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
         do csd = istart, iend
            psd = pstart(csd)
            do xjmp = xjmpstart,xjmpend
              !--------- FETCH MATRIX ELEMENT...............................
              xme = n2b_me(xjmp)
              !--------- GET PHASE.........................................
              xme = xme*n2b_phase(xjmp)
              statei = n2b_1sd(xjmp) + psd
              statef = n2b_2sd(xjmp) + psd
              !---------- GET INITIAL, FINAL SDs and place in basis..............
              ! if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
              if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) then
                 print *, "Wrong exe NN1 ", iproc, ibundle, statef, n2b_2sd(xjmp), mythread, &
				   mybasisstart(mythread), mybasisstop(mythread)
                 stop 11
              end if
              !vecout(statef) = vecout(statef) +  xme*vecin(statei)
              locationi = (statei-1)*dimblock
              locationf = (statef-1)*dimblock
      
               do ivec = 1,dimblock
                    vecout(ivec+locationf) = vecout(ivec+locationf) + xme*vecin(ivec+locationi)
               end do
            end do  !
         end do  ! xjmp
      end if
   end do ! ibundle
!$omp end parallel
   return
 end subroutine applyhNNbundled_block_omp

 
subroutine applyhPNbundled_block_omp (vchar,hchar,startbundle,endbundle )

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
   use contigpointervectors
   use localblocks
   use sectors

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

   integer(kind=8) :: bleft, bright
   integer :: fs
   integer :: ivec
   integer(kind=8) :: locationi,locationf

   if (vchar /= 'n') then
      print *, "Vchar PN can only be n not ", vchar
      stop 4
   end if
   
   vecin  => blockvec_1
   vecout => blockvec_2

! schedule is dynamic because runtime per step is variable
!$omp parallel private(ibundle, fs, ivec, bleft, bright)           &
!$omp     private(mythread)                          &
!$omp     private(pjmpstart, pjmpend, njmpstart, njmpend)       &
!$omp     private(numpthreads, numnthreads)                     &
!$omp     private(startp_thread, npjmps_thread)                 &
!$omp     private(pjmp, njmp, psdi, psdf)                       &
!$omp     private(a, b, c, d)                                   &
!$omp     private(locationi,locationf)&
!$omp     private(phasep, phasen, coplabel, doplabel)           &
!$omp     private(xme, nsdi, nsdf, statei, statef)              

   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'PN')cycle
      if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
	  
      mythread = omp_get_thread_num()
      
      fs = opbundle(ibundle)%fsector
      bleft = xsd(1)%sector(fs)%basisstart
      bright = xsd(1)%sector(fs)%basisend
!!!      if (bleft > mybasisstop(mythread) .or. bright < mybasisstart(mythread)) cycle
      
      !...... EXTRACT INFORMATION FROM OPBUNDLE ........
      !
      pjmpstart = opbundle(ibundle)%pxstart
      pjmpend   = opbundle(ibundle)%pxend
      njmpstart = opbundle(ibundle)%nxstart
      njmpend   = opbundle(ibundle)%nxend
      numpthreads = opbundle(ibundle)%numpthreads
      numnthreads = opbundle(ibundle)%numnthreads
      !
      if (hchar /= 'b') then

         !  starting position for proton 1-body jumps for this thread
         ! bit of a hack here to keep this routine similar to _orig version
         ! we take the complete thread range and ignore the divisions
         startp_thread = opbundle(ibundle)%startp_thread(0)     
         npjmps_thread = opbundle(ibundle)%startp_thread(numpthreads) - startp_thread

         ! KSM:  start/stop set up so that each proton final state appears on only one thread
         ! KSM:  prevents collison over update of voutp(statef) below
         !---------   Forward direction ------------------------------
         do pjmp = pjmpstart, pjmpend
            psdi = p1b_isd(pjmp)       ! initial proton slater determinant
            psdf = p1b_fsd(pjmp)       ! final proton SD
            phasep = p1b_phase(pjmp)   ! phase of proton jumps
            a = p1b_cop(pjmp)     ! KSM: Proton 1-body creation label
            c = p1b_dop(pjmp)     ! KSM: Proton 1-body destruction label
            !--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
            do njmp = njmpstart,njmpend
               nsdf = n1b_fsd(njmp)
               statef = nsdf + psdf                ! final state in combined basis
               if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
               
               !----------- FIND MATRIX ELEMTN --------------------------------------------
               b = n1b_cop(njmp)  ! KSM: Neutron 1-body creation label
               d = n1b_dop(njmp)  ! KSM: Neutron 1-body destruction label
               phasen = n1b_phase(njmp)
               coplabel = cpnpair(b,a)
               doplabel = dpnpair(d,c)
               xme = hmatpn(coplabel + doplabel)   ! get matrix element
               xme = xme*phasep*phasen             ! multiply matrix element by jump phases
               nsdi = n1b_isd(njmp)
               statei = nsdi + psdi                ! initial state in combined basis
               locationi = (statei-1)*dimblock
               locationf = (statef-1)*dimblock
           
               do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + xme*vecin(ivec+locationi)
               end do
               !vecout(statef) = vecout(statef) + xme*vecin(statei)
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
         
         do njmp = njmpstart,njmpend     
               nsdi = n1b_isd(njmp)
               nsdf = n1b_fsd(njmp)
               b  = n1b_cop(njmp)
               d  = n1b_dop(njmp)
               phasen = n1b_phase(njmp)
               do pjmp = pjmpstart,pjmpend
                  psdf = p1b_fsd(pjmp)       ! final proton SD
                  psdi = p1b_isd(pjmp)       ! initial proton slater determinant
                  statef = nsdf + psdf                  ! final state in combined basis
                  statei = nsdi + psdi                  ! initial state in combined basis
                  if (statei < mybasisstart(mythread) .or. statei > mybasisstop(mythread)) cycle

                  phasep = p1b_phase(pjmp)   ! phase of proton jumps
                  a  = p1b_cop(pjmp) 
                  c  = p1b_dop(pjmp)
                  !--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
                  !----------- FIND MATRIX ELEMENT -------------------------------------------
                  coplabel = cpnpair(b,a)
                  doplabel = dpnpair(d,c)	   
                  xme = hmatpn(coplabel + doplabel)     ! get matrix element
                  xme = xme*phasep*phasen               ! multiply matrix element by jump phases
                  locationi = (statei-1)*dimblock
                  locationf = (statef-1)*dimblock

                  do ivec = 1,dimblock
                    vecout(ivec+locationi) =vecout(ivec+locationi) + xme*vecin(ivec+locationf)
                  end do
                  !vecout(statei) = vecout(statei) + xme*vecin(statef)
               end do  ! pjmp
         end do  ! njmp
            
      end if
   end do  ! ibundle
   !$omp end parallel
   
   return
 end subroutine applyhPNbundled_block_omp

! Alternate version that threads across bundles, not slices
subroutine applySPEbundled_block_omp(vchar,startbundle,endbundle )
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
   use contigpointervectors
   use localblocks
   implicit none

   integer :: ibundle,startbundle,endbundle
   character(1) :: vchar
!------------------------------------------------------------

   integer(kind=8) nsdstart,nsdend,psdstart,psdend
   integer(kind=8) xjmp,xjmpstart,xjmpend

   integer(kind=8) :: statei, statef
!-------- OpenMP functions ---------------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: num_threads
   integer :: mythread
   integer(kind=basis_prec) :: vs

   real(4)    ::  pspe,nspe
   integer(kind=8) :: ip,in

   integer :: fs
   integer(kind=8) :: bleft, bright
   integer :: ivec
   integer(kind=8) :: locationi,locationf
   

   if (vchar /= 'n') then
      print *, "Vchar NN can only be n not ", vchar
      stop 4
   end if
   
   if(vchar == 'n')then
      vecin  => blockvec_1
      vecout => blockvec_2
   else
      vecin  => blockvec_2
      vecout => blockvec_1
   end if

! note: first private says that each such private copy is initialized with 
! the value before the parallel pragma

!$omp parallel private(ibundle,fs, bleft, bright, ip,in,pspe,nspe,statei, statef) & 
!$omp  private(psdstart, psdend, nsdstart, nsdend)        &
!$omp  private(locationi,locationf)&
!$omp  private(mythread)                              
   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'SPE')cycle
!	  if(opbundle(ibundle)%annexed)cycle
	  
      mythread = omp_get_thread_num()
      if (threadStart(ibundle, mythread) > threadStop(ibundle, mythread)) cycle

      fs = opbundle(ibundle)%fsector
      bleft = xsd(1)%sector(fs)%basisstart
      bright = xsd(1)%sector(fs)%basisend
!!!      if (bleft > mybasisstop(mythread) .or. bright < mybasisstart(mythread)) cycle

      psdstart = opbundle(ibundle)%pxstart
      psdend   = opbundle(ibundle)%pxend
      nsdstart = opbundle(ibundle)%nxstart
      nsdend   = opbundle(ibundle)%nxend

!.... LOOP OVER PROTON SDs in that sector..........
!      do ip = psdstart,psdend
      do ip = threadStart(ibundle, mythread), threadStop(ibundle, mythread)
           pspe = pspe_h(ip)   ! the proton contribution to the single-particle energies
           statei = pstart(ip)
           statef = pstart(ip) + nstart(nsdstart)
           locationi = (statei-1)*dimblock
           
!......... LOOP OVER NEUTRON SDS in sector jsc........................
           do in = nsdstart,nsdend
              nspe = nspe_h(in)  ! neutron contributions to s.p.e.
              !vecout(statef) = vecout(statef) + vecin(statef)*( pspe + nspe )  ! add spes

              locationf = (statef-1)*dimblock
               do ivec = 1,dimblock
                    vecout(ivec+locationf) =vecout(ivec+locationf) + vecin(ivec+locationf)*( pspe + nspe )
               enddo
               statef = statef+1   ! neutron SDs are contiguous
            end do  !in
      end do  !ip
   end do ! ibundle
!$omp end parallel 
   return
end subroutine applySPEbundled_block_omp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module applyhamblockomp


