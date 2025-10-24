!=====================================================================
!
!   BDENSLIB3.f90
!
!  reorganization in 7.6.7:
!  here are the routines for actually carrying out operations

!====================================================================
!
! based upon applyhPPbundled_g (as of 7.3.7)
! but with order of loops switched
! no OpenMP; should be fast enough to be not needed 
!
subroutine applyP1Bdenbundled_g (startbundle,endbundle )

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
  use densities
  use bmpi_mod
  implicit none

  integer :: ibundle
  character(1) :: hchar,vchar
  integer :: startbundle,endbundle
! --------- NOTE: basestart, basestop stored in module fragments

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=obs_prec)   xme,xmetr
  integer(kind=8) :: statei, statef,nsd
  integer(kind=8) :: statei0,statef0
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(4) :: mythread
  integer :: a,b

  if(startbundle == 0 .or. endbundle==0)return
  
  call clocker('p1b','sta')

  do ibundle = startbundle,endbundle
!	  print*,ibundle, opbundle(ibundle)%optype
     if(opbundle(ibundle)%optype /= 'P1B')cycle
  hchar = opbundle(ibundle)%hchar
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
  csdstart = opbundle(ibundle)%nxstart
  csdend   = opbundle(ibundle)%nxend
  xjmpstart = opbundle(ibundle)%pxstart
  xjmpend   = opbundle(ibundle)%pxend
  cstride   = opbundle(ibundle)%cstride

  ncstates = (csdend +cstride -csdstart)/cstride

!--------- OUTER LOOP OVER PROTON JUMPS ---------

  istart = 1
  iend = ncstates
  csd_index = csdstart + (istart - 1)*cstride - cstride
  if(istart <= iend)then
  csd_index = nstart(csd_index+1)-1

! KSM: get start point for neutron slater det.   This is the beginning
! KSM: of a diagonal block.
     select case (hchar)  ! choose whether 'forward', 'backwards' or 'hermitian '

     case ('f','h')
!--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............
! OMP ADDED 7.9.5; 
! OMP disabled 7.9.7; appears to cause problems in some W-truncated cases see note Sept 7, 2020
! however reinstated 7.11.0; appears to work
!$OMP parallel do private(a,b,statei,statef,statei0,statef0,csd,xme), reduction(+:p1bopme),firstprivate(ncstates,csd_index)
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
        a    = p1b_cop(xjmp) 
        b    = p1b_dop(xjmp)
!--------- GET PHASE
        xme = 0.0
        statei0 = p1b_isd(xjmp)+csd_index   ! CWJ initial PSD + initial NSD
        statef0 = p1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial NSD
!------------- LOOP OVER CONJUGATE SDs ------------------------
!  In this case, OpenMP seems to  just slow things down
!!$OMP parallel do private(statei,statef), reduction(+:xme)
           do csd = 1,ncstates
              statei = statei0+cstride*(csd)
              statef = statef0+cstride*(csd)
              xme = xme + real(vec1(statei),obs_prec)*real(vec2(statef),obs_prec)
          enddo  ! csd
          p1bopme(a,b) = p1bopme(a,b)+ xme*p1b_phase(xjmp)
    end do ! xjmp

    case('b')
   !--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............
!... ALSO DISABLED in 7.9.7
!$OMP parallel do private(a,b,statei,statef,statei0,statef0,csd,xme), reduction(+:p1bopme),firstprivate(ncstates,csd_index)
   
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
        a    = p1b_cop(xjmp) 
        b    = p1b_dop(xjmp)
!--------- GET PHASE
        xme = 0.0
        statei0 = p1b_isd(xjmp)+csd_index   ! CWJ initial PSD + initial NSD
        statef0 = p1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial NSD
!------------- LOOP OVER CONJUGATE SDs ------------------------
!!$OMP parallel do private(statei,statef), reduction(+:xme)
        do csd = 1,ncstates
              statei = statei0+cstride*(csd)
              statef = statef0+cstride*csd			
              xme = xme + real(vec1(statef),obs_prec)*real(vec2(statei),obs_prec)

        enddo  ! csd
        p1bopme(b,a) = p1bopme(b,a)+ xme*p1b_phase(xjmp)
    end do ! xjmp

!   case ('h')
!--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............
!     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
!        a    = p1b_cop(xjmp) 
!        b    = p1b_dop(xjmp)
!--------- GET PHASE
!        xme = 0.0
!        xmetr = 0.0
!        statei = p1b_isd(xjmp)+csd_index   ! CWJ initial PSD + initial NSD
!        statef = p1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial NSD
!           do csd = 1,ncstates
!              statei = statei+cstride
!              statef = statef+cstride
!              xme = xme + real(vec1(statei),obs_prec)*real(vec2(statef),obs_prec)
!           enddo  ! csd
!           p1bopme(a,b) = p1bopme(a,b)+ xme*p1b_phase(xjmp)
!    end do ! xjmp
	
    end select
  end if
  end do ! ibundle
  
  call clocker('p1b','end')
  
  return
end subroutine applyP1Bdenbundled_g
!===================================================================
!
! based upon applyhNNbundled_g (as of 7.3.7)
! but with order of loops switched
! OpenMP added in 7.5.5 for cases with many states
!
subroutine applyN1Bdenbundled_g (startbundle,endbundle )

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
  use densities
  use bmpi_mod
  implicit none

  integer :: ibundle
  character(1) :: hchar,vchar
  integer :: startbundle,endbundle
! --------- NOTE: basestart, basestop stored in module fragments

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index,pstride
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=obs_prec)   xme,xmetr
  integer(kind=8) :: statei, statef,nsd,statei0,statef0
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(4) :: mythread
  integer :: a,b

  if(startbundle == 0 .or. endbundle==0)return
  call clocker('n1b','sta')

  pstride = 0

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'N1B')cycle
  hchar = opbundle(ibundle)%hchar
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
  csdstart = opbundle(ibundle)%pxstart
  csdend   = opbundle(ibundle)%pxend
  xjmpstart = opbundle(ibundle)%nxstart
  xjmpend   = opbundle(ibundle)%nxend
  cstride   = opbundle(ibundle)%cstride
!  print*,ibundle,csdstart,cstride
  ncstates = (csdend +cstride -csdstart)/cstride

  if(ncstates > 1)pstride   = pstart(csdstart+cstride)-pstart(csdstart)


  istart = 1
  iend = ncstates
  if(istart <= iend)then
  csd_index = csdstart
  csd_index = pstart(csd_index)
  
!!! =========== OPENMP added in 7.5.5=============

! firstprivate gives each thread its own copy, but initializes it
!    better than shared for read-only vars
! private gives each thread its own copy
! KSM: get start point for neutron slater det.   This is the beginning
! KSM: of a diagonal block.
     select case (hchar)  ! choose whether 'forward', 'backwards' or 'hermitian '
     case ('f','h')
!--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............
!--------- OUTER LOOP OVER NEUTRON JUMPS ---------
! disabled in 7.9.7 to avoid crashes
!$OMP parallel do private(a,b,statei,statef,statei0,statef0,csd,xme), reduction(+:n1bopme),firstprivate(csd_index,ncstates)
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
        a    = n1b_cop(xjmp) 
        b    = n1b_dop(xjmp)
!--------- GET PHASE
        xme = 0.0
        statei0 = n1b_isd(xjmp)+csd_index   ! CWJ initial NSD + initial PSD
        statef0 = n1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial PSD

!------------- LOOP OVER CONJUGATE PROTON SDs ------------------------
!!$OMP parallel do private(statei,statef), reduction(+:xme)

        do csd = 1,ncstates
            statei = statei0+pstride*(csd-1)
            statef = statef0+pstride*(csd-1)
              xme = xme + real(vec1(statei),obs_prec)*real(vec2(statef),obs_prec)

        enddo  ! csd
        n1bopme(a,b) = n1bopme(a,b)+ xme*n1b_phase(xjmp)
    end do ! xjmp

    case('b')
!--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............
! OMP ADDED 7.9.5
! OMP disabled 7.9.7 to address bug in some cases  see note Sept 7, 2020
! however reinstated in 7.11.0; not clear if bug is fixed...
!$OMP parallel do private(a,b,statei,statef,statei0,statef0,csd,xme), reduction(+:n1bopme),firstprivate(csd_index,ncstates)
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
        a    = n1b_cop(xjmp) 
        b    = n1b_dop(xjmp)
!--------- GET PHASE
        xme = 0.0
        statei0 = n1b_isd(xjmp)+csd_index   ! CWJ initial PSD + initial NSD
        statef0 = n1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial NSD
!------------- LOOP OVER CONJUGATE SDs ------------------------
!!$OMP parallel do private(statei,statef), reduction(+:xme)
        do csd = 1,ncstates
            statei = statei0+pstride*(csd-1)
            statef = statef0+pstride*(csd-1)
              xme = xme + real(vec1(statef),obs_prec)*real(vec2(statei),obs_prec)

        enddo  ! csd
        n1bopme(b,a) = n1bopme(b,a)+ xme*n1b_phase(xjmp)
		
    end do ! xjmp

!   case ('h')
!--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............

!     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
!        a    = n1b_cop(xjmp) 
!        b    = n1b_dop(xjmp)
!--------- GET PHASE
!        xme = 0.0
!        statei = n1b_isd(xjmp)+csd_index   ! CWJ initial PSD + initial NSD
!        statef = n1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial NSD
!------------- LOOP OVER CONJUGATE SDs ------------------------

!           do csd = 1,ncstates
!              xme = xme + real(vec1(statei),obs_prec)*real(vec2(statef),obs_prec)
!              statei = statei+pstride
!              statef = statef+pstride
!           enddo  ! csd
!           n1bopme(a,b) = n1bopme(a,b)+ xme*n1b_phase(xjmp)
!    end do ! xjmp
	
    end select

  end if
  end do ! ibundle
  
  call clocker('n1b','end')
  
  return
end subroutine applyN1Bdenbundled_g
!====================================================================
!
! based upon applyP1bopbundled_g (as of 7.367)
! but with order of loops switched
! no OpenMP; should be fast enough to be not needed 
! to replace routine below
!
subroutine applyP1Bopbundled_g (startbundle,endbundle )

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
  use densities
  use bmpi_mod
  use onebodypot
  implicit none

  integer :: ibundle
  character(1) :: hchar,vchar
  integer :: startbundle,endbundle
! --------- NOTE: basestart, basestop stored in module fragments

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=obs_prec)   xme,xmetr
  integer(kind=8) :: statei, statef,nsd
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(4) :: mythread
  integer :: a,b
  if(startbundle == 0 .or. endbundle==0)return

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'P1B')cycle
  hchar = opbundle(ibundle)%hchar  
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
  csdstart = opbundle(ibundle)%nxstart
  csdend   = opbundle(ibundle)%nxend
  xjmpstart = opbundle(ibundle)%pxstart
  xjmpend   = opbundle(ibundle)%pxend
  cstride   = opbundle(ibundle)%cstride

  ncstates = (csdend +cstride -csdstart)/cstride

!--------- OUTER LOOP OVER PROTON JUMPS ---------

  istart = 1
  iend = ncstates
  csd_index = csdstart + (istart - 1)*cstride - cstride
  if(istart <= iend)then
  csd_index = nstart(csd_index+1)-1

! KSM: get start point for neutron slater det.   This is the beginning
! KSM: of a diagonal block.
     select case (hchar)  ! choose whether 'forward', 'backwards' or 'hermitian '

     case ('f','h')
!--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............

     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
        a    = p1b_cop(xjmp) 
        b    = p1b_dop(xjmp)
!--------- GET PHASE

        xme =  pPOT_h(a,b)

   !--------- GET PHASE
        xme = xme*p1b_phase(xjmp)
!        xmetr = xmetr*p1b_phase(pjmp)		
		
        statei = p1b_isd(xjmp)+csd_index   ! CWJ initial PSD + initial NSD
        statef = p1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial NSD
!------------- LOOP OVER CONJUGATE SDs ------------------------
           do csd = 1,ncstates
              statei = statei+cstride
              statef = statef+cstride
			  vec2(statef)=vec2(statef)+ xme*vec1(statei)
          enddo  ! csd
    end do ! xjmp

    case('b')
   !--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
        a    = p1b_cop(xjmp) 
        b    = p1b_dop(xjmp)
!--------- GET PHASE
        xme = pPot_h(b,a)*p1b_phase(xjmp)
        statei = p1b_isd(xjmp)+csd_index   ! CWJ initial PSD + initial NSD
        statef = p1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial NSD
!------------- LOOP OVER CONJUGATE SDs ------------------------
        do csd = 1,ncstates
              statei = statei+cstride
              statef = statef+cstride			  
			  vec2(statei)= vec2(statei)+xme*vec1(statef)
        enddo  ! csd
    end do ! xjmp
	
    end select
  end if
  end do ! ibundle
  return
end subroutine applyP1Bopbundled_g
!===================================================================
!
! based upon applyN1Bopbundled_g (as of 7.6.7)
! but with order of loops switched
!
! to replace routine below
!
subroutine applyN1Bopbundled_g (startbundle,endbundle )

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
  use densities
  use bmpi_mod
  use onebodypot
  use jumplimits

  implicit none

  integer :: ibundle
  character(1) :: hchar,vchar
  integer :: startbundle,endbundle
! --------- NOTE: basestart, basestop stored in module fragments

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index,pstride
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=obs_prec)   xme,xmetr
  integer(kind=8) :: statei, statef,nsd
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(4) :: mythread
  integer :: a,b
  integer :: ierr

  pstride = 0
  if(startbundle == 0 .or. endbundle==0)return

  
  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'N1B')cycle
  hchar = opbundle(ibundle)%hchar
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
  csdstart = opbundle(ibundle)%pxstart
  csdend   = opbundle(ibundle)%pxend
  xjmpstart = opbundle(ibundle)%nxstart
  xjmpend   = opbundle(ibundle)%nxend
  cstride   = opbundle(ibundle)%cstride
  ncstates = (csdend +cstride -csdstart)/cstride

  if(ncstates > 1)pstride   = pstart(csdstart+cstride)-pstart(csdstart)


  istart = 1
  iend = ncstates
  if(istart <= iend)then
  csd_index = csdstart
  csd_index = pstart(csd_index)
  
!!! =========== OPENMP added in 7.5.5=============

! firstprivate gives each thread its own copy, but initializes it
!    better than shared for read-only vars
! private gives each thread its own copy
! KSM: get start point for neutron slater det.   This is the beginning
! KSM: of a diagonal block.
!print*,' (X2 )',iproc,ibundle,hchar


     select case (hchar)  ! choose whether 'forward', 'backwards' or 'hermitian '
     case ('f','h')
!--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............
!--------- OUTER LOOP OVER NEUTRON JUMPS ---------

!print*,' (X2 )',iproc,ibundle
!call BMPI_BARRIER(icomm,ierr)

     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
        a    = n1b_cop(xjmp) 
        b    = n1b_dop(xjmp)
!--------- GET PHASE
        xme =  nPOT_h(a,b) * n1b_phase(xjmp)
        statei = n1b_isd(xjmp)+csd_index   ! CWJ initial NSD + initial PSD
        statef = n1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial PSD
!------------- LOOP OVER CONJUGATE PROTON SDs ------------------------
        do csd = 1,ncstates
			  vec2(statef)=vec2(statef)+xme*(vec1(statei))
!			  if(statef==28)print*,iproc,' (A n)',xme,ibundle,xjmp,vec1(statei),vec2(statef)
			  
              statei = statei+pstride
              statef = statef+pstride
        enddo  ! csd
    end do ! xjmp

    case('b')
!--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............
!print*,' (X3 )',iproc,ibundle

     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
        a    = n1b_cop(xjmp) 
        b    = n1b_dop(xjmp)
!--------- GET PHASE
        xme =  nPOT_h(b,a) * n1b_phase(xjmp)
        statei = n1b_isd(xjmp)+csd_index   ! CWJ initial PSD + initial NSD
        statef = n1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial NSD
!------------- LOOP OVER CONJUGATE SDs ------------------------
        do csd = 1,ncstates
		      vec2(statei)=vec2(statei)+xme*(vec1(statef))
!			  if(statei==28)print*,iproc,' (B n)',xme,statef,vec1(statef),vec2(statei)
			
!              xme = xme + real(vec1(statef),obs_prec)*real(vec2(statei),obs_prec)
              statei = statei+pstride
              statef = statef+pstride
        enddo  ! csd
!        n1bopme(b,a) = n1bopme(b,a)+ xme*n1b_phase(xjmp)
    end do ! xjmp

!   case ('h')
!--------- LOOP OVER 1-BODY JUMPS IN THIS SECTOR JUMPS.............

!     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT
!        a    = n1b_cop(xjmp) 
!        b    = n1b_dop(xjmp)
!--------- GET PHASE
!        xme = 0.0
!        statei = n1b_isd(xjmp)+csd_index   ! CWJ initial PSD + initial NSD
!        statef = n1b_fsd(xjmp)+csd_index   ! CWJ final PSD + initial NSD
!------------- LOOP OVER CONJUGATE SDs ------------------------

!           do csd = 1,ncstates
!              xme = xme + real(vec1(statei),obs_prec)*real(vec2(statef),obs_prec)
!              statei = statei+pstride
!              statef = statef+pstride
!           enddo  ! csd
!           n1bopme(a,b) = n1bopme(a,b)+ xme*n1b_phase(xjmp)
!    end do ! xjmp
	
    end select

  end if
  end do ! ibundle
  return
end subroutine applyN1Bopbundled_g
!======================================================
