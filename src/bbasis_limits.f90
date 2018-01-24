!===========================================================================
!  BIGSTICK configuration-interaction shell-model code
!
!===========================================================================
!
!  bbasislib5.f90
!
!  fifth set of routines to generate the basis:  sets up for sectors of SDs
!  Determines the quantum # limits for SDs generated by combining blocks of haikus
!  sets up construct limxsd
!
!===========================================================================
!  SUBROUTINES:
!     mastersdlimits: master routine for limits on q#s for SDs
!     sdparlimits:  finds limits on parity of SDs for given species
!     reconcile_sdpar: compares parity limits for 2 species
! 
!     sdjzlimits:  finds limits on jzs for given species
!     reconcile_sdjz: compares jz limits for 2 species
! 
!     sdwlimits:  finds limits on W for given species
!     reconcile_sdw: compares w limits for 2 species
! 
!
!===========================================================================

!---------------------------------------------------------------------------
! limits on quantum numbers useful for constructing SDs, basis
!---------------------------------------------------------------------------
module slaterlimits

  implicit none
  
!--------------------- NEW DERIVED TYPES TO REPLACE OLD---------------------

  type limsd3
     integer          :: nw
     integer, pointer :: wlist(:)
     integer, pointer :: wmap(:)
     integer, pointer :: nhmin(:),nhmax(:)  
  end type limsd3

  type limsd2
     integer jzmin,jzmax
     type (limsd3), pointer :: jz(:)
  end type limsd2
  
  type limsd1
     integer       :: parmin,parmax
     type (limsd2) :: par(2)
  end type limsd1
  
  type (limsd1) :: limxsd(2)

  contains
	  
!      subroutine mastersdlimits
!
!  SUBROUTINES CALLED:
!   sdparlimits -- computes parity limits for SDs for each species
!   reconcile_sdpar -- looks at combined parities and changes limits if needed
!
!===========================================================================
subroutine mastersdlimits
  use system_parameters
  use sporbit
  use haiku_info
  use blocks

  implicit none
  integer :: it
  integer :: jzp,jzn
  integer :: par,parp,parn
  integer minwx,maxwx

  call sdparlimits(1)
  call sdparlimits(2)

  call reconcile_sdpar

  do parp = limxsd(1)%parmin,limxsd(1)%parmax

      call sdjzlimits(1,parp)

      parn = parmult(iparity,parp)
      call sdjzlimits(2,parn)

      call reconcile_sdjz(parp)

      do jzp = limxsd(1)%par(parp)%jzmin, limxsd(1)%par(parp)%jzmax,2
        jzn = Jz - jzp
        call sdwlimits(1,parp,jzp,minwx,maxwx,.false.)
        call sdwlimits(1,parp,jzp,minwx,maxwx,.true.)
        call sdwlimits(2,parn,jzn,minwx,maxwx,.false.)
        call sdwlimits(2,parn,jzn,minwx,maxwx,.true.)
        call reconcile_sdw(parp,jzp)

      enddo  ! jzp

  enddo ! par

  return
end subroutine mastersdlimits


!===========================================================================
!      subroutine sdparlimits
!
!  finds, for a given species, what are the min, max of parity 
!  
!  INPUT: it = species
!
!===========================================================================
subroutine sdparlimits(it)

  use system_parameters
  use sporbit
  use haiku_info
  use blocks
  use butil_mod

  implicit none
  integer :: it
  integer :: N
  integer :: nh,nhc
  integer :: parmin,parmax,parh,parhc
  integer :: parhmin,parhmax,parhcmin,parhcmax
  integer :: par

  if ( it < 0 ) then  ! error trap
     print*,' Idiot! Variable it should be > 0 '
     stop
  endif

  N = np(it)

  parmin = 2
  parmax = 1

  do nh = blockmap(it)%nhmin,blockmap(it)%nhmax
   nhc = N - nh
   if(nhc > blockmap(-it)%nhmax .or. nhc < blockmap(-it)%nhmin)cycle
   parhmin = blockmap(it)%nh(nh)%parmin
   parhmax = blockmap(it)%nh(nh)%parmax
   parhcmin = blockmap(-it)%nh(nhc)%parmin
   parhcmax = blockmap(-it)%nh(nhc)%parmax

   if(parhmin /= parhmax .or. parhcmin /= parhcmax)then
     parmin = 1
     parmax = 2
   else
     par = parmult(parhmin,parhcmin)
     parmin = bmin(parmin,par)
     parmax = bmax(parmax,par)
   endif
  enddo ! nh

  limxsd(it)%parmin = parmin
  limxsd(it)%parmax = parmax
  return
end subroutine sdparlimits
!=========================================================================
!  subroutine reconcile_sdpar
!  reconciles parity limits between protons and neutrons
!
subroutine reconcile_sdpar
  use system_parameters
  use sporbit
  use haiku_info
  use blocks

  implicit none

  integer par
  integer partmp

  par = iparity
! IF BOTH SPECIES ADMIT BOTH PARITIES NO PROBLEM
  if(limxsd(1)%parmin /= limxsd(1)%parmax .and. limxsd(2)%parmin/=limxsd(2)%parmax)return

! FIND OUT WHICH ONE IS LIMITED

  if(limxsd(1)%parmin == limxsd(1)%parmax)then
     partmp = parmult(par, limxsd(1)%parmin)
!........... ERROR TRAP .....................
     if( partmp < limxsd(2)%parmin .or. partmp > limxsd(2)%parmax)then
        print*,' Cannot make requested parity (1)'
        stop
     endif
     limxsd(2)%parmin = partmp
     limxsd(2)%parmax = partmp

  else
     partmp = parmult(par, limxsd(2)%parmin)
!........... ERROR TRAP .....................
     if( partmp < limxsd(1)%parmin .or. partmp > limxsd(1)%parmax)then
        print*,' Cannot make requested parity (2) '
        stop
     endif
     limxsd(1)%parmin = partmp
     limxsd(1)%parmax = partmp
  endif
  return

end subroutine reconcile_sdpar

!=========================================================================
!subroutine sdjzlimits
subroutine sdjzlimits(it,parx)

  use system_parameters
  use sporbit
  use haiku_info
  use blocks
  use butil_mod

  implicit none
  integer :: it
  integer :: N
  integer :: nh,nhc
  integer :: jzmax,jzmin
  integer :: jzh,jzhc
  integer :: parmin,parmax,parh,parhc

  integer :: parx
  integer :: aerr

  if ( it < 0 ) then  ! error trap
     print*,' Idiot! Variable it should be > 0 '
     stop
  endif

  N = np(it)

  jzmin = 1000
  jzmax = -1000
!  print*,' min max ',blockmap(it)%nhmin,blockmap(it)%nhmax
  do nh = blockmap(it)%nhmin,blockmap(it)%nhmax
     nhc = N - nh
     if(nhc < blockmap(-it)%nhmin .or. nhc > blockmap(-it)%nhmax)cycle
     do parh = blockmap(it)%nh(nh)%parmin,blockmap(it)%nh(nh)%parmax
       parhc = parmult(parx,parh)
       if(parhc > blockmap(-it)%nh(nhc)%parmax .or. parhc < blockmap(-it)%nh(nhc)%parmin)cycle
       do jzh = blockmap(it)%nh(nh)%par(parh)%jzmin,blockmap(it)%nh(nh)%par(parh)%jzmax,2
!   print*,it,nhc,parhc, blockmap(-it)%nh(nhc)%par(parhc)%jzmin, blockmap(-it)%nh(nhc)%par(parhc)%jzmax
          do jzhc = blockmap(-it)%nh(nhc)%par(parhc)%jzmin, blockmap(-it)%nh(nhc)%par(parhc)%jzmax, 2
            jzmin = bmin(jzmin,jzh+jzhc)
            jzmax = bmax(jzmax,jzh+jzhc)
!            print*,it,nh,nhc,parh,parhc,jzh,jzhc
          enddo ! jzhc
       enddo ! jzh
     enddo ! parh
  enddo ! nh

  limxsd(it)%par(parx)%jzmin = jzmin
  limxsd(it)%par(parx)%jzmax = jzmax
  if(jzmin <= jzmax)then
!    print*,jzmin,jzmax
    if(associated( limxsd(it)%par(parx)%jz ) )nullify( limxsd(it)%par(parx)%jz)
    allocate(limxsd(it)%par(parx)%jz(jzmin:jzmax), stat=aerr)
    if(aerr /= 0) call memerror("sdjzlimits")
  endif
  return
end subroutine sdjzlimits

!=========================================================================
!subroutine reconcile_sdjz

subroutine reconcile_sdjz(parp)
   use system_parameters
   use haiku_info
   use blocks
   use butil_mod
   use sporbit
   implicit none

   integer parp
   integer :: parn
   integer :: par


   parn = parmult(parp,iparity)

   limxsd(1)%par(parp)%jzmax = bmin( limxsd(1)%par(parp)%jzmax, jz - limxsd(2)%par(parn)%jzmin)
   limxsd(2)%par(parn)%jzmax = bmin( limxsd(2)%par(parn)%jzmax, jz - limxsd(1)%par(parp)%jzmin)

   limxsd(1)%par(parp)%jzmin = bmax( limxsd(1)%par(parp)%jzmin, jz - limxsd(2)%par(parn)%jzmax)
   limxsd(2)%par(parn)%jzmin = bmax( limxsd(2)%par(parn)%jzmin, jz - limxsd(1)%par(parp)%jzmax)

   return
end subroutine reconcile_sdjz

!=========================================================================
!subroutine sdwlimits

subroutine sdwlimits(it,parx,jzx,minwx,maxwx,fill)

use system_parameters
use W_info
use haiku_info
use blocks
use butil_mod
use sporbit

implicit none
integer :: it
integer :: parx,jzx
integer :: minwx,maxwx
logical :: fill
!..............................
integer N
integer :: nh,nhc
integer :: parh,parhc
integer :: jzh,jzhc
!.......................
integer :: iwh,iwhc,w,wh,whc
integer :: nw

logical, allocatable :: warray(:)
integer :: aerr

  if ( it < 0 ) then  ! error trap
     print*,' Idiot! Variable it should be > 0 '
     stop
  endif

  N = np(it)

  if(.not.fill)then
    minwx = 1000
    maxwx = 0
	 allocate(warray(1:1), stat=aerr) ! to suppress uninitialized msg from compiler
    if(aerr /= 0) then
       call memerror("sdwlimits 0")
       stop 5
    end if
  else
    allocate(warray(minwx:maxwx), stat=aerr)
    if(aerr /= 0) then
       call memerror("sdwlimits 1")
       stop 5
    end if
    warray(:) = .false.
  endif

  do nh = blockmap(it)%nhmin,blockmap(it)%nhmax
     nhc = N - nh
     if(nhc < blockmap(-it)%nhmin .or. nhc > blockmap(-it)%nhmax)cycle

     do parh = blockmap(it)%nh(nh)%parmin,blockmap(it)%nh(nh)%parmax
       parhc = parmult(parx,parh)
       if(parhc > blockmap(-it)%nh(nhc)%parmax .or. parhc < blockmap(-it)%nh(nhc)%parmin)cycle

       do jzh = blockmap(it)%nh(nh)%par(parh)%jzmin,blockmap(it)%nh(nh)%par(parh)%jzmax,2
         jzhc = jzx - jzh
         if(jzhc > blockmap(-it)%nh(nhc)%par(parhc)%jzmax .or. jzhc < blockmap(-it)%nh(nhc)%par(parhc)%jzmin)cycle
         do iwh = 1,blockmap(it)%nh(nh)%par(parh)%jz(jzh)%nw
            wh = blockmap(it)%nh(nh)%par(parh)%jz(jzh)%wlist(iwh) 
            do iwhc = 1,blockmap(-it)%nh(nhc)%par(parhc)%jz(jzhc)%nw
               whc = blockmap(-it)%nh(nhc)%par(parhc)%jz(jzhc)%wlist(iwhc) 
               w= wh+whc
               if(w > maxW(it))cycle
               if(.not.fill)then
                 minwx = bmin(minwx,w)
                 maxwx = bmax(maxwx,w)
               else
                 warray(w) = .true.
               endif

            enddo ! iwhc
         enddo ! iwh

       enddo !jzh
     enddo  ! parh
  enddo ! nh

  if(fill)then
    nw = 0
    do w = minwx,maxwx
      if(warray(w))nw = nw+1
    enddo
    limxsd(it)%par(parx)%jz(jzx)%nw = nw
    if( associated( limxsd(it)%par(parx)%jz(jzx)%wlist) ) & 
              nullify( limxsd(it)%par(parx)%jz(jzx)%wlist )
    allocate(limxsd(it)%par(parx)%jz(jzx)%wlist(nw), stat=aerr)
    if(aerr /= 0) call memerror("sdwlimits 2")
    nw = 0
    do w = minwx,maxwx
      if(warray(w))then
         nw = nw+1
         limxsd(it)%par(parx)%jz(jzx)%wlist(nw)=w
      endif
    enddo
  endif
  deallocate(warray)
  
return
end subroutine sdwlimits

!=========================================================================
! subroutine reconcile_wlimits

subroutine reconcile_sdw(parp,jzp)
use system_parameters
use haiku_info
use blocks
use W_info
use system_parameters
use sporbit

implicit none
integer :: parp,parn
integer :: jzp,jzn
integer :: nwp,nwn
integer :: wp,wn
logical :: okay


parn = parmult(parp, iparity)
jzn = jz - jzp
!print*,jzp,jzn
nwp = limxsd(1)%par(parp)%jz(jzp)%nw
nwn = limxsd(2)%par(parn)%jz(jzn)%nw

if(nwp == 0)limxsd(2)%par(parn)%jz(jzn)%nw = 0
if(nwn == 0)limxsd(1)%par(parp)%jz(jzp)%nw = 0

if(nwp == 0 .or. nwn==0)return

wp = limxsd(1)%par(parp)%jz(jzp)%wlist(1)
wn = limxsd(2)%par(parn)%jz(jzn)%wlist(1)

if(wp+wn > maxWtot)then
   limxsd(1)%par(parp)%jz(jzp)%nw = 0
   limxsd(2)%par(parn)%jz(jzn)%nw = 0

   return
!   print*,' We have a W problem ',parp,jzp,wp,wn
!   stop
endif

wp = limxsd(1)%par(parp)%jz(jzp)%wlist(1)
wn = limxsd(2)%par(parn)%jz(jzn)%wlist(nwn)

if(wp+wn > maxWtot)then
  okay = .false.

  do while(.not.okay)
     nwn = nwn-1
     wn = limxsd(2)%par(parn)%jz(jzn)%wlist(nwn)
     if(wp+wn <= maxWtot)okay = .true.
  enddo
  limxsd(2)%par(parn)%jz(jzn)%nw = nwn
endif


wp = limxsd(1)%par(parp)%jz(jzp)%wlist(nwp)
wn = limxsd(2)%par(parn)%jz(jzn)%wlist(1)

if(wp+wn > maxWtot)then
  okay = .false.

  do while(.not.okay .and. nwp > 1)
     nwp = nwp-1
     wp = limxsd(1)%par(parp)%jz(jzp)%wlist(nwp)
     if(wp+wn <= maxWtot)okay = .true.
  enddo
  limxsd(1)%par(parp)%jz(jzp)%nw = nwp
endif

return
end subroutine reconcile_sdw

!=========================================================================
end module slaterlimits