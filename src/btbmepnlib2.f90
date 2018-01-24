!==========================================================================
! BTBMEPNLIB2.f90
!
! second set of routines for handling pn, ppn, and pnn
! intended specifically for parallel MPI applications
! when Lanczos vectors are fragmented
!
! initialized July 2014 by CWJ @ SDSU
!
! Basic idea: between fragments, identify allowed Delta M, Delta Pi, and Delta W
! This will limit the matrix elements which can be stored
! 
!==========================================================================
!
!  routine to determine change in quantum numbers between sectors
!
!  INPUT:
!     it:  species
!     is,fs: initial, final sector
!
!  OUTPUT:
!     dM = change in M
!     dPi= change in parity (pi)
!     dW = change in W (weighting for truncation)
!
subroutine deltaqsector(it,is,fs,dM,dPi,dW)
   use sectors
   implicit none
!..... INPUT......
   integer :: it  ! species
   integer :: is, if ! initial, final sectors

!..... OUTPUT......
   integer dM, dPi, dW  ! change in quantum numbers
!...... FUNCTIONS CALLED....
!! interface --   integer parmult

   dM = (xsd(it)%sector(fs)%jzX-xsd(it)%sector(fs)%jzX)/2
   dW =  xsd(it)%sector(fs)%wX -xsd(it)%sector(fs)%wX
   dPi=  parmult (  xsd(it)%sector(fs)%parX, xsd(it)%sector(is)%parX)

   return
end subroutine deltaqsector

!==========================================================================
!
! routine to aggregate changes in quantum numbers between fragments
! also applies when not fragmented -- i.e., the entire space
!
! INPUT: 
!   ifrag, ffrag:  indices for initial, final fragments
!
! SUBROUTINES CALLED:
!   deltaqsector
!

subroutine alldeltaqfragment(ifrag,ffrag)
   use fragments
   use sectors
   use interaction

   implicit none
!.... INPUT.......
   integer :: ifrag, ffrag  ! indices of initial and final fragments
!.....INTERNAL.......

   integer :: isp,fsp, isn, fsn   ! initial, final sectors for protons, neutrons



   return
end subroutine alldeltaqfragment
!==========================================================================
!
!
