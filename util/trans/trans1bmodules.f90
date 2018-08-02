!===========================================================================
! MODULES
!===========================================================================

!...........................................................................
module menu_choices
  implicit none
  character(len=2)  :: menu_char
  integer ichoice
end module menu_choices


!---------------------------------------------------------------------------
!  MODULE VERBOSITY
!---------------------------------------------------------------------------
!  Various flags for printing out information.
!  This allows one to easily print out or suppress information.
!
!  NOTE: verbose_X  means write out information to terminal
!        chatty_X  means write out information as being computed (rare)
!        print_X    means write to file
!        X_file     is unit number for file
!---------------------------------------------------------------------------
module verbosity
  implicit none

  logical :: print4modelinfo   = .false.  ! prints out information needed for modeling
  integer :: modelinfo_file    =  44
  logical :: verbose_winnow    = .false.
  logical :: print_groups      = .false.
  integer :: group_file        = 37
  logical :: print_hsps        = .false.
  integer :: hsps_file         = 37
  logical :: verbose_templates =  .false. 
  logical :: verbose_haikus    = .false.
  logical :: print_haikus      = .false.  
  integer :: haiku_file        =  37
  logical :: verbose_blocks    = .false.
  logical :: verbose_hops      = .false.
  logical :: print_hops        = .false.
  integer :: hop_file          =  38 
  logical :: verbose_sectors   = .false.
  logical :: print_sectors     = .false.
  integer :: sector_file       = 37  
  logical :: chatty_descent    = .false.
  logical :: chatty_links      = .false.

  logical :: print_genes       = .false.
  integer :: gene_file         = 38
  logical :: print_basis       = .false.
  integer :: basis_file        = 31

  logical :: print_matrix      = .false.
  integer :: matrix_file       = 55
  logical :: verbose_orthog    = .false.

  logical :: print_block_ops   = .false.
  logical :: print_node_load   = .false.

  logical :: print_sectorjumps = .false. !.true.    ! added 5/2010 CWJ
  integer :: sectorjump_file   = 39

end module verbosity

!.....................................................................


!--------------------------------------------------------------------------
!  single-particle ORBIT information
!  key derived variable: orbqn
!---------------------------------------------------------------------------
module sporbit
   implicit none
   integer :: numorb(2)  ! # of s.p. orbits
   integer :: numorbmax	        ! max of numorb
   logical :: isoflag           ! in isospin formalism
   logical :: spinless          ! for "spinless" fermions
   logical :: allsameparity     ! flag for all sp orbits same parity
   logical :: allsameW          ! flag for all sp orbits same W
!------------ CREATE A DEFINED TYPE-----------------------------------------
   type orb
      integer :: nr            ! radial quantum number
      integer :: j             ! 2 x j
      integer :: l             ! L
      integer :: par           ! parity = +/- 1
      integer :: w             ! excitation 
   end type orb
   type (orb),allocatable :: orbqn(:,:) ! orbqn(species,iorb)
end module sporbit
!----------------------------------------------------
module phonon1b
   implicit none

   integer Jtrans, Ttrans  ! J and T of transition
   real, allocatable :: t1bme(:,:),p1bme(:,:),n1bme(:,:)
   logical :: pnformal,xpn

end module phonon1b




