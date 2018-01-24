!===============================================================
!
!    file B3MODULES.F90
!
!    modules used for 3-body interactions in BIGSTICK
!
!    started Sept 2009 by CWJ @ SDSU
!
!================================================================ 
!
!  flags for controlling 3-body forces
!
  module flags3body
  implicit none

  logical :: threebodycheck = .false.  ! flag whether or not even to ask 
                                                  ! about 3-body
  logical threebody

  end module flags3body

!================================================================ 
!
!  module of arrays for 3-body jumps
!
  module jump3body
  use jumpdef
  use precisions
  implicit none
!---- MADE EXPLICIT IN 7.7.4 -- memory storage per jump
  integer(4) :: bytesper3Bjump = 25    ! 2 x 8 isd/fsd + 8 m.e. label + 1 phase

  type (jumpsect),target :: x3bjump(2)
  integer(kind=8) :: totn3bjumps, totp3bjumps
  integer(4), allocatable,target:: n2b_cop(:), n2b_dop(:)
  integer(4), allocatable,target :: p2b_cop(:), p2b_dop(:)

  integer(kind=basis_prec),allocatable, target :: n3b_isd(:),n3b_fsd(:),p3b_isd(:),p3b_fsd(:)
  integer(kind=4),allocatable, target :: n3b_op(:),p3b_op(:)
  real(4), allocatable, target :: n3b_me(:), p3b_me(:)

  integer(kind=1), allocatable, target :: n3b_phase(:),p3b_phase(:)

  type (jumpsect), target :: PPNjump(2), PNNjump(2)

  end module jump3body
!==============================================================

  module tripdef

  implicit none
!  information on 3-body triplets | ia, ib, ic: JT >
!
  type triplet_qn
    integer :: M
    integer :: par
    integer :: W
    integer :: indx
    integer :: ia,ib,ic
  end type triplet_qn

!      type pairinfo
!        type (pair_qn),pointer :: pair(:)
!        integer, pointer       :: meref(:,:), mestart(:,:)
!      end type pairinfo

   type tripinfo
        integer :: jzstart,jzend
        type (triplet_qn),pointer :: trip(:)
        integer, pointer       :: meref(:,:), mestart(:,:)
   end type tripinfo

  end module tripdef

!====================================================================
!
!  module of arrays for 3-body interactions
!
  module interactions3body
  use pairdef
  use tripdef
  implicit none

  integer(kind=8),target :: nmatppp,nmatppn, nmatpnn,nmatnnn

  real, allocatable, target :: hmatPPP(:),hmatNNN(:)
  real, allocatable, target :: hmatPPN(:),hmatPNN(:)

  type (tripinfo) :: XXX3(2)  
  integer :: nXXX(2)
  integer, allocatable, target :: mapPPP(:),mapNNN(:)
  
  integer, allocatable, target :: mapPPN(:,:), mapPNN(:,:)
  integer :: nXXY(2)
  type (tripinfo) :: XXY(2)

  integer, allocatable,target :: cppntriplet(:,:), dppntriplet(:,:)
  integer, allocatable,target :: cpnntriplet(:,:), dpnntriplet(:,:)
!............. IN REVERSE ORDER TO MINIMIZE CACHE CALLS
  integer, allocatable,target :: cppntripletC(:,:), dppntripletC(:,:)
  integer, allocatable,target :: cpnntripletC(:,:), dpnntripletC(:,:)

  integer(1), allocatable, target :: phaseppn(:,:), phasepnn(:,:)

  real :: scale3body


  end module interactions3body
!=====================================================================

