 module localblocks
   use precisions
   implicit none
!---------------------FOR BLOCK LANCZOS ROUTINE------------------------------
! dimblock: is used to keep track of the size of our blocks of lanczos vectors
! nvec: serves as our main index through row space of block_1, and block_2
! block_1, block_2: analog of vec1 and vec2 in standard lanczos
! We define  block_1(dimblock, v1s:v1e) and block2(dimblock, v2s:v2e)

! RMZ - 8/18

integer :: dimblock !number of vectors in the columnspace of a block
integer(4) :: maxblockiter
integer :: block_flag ! used for lanczos_p to switch over to block lanczos routine
integer :: dimlanc ! dim of lanczos matrix
character :: piv_flag
real(kind=lanc_prec), allocatable, target :: block_1(:,:), block_2(:,:) !block* is taken by fortran intrinsic
! ADDED  BY CWJ IN 7.8.8, 11/2018 to allow reduction after Hmult
! block_i(ivec,istate) = blockvec_i( ivec + (istate-1)*dimblock)
! allocate blockvec_1(1+dimblock*(v1s-1):dimblock*v1e), etc
!
real(kind=lanc_prec), allocatable, target :: blockvec_1(:),blockvec_2(:)

! VARIABLES FOR BLOCK STRENGTH FUNCTIONS -- added 7.9.4 ----

integer, allocatable :: pivot_choice(:)
real(kind=8), allocatable :: dnormblock(:)
end module localblocks

