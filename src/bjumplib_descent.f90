!
! "descendents" help to set up construction of jumps
!
! Created in subroutine count_create_descendents in bjumplib1.f90
!
!

module descendents

implicit none

integer completed1bodyflag       ! this is a flag to signal completion of 1-body

type gjblock
   integer(8) :: ndescendents
   integer, pointer :: descent(:,:)   ! stores a list of block hops
                                      !  if > 0 then right blockhop, if < 0 then left
   integer :: ilblock, irblock  ! initial blocks before the descent
   integer,pointer :: flblock(:), frblock(:)  ! final blocks after the descent
end type gjblock

type gjsector
   integer nblocks  ! for convenience
   type (gjblock), pointer :: block(:)
   integer(8) :: totdescendents   ! ADDED 7.5.9: total number of descendents in this sector
                                   ! used (possibly) to estimate scaling and to distribute work  
end type gjsector

type groupjump
   integer nsectors  ! for convenience
   type (gjsector), pointer :: sector(:)
end type groupjump

! type (groupjump) :: descent(2)

contains

!====================================================================
!
! SUBROUTINES CALLED
!  blockdescent
!
subroutine sectordescent(it,is,nbody,sdescent)

use verbosity
use sectors
implicit none

integer it  ! species
integer nbody !  # of descendents to get
type (gjsector) :: sdescent
integer is
!...................................
integer iblock
integer(8) :: totaldescendents
integer :: aerr

sdescent%nblocks = xsd(it)%sector(is)%nhblocks
sdescent%totdescendents = 0

allocate(sdescent%block(  xsd(it)%sector(is)%nhblocks ), stat=aerr)
if(aerr /= 0) call memerror("sectordescent")

totaldescendents = 0
do iblock = 1,sdescent%nblocks
    if(chatty_descent)print*,' SECTOR ',is,', BLOCK ',iblock
    sdescent%block(iblock)%ilblock = xsd(it)%sector(is)%lhblock(iblock)
    sdescent%block(iblock)%irblock = xsd(it)%sector(is)%rhblock(iblock)
    call blockdescent(it,nbody,sdescent%block(iblock),totaldescendents)
enddo ! iblock
sdescent%totdescendents = totaldescendents

return
end subroutine sectordescent

!====================================================================
!
! CALLED BY:  
!     sectordescent
!
! SUBROUTINES CALLED:
!   count_create_descendents
!
subroutine blockdescent(it,nbody,bdescent,totaldescendents)

use verbosity
implicit none

integer :: it
integer :: iblock,nbody
type (gjblock) :: bdescent
integer(8) :: totaldescendents
!.........................................

integer :: igen
logical :: create
integer(8) :: ndescendents

ndescendents = 1
bdescent%ndescendents = 1

do igen = 1,nbody
    create = .false.
    if(chatty_descent)print*,' count '
    call count_create_descendents(it,nbody,igen,create,bdescent,ndescendents)
    create = .true.
    if(chatty_descent)print*,' create ',ndescendents

    call count_create_descendents(it,nbody,igen,create,bdescent,ndescendents)
    totaldescendents = totaldescendents + ndescendents
enddo  ! igen

return
end subroutine blockdescent

!====================================================================
!
!  NB: One can speed up by factor of x2 by using linked lists. This, however,
!  makes things rather complicated.
!
!  CALLED BY: 
!     blockdescent
!
! SUBROUTINES CALLED:
!   applyblockhop
!   nextblockgen
!
subroutine count_create_descendents(it,ngen,gen,create,bdescent,ndescendents)

use verbosity
use spstate
implicit none

integer :: it
integer :: ngen   ! total # of generations, i.e., # of creation/annihilation operators
integer :: gen    ! which generation we are creating
logical :: create
type (gjblock) :: bdescent

!................................................................

integer(8) :: ndescendents
integer, allocatable :: tmp(:,:)
integer(8) :: idescend
integer ::  lblock,rblock
integer :: lblock0,rblock0
integer :: igen
logical :: success
integer :: ihop,gop,g
integer :: aerr

allocate(tmp(1,1), stat=aerr)  ! prevent compiler warnings about size being uninitialized
if(aerr /= 0) call memerror("count_create_descendents 0")

if(gen > 1)then
   if(allocated(tmp))deallocate(tmp)
   allocate(tmp(ngen,bdescent%ndescendents), stat=aerr)
   if(aerr /= 0) call memerror("count_create_descendents 1")
   tmp = bdescent%descent
endif
if(create)then
    if(associated(bdescent%descent))nullify(bdescent%descent)
    allocate(bdescent%descent(ngen,ndescendents), stat=aerr)
    if(aerr /= 0) call memerror("count_create_descendents 2")
    if(gen==ngen) then
       allocate(bdescent%flblock(ndescendents),bdescent%frblock(ndescendents), stat=aerr)
       if(aerr /= 0) call memerror("count_create_descendents 3")
    end if
endif

ndescendents = 0
do idescend = 1,bdescent%ndescendents
   lblock = bdescent%ilblock
   rblock = bdescent%irblock
   if(gen > 1 )then
       do igen = 1,gen-1
          ihop = tmp(igen,idescend) 
          if(ihop ==0)then
             stop
          endif
          call applyblockhop(it,ihop,'d',lblock,rblock,gop)
       enddo
   else
       gop = ngroups(it)
   endif

!---------------- NEXT, LOOP OVER ALLOWED GROUPS
!                 including for both right ( g > 0) and left ( g< 0) 
   do g = gop,-ngroups(-it),-1
       if(g == 0)cycle  ! can't use this
       lblock0 = lblock
       rblock0 = rblock
       call nextblockgen(it,g,'d',lblock0,rblock0,.true.,success,ihop)
       if(ihop == 0 .and. success)then
          stop
       endif
       if(success)then
           ndescendents = ndescendents + 1

           if(create)then
             if(gen > 1)then
               do igen = 1,gen-1
                 bdescent%descent(igen,ndescendents) = tmp(igen,idescend)
               enddo
             endif
             bdescent%descent(gen,ndescendents) = ihop
!------------------- SAVE FINAL BLOCKS FOR LATER
             if(gen==ngen)then
               bdescent%frblock(ndescendents) = rblock0
               bdescent%flblock(ndescendents) = lblock0
             endif
           endif  ! end if create

       endif  ! end if success

   enddo  ! g
enddo  ! idescend

if(create)bdescent%ndescendents = ndescendents
return
end subroutine count_create_descendents


!==========================================================================
!
!  subroutine nextblockgen
!
!  given initial blocks lblock and rblock, takes group operator label gop
!  and hops to the next block
!  if gop > 0, rblock hops; if gop < 0, lblock hops
!  if optype == 'c' then creation hop, if 'd' then destruction hop
!
! INPUT: 
!   it = species 
!   gop = group of operator; if < 0, then acts on lblock, if > 0 acts on right block
!   optype = 'c' creation 'd' destruction
!   lblock, rblock = left and right blocks
!   replace = flag to tell if to replace lblock, rblock 
!   
! OUTPUT: modified left (if gop < 0) or right (gop > 0) block
!   success = flag to signal a hop was found
!   ihop = labels which blockhop
!
! NB: FOR LARGE SYSTEMS THIS SUBROUTINE IS CALLED A LOT
!  to speed up: create arrays in blockhops that tell you where to look, to avoid searching
!
!  CALLED BY: 
!    count_create_descendents
!
subroutine nextblockgen(it,gop,optype,lblock,rblock,replace,success,ihop)
   use verbosity
   use hoppy
   implicit none

   integer :: it
   integer :: gop
   character :: optype*1
   integer :: lblock,rblock
   logical :: replace
   logical :: success
   integer :: ihp
   integer :: ihop

   type (grhoplink), pointer :: current

   success = .true.
   if(gop > 0)then

     select case (optype)

     case ('d','D')
   !------------- SEARCH THROUGH BLOCK HOPS TO MATCH GOP

       current => bunny(it)%dgrhoplist(gop)%start
       ihp = current%bhop
       do while(ihp /= 0)
          if(chatty_descent)print*,gop,ihp,rblock, bunny(it)%dblockhop(ihp)%iblock

          if(rblock == bunny(it)%dblockhop(ihp)%iblock)then
              if(replace)rblock = bunny(it)%dblockhop(ihp)%fblock
              ihop = ihp
               return
          endif
          current => current%next
          ihp = current%bhop
       enddo
       success = .false.
       return

     case ('c','C')
   !------------- SEARCH THROUGH BLOCK HOPS TO MATCH GOP
       do ihp = 1,bunny(it)%ncbhops
         if(bunny(it)%cblockhop(ihp)%iblock==rblock .and. bunny(it)%cblockhop(ihp)%group == gop)then
           
           if(replace)rblock = bunny(it)%cblockhop(ihp)%fblock
           ihop = ihp
           return
         endif
       enddo  !ibh

     case default
       print*,' error in optype ',optype
       stop
     end select

   else
     select case (optype)

     case ('d','D')
   !------------- SEARCH THROUGH BLOCK HOPS TO MATCH GOP

       current => bunny(-it)%dgrhoplist(-gop)%start
       ihp = current%bhop
       do while(ihp /= 0)
           if(lblock == bunny(-it)%dblockhop(ihp)%iblock)then
   !           write(39,*)' block hop ',-it,-gop,ihp,lblock
              if(replace)lblock = bunny(-it)%dblockhop(ihp)%fblock
              ihop = -ihp
               return
           endif
           current => current%next
           ihp = current%bhop
       enddo
       success = .false.
       return


     case ('c','C')

   !------------- SEARCH THROUGH BLOCK HOPS TO MATCH GOP
       do ihp = 1,bunny(-it)%ncbhops
         if(bunny(-it)%cblockhop(ihp)%iblock==lblock.and.bunny(-it)%cblockhop(ihp)%group == -gop)then
           if(replace)lblock = bunny(-it)%cblockhop(ihp)%fblock
           ihop = -ihp
           return
         endif
       enddo  !ibh

     case default
       print*,' error in optype ',optype
       stop
     end select

   endif
   success = .false.

   return
end subroutine nextblockgen
!==========================================================================
!
!  subroutine applyblockhop
!
!  given initial blocks lblock and rblock, hops to the next block
!  if gop > 0, rblock hops; if gop < 0, lblock hops
!  if optype == 'c' then creation hop, if 'd' then destruction hop
!
! INPUT: 
!   it = species 
!   ihop = which label of hop
!   optype = 'c' creation 'd' destruction
!   lblock, rblock = left and right blocks
! OUTPUT: 
!   modified left (if gop < 0) or right (gop > 0) block
!   gop = group of operator; if < 0, then acts on lblock, if > 0 acts on right block
!
!  CALLED BY:
!   count_create_descendents
!
subroutine applyblockhop(it,ihop,optype,lblock,rblock,gop)
   use hoppy
   implicit none
   integer :: it
   integer :: ihop,gop
   character :: optype*1
   integer :: lblock,rblock

   if(ihop ==0)then
      print*,' problem '
      stop
   endif
   if(ihop > 0)then

     select case (optype)

     case ('d','D')
       rblock = bunny(it)%dblockhop(ihop)%fblock
       gop    = bunny(it)%dblockhop(ihop)%group

     case ('c','C')
       rblock = bunny(it)%cblockhop(ihop)%fblock
       gop    = bunny(it)%cblockhop(ihop)%group


     case default
       print*,' error in optype ',optype
       stop
     end select

   else
     select case (optype)

     case ('d','D')
       lblock = bunny(-it)%dblockhop(-ihop)%fblock
       gop    =-bunny(-it)%dblockhop(-ihop)%group

     case ('c','C')

       lblock = bunny(-it)%cblockhop(-ihop)%fblock
       gop    =-bunny(-it)%cblockhop(-ihop)%group

     case default
       print*,' error in optype ',optype
       stop
     end select

   endif
   return
end subroutine applyblockhop

end module descendents
