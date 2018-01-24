!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK configuration-interaction shell-model code
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!---------------------------------------------------------------------------
! module for organizing the hops 
!
! "hops" are fundamental operations, essentially the action of 
! fermion creation/annihilation operators.
!
! "jumps" organize the action of N-body operators between Slater determinants
! of the same species. Jumps are constructed from hops.
!
!  In parallel, but at a lower level, hops are the action of creation/annihilation
!  operators on half-Slater determinants or "haikus".
!
!  Hops are constructed in bhopslib.f90
!
!---------------------------------------------------------------------------
      module hoppy   
!
!  by CWJ 10/2008
!
      implicit none

!--------- DERIVED TYPE FOR LINKED LISTS

      type grhoplink
         type (grhoplink), pointer :: next
         integer :: bhop
      end type grhoplink

      type grhopstart
         type (grhoplink), pointer :: cur,start
      end type grhopstart
!.............................
      type basehop
       integer :: iaddr ! initial address
       integer :: faddr ! final address
       integer :: phase ! phase
       integer :: op    ! operator
      end type basehop
      
      type oplist
       integer                      :: group    ! which group the operators are in
       integer                      :: iblock    ! initial block
       integer                      :: fblock    ! final block
       integer                      :: nhops   ! # of hops in this block hop
       type(basehop), pointer :: hop(:) ! list of hops that can take
                                                 ! you from the initial to the final block
      end type oplist

      type hoplist
      ! anniliation blockhops
       integer                   :: ndjumps      ! TOTAL number of annihilation hops
       integer                   :: ndbhops      ! number of annihilation blockhops
       type(grhopstart), pointer :: dgrhoplist(:)

       type(oplist), pointer :: dblockhop(:)! annihilation blockhops

      ! creation blockhops
       integer                   :: ncjumps      ! TOTAL number of creation hops
       integer                   :: ncbhops      ! number of creation blockhops
       type(grhopstart), pointer :: cgrhoplist(:)
       type(oplist), pointer :: cblockhop(:)! annihilation blockhops

      end type hoplist

      type(hoplist) :: bunny(-2:2)

	  contains
		  
!===========================================================================
! This set of routines is used to calculate the annihilation and
! creation operators on the haikus
!===========================================================================

!=====================================================================
!
!  hopmaker
! ... 'hops' are fundamental creation/annihilation operations between haikus
! ... n-body 'jumps' are constructed from hops
!
!  separated out in 7.4.4
!...........................................................................
!
!  SUBROUTINES CALLED
!     clocker
!     bunnymaster
!     print_out_hops
!     linkupgrouphops
!
!  called by
!     hoopmaster
!     density1b_from_oldwfn
!     applicator1b
!
!

subroutine hopmaker
   use timing
   use verbosity
   use nodeinfo
   implicit none

   integer :: ith, its	!ok8


   call clocker('hop','sta')

!$omp parallel do private(its,ith) SCHEDULE(STATIC)
   do its = -1,2
      ith = its
      if (ith < 1)ith = ith-1
      call bunnymaster(ith)
   end do 
!$omp end parallel do

   if ( print_hops .and. iproc==0) then 
      call print_out_hops(1)
      call print_out_hops(-1)
   end if

   call linkupgrouphops(1)
   call linkupgrouphops(-1)
   call linkupgrouphops(2)
   call linkupgrouphops(-2)
   call clocker('hop','end')
   return
end subroutine hopmaker


!===========================================================================
! This subroutine converts a haiku from binary word into an array of 1s 
! and 0s 
!
! INPUT:
! (1) ith  --> species/handness
! (2) hsd  --> haiku
!
! OUTPUT:
! (1) iarray --> array filled with 1s and 0s
!===========================================================================
subroutine convertbin(ith,hsd,iarray)
  use haiku_info
  use spstate
  use bitstuff

  implicit none
  integer :: ith
  integer :: hsd(nword(abs(ith)))
  integer :: n
  integer :: iarray(nhspsmax)
  
  integer :: imask
  integer :: i, j
  integer :: iword
  integer :: jmax
  integer :: k
  integer :: hsd0

  iarray(:) = 0
  k = 0
  do iword = 1, nword(abs(ith))
     if( iword == nword(abs(ith)) ) then
        jmax = nhspsmax
     else
        jmax = max_bit_word
     endif
     imask = 1
     hsd0 = hsd(iword)
     do j = 1, jmax
        k = k+1
        if ( iand(imask,hsd0) /= 0 ) iarray(k) = 1
        imask = ishft(imask,1)
     end do
  end do
  return
end subroutine convertbin

!===========================================================================
!  This subroutine returns phase when eliminated i'th particle
!
! INPUT:
! (1) ith  --> species/handness
! (2) hsd  --> a pointer of dimension 1 with the haiku
! (3) i    --> location of the particle in the haiku counted from 1... 
!
! OUTPUT:  
! (1) phase --> phase    
!
! EXAMPLE: 
! If the haiku is 0101100100011 and i = 8 the phase is 
! -1 because jump over 3 occupied bits
!===========================================================================
subroutine getphase(ith,hsd,i,phase)
  use haiku_info
  use bitstuff
  use spstate
  implicit none
  
  integer, pointer ::  hsd(:)

  integer :: phase
  integer :: ith
  integer :: i
  integer :: imask
  integer :: j
  integer :: iword
  integer :: jmax
  integer :: k
  integer :: icount

  k = 0
  phase = 1
  do iword = 1, nword(abs(ith))
     if ( iword == nword(abs(ith)) ) then
        jmax = nhspsmax
     else
        jmax = max_bit_word
     end if
     imask = 1
     do j = 1, jmax
        k = k + 1
        if ( k == i ) return
        if ( iand(imask,hsd(iword)) /= 0 ) phase =-phase
        imask = ishft(imask,1)
     end do
  end do
  print*, ' should not be here in getphase '
  print*, i
  print*, hsd
  stop
end subroutine getphase

!===========================================================================
subroutine search_haiku(ith,hsdin,nhsd,hsdlist,ilist)
  use haiku_info

  implicit none
!-----------Input-----------------------------------------------------------
  integer :: ith
  integer :: nhsd

  integer, pointer :: hsdin(:)
  integer, pointer :: hsdlist(:,:)
  integer mid,nword0
!-----------Output----------------------------------------------------------
  integer :: ilist

!----------intermediate-----------------------------------------
  integer iword
  nword0 = nword(abs(ith))

  do ilist = 1, nhsd
     iword = 1
     do while ( hsdin(iword) == hsdlist(iword,ilist) )
        if ( iword == nword0 ) return
        iword = iword +1
     end do
  end do

  print*, ' did not find haiku in list! '
  print* ,hsdin,nhsd,nword(abs(ith))
  do ilist = 1, nhsd
     print*, hsdlist(:,ilist)
  end do
  stop
end subroutine search_haiku

!-----------------------------------------------------------
!================================================================
! This subroutine returns the block attributes  for a given set of 
! quantum numbers: nh, jzh, parh, and wh
!
! INPUT:
! (1) nh   --> number of particles in a haiku
! (2) jzh  --> values of (2x)jz
! (3) parh --> parity
! (4) wh   --> W
!
! OUTPUT:
! (1) id   --> block id /This number is assigned to identify the blocks/
! (2) nhsd --> number of haikus in the block
! (3) hsd  --> a pointer of dimension 2 with the haikus in the block
!===========================================================================
subroutine hblock_find(ith,nh,jzh,parh,wh,id)
  use blocks
  implicit none

  integer :: ith
  integer :: nblocks
  integer :: nh,jzh,parh,wh
  integer :: id
  integer :: i
  integer :: nhblocks
  logical :: success


  success = .false.

  nhblocks = hblock(ith)%nhblocks
  do i = 1, nhblocks
     if ( nh   == hblock(ith)%list(i)%nh   .and. &
          jzh  == hblock(ith)%list(i)%jzh  .and. &   
          parh == hblock(ith)%list(i)%parh .and. &
          wh   == hblock(ith)%list(i)%wh ) then
        id   = hblock(ith)%list(i)%blockid
 
        success = .true.
        return
     end if
  end do
  
  if ( .not. success  ) then
     write (6,*) 'Block with this set of quantum numbers:'
     write (6,*) ith,nh, jzh, parh, wh
     write (6,*) 'does not exist!'
     stop
  end if

end subroutine hblock_find

!-----------------------------------------------------------
!subroutine bunnymaster
!
!  master routine for generating bunny hops 
!  CWJ 10/2008
!

subroutine bunnymaster(ith)  

  implicit none

  integer ith

  call shiva(ith)
  call brahma(ith)
  call sorthops(ith)

  return
end subroutine bunnymaster

!-----------------------------------------------------------------
!===========================================================================
! Calculates all destruction operators
! For each haiku with nh particles there are exactly
! nh destruction operators
!
! routines called: 
!  convertbin
!  parmult
!===========================================================================
subroutine shiva(ith)

  use system_parameters
  use verbosity
  use haiku_info
  use haikus
  use blocks
  use spstate
  use sporbit
  use bitstuff
  use nodeinfo
  use bmpi_mod
  use basis
  implicit none
!  include 'binterfaces.inc'

  integer :: ierr

  integer :: ith, it
  integer :: nh, parh, wh, iw, jzh
  integer :: nhblocks
  integer :: nhaiku
  integer :: iblock, fblock
  integer :: i, j, g,m
  integer :: n
  integer :: iop, djz, dw, dpar
  integer :: jzf
  integer :: parf,wf
  integer :: fw
  integer :: ntmp
  integer :: blockid
  integer :: ihaiku, iadd, nadd
  integer :: iword
  integer :: imask
  integer :: phase
  integer :: fadd
  integer(8) :: ndjumps
  integer :: ndbhops
  integer :: nhops
  
  integer, pointer     :: hsd(:,:)
  integer, pointer     :: htmp(:,:)
  integer, pointer     :: hsdi(:)
  integer, allocatable :: occcom(:)
  integer, allocatable :: iarray(:)
  
  logical, allocatable :: occgroup(:)
  integer, allocatable :: groupmap(:),fmap(:)
  
  integer :: extraphase 
  integer :: aerr

  it = abs(ith) ! nspecies

  nhblocks = hblock(ith)%nhblocks
  
  allocate ( iarray(nhspsmax) , stat=aerr)
  if(aerr /= 0) then
     call memerror("shiva 1")
     stop 5
  end if
  allocate ( occcom(nword(it)) , stat=aerr)
  if(aerr /= 0) then
     call memerror("shiva 2")
     stop 5
  end if
  allocate ( occgroup(ngroups(ith)), stat=aerr)
  if(aerr /= 0) then
     call memerror("shiva 3")
     stop 5
  end if
  allocate ( groupmap(nhspsmax) , fmap(nhspsmax) , stat=aerr)
  if(aerr /= 0) then
     call memerror("shiva 4")
     stop 5
  end if

  if ( iproc == 0 .and. verbose_hops) print*,ngroups(ith),' groups '

!  max_addr = 0
!  max_hops = 0
!------------------- FIRST DETERMINE THE NUMBER OF BLOCK HOPS

  ndbhops = 0
  ndjumps = 0

  do iblock = 1, nhblocks  ! loop over all blocks
     nh      = hblock(ith)%list(iblock)%nh
     if ( nh > maxNh(ith) .or. nh == 0 ) cycle
     nhaiku  = hblock(ith)%list(iblock)%nhsd
     hsd    => hblock(ith)%list(iblock)%hsd(:,:)
     occcom  = 0
!----search through all haikus and determine common occupied states --------
     do i = 1, nhaiku
        occcom(:) = ior(occcom(:),hsd(:,i))
     end do
     call convertbin(ith,occcom,iarray)

!--------- DETERMINE OCCUPIED GROUPS ---------------------------

     occgroup = .false.
     do i = 1,nhspsmax
        if(iarray(i)/=0)occgroup( hspsqn(ith,i)%group) = .true.
     enddo ! i
     
     do i = 1,ngroups(ith)
        if(occgroup(i))ndbhops = ndbhops + 1
     end do !i
  enddo ! iblock

  bunny(ith)%ndbhops = ndbhops

  allocate( bunny(ith)%dblockhop(ndbhops), stat=aerr)
  if(aerr /= 0) call memerror("shiva 10")


  if( iproc == 0  .and. verbose_hops) print*,ndbhops,' block hops '
!---------- LOOP OVER BLOCK HOPS AND CREATE HOPS ------------

  n = 0
  do iblock = 1, nhblocks  ! loop over all blocks

     nh      = hblock(ith)%list(iblock)%nh
     if ( nh > maxNh(ith) .or. nh == 0 ) cycle

     nhaiku  = hblock(ith)%list(iblock)%nhsd
     hsd    => hblock(ith)%list(iblock)%hsd(:,:)
     occcom  = 0
!----search through all haikus and determine common occupied states --------
     do i = 1, nhaiku
        occcom(:) = ior(occcom(:),hsd(:,i))
     end do
     call convertbin(ith,occcom,iarray)

!--------- DETERMINE OCCUPIED GROUPS ---------------------------

     occgroup = .false.
     do i = 1,nhspsmax
        if(iarray(i)/=0)occgroup( hspsqn(ith,i)%group) = .true.
     end do ! i
     
     groupmap = 0
     do g = 1,ngroups(ith)
        if(occgroup(g))then
           n = n+1
           bunny(ith)%dblockhop(n)%iblock = iblock
           bunny(ith)%dblockhop(n)%group = g
           bunny(ith)%dblockhop(n)%nhops = 0
!--------------- MAKE GROUPMAP
!              maps groups from this iblock to a particular blockhop

           groupmap(g) = n
        endif
     enddo !i 

!----LOOP OVER HAIKUS AND COUNT UP HOPS  --------
     do i = 1, nhaiku
        call convertbin(ith,hsd(:,i),iarray)
        do j = 1,nhspsmax
           if (iarray(j) == 1 ) then
              g = hspsqn(ith,j)%group
              if ( g == 0 ) then
                 if ( iproc == 0 )print*,' error in group '
                 call BMPI_ABORT(icomm,101,ierr)
                 stop
              endif
              m = groupmap(g)
              if ( m == 0 .and. iproc == 0 ) then
                 print*,' problem with finding group '
                 print*,g,n,j,iblock
                 print*,groupmap
              end if
              bunny(ith)%dblockhop(m)%nhops = bunny(ith)%dblockhop(m)%nhops+1
           end if
        end do !j 
     end do  !i 

!--------------- ALLOCATE SPACE FOR HOPS ------------------------
     do g = 1,ngroups(ith)
        m = groupmap(g)
        if(m==0)cycle  ! this group not used
        nhops = bunny(ith)%dblockhop(m)%nhops
        allocate( bunny(ith)%dblockhop(m)%hop(nhops), stat=aerr)
        if(aerr /= 0) call memerror("shiva 20")
        bunny(ith)%dblockhop(m)%nhops = 0   ! this will count back up again later
     enddo
!--------------- MAKE LIST OF FINAL BLOCKS AS FUNCTION OF GROUP -----------
!             GET INITIAL QUANTUM NUMBERS
     wh   = hblock(ith)%list(iblock)%wh
     parh = hblock(ith)%list(iblock)%parh
     iw   = hblock(ith)%list(iblock)%iwh
     jzh  = hblock(ith)%list(iblock)%jzh
     
!--------- LOOP OVER FINAL GROUPS -------------------
     fmap  = 0

     do g = 1,ngroups(ith)
        m = groupmap(g)
        if(m==0)cycle  ! this group not used
        djz  = -group(ith)%jz(g)
        dw  = -group(ith)%w(g)
        dpar = group(ith)%par(g)
        dpar = (3-dpar)/2
        jzf = jzh + djz
!--------------------- HISTORIC ERROR TRAP ---------------------

        if ( (ith > 0 .and. jzf < 0) .or. (ith < 0 .and. jzf > 0) ) then
           if ( iproc == 0 ) then
              write(6,*) ' jzf problem ',ith, jzh, djz, jzf
              write(6,*)iop,hspsqn(ith,iop)%m
           end if
           call BMPI_ABORT(icomm,101,ierr)
           stop
        end if

        parf = parmult(parh,dpar)
        wf   = wh +dw
        if(wf < 0)cycle
!

        if ( nh == 1 .and. parf == 2 ) then
           if ( iproc==0 ) print*,' problem '
           call BMPI_ABORT(icomm,101,ierr)
           stop
        end if
        call hblock_find(ith,nh-1,jzf,parf,wf,fblock)
        fmap(g) = fblock

     end do ! g

!---------------- LOOP OVER HAIKUS AND STORE HOPS ----------------------------------
     do iop = 1,nhsps(ith)
        g = hspsqn(ith,iop)%group
        if(g==0)cycle
        m = groupmap(g)
        if(m==0)cycle
!-------------- FIND FINAL BLOCK
        fblock = fmap(g)
        if ( fblock == 0 ) then
           if ( iproc == 0 )print*,' Error finding fblock ',g
           call BMPI_ABORT(icomm,101,ierr)
           stop
        endif
        bunny(ith)%dblockhop(m)%fblock = fblock
        ntmp = hblock(ith)%list(fblock)%nhsd
        htmp => hblock(ith)%list(fblock)%hsd(:,:)
        iword = (iop-1)/max_bit_word+1
        j = mod(iop-1,max_bit_word)+1 
        imask = ishft(1,j-1)

!--------------------------------------------------------------------
!     CRITICAL EXTRA PHASE
!     choose SD looks like ( left haiku creation ) ( right haiku creation ) | 0 >
!     therefore, right haikus pick up an EXTRA PHASE from operators commuting
!     past the left haiku
!        extraphase = (-1)**(np(it) - nh )

        do iadd = 1, nhaiku
           if ( iand(imask,hsd(iword,iadd)) /=0 ) then
              bunny(ith)%dblockhop(m)%nhops = bunny(ith)%dblockhop(m)%nhops+1
              nhops = bunny(ith)%dblockhop(m)%nhops
!---------------
              hsdi => hsd(:,iadd)
              call getphase(ith,hsdi,iop,phase)
!--------------------------------------------------------------------
!     CRITICAL EXTRA PHASE
!     choose SD looks like ( left haiku creation ) ( right haiku creation ) | 0 >
!     therefore, right haikus pick up an EXTRA PHASE from operators commuting
!     past the left haiku
!              if(ith > 0)then
!                phase = phase*extraphase
!              endif
!----Remove that particle---------------------------------------------------
              hsd(iword,iadd) = hsd(iword,iadd)-imask
!----Then find final address------------------------------------------------   
! need to fetch itmp or rewrite the following call
              call search_haiku(ith,hsdi,ntmp,htmp,fadd)

!              max_addr = bmax(max_addr,fadd)
              bunny(ith)%dblockhop(m)%hop(nhops)%iaddr = iadd
              bunny(ith)%dblockhop(m)%hop(nhops)%faddr = fadd
              bunny(ith)%dblockhop(m)%hop(nhops)%phase = phase
              bunny(ith)%dblockhop(m)%hop(nhops)%op = iop

              ndjumps  = ndjumps + 1
!----Replace particle-------------------------------------------------------
              hsd(iword,iadd) = hsd(iword,iadd)+imask


           endif
        enddo !iadd 
     end do  !iop 
     
  enddo ! iblock

  if( verbose_hops .and. iproc == 0 ) then
     print*,' '
     print*,' there are ',ndjumps,' destruction hops '
  end if

  deallocate ( occcom, iarray, occgroup, groupmap, fmap ) ! Free memory
  return
end subroutine shiva
!===========================================================================
!subroutine brahma

subroutine brahma(ith)

  use verbosity
  use haiku_info
  use haikus
  use blocks
  use spstate
  use nodeinfo
  use bmpi_mod
  implicit none

  integer(4) :: ierr

  integer(4) :: ith, it
  integer(4) :: nhblocks
  integer(4) :: ib,ih
  integer(4) :: nbhops
  integer(4) :: nhops
  integer :: aerr

  it = abs(ith) ! nspecies
  nhblocks = hblock(ith)%nhblocks
  nbhops = bunny(ith)%ndbhops

  bunny(ith)%ncbhops = nbhops
  allocate(bunny(ith)%cblockhop(nbhops), stat=aerr) 
  if(aerr /= 0) call memerror("brahma 1")

  do ib = 1,nbhops
     bunny(ith)%cblockhop(ib)%group = bunny(ith)%dblockhop(ib)%group
     bunny(ith)%cblockhop(ib)%iblock = bunny(ith)%dblockhop(ib)%fblock
     bunny(ith)%cblockhop(ib)%fblock = bunny(ith)%dblockhop(ib)%iblock
     nhops = bunny(ith)%dblockhop(ib)%nhops

     bunny(ith)%cblockhop(ib)%nhops = nhops
     allocate(bunny(ith)%cblockhop(ib)%hop(nhops), stat=aerr)
     if(aerr /= 0) call memerror("brahma 2")
     do ih = 1,nhops
       bunny(ith)%cblockhop(ib)%hop(ih)%iaddr = bunny(ith)%dblockhop(ib)%hop(ih)%faddr
       bunny(ith)%cblockhop(ib)%hop(ih)%faddr = bunny(ith)%dblockhop(ib)%hop(ih)%iaddr
       bunny(ith)%cblockhop(ib)%hop(ih)%phase = bunny(ith)%dblockhop(ib)%hop(ih)%phase
       bunny(ith)%cblockhop(ib)%hop(ih)%op = bunny(ith)%dblockhop(ib)%hop(ih)%op

     enddo ! ih
  enddo  ! ib

  return
end subroutine brahma
!===========================================================================
! subroutine linkupgrouphops
!
!  creates linked lists of blockhops that have the same group operator
!  this reduces searching later on -- see routine nextblockgen in bjumplib1
!
!
subroutine linkupgrouphops(ith)

use spstate
implicit none

integer ith
integer ihop
type (grhoplink), pointer :: current
integer g

logical verbosen
integer :: aerr

verbosen = .false.

allocate(bunny(ith)%cgrhoplist( ngroups(ith)), stat=aerr)
if(aerr /= 0) call memerror("linkupgrouphops 1")
allocate(bunny(ith)%dgrhoplist( ngroups(ith)), stat=aerr)
if(aerr /= 0) call memerror("linkupgrouphops 2")

do g = 1,ngroups(ith)
  allocate( bunny(ith)%cgrhoplist(g)%start, stat=aerr)
  if(aerr /= 0) call memerror("linkupgrouphops 10")
  nullify( bunny(ith)%cgrhoplist(g)%start%next)
  bunny(ith)%cgrhoplist(g)%start%bhop = 0  ! initialize
  bunny(ith)%cgrhoplist(g)%cur => bunny(ith)%cgrhoplist(g)%start

  allocate( bunny(ith)%dgrhoplist(g)%start, stat=aerr)
  if(aerr /= 0) call memerror("linkupgrouphops 11")
  nullify( bunny(ith)%dgrhoplist(g)%start%next)
  bunny(ith)%dgrhoplist(g)%start%bhop = 0  ! initialize
  bunny(ith)%dgrhoplist(g)%cur => bunny(ith)%dgrhoplist(g)%start
enddo  !g
!----------------------------------- CREATION BLOCK HOPS-----------------------

do ihop = 1,bunny(ith)%ncbhops
  g = bunny(ith)%cblockhop(ihop)%group
  bunny(ith)%cgrhoplist(g)%cur%bhop = ihop
  allocate(bunny(ith)%cgrhoplist(g)%cur%next, stat=aerr)
  if(aerr /= 0) call memerror("linkupgrouphops 20")
  nullify( bunny(ith)%cgrhoplist(g)%cur%next%next)
  bunny(ith)%cgrhoplist(g)%cur => bunny(ith)%cgrhoplist(g)%cur%next
  bunny(ith)%cgrhoplist(g)%cur%bhop = 0
enddo  ! ihop


!----------------------------------- DESTRUCTION BLOCK HOPS-----------------------

do ihop = 1,bunny(ith)%ndbhops
  g = bunny(ith)%dblockhop(ihop)%group
  bunny(ith)%dgrhoplist(g)%cur%bhop = ihop
  allocate(bunny(ith)%dgrhoplist(g)%cur%next, stat=aerr)
  if(aerr /= 0) call memerror("linkupgrouphops 30")
  nullify( bunny(ith)%dgrhoplist(g)%cur%next%next)
  bunny(ith)%dgrhoplist(g)%cur => bunny(ith)%dgrhoplist(g)%cur%next
  bunny(ith)%dgrhoplist(g)%cur%bhop = 0
enddo  ! ihop

do g = 1,ngroups(ith)
  bunny(ith)%cgrhoplist(g)%cur => bunny(ith)%cgrhoplist(g)%start
  bunny(ith)%dgrhoplist(g)%cur => bunny(ith)%dgrhoplist(g)%start
enddo
!------------------ SAMPLE USE ------------------------

if(verbosen)then
print*,' results '
do g = 1,ngroups(ith)
   current => bunny(ith)%cgrhoplist(g)%start
   ihop = 1
   do while(ihop == 1)
   print*,ith,g,current%bhop
   if(associated(current%next))then
      current => current%next
   else
      ihop = 2
   endif
   enddo
enddo  ! g
endif
!-----------------------------------------------------
return
end subroutine linkupgrouphops
!===========================================================================

 subroutine print_out_hops(ith)

  use haiku_info
  use haikus
  use blocks
  use spstate
  use verbosity
  implicit none

  integer ith, it
  integer nhblocks
  integer ib,ih
  integer nbhops
  integer nhops,ihop

  if (.not.print_hops)return

  it = abs(ith) ! nspecies
  nbhops = bunny(ith)%ndbhops
  
  write(hop_file,*)' '
  write(hop_file,*)' HOPS BY BLOCK ',ith
  do ib = 1,nbhops
     write(hop_file,901)ib, bunny(ith)%dblockhop(ib)%iblock,bunny(ith)%dblockhop(ib)%fblock, &
                             bunny(ith)%dblockhop(ib)%group
     nhops = bunny(ith)%dblockhop(ib)%nhops
901  format(i3,' block : ',i3, '->',i3,' by group ',i3)
     do ihop = 1,nhops
        write(hop_file,1001)bunny(ith)%dblockhop(ib)%hop(ihop)%iaddr, & 
                             bunny(ith)%dblockhop(ib)%hop(ihop)%faddr,& 
                             bunny(ith)%dblockhop(ib)%hop(ihop)%op,   &
                             bunny(ith)%dblockhop(ib)%hop(ihop)%phase
     enddo
1001 format('address : ',i4,' -> ',i4,' by op ',i4,' phase = ',i2)

  enddo  ! ib

 return
 end subroutine print_out_hops

!End of BHOPSLIB.
!===========================================================================
!
! routines called:
!   hippity
!
subroutine sorthops(ith)

implicit none
!include 'binterfaces.inc'
integer ith
type(basehop), pointer :: hop(:)
integer ibh,j
integer nhops

do ibh = 1, bunny(ith)%ndbhops
   hop => bunny(ith)%dblockhop(ibh)%hop
   nhops =bunny(ith)%dblockhop(ibh)%nhops
   call hippity(hop,nhops)

end do  ! ibh
do ibh = 1, bunny(ith)%ncbhops
   hop => bunny(ith)%cblockhop(ibh)%hop
   nhops =bunny(ith)%cblockhop(ibh)%nhops
   call hippity(hop,nhops)

end do  ! ibh
return
end subroutine sorthops
!=======================================================================
subroutine hippity(hop,nhops)

 implicit none
 type (basehop), pointer :: hop(:)
 type (basehop) :: tmphop
 integer nhops
  integer i,ir,j,l
  integer iadd

  l = nhops/2+1
  ir = nhops

  do while(ir > 0)
    if(l > 1)then
        l = l-1
        tmphop = hop(l)
    else
        iadd = hop(ir)%iaddr
        tmphop = hop(ir)
        hop(ir) = hop(1)
        ir = ir -1
        if(ir == 1)then
            hop(1) = tmphop
            return
        end if
    end if    
    i = l
    j = l+l
    do while(j <= ir)
         if(j < ir) then
             if(hop(j)%iaddr < hop(j+1)%iaddr) j = j+1
         end if
         if( tmphop%iaddr < hop(j)%iaddr) then
            hop(i) = hop(j)
            i = j
            j = j + j
         else
            j = ir +1
         end if

    end do
    hop(i) = tmphop
  end do

 
  return
end subroutine hippity
!=======================================================================

end module hoppy

