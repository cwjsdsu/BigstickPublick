!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  bjumplib2
!
!  subroutines:
!  MATCHMAKER: combines two streams of "descendents" to make "geneologies"
!  CHAINGANG: 
!  GETGENPHASE:
!  FORGECHAINS:
!  ADD_A_LINK
!  DIAGLINK
!
!====================================================================

!====================================================================
! subroutine matchmaker
!
!  Combines two streams of descendents to make a geneology
!
! INPUT
!  it = which species
!  igen = # of generations in initial descent
!  is = initial sector
!  isdescent = descents from initial sector
!  fgen = # of generations in final descent
!  fs = final sector
!  fsdescent = descents from final sector
!  ordered  = logical flag, only used for diagonal geneology (same initial, starting sector)
!             makes sure the genologies aren't repeated
! OUTPUT
!
!  geneology    : linked list of geneologies starting from each initial haiku block pairs 
!  ngeneology    # of geneologies
!
!  NOTES: this routine is slow for very large problems. To speed things up, find out 
!  a priori which sectors can even talk to one another. 
!

subroutine matchmaker(it,igen,is,isdescent,fgen,fs,fsdescent,ordered,ngeneologies)

use sectors
use descendents
use geneologies
implicit none

integer :: it
integer igen,fgen  ! # of generations in initial, final descents
integer :: is,fs
type (gjsector) :: isdescent,fsdescent
logical          :: ordered
integer(8) :: ngeneologies
!.......................................................

type (genlink), pointer ::  curgen
integer :: nblocki,nblockf
integer :: iblocki,iblockf
integer(8) :: idescend,fdescend
integer :: lblock,rblock
integer :: fstop
integer(8) :: descentstop
 integer :: i
 integer :: aerr

 nblocki = xsd(it)%sector(is)%nhblocks
 nblockf = xsd(it)%sector(fs)%nhblocks

 ngeneologies = 0
 curgen => geneology
 do iblocki = 1,nblocki
    if(ordered .and. is == fs)then
      fstop = iblocki
    else
      fstop = nblockf
    endif
    do idescend = 1,isdescent%block(iblocki)%ndescendents
      lblock = isdescent%block(iblocki)%flblock(idescend)
      rblock = isdescent%block(iblocki)%frblock(idescend)
      
      do iblockf = 1,fstop
         descentstop = fsdescent%block(iblockf)%ndescendents
         do fdescend = 1,descentstop
             if(lblock /= fsdescent%block(iblockf)%flblock(fdescend))cycle
             if(rblock /= fsdescent%block(iblockf)%frblock(fdescend))cycle
             ngeneologies = ngeneologies +1 

             curgen%iblock  = iblocki
             curgen%fblock  = iblockf

             curgen%ilblock = isdescent%block(iblocki)%ilblock
             curgen%irblock = isdescent%block(iblocki)%irblock
             curgen%flblock = fsdescent%block(iblockf)%ilblock
             curgen%frblock = fsdescent%block(iblockf)%irblock

             allocate(curgen%gen(igen+fgen), stat=aerr)
             if(aerr /= 0) call memerror("matchmaker 1")
             do i = 1,igen
               curgen%gen(i) = isdescent%block(iblocki)%descent(i,idescend)
             enddo  !i
             do i = 1,fgen,1
               curgen%gen(i+igen) = fsdescent%block(iblockf)%descent(fgen+1-i,fdescend)
             enddo  !i
!------------------------ SET UP FOR NEXT
             if(.not.associated(curgen%next))then
             allocate(curgen%next, stat=aerr)
             if(aerr /= 0) call memerror("matchmaker 2")
             nullify(curgen%next%next)
             end if
             curgen => curgen%next
         enddo ! fdescend
      enddo  ! iblockf
    enddo ! idescend
 enddo ! iblocki
 curgen => geneology
 return

end subroutine matchmaker


!==============================================================
module chains

implicit none

integer :: maxgen = 6
type chlink
   integer :: iadd
   integer :: fadd
   integer :: phase
   integer :: oplist(6)  ! must = maxgen
   logical :: sterile ! 
end type chlink

type chainbase
   type (chainbase), pointer :: next
   type (chlink):: chain
   logical      :: diag
!-------------- NOTE: MIGHT BE ABLE TO GET RID OF THE FOLLOWING 
! by referring to iblock, fblock up in handchain; but will take some revision
   integer      :: iblock,fblock
end type chainbase

type handchain
   type (handchain), pointer :: next
   type (chainbase) :: rchain,lchain
   integer         :: iblock,fblock
   integer(kind=1) :: genphase
end type handchain

type (chainbase), target :: Rchain,Lchain

type (handchain), target :: drchain

integer :: currentID
integer :: numberchains

contains
	
!=====================================================================
!
! fundamental routine to create "chains," an intermediate step in creating jumps.
! Chains are constructed from geneologies, then forged together to make jumps
!
! Because it is difficult to count ahead of time the # of links in a chain, 
! and because we want to reuse them (in case of memory leaks),
! we used linked lists to construct them. 
!
! INPUT:
!  it = species
!  ngen = total # of generations (1,2, or 3 for 1,2,3-body ops)
!  ngeneologoes: # of geneologies
!  optype = 'd' or 'c' for destruction/creation operators
!
! MODIFIED:
!  geneology = pointer to a geneology; links are added to it
!
!  SUBROUTINES CALLED:
!    forgechains
!    getgenphase
!

  subroutine chaingang(it,ngen,curgen,optype)

  use verbosity
  use sectors
  use geneologies
!  use chains
  implicit none

  integer it
!  integer is,fs
  integer ngen
  integer(8) :: ngeneologies
  type (genlink), pointer :: curgen
  character, pointer :: optype(:)*1
!  type (handchain), allocatable :: mrchain(:)
  type (handchain),pointer :: mrchain
  integer handsign
  integer ig
  integer(kind =1) :: genphase

  mrchain => drchain

  mrchain%iblock = curgen%iblock
  mrchain%fblock = curgen%fblock

!------------------ COMPUTE "GENPHASE"-----------
! this is the phase the right haikus pick up after operators 
! commute past the left haiku
  call getgenphase(it,curgen%ilblock,ngen,curgen%gen,genphase)
  mrchain%genphase = genphase
!------------------ DO LEFT--------------
  handsign = -1
  call forgechains(it,ngen,curgen,optype,handsign,mrchain%lchain)
!                     RIGHT
  handsign = +1
  call forgechains(it,ngen,curgen,optype,handsign,mrchain%rchain)

  curgen => curgen%next
  return
  end subroutine chaingang
!=====================================================================
!
! subroutine getgenphase
!
! a subtle but important phase
! when one jumps over a left-haiku, one can pick up a phase if 
! there an odd number of particles in that left haiku. 
! This is compounded by the fact that particles are created and 
! destroyed in a geneology.
! Computes the phase from anticommuting past a left haiku
! even as particles are dynamically removed from that haiku
!
! input:
!   it: species
!   lblock: which (left) haiku block
!   ngen = total # of generations
!
! output:
!   genphase: phase picked up by anticommuting past particles in 
!           left haiku
!
! called by: CHAINGANG

  subroutine getgenphase(it,lblock,ngen,gen,genphase)

  use blocks
  implicit none
  integer it
  integer lblock
  integer ngen
  integer gen(ngen)
  integer(kind =1),intent(out) :: genphase

  integer nhl
  integer(kind = 1) :: lphase
  integer igen
  
  genphase = 1
  nhl = hblock(-it)%list(lblock)%nh  ! # of particles in left haiku
  lphase = int((-1)**nhl,1)
  do igen = 1,ngen
    if(gen(igen) > 0)then  ! NB this was erroneously reversed in original version; fixed August 2009 CWJ
       genphase = genphase*lphase
    else
       lphase = -lphase
    endif
  enddo
  return
  end subroutine getgenphase
!====================================================================
!subroutine forgechains
!
!  Runs down geneologies and forges chains to create jumps
!
!  Chains are created recursively
!
! INPUT:
!  it: species
!  ngen: # of generations (typically 1,2,or 3 for 1,2 or 3-body ops)
!  ageneology: working geneology
!  optype = 'd' or 'c' for creation or destruction operator 
!  handsign = + or - for right or left haikus
! MODIFIED: 
!  hchain: chain to which a link is added
!
! subroutines called:
!	ADD_A_LINK
!	DIAGLINK
!
! called by: CHAINGANG
!
subroutine forgechains(it,ngen,ageneology,optype,handsign,hchain)

use hoppy
use blocks
!use chains
use geneologies

implicit none

integer :: it  ! species
integer :: ngen
type (genlink) :: ageneology
character, pointer :: optype(:)*1
integer handsign
type (chainbase) :: hchain

!.................................
integer igen
integer :: nhgen,ihgen
  integer bhop
logical lookback
integer i,looksign
integer group1,group2

!-------------- compute # of handed generations
    if(handsign < 0)then
      hchain%iblock = ageneology%ilblock

      hchain%fblock = ageneology%flblock

    else
      hchain%iblock = ageneology%irblock
      hchain%fblock = ageneology%frblock

    endif

  nhgen = 0
  do igen = 1,ngen
     if(ageneology%gen(igen)*handsign >0) nhgen = nhgen +1
  enddo ! igen
!-------------------- IF NO OPERATORS THEN ALL LINKS ARE 'DIAGONAL '
  if(nhgen ==0)then
     call diaglink(it*handsign,ageneology%ilblock, ageneology%irblock, hchain)
     return
  endif
!------------------ RECURSIVELY SET UP LINKS -----------------
  ihgen = 0
  hchain%diag = .false.
  do igen = 1,ngen
     if(ageneology%gen(igen)*handsign < 0) cycle
     if(optype(igen) == 'c')then
        looksign = 1 
     else
        looksign = -1 
     endif
     ihgen = ihgen +1 
     bhop = abs(ageneology%gen(igen))
!----------- SET UP LOOKBACK TO GET ORDERING CORRECT ----------- 
     lookback = .false.
     if(igen > 1)then
       if(ageneology%gen(igen)*ageneology%gen(igen-1) > 0 )then
       group1 = bunny(it*handsign)%dblockhop( abs(ageneology%gen(igen)))%group
       group2 = bunny(it*handsign)%dblockhop( abs(ageneology%gen(igen-1)))%group

       if(optype(igen) == optype(igen-1) .and. group1 == group2 )then
           lookback = .true.
       endif
     endif
     endif
!------------------- ADD NEW LINKS -----------------
     call add_a_link(it*handsign,igen,nhgen,ihgen,bhop,optype,lookback,looksign,hchain)
  enddo !igen

  return
end subroutine forgechains

!====================================================================
!
!    subroutine add_a_link
!
! This routine recursively adds links to hchain
!
!
  subroutine add_a_link(ith,igen,nhgen,ihgen,bhop,optype,lookback,looksign,hchain)

  use verbosity
!  use chains
  use geneologies
  use hoppy

  implicit none
  integer :: ith
  integer :: igen
  integer :: nhgen
  integer :: ihgen
  integer :: bhop
  integer :: looksign
  logical :: lookback
  character, pointer :: optype(:)*1
  type (chainbase),target :: hchain

  type (chainbase), pointer :: curch,nextch
  integer :: nhops
  integer :: iadd
  integer :: fadd
  integer :: phase
  integer :: ops(maxgen)
  type(basehop), pointer :: hop(:)
  integer :: ihop
  integer :: nlink, ilink
  logical :: first
  integer :: nsterile
  integer :: oldop,i,j
  logical ::  test
  integer :: itest
  integer :: aerr

  curch =>hchain

!-------------- COUNT UP # OF LINKS SO FAR ----------
!---------------AND SET nextch TO END OF LINKS SO FAR
  if(ihgen == 1)then
    nextch => hchain
    nlink = 0
  else
!------------ loop through to find the end
    nextch => hchain
    nlink = 0
    do while(nextch%chain%iadd /= 0)
      nlink = nlink +1
      nextch => nextch%next
    enddo

  endif

  select case (optype(igen))

     case ('d', 'D')
        nhops = bunny(ith)%dblockhop(bhop)%nhops
        hop => bunny(ith)%dblockhop(bhop)%hop
     case ('c','C')
        nhops = bunny(ith)%cblockhop(bhop)%nhops
        hop => bunny(ith)%cblockhop(bhop)%hop
     case default
        print*,' ooh  noo ',optype,igen
        stop
  end select

!----------------------------------- INITIALIZE IF FIRST LINK
  if(nlink == 0)then
!--------------- ERROR TRAP ----------

     if(ihgen /= 1)then  
        print*,' some problem, I expect ihgen = 1 ',ihgen
        print*,igen,ihgen,nhgen
        stop
     endif

!--------------- END ERROR TRAP ------------
     do ihop = 1,nhops

        nextch%chain%sterile = .false.
        nextch%chain%iadd = hop(ihop)%iaddr
        nextch%chain%fadd = hop(ihop)%faddr
        nextch%chain%phase= hop(ihop)%phase

        nextch%chain%oplist = 0
        nextch%chain%oplist(igen) = hop(ihop)%op

        if(.not.associated(nextch%next))then
          allocate(nextch%next, stat=aerr)
          if(aerr /= 0) call memerror("add_a_link 1")
          nullify(nextch%next%next)
        endif


        nextch=>nextch%next
        nextch%chain%iadd = 0

     enddo

     return
  endif

!------------------------------ AFTER FIRST LINKS ------------------
  do ilink = 1,nlink

     if(curch%chain%sterile)then
        curch=>curch%next
                   cycle
     endif

! -------- SET UP TEMPORARY LINK -------------------
!     if(curch%chain%sterile)cycle
     ops    = curch%chain%oplist
     iadd   = curch%chain%iadd
     fadd   = curch%chain%fadd
     phase  = curch%chain%phase

!---------------- CHANGE THIS DEPENDING ON OPTYPE
     oldop = 0

     if(lookback)then
      oldop = ops(igen -1)

     endif

     first = .true.
!--------- NOW FIND ALL THE SUBSEQENT HOPS ---------
!          SEARCH THROUGH AND FIND INITIAL ADDRESSES that equal = FINAL ADDRESS
     do ihop = 1,nhops
!        if( fadd < hop(ihop)%iaddr)exit

        if(fadd == hop(ihop)%iaddr)then

!---------------- CHECK HERE THAT WE HAVE PROPER ORDERING FOR OPERATORS FROM SAME GROUP
            if(lookback .and. looksign*oldop > looksign*hop(ihop)%op)cycle

!------------ WE HAVE A WINNER; ADD ANOTHER LINK
           if(first)then  ! replace the current link

              curch%chain%sterile = .false.

              curch%chain%fadd = hop(ihop)%faddr
              curch%chain%phase= phase*hop(ihop)%phase
              curch%chain%oplist(igen) = hop(ihop)%op
              first = .false.
           else            ! jump down and add a new link

              nextch%chain%sterile = .false.

              nextch%chain%iadd = iadd
              nextch%chain%fadd = hop(ihop)%faddr
              nextch%chain%phase= phase*hop(ihop)%phase

              nextch%chain%oplist = ops
              nextch%chain%oplist(igen) = hop(ihop)%op
              if(.not.associated(nextch%next))then
                allocate(nextch%next, stat=aerr)
                if(aerr /= 0) call memerror("add_a_link 2")
                nullify(nextch%next%next)
              endif
              nextch=>nextch%next
              nextch%chain%iadd = 0
          endif

        endif
     enddo ! ihop
    if(first)then  ! never found a hop; make this link sterile 

     curch%chain%sterile = .true.

    endif
     curch => curch%next


  enddo     ! ilink
!--------------- CHECK
!    nextch => hchain
!    nlink = 0
!    nsterile = 0
!    do while(nextch%chain%iadd /= 0)
      
!      if(nextch%chain%sterile)then
!        nsterile = nsterile + 1
!      else
!        nlink = nlink +1
!      endif
!      nextch => nextch%next
!    enddo
!    print*,nlink,' LINKS <-- ,',nsterile,' sterile '
  return 
  end subroutine add_a_link
!=====================================================================
! subroutine diaglink
!
!  INPUT:
!	ith: species/handness
!	lblock,rblock: label of left and right haiku blocks
!		only use lblock OR rblock depending if ith < or > 0
!	hchain: chain to be added to
!
! called by: FORGECHAINS
!
  subroutine diaglink(ith,lblock,rblock,hchain)

  use blocks
!  use chains
  implicit none

  integer :: ith
  integer :: lblock,rblock
  type (chainbase),target :: hchain
  type (chainbase), pointer :: curch

  integer :: nadd
  integer :: iadd
  integer :: block
  integer :: aerr

  if(ith < 0)then
     block = lblock
  else
     block = rblock
  endif

  nadd = hblock(ith)%list(block)%nhsd

  curch => hchain
  curch%chain%iadd = 0
  do iadd = 1,nadd
    curch%chain%iadd = iadd
    curch%chain%fadd = iadd
    curch%chain%sterile = .false.
    curch%chain%phase = 1
    curch%chain%oplist = 0
    if(.not.associated(curch%next))then
      allocate(curch%next, stat=aerr)
      if(aerr /= 0) call memerror("diaglink")
      nullify(curch%next%next)
    endif
    curch => curch%next
    curch%chain%iadd = 0

  enddo
  return
  end subroutine diaglink
!=====================================================================
end module chains

