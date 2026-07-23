!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK configuration-interaction shell-model code
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  bbasislib4.f
!
!  fourth set of routines to generate the basis: BLOCKS of haikus
!
!  SUBROUTINES:
!    blockmaster: master calling routine to generate blocks
!       calls:
!               count_blocks
!               fill_blocks!
!    count_blocks: counts up available blocks of haikus
!    fill_blocks:  fills up available blocks
!    (printblockhaikus: for testing purposes only, prints out haikus by blocks)
!
!====================================================================
!  GOALS: to construct and fill the blocks listed in hblock
!  which has the quantum #s and lists of haikus; this is done in fill_blocks.
!
!  Towards this goal, routine count_blocks constructs blockmap
!  which, for given quantum numbers, points towards the appropriate blocks
!
!===============================================================
!
!  NOTE: These routines assume I have sorted the haikus by parity and jz.
!  I did that because I was hoping to use pointers to the haiku arrays;
!  however I did not succeed and just reassigned the haikus.
!  In that case one does not need to sort, although the algorithm should
!  be revised. In all but the largest cases the sorting does not slow
!  the code down by much.
!
!---------------------------------------------------------------------------
!  blocks of haikus, organized by nh, jz,etc.
!  key derived variable: hblock
!---------------------------------------------------------------------------
module blocks
  implicit none

  type baseHaikuBlock
     integer :: blockid ! block id
     integer :: nh
     integer :: jzh
     integer :: parh
     integer :: iwh
     integer :: wh
     integer :: nhsd
	 logical :: used     ! added in 7.7.2, to determine if blocks of haikus go unused
     integer, pointer :: hsd(:,:)
  end type baseHaikuBlock

  type speciesblocks
     type (baseHaikuBlock), pointer :: list(:)
     integer :: nhblocks
  end type speciesblocks

  type (speciesblocks) :: hblock(-2:2)

!--------------- A DERIVED VARIABLE TO FIND BLOCKS--------------------------
!
!  Note -- a particular ordering for the quantum numbers of the blocks 
!          has been chosen here:  nh ->  par -> jz -> W
!          This follows the ordering of the SDs in the proton/neutron basis
!          but is different from that after sorting the haikus
!          (nh -> w -> par -> jz )
!   
  type block4
     integer          :: nw
     integer          :: wmin,wmax
     integer, pointer :: wlist(:)
     integer, pointer :: wmap(:)
     integer, pointer :: hblock(:) 
  end type block4

  type block3
     integer                :: jzmax,jzmin
     type (block4), pointer :: jz(:)
  end type block3

  type block2
     integer       :: parmin,parmax
     type (block3) :: par(2)  ! assume both parities; won't take up much space if wrong
  end type block2
      
  type block1
     integer                :: nhmax,nhmin
     type (block2), pointer :: nh(:)
  end type block1
  
  type (block1) :: blockmap(-2:2)
  
  contains
	  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine blockmaster
!
!  INPUT: ith = species/handness
!
!  SUBROUTINES CALLED:
!     count_blocks
!     fill_blocks
!
subroutine blockmaster(ith)

use system_parameters
use haiku_info
use verbosity
use butil_mod

implicit none
integer ith
!................................
integer nh
integer nblocks
integer nhaiku
integer sumhaiku
integer maxhaiku,maxh
integer :: iblock
integer :: aerr

blockmap(ith)%nhmin = Np(abs(ith))-maxNh(-ith)
blockmap(ith)%nhmax = maxNh(ith)

allocate(blockmap(ith)%nh(0:maxNh(ith)), stat=aerr)
if(aerr /= 0) call memerror("blockmaster 1")
nblocks = 0
do nh = 0,maxNh(ith)
   call count_blocks(ith,nh,nblocks)
enddo  ! nh      
hblock(ith)%nhblocks = nblocks
allocate(hblock(ith)%list(nblocks), stat=aerr)
if(aerr /= 0) call memerror("blockmaster 2")

if(verbose_blocks)then
   print*,' There are ',nblocks,' blocks '
endif
nblocks = 0
sumhaiku = 0
maxhaiku = 0
do nh = 0,maxNh(ith)
   call fill_blocks(ith,nh,nblocks,maxh,nhaiku)
   maxhaiku = bmax(maxhaiku,maxh)
   sumhaiku = sumhaiku + nhaiku
enddo  ! nh      

if(verbose_blocks)then
   print*,' Max # of haiku in a block is ',maxhaiku
   print*,' Check: # of haikus = ',sumhaiku
endif

!------Assigns ids to the hblocks------------
do iblock = 1, nblocks
   hblock(ith)%list(iblock)%blockid = iblock
end do

return
end subroutine blockmaster
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine count_blocks
!
!   count up and allocate blocks of haikus
!
!  Note -- a particular ordering for the quantum numbers of the blocks 
!   has been chosen here:  nh -> par -> jz -> W
!   This follows the ordering of the SDs in the proton/neutron basis
!   but is different from that after sorting the haikus
!   (nh -> w -> par -> jz )
!
!  FURTHER NOTE: this routine does not require the haikus to be sorted
!
!  KEY RESULT: TO COUNT UP BLOCKS AND CONSTRUCT blockmap
!
!  INPUT:
!   ith = species/handness
!   nh = # of particles
!
! OUTPUT:
!   nb = # of blocks for this nh
!

subroutine count_blocks(ith,nh,nb)

use sporbit
use haikus
use butil_mod

implicit none
integer ith
integer nh
!-----------------------------
integer nb
!---------------------------
integer ih
integer jzmin,jzmax,jz
integer iw,w
integer parmin,parmax,ip
integer wmax,wmin,wmx,wmn
integer n,in
integer m
integer parh
integer :: aerr

logical, allocatable :: wallow(:)

!------------- GO THROUGH AND DETERMINE WLIMITS
wmin = 1000
wmax = 0
do iw = 1,haiku(ith)%of_nh(nh)%nw
   w = haiku(ith)%of_nh(nh)%wlist(iw)
   wmax = bmax(w,wmax)
   wmin = bmin(w,wmin)
enddo
allocate(wallow(wmin:wmax), stat=aerr)
if(aerr /= 0) then
   call memerror("count_blocks 1")
   stop 5
end if

!----------------PARITY ----------------------------
!
!  assume that for any jz, it comes in both parities;
!  this is probably a good assumption, except for large jzs
if(allsameparity .or. nh == 0)then
           parmin = 1
           parmax = 1
else
!-------------- need to loop through and determine min,max parity
   parmin = 2
   parmax = 1
   do iw = 1,haiku(ith)%of_nh(nh)%nw
      n = haiku(ith)%of_nh(nh)%of_w(iw)%nhaiku
      do in = 1,n
          parmin = bmin(parmin,haiku(ith)%of_nh(nh)%of_w(iw)%par(in))
          parmax = bmax(parmax,haiku(ith)%of_nh(nh)%of_w(iw)%par(in))
      enddo
   enddo  !iw
endif
blockmap(ith)%nh(nh)%parmin = parmin
blockmap(ith)%nh(nh)%parmax = parmax

!         write(26,*)ith,nh,jz,' par ',parmin,parmax

!--------------- JZ ---------------------------------
!:         DETERMINE MIN, MAX JZ

do ip = parmin,parmax
        jzmin = 1000
        jzmax = -1000
        do iw = 1,haiku(ith)%of_nh(nh)%nw
           n = haiku(ith)%of_nh(nh)%of_w(iw)%nhaiku
           do in = 1,n
             if(ip==haiku(ith)%of_nh(nh)%of_w(iw)%par(in))then
               jzmin = bmin(jzmin,haiku(ith)%of_nh(nh)%of_w(iw)%jz(in))
               jzmax = bmax(jzmax,haiku(ith)%of_nh(nh)%of_w(iw)%jz(in))
             endif
           enddo
        enddo  !iw
        if(jzmin <= jzmax)then
          blockmap(ith)%nh(nh)%par(ip)%jzmin = jzmin
          blockmap(ith)%nh(nh)%par(ip)%jzmax = jzmax
          allocate(blockmap(ith)%nh(nh)%par(ip)%jz(jzmin:jzmax), stat=aerr)
         if(aerr /= 0) call memerror("count_blocks 10")
        else  ! there is a problem
          print*,' There is a problem, cannot find jzmin,max ',ip
          print*,jzmin,jzmax,nh
          print*,ith
          print*,haiku(ith)%of_nh(nh)%nw,haiku(ith)%of_nh(nh)%of_w(1)%nhaiku
          print*,haiku(ith)%of_nh(nh)%of_w(1)%jz
          print*,ip,parmin,parmax
          print*,haiku(ith)%of_nh(nh)%of_w(1)%par

          stop
        endif
      enddo ! par

      do ip = parmin,parmax
        jzmin = blockmap(ith)%nh(nh)%par(ip)%jzmin
        jzmax = blockmap(ith)%nh(nh)%par(ip)%jzmax
        do jz = jzmin,jzmax,2
!----------------W --------------------------------
!----------- FIND "OCCUPIED" W FOR FIXED JZ,PARITY
          wallow = .false.
          do iw = 1,haiku(ith)%of_nh(nh)%nw

            w = haiku(ith)%of_nh(nh)%wlist(iw)
            n = haiku(ith)%of_nh(nh)%of_w(iw)%nhaiku

            do ih = 1,n
               m = haiku(ith)%of_nh(nh)%of_w(iw)%jz(ih)
               parh=haiku(ith)%of_nh(nh)%of_w(iw)%par(ih)
               if(m == jz .and. ip == parh)then
                 wallow(w) = .true.
               endif
            enddo
         enddo  !iw
!--------------- SET UP W ARRAYS ---------------------
         wmn = 1000
         wmx = 0
         n= 0
         do w = wmin,wmax
          if(wallow(w))then
            n = n+1
            wmx = bmax(wmx,w)
            wmn = bmin(wmn,w)
          endif
         enddo  ! w
         blockmap(ith)%nh(nh)%par(ip)%jz(jz)%nw = n
         if(n >0)then
            allocate(blockmap(ith)%nh(nh)%par(ip)%jz(jz)%wlist(n), stat=aerr)
            if(aerr /= 0) call memerror("count_blocks 10")
            allocate(blockmap(ith)%nh(nh)%par(ip)%jz(jz)%hblock(n), stat=aerr)
            if(aerr /= 0) call memerror("count_blocks 11")
            allocate(blockmap(ith)%nh(nh)%par(ip)%jz(jz)%wmap(wmn:wmx), stat=aerr)
            if(aerr /= 0) call memerror("count_blocks 12")
         endif
         n = 0
         do w = wmn,wmx
          if(wallow(w))then
            n = n+1
            blockmap(ith)%nh(nh)%par(ip)%jz(jz)%wlist(n) = w
            blockmap(ith)%nh(nh)%par(ip)%jz(jz)%wmap(w) = n
            nb = nb +1
            blockmap(ith)%nh(nh)%par(ip)%jz(jz)%hblock(n)= nb
!            write(24,*)ith,nh,jz,ip,n,nb
          endif

         enddo  ! w
         enddo  ! jz
!--------------------------------------------------
enddo  ! ip

deallocate(wallow)
return
end subroutine count_blocks

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine fill_blocks
!
!  USING TYPE blockmap, CONSTRUCTS hblock
!
!  INPUT:
!   ith = species/handness
!   nh = # of particles
!   nb = # of blocks so far (can add to them)
!
!  OUTPUT
!    maxhaiku = max # of haikus for a block for this ith,nh
!    sumhaiku = # of haikus for this ith/nh
!

subroutine fill_blocks_OLD(ith,nh,nb,maxhaiku,sumhaiku)

use sporbit
use haiku_info
use haikus
use verbosity
use butil_mod
implicit none
integer ith
integer nh
!-----------------------------
integer nb
integer maxhaiku,sumhaiku
!---------------------------
integer ih
integer jzmin,jzmax,jz,jzh
integer iw,w,iwh
integer parmin,parmax,ip
integer wmax,wmin,wmx,wmn
integer nw
integer m
integer parh
integer nhaiku
integer nhsd
integer jzstart,jzstop
integer it
integer :: aerr

jzstop = -1000000  ! KSM - initialize for -Wuninitialized

!................................BEGIN..................

it = abs(ith)

parmin = blockmap(ith)%nh(nh)%parmin
parmax = blockmap(ith)%nh(nh)%parmax

maxhaiku = 0
sumhaiku = 0
do ip = parmin,parmax
   jzmin = blockmap(ith)%nh(nh)%par(ip)%jzmin
   jzmax = blockmap(ith)%nh(nh)%par(ip)%jzmax 
   do jzh = jzmin,jzmax,2
      jz = jzh
      nw = blockmap(ith)%nh(nh)%par(ip)%jz(jzh)%nw
      do iw = 1,nw
          w = blockmap(ith)%nh(nh)%par(ip)%jz(jzh)%wlist(iw)
          iwh = haiku(ith)%of_nh(nh)%wmap(w)
          nhaiku = haiku(ith)%of_nh(nh)%of_w(iwh)%nhaiku
!-------------------- find start, stop of these haikus --------
!  NOTE: Here I am using the fact that the haikus have been sorted
!  however it is possible to not sort them and just assign the values
!  but this possible revision can wait for later and may not be needed
!
          jzstart = 0

          do ih = 1,nhaiku
             parh = haiku(ith)%of_nh(nh)%of_w(iwh)%par(ih)
             if(parh /= ip)cycle
             m    = haiku(ith)%of_nh(nh)%of_w(iwh)%jz(ih)
             if(m == jzh)then
                 jzstart = ih
                 exit
             endif
          enddo

          if(jzstart == 0)then  ! none fit
              jzstop = 0
              nhsd = 0
              cycle
          endif
          do ih = jzstart,nhaiku
             parh = haiku(ith)%of_nh(nh)%of_w(iwh)%par(ih)
             if(parh /= ip)cycle
             m    = haiku(ith)%of_nh(nh)%of_w(iwh)%jz(ih)
             if(m  > jzh)exit
             jzstop = ih
          enddo

          nhsd = jzstop - jzstart + 1
          nb = nb + 1
          hblock(ith)%list(nb)%nh = nh
          hblock(ith)%list(nb)%wh = w
          hblock(ith)%list(nb)%jzh = jzh
          hblock(ith)%list(nb)%parh = ip
          hblock(ith)%list(nb)%nhsd = nhsd
!--------------------
! ONE CONCERN I HAVE IS THAT I AM DUPLICATING THE STORAGE OF THE HAIKUS
! this may or may not be an issue later on, but one we should be aware of
!  Merely deallocating the original set of haikus may also not be effective 
!  because of "memory leaks" where deallocated memory ends up unusable
!            
          allocate(hblock(ith)%list(nb)%hsd(nword(it),nhsd), stat=aerr)
          if(aerr /= 0) call memerror("fill_blocks_OLD")
          maxhaiku = bmax(maxhaiku,nhsd)
          sumhaiku = sumhaiku+nhsd
          do ih = 1,nhsd
            hblock(ith)%list(nb)%hsd(:,ih) =           haiku(ith)%of_nh(nh)%of_w(iwh)%hsd(:,jzstart+ih-1)
          enddo
        enddo  ! iw
    enddo ! jzh
enddo ! ip
return
end subroutine fill_blocks_OLD
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine fill_blocks
!
!  USING TYPE blockmap, CONSTRUCTS hblock
!
!  REVISED March 2013 to avoid sorting
!
!  INPUT:
!   ith = species/handness
!   nh = # of particles
!   nb = # of blocks so far (can add to them)
!
!  OUTPUT
!    maxhaiku = max # of haikus for a block for this ith,nh
!    sumhaiku = # of haikus for this ith/nh
!

subroutine fill_blocks(ith,nh,nb,maxhaiku,sumhaiku)

use sporbit
use haiku_info
use haikus
use verbosity
use butil_mod
implicit none
integer ith
integer nh
!-----------------------------
integer nb
integer maxhaiku,sumhaiku
!---------------------------
integer ih
integer jzmin,jzmax,jz,jzh
integer iw,w,iwh
integer parmin,parmax,ip
integer wmax,wmin,wmx,wmn
integer nw
integer m
integer parh
integer nhaiku
integer nhsd
integer jzstart,jzstop
integer it
integer :: aerr

!................................BEGIN..................

it = abs(ith)

parmin = blockmap(ith)%nh(nh)%parmin
parmax = blockmap(ith)%nh(nh)%parmax

maxhaiku = 0
sumhaiku = 0
do ip = parmin,parmax
   jzmin = blockmap(ith)%nh(nh)%par(ip)%jzmin
   jzmax = blockmap(ith)%nh(nh)%par(ip)%jzmax 
   do jzh = jzmin,jzmax,2
      jz = jzh
      nw = blockmap(ith)%nh(nh)%par(ip)%jz(jzh)%nw
      do iw = 1,nw
          w = blockmap(ith)%nh(nh)%par(ip)%jz(jzh)%wlist(iw)
          iwh = haiku(ith)%of_nh(nh)%wmap(w)
          nhaiku = haiku(ith)%of_nh(nh)%of_w(iwh)%nhaiku
!-------------- DO NOT ASSUME HAIKUS ARE SORTED -------------------
!               COUNT UP # OF HAIKUS WITH PARITY ip, JZ m
          nhsd = 0
          do ih = 1,nhaiku
             parh = haiku(ith)%of_nh(nh)%of_w(iwh)%par(ih)
             m    = haiku(ith)%of_nh(nh)%of_w(iwh)%jz(ih)
             if(parh == ip .and. m==jzh) nhsd = nhsd + 1
          enddo

          nb = nb + 1
          hblock(ith)%list(nb)%nh = nh
          hblock(ith)%list(nb)%wh = w
          hblock(ith)%list(nb)%jzh = jzh
          hblock(ith)%list(nb)%parh = ip
          hblock(ith)%list(nb)%nhsd = nhsd
          hblock(ith)%list(nb)%used = .false.
		  
!--------------------
! ONE CONCERN I HAVE IS THAT I AM DUPLICATING THE STORAGE OF THE HAIKUS
! this may or may not be an issue later on, but one we should be aware of
!  Merely deallocating the original set of haikus may also not be effective 
!  because of "memory leaks" where deallocated memory ends up unusable
!           
          allocate(hblock(ith)%list(nb)%hsd(nword(it),nhsd), stat=aerr)
          if(aerr /= 0) call memerror("fill_blocks")
          maxhaiku = bmax(maxhaiku,nhsd)
          sumhaiku = sumhaiku+nhsd
          nhsd = 0
          do ih = 1,nhaiku
             parh = haiku(ith)%of_nh(nh)%of_w(iwh)%par(ih)
             m    = haiku(ith)%of_nh(nh)%of_w(iwh)%jz(ih)
             if(parh == ip .and. m==jzh)then
		nhsd = nhsd + 1
                hblock(ith)%list(nb)%hsd(:,nhsd) = haiku(ith)%of_nh(nh)%of_w(iwh)%hsd(:,ih)
             end if
          enddo

        enddo  ! iw
    enddo ! jzh
enddo ! ip
return
end subroutine fill_blocks

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! for testing only
!
! CALLS:
!  convert_haiku_array
!
subroutine printblockhaikus(ith)

use haiku_info
use haikus
use verbosity
implicit none
!include 'binterfaces.inc'

integer ith
integer nh
integer ib
integer wh,iwh,nw
integer jzh,parh
integer nhaiku,ih
      
integer,allocatable :: occ(:)
integer,pointer  :: hsd(:)
integer :: aerr

!integer :: nh_of_occ,wh_of_occ
if(.not.print_haikus)return

allocate(occ(nhsps(abs(ith))), stat=aerr)
if(aerr /= 0) call memerror("printblockhaikus")
write(haiku_file,*)' '
write(haiku_file,*)' HAIKUS ',ith
do ib = 1,hblock(ith)%nhblocks
   nh = hblock(ith)%list(ib)%nh
   wh = hblock(ith)%list(ib)%wh
   parh = hblock(ith)%list(ib)%parh
   jzh = hblock(ith)%list(ib)%jzh
   write(haiku_file,101)ib,nh,jzh,parh,wh
101    format('block ',i3,' : ',3i3)
   do ih = 1,hblock(ith)%list(ib)%nhsd
      hsd => hblock(ith)%list(ib)%hsd(:,ih)

      call convert_haiku_array(ith,hsd,occ)
      write(haiku_file,103)ih,occ
103      format(6x,i3,':',3x,60i1)
   enddo

enddo ! ih

deallocate(occ)
return
end subroutine printblockhaikus

!---------------------------------------
!
!  count up how many used and unused haikus there are
!  FOR DIAGNOSIS PURPOSES ONLY -- added 7.7.2
!
!  INITIAL FINDINGS: for large W cuts, between 75% and 90-95% of 
!  haikus go unused. 
!  Solution strategy: 
!  (1) calculate a "skeleton," that is, a list of haikus of every possible nh and jz etc
!  (2) Find limits for basis and thus put a limit (max/ Jz on each template)
!  (3) With this, one can count ahead of time how many haikus in a block; this might
!      mitigate the need to sort and re-store haikus as currently done
!

subroutine countusedhaikus
	
	use nodeinfo
	
	implicit none
	integer(8) :: allhaikus,usedhaikus
	integer :: it,hsign,nh,jzh,parh,wh
	integer :: ith
	integer :: iblock
	
	if(iproc/=0)return
	
	
	allhaikus = 0
	usedhaikus = 0
	
	do it = 1,2
		do hsign = 0,1
			ith = it*(-1)**hsign
			do iblock = 1,hblock(ith)%nhblocks
				allhaikus = allhaikus + hblock(ith)%list(iblock)%nhsd
				if( hblock(ith)%list(iblock)%used)usedhaikus = usedhaikus + hblock(ith)%list(iblock)%nhsd
				
			end do
			
			
		end do
		
		
	end do
	print*,usedhaikus,' out of ',allhaikus,' haikus used '
	
	
	return
	
	
end subroutine countusedhaikus

!----------------------------------------

end module blocks

