!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC.true.
!  BIGSTICK CI shell-model code
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  MAKING JUMPS BETWEEN SDs
!
!  JUMPS ARE ORGANIZED BETWEEN INITIAL AND FINAL SECTORS 
!  (a "sector" is a set of SDs that have the same quantum #s)
!  
!  WITHIN SECTORS, SETS OF INITIAL AND FINAL BLOCKS ARE FOUND 
!  Once initial/final blocks are found, hops from block to blocks are used to generate the jumps
!
!  NOTE that both of these processes can and should be distributed in parallel
!
!=======================================================================
!
!  HOW TO MANUFACTURE JUMPS FROM HOPS
!
!  A little review:
! 
!  We organize haikus into blocks (labeled by quantum numbers nh, jz, parity, W) and combine
!  left and right blocks into sectors of SDs (also labeled by combined quantum numbers)
!  Every SD is a combination of a left haiku and a right haiku. Within a sector there are
!  several conjugate sets of left and right haiku blocks.
!
!  For example, suppose we have a sector of proton SDs  with Jz = +3. (Ignore parity, W for now)
!  Let Z = 4. Then we might have
!
!     left block                           right block
!     nh = 0   jz = 0                       nh = 4  jz = 3
!     nh = 1   jz = -1/2                    nh = 3  jz = 7/2
!     nh = 1   jz = -3/2                    nh = 3  jz = 9/2
!     nh = 2   jz = -1                      nh = 2  jz = 4
!     nh = 2   jz = -2                      nh = 2  jz = 5
!       and so on
!
!  Single-particle states are organized by groups (quantum numbers m, parity, W, not J)
!
!  Hops are the actions of creation / annihilation operators on haikus. These are 
!  organized by quantum number into block hops between blocks. Because of this. 
!  all hops within a block hop have operators within the same group.
!  So all the hops between (say) nh=3, jz =9/2 and nh =2 jz = 2 would be mediated by 
!  all single-particle destruction operators with jz = 5/2 (which all belong to the same group) 
!
!  Now to jumps: Jumps are the actions of several hops, that is, a^+a. a^+ a^+ a a, etc.
!
!  Jumps are organized between sectors of SDs. Thus they too are organized by quantum numbers.
!  Jumps are created in several steps. 
!
!  STEP 1: SPAWN DESCENDENTS
!    A DESCENT is a list of operator "groups" that move us from one set of haiku blocks
!  to another. These are lists of the group labels. The convention is, if the group label
!  is < 0, then it applies to the left block; if > 0 then applies to the right block.
!     We do only destruction descents; for N-body jumps, we apply N hops to get a descent;
!  the descent is a list of N group labels. To prevent  double counting they are ordered
!  in decreasing value. However one can in principle have a group label repeated (because you
!  can have different operators from the same group). 
!     THIS IS DONE recursively by action of the routine SPAWNDESCENT starting from a sector
!  For a N-body jump, each descent is an array of length N.

!   
!  STEP 2: MATCH DESCENTS AND MAKE GENEOLOGIES
!     For N-body jumps between two sectors, we find all the groupjumps from N applications 
!  of destructions operators, starting from each block in each sector. By applying the 
!  group operators, we find out the final left and right blocks at the end of each 
!  groupjump. We then search to match up left and right haiku blocks. Because there are 
!  far fewer blocks than haikus (and far fewer haikus than SDs) this is more efficient 
!  than searching SDs. 
!     When two groupjumps are matched, we then construct a "geneology." A geneology is 
!  just like a descent, although includes both destruction and creation group operations.
!  We find all geneologies between two sectors. The geneology is stored just like a groupjump.
!  and takes us from some initial pair of haiku blocks to a final pair of haiku blocks.
!     This is done by routine MATCHMAKER
!
!  STEP 3: MAKE LEFT AND RIGHT CHAINS
!    Each geneology is a set of group operations that take us from an initial pair of 
!  haiku blocks to final blocks. We now split those geneologies into chains of haiku hops.
!  because the geneologies factorize naturally. 
!    Chains are *almost* like jumps, in that they consist of:
!    initial sd, final sd, phase, and operator info
!   This is done by routine **** which class routine ADD_A_LINK
!
!  STEP 4: COMBINE CHAINS INTO JUMPS
!    Done by routine WELD. Jumps factorize, which means given a geneology, 
!  we can multiply all the left chains for that geneology by all the right chains 
!  and get out all the jumps.
!    
!
!================================================================================
!
!  "GENEOLOGIES"  -- the convention
!   a geneology starts from a specific haiku block
!   each haiku block has some number of geneologies descended from it (ngeneologies, nspawn, etc)
!   each geneology is specified by a list of indices of the haiku blocks
!
!   
!=====================================================================
!
!  BJUMPLIB1.f90: master jump generation routines 
! + routines to find 1-body block jumps between sectors
!
!==========================================================================
!
!
!  NOTE: One way to speed up is to call matchmaker fewer times. We can do this by 
!  checking quantum numbers; not all sectors can connect by a 1-body operator, etc.
!
!==============================================================
module jump_mod
   
contains
!=====================================================================
!
!  subroutine hoopmaster
!
!  routines to prepare for jumps
!  CALLED BY: bigstick_main program
!
!  SUBROUTINES CALLED:
!	clocker 	: internal timing
!	bunnymaster	: creation/annihilation "hops"
!	print_out_hops  : (optional)
!       linkupgrouphops 
! 	threebodyquestion : inquires if 3-body forces will be used
!	master_q_ntuple :  used to restrict quantum numbers in XX, XXX jumps
!  	master_cross_ntuple : restrict quantum numbers in XY, XXY jumps
!	setlinkID	: mostly obsolete
!	prepare_to_uncoupleXXtbme	
!     	prepare_to_uncouplePNtbme
!       threebodysetup
! 	master3bodyjumps
!       masterfindconjumps3b
!       set1bsectorjumps
!       set2bsectorjumps
!       masterXXYconjugatejumps
! 	master2bodyjumps
! 	master1bodyjumps
!	set2bsectorjumps
!       masterfindconjumps2b
!	masterfindconjumps1b

subroutine hoopmaster(hermflag,whermflagp,whermflagn,change1bodyqs,change2bodyqs,change3bodyqs)
   use system_parameters
   use io
   use menu_choices
   use verbosity
   use flags3body
   use nodeinfo
   use timing
   use btbme_mod
   use hoppy
   use ntuple_info
   
   implicit none
!----------- FLAGS FOR RESTRICTION ON CHANGING q#s in jumps-----------------
  logical    :: change1bodyqs, change2bodyqs, change3bodyqs
  logical    :: hermflag, whermflagp, whermflagn
  integer :: ith, its
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!--------------------- CREATION OF INTERACTION / OPERATORS --------------
! 1. create 'hops' (fundamental creation/annihilation)
! 2. read in interaction matrix elements
! 3. create n-body 'jumps'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  jumptime = 0.d0

  call hopmaker

!--------------------------- FIRST ROUND OF CREATING JUMPS --------

!----------------- RESTRICTIONS
! In order to save space, one needs to restrict the uncoupled matrix elements 
! created.  In order to do this, one needs the hops beforehand
  call master_q_ntuple(2,1)
  call master_q_ntuple(2,2)
  call master_cross_ntuple(1,1)   ! DON'T GET RID OF THIS... YET

!--------------- JUMPS------------------------------
!  NOTE: It is important to create 2-body jumps before 1-body
!  This is because certain linked lists are reused to plug memory leaks;
!  and an array contained in 2-body lists must be larger than for 1-body 
!
  call setlinkID  ! for debugging
  if(iproc==0)print*,' '

!----------------- COUNT UP JUMPS/OPERATIONS -------------
!........ SET UP NEEDED ARRAYS.....
  if(iproc == 0)then
     print*,' '
     print*,' .... Setting up arrays for matrix elements ... '
     print*,' '
  end if
  call clocker('mes','sta')

!........ THIS WILL HAVE TO BE PUT ELSEWHERE....

  call prepare_to_uncoupleXXtbme(1)
  call prepare_to_uncoupleXXtbme(2)
!  call prepare_to_uncouplePNtbme
  call clocker('mes','end')

  if(iproc == 0)then
     print*,' '
     print*,' .... Counting jumps ... '
     if(verbose_jumps) print*,' Information about jumps: '
     print*,' '
  end if
  call clocker('jmc','sta')
  if(iproc==0 .and. print_sectorjumps)open(unit=sectorjump_file,file='sectorjumps.bigstick',status='unknown')
!................... 3-body...........................
  if ( threebody ) then
        change2bodyqs = .true.
        call threebodysetup
!..................... PPP/NNN............................
        if ( iproc == 0 ) print*,' 3-body jumps '		
        call master3bodyjumps(1,.false. , change3bodyqs )		
        call master3bodyjumps(2,.false. , change3bodyqs )
        call masterfindconjumps3b(1)
        call masterfindconjumps3b(2)

!.......................PPN/PNN ........................
!........... set up sector jumps ................
!......... NB: I do not use hermiticity for 1-body jumps, but automatically do for 2-body jump

!.....  PAY ATTENTION TO HERMITICITY!!!

        if( np(1)*(np(2) -1 ) > 0)then   ! PNN
           call set1bsectorjumps(1, .false. , .true. , .false. ,   .true. )
           call set1bsectorjumps(1, .true. , .true. , .false. ,   .true. )
           call set2bsectorjumps(2,.false., .true., .true.,   .true. , .true., 1)
           call set2bsectorjumps(2,.true., .true., .true.,   .true. ,.true., 1)
           call masterXXYconjugatejumps(2)

        endif

        if( np(2)*(np(1) -1 ) > 0)then   ! PPN
           call set1bsectorjumps(2, .false. , .true. , .false. ,   .true. )
           call set1bsectorjumps(2, .true. , .true. , .false. ,   .true. )

           call set2bsectorjumps(1,.false.,  .true., .true.,  .true. ,.true., 1)
           call set2bsectorjumps(1,.true.,  .true., .true.,  .true. ,.true., 1)

           call masterXXYconjugatejumps(1)
        endif

        if(iproc==0)print*,' 2-body jumps for PPN/PNN '

        call master2bodyjumps(1,.false. )
        call master2bodyjumps(2,.false.)

        if ( iproc == 0 ) then
           print*,' '
           print*,' 1-body jumps for PPN/PNN '
        end if
        call master1bodyjumps(1,.false. )
        call master1bodyjumps(2,.false. )                  !  changed 12/16/2012 by weo to remove 2nd logical

     endif
!.....................................................
! --- END OF 3-BODY OPERATION

     if(.not.threebody .and. menu_char /= 'dx'.and. menu_char /= 'o')then     

     if(iproc==0)print*,' 2-body jumps '

     if(np(1) > 1)then

         call set2bsectorjumps(1,.false.,  hermflag, .true. ,change2bodyqs, .true., 2 )
         call set2bsectorjumps(1,.true.,  hermflag, .true. ,change2bodyqs,.true., 2 )
         call masterfindconjumps2b(1)
        call master2bodyjumps(1,.false.  )

     endif

     if(np(2) > 1)then
         call set2bsectorjumps(2,.false.,  hermflag, .true. ,change2bodyqs,.true., 2 )
         call set2bsectorjumps(2,.true.,  hermflag, .true. ,change2bodyqs,.true., 2 )
       call masterfindconjumps2b(2)
       call master2bodyjumps(2,.false. )
     endif

!     call Analyze3bodyJumps(1)
!     call Analyze3bodyJumps(2)
     if ( iproc == 0 ) then
        print*,' '
        print*,' 1-body jumps '
     end if

     if( np(1)*np(2) > 0)then
!................ FIND SECTORS CONNECTED BY 1-BODY JUMPS.........
! NB: logical flags are: fill; hermflag; whermflag 
!
        call set1bsectorjumps(1, .false. , hermflag ,whermflagp , change1bodyqs )
        call set1bsectorjumps(1, .true.  , hermflag, whermflagp , change1bodyqs )

! NOTE that I use "hermiticity" on 1-body jumps but not "W-hermiticity"

        call set1bsectorjumps(2, .false. , hermflag, whermflagn , change1bodyqs )
        call set1bsectorjumps(2, .true.  , hermflag, whermflagn, change1bodyqs )
        call masterfindconjumps1b(1,change1bodyqs)
        call masterfindconjumps1b(2,change1bodyqs)  ! not really needed

        call master1bodyjumps(1,.false.)
        call master1bodyjumps(2,.false. )
!--------------- find conjugate jumps
        call masterfindconjumps1b(1,change1bodyqs)

        call masterfindconjumps1b(2,change1bodyqs)  ! not really needed

     endif     

     endif    ! if .not. threebody
     if(iproc == 0)print*,'  '
  call clocker('jmc','end')

     return

end subroutine hoopmaster

!=====================================================================
!
!  subroutine jumpmaster
!
!  routines to create jumps
!  CALLED BY: bigstick_main program
!
!  SUBROUTINES CALLED:
!	master3bodyjumps
!	master2bodyjumps
! 	master1bodyjumps
!	write_out_jumps	: (optional)
!
subroutine jumpmaster(hermflag,whermflagp,whermflagn,change1bodyqs,change2bodyqs,change3bodyqs)

   use system_parameters
   use io
   use menu_choices
   use verbosity
   use flags3body
   use nodeinfo
   use jumplimits
   implicit none
!----------- FLAGS FOR RESTRICTION ON CHANGING q#s in jumps-----------------
   logical    :: change1bodyqs, change2bodyqs, change3bodyqs
   logical    :: hermflag, whermflagp, whermflagn

   if(iproc == 0)then
      print*,' '
      print*,' .... Building jumps ... '
      print*,' '
   end if

!----------------------- GENERATE JUMPS------------------
!................... 3-body...........................
   if ( threebody ) then
        change2bodyqs = .true.
!..................... PPP/NNN............................
        if ( iproc == 0 ) print*,' 3-body jumps '
        if(ham_readin)call master3bodyjumps(1,.true., change3bodyqs)
        if(ham_readin)call master3bodyjumps(2,.true., change3bodyqs)

!.......................PPN/PNN ........................
!........... set up sector jumps ................
!......... NB: I do not use hermiticity for 1-body jumps, but automatically do for 2-body jump

!.....  PAY ATTENTION TO HERMITICITY!!!


        if(iproc==0)print*,' 2-body jumps for PPN/PNN '
        if(ham_readin)call master2bodyjumps(1, .true. )		
        if(ham_readin)call master2bodyjumps(2, .true. )
		

        if ( iproc == 0 ) then
           print*,' '
           print*,' 1-body jumps for PPN/PNN '
        end if
        if(ham_readin)call master1bodyjumps(1,.true. )
        if(ham_readin)call master1bodyjumps(2,.true.)
       
     endif

!.....................................................
! --- END OF 3-BODY OPERATION

   if(.not.threebody .and. menu_char /= 'dx'.and. menu_char /= 'o')then     
      if(iproc==0)print*,' 2-body jumps '
      if(np(1) > 1)then
         if(ham_readin)call master2bodyjumps(1,.true. )
      endif
      if(np(2) > 1)then
         if(ham_readin)call master2bodyjumps(2,.true.)
      endif

      if ( iproc == 0 ) then
         print*,' '
         print*,' 1-body jumps '
      end if
      if( np(1)*np(2) > 0)then
         if(ham_readin)call master1bodyjumps(1,.true.)
         if(ham_readin)call master1bodyjumps(2,.true. )
      endif     

   endif    ! if .not. threebody

   if(print_jumps)  then
      call write_out_jumps(1,1,'o')  ! open a file
      call write_out_jumps(2,1,'x')
      call write_out_jumps(1,2,'x')
      call write_out_jumps(2,2,'c')  ! close file when finished
   end if

   return
end subroutine jumpmaster

!=====================================================================


subroutine master1bodyjumps(it,fill)

  use system_parameters
  use verbosity
  use sectors
  use descendents
  use geneologies
  use chains
  use welder
  use jumpNbody
  use nodeinfo
  use jumplimits
  use flagger
  use timing
  use bmpi_mod

  implicit none
 
! KSM TEST  include 'binterfaces.inc'
  integer :: ierr

  integer :: it  ! species
  type (groupjump) :: descent
  type (genlink), pointer :: curgen
  integer(8) :: ngeneologies,igen
  integer :: is,fs
  integer :: nbody
  character, pointer :: optype(:)*1
  integer :: isjmp
  integer(8) :: njumps,nlinks,nzero
  integer(8) ::  totjumps,totlinks
  logical :: fill
  logical :: alreadyallocated   

  integer(8), allocatable :: jumpcount(:),jumpreduce(:)   
                                      ! used to pass jump information in MPI
  integer(8):: minjump,maxjump,startjump,endjump,shiftjump
  integer :: iprocs
  integer :: aerr

  
!--------- KSM - initialize for -Wuninitialized
! set to values that would cause problems if used
  minjump = 100000000
  maxjump = -10000000

!--------------------------------- SET UP DESCENDENTS----------------
!       NOTE FOR FUTURE PARALLELIZATION 
!       CAN DISTRIBUTE -- ONLY GET DESCENDENTS FOR RELEVANT NODES
!

!  if(np(1)*np(2)==0)return  ! no proton-neutron interaction
  nbody = 1

  if(associated(optype))nullify(optype)
  if(.not.associated(optype)) then
     allocate(optype(2*nbody), stat=aerr)
     if(aerr /= 0) call memerror("master1bodyjumps 1")
  end if

  optype(1) = 'd'
  optype(2) = 'c'

!-------------------- FIGURE OUT WHICH SECTORS WE WANT TO GO TO
! NB: taken outside subroutine May 2010 by CWJ
!     in order to expedite 3-body interactions

  descent%nsectors = nsectors(it)
  allocate(descent%sector(nsectors(it)), stat=aerr)
  if(aerr /= 0) call memerror("master1bodyjumps 2")
  call clocker('des','sta')
!$OMP PARALLEL DO SCHEDULE(STATIC)
  do is = 1,nsectors(it)
     call sectordescent(it,is,nbody,descent%sector(is))
  enddo ! is
!$OMP END PARALLEL DO

  call clocker('des','end')
  
  if(.not.fill) then
     jumptimestart = BMPI_WTIME()  ! MPI clock
  end if
!==========================================================================

  totjumps = 0
  if(.not.fill )then
      allocate( jumpcount( x1bjump(it)%nsectjumps ), jumpreduce(x1bjump(it)%nsectjumps) , stat=aerr)
      if(aerr /= 0) call memerror("master1bodyjumps 3")
      jumpcount(:) = 0
      jumpreduce(:) = 0
  end if
  if(fill)then
     if(restrictjumps .and. nproc > 1 )then
        if(it ==1)then
           if(.not.makep1bjumps(iproc))return
           minjump = startp1bjumps(iproc)
           maxjump =  stopp1bjumps(iproc)
        end if
        if( it ==2)then
!          if(iproc ==0)print*,' NEUTRON 1B ',maken1bjumps(iproc), startn1bjumps(iproc), stopn1bjumps(iproc)
           if(.not.maken1bjumps(iproc))return
           minjump = startn1bjumps(iproc)
           maxjump =  stopn1bjumps(iproc)
        end if

     else
		 if(iproc==0 .and. nproc > 1)print*,' NOT RESTRICTING JUMP STORAGE '
        minjump = 1
        if(it ==1)then
          maxjump = totp1bjumps
        end if
        if(it ==2)then
         maxjump = totn1bjumps
        end if

     end if 

     if(it == 1)then
        if(allocated(p1b_isd))then
           if(size(p1b_isd) < totp1bjumps)then
              if ( iproc == 0 ) then
                 print*,' WARNING RESIZING proton 1b jumps  '
                 print*,totp1bjumps,' jumps increased from ',size(p1b_isd)
			 end if
             deallocate(p1b_isd, p1b_fsd,p1b_cop,p1b_dop,p1b_phase)
             allocate(p1b_isd( totp1bjumps), stat=aerr)
             if(aerr /= 0) call memerror("master1bodyjumps 10")
             allocate(p1b_fsd( totp1bjumps), stat=aerr)
             if(aerr /= 0) call memerror("master1bodyjumps 11")
             allocate(p1b_cop( totp1bjumps), stat=aerr)
             if(aerr /= 0) call memerror("master1bodyjumps 12")
             allocate(p1b_dop( totp1bjumps), stat=aerr)
             if(aerr /= 0) call memerror("master1bodyjumps 13")
             allocate(p1b_phase( totp1bjumps), stat=aerr)
             if(aerr /= 0) call memerror("master1bodyjumps 14")
             if(print_jumps)then
                 deallocate(p1b_isd0,p1b_fsd0)
                 allocate(p1b_isd0(totp1bjumps) , stat=aerr)
                 if(aerr /= 0) call memerror("master1bodyjumps 15")
                 allocate(p1b_fsd0(totp1bjumps) , stat=aerr)
                 if(aerr /= 0) call memerror("master1bodyjumps 16")
             endif 
           end if
        else
           allocate(p1b_isd( minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 20")
           allocate(p1b_fsd( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 21")
           allocate(p1b_cop( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 22")
           allocate(p1b_dop( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 23")
           allocate(p1b_phase(minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 25")

           if(print_jumps)then
              allocate(p1b_isd0(totp1bjumps) , stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 26")
              allocate(p1b_fsd0(totp1bjumps) , stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 27")
           endif
        endif
        
     else
        if(allocated(n1b_isd))then
           if ( size(n1b_isd) < totn1bjumps ) then
              if ( iproc == 0 ) then
                 print*,' WARNING RESIZING neutron 1b jumps '
                 print*,totn1bjumps,' jumps increased from ',size(n1b_isd)
              end if
              deallocate(n1b_isd, n1b_fsd,n1b_cop,n1b_dop,n1b_phase)
              allocate(n1b_isd( minjump:maxjump), stat=aerr)
              if(aerr /= 0) call memerror("master1bodyjumps 30")
              allocate(n1b_fsd( minjump:maxjump), stat=aerr)
              if(aerr /= 0) call memerror("master1bodyjumps 31")
              allocate(n1b_cop( minjump:maxjump), stat=aerr)
              if(aerr /= 0) call memerror("master1bodyjumps 32")
              allocate(n1b_dop( minjump:maxjump), stat=aerr)
              if(aerr /= 0) call memerror("master1bodyjumps 33")
              allocate(n1b_phase( minjump:maxjump), stat=aerr)
              if(aerr /= 0) call memerror("master1bodyjumps 34")
              if(print_jumps)then
                 deallocate(n1b_isd0,n1b_fsd0)
                 allocate(n1b_isd0(totn1bjumps) , stat=aerr)
                 if(aerr /= 0) call memerror("master1bodyjumps 35")
                 allocate(n1b_fsd0(totn1bjumps) , stat=aerr)
                 if(aerr /= 0) call memerror("master1bodyjumps 36")
              endif 
           end if
        else  ! not allocated
			if(startn1bjumps(iproc) /= minjump .or. stopn1bjumps(iproc)/=maxjump)then
				print*,' problem with unmmatched jump limits ',iproc,minjump,maxjump,startn1bjumps(iproc),stopn1bjumps(iproc)
			end if
			
           allocate(n1b_isd( minjump:maxjump), stat=aerr)
           if(aerr /= 0) then
			   print*,iproc,' n1b_isd ',minjump,maxjump
			   if(iproc==189)print*,' target proc ',startn1bjumps(iproc),stopn1bjumps(iproc)
			    call memerror("master1bodyjumps 40")
			end if
           allocate(n1b_fsd( minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 41")
           allocate(n1b_cop( minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 42")
           allocate(n1b_dop( minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 43")
           allocate(n1b_phase( minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master1bodyjumps 44")
           if(print_jumps)then
              allocate(n1b_isd0(totn1bjumps) , stat=aerr)
              if(aerr /= 0) call memerror("master1bodyjumps 45")
              allocate(n1b_fsd0(totn1bjumps) , stat=aerr)
              if(aerr /= 0) call memerror("master1bodyjumps 46")
           endif
        endif
     endif

  endif

  if(fill .and. restrictjumps .and. compactjumpstorage) call surveysectorjumps(it,1)
  do isjmp = 1,x1bjump(it)%nsectjumps
     if(.not.fill .and. nproc > 1  & 
           .and. mod(isjmp,nproc) /= iproc .and. MPIjumps)cycle  ! for counting up jumps

!---- NEED TO REWRITE---
!     if( fill  .and. restrictjumps .and. & 
!             x1bjump(it)%sjmp(isjmp)%nstart+x1bjump(it)%sjmp(isjmp)%njumps < minjump )cycle  
!	 if(fill .and. restrictjumps .and. compactjumpstorage)then
!		 startjump = x1bjump(it)%sjmp(isjmp)%nstart+1
!		 if(it==1)then
!			 call fetchjumpshift('P1B',startjump)
!		 else
!			 call fetchjumpshift('N1B',startjump)
!		 end if
!		 x1bjump(it)%sjmp(isjmp)%nstart = startjump-1
!		 if (x1bjump(it)%sjmp(isjmp)%nstart< 0)cycle  ! nstart < 0 signals to skip these jumps
!	 end if

!........ DEFAULT CHOICES......
     shiftjump = 0
	 startjump = 	 x1bjump(it)%sjmp(isjmp)%nstart+1
	 endjump   =     x1bjump(it)%sjmp(isjmp)%nstart+x1bjump(it)%sjmp(isjmp)%njumps
	 if(fill .and. restrictjumps .and. compactjumpstorage)then  ! SET LIMITS ON STORAGE, SHIFT IN STORAGE
		 if(.not.x1bjump(it)%sjmp(isjmp)%containsexons)cycle
		 startjump = x1bjump(it)%sjmp(isjmp)%exstart
		 endjump   = x1bjump(it)%sjmp(isjmp)%exstop
		 shiftjump = x1bjump(it)%sjmp(isjmp)%exshift
		 
	 end if
     if( fill .and. restrictjumps .and. x1bjump(it)%sjmp(isjmp)%nstart-shiftjump > maxjump )cycle  
	 	 
!.......... SKIP IF NO CONJUGATE JUMPS...........
     if(.not.fill ) x1bjump(it)%sjmp(isjmp)%njumps = 0

     if( x1bjump(it)%csjmp(isjmp)%ncjmps ==0) then
         x1bjump(it)%sjmp(isjmp)%njumps = 0
         x1bjump(it)%sjmp(isjmp)%nstart = 0 ! totjumps
         cycle
     endif
     is = x1bjump(it)%isector(isjmp)
     fs = x1bjump(it)%fsector(isjmp)
     if(.not.fill)then
         x1bjump(it)%sjmp(isjmp)%nstart = 0  !totjumps
         x1bjump(it)%sjmp(isjmp)%njumps = 0
     endif

     call matchmaker(it,nbody,is,descent%sector(is),nbody,fs,descent%sector(fs),.false.,  ngeneologies)
     if(ngeneologies == 0)cycle

     curgen =>geneology
     njumps =   x1bjump(it)%sjmp(isjmp)%nstart
!	   if(fill .and. it==1 .and. iproc==1)print*,' going into Weld',isjmp,startjump,endjump,shiftjump
     do igen = 1,ngeneologies
        call chaingang(it,nbody*2,curgen,optype)
        call weld(it,is,fs,nbody,isjmp,.false., fill , x1bjump(it)%sjmp(isjmp),startjump,endjump,shiftjump,njumps)
     end do  ! igen
     njumps = njumps -x1bjump(it)%sjmp(isjmp)%nstart
     if(.not.fill) then
            x1bjump(it)%sjmp(isjmp)%njumps = njumps
            jumpcount(isjmp) = njumps
     end if

  enddo ! isjmps

  if(.not.fill ) then
     jumptime= BMPI_WTIME()-jumptimestart  ! MPI clock
!	 print*,iproc,jumptime
	 
  end if

  if(.not.fill .and. nproc > 1 .and. MPIjumps)then   ! broadcast jumps via MPI
     call BMPI_ALLREDUCE(jumpcount,jumpreduce,x1bjump(it)%nsectjumps, & 
          MPI_SUM,icomm,ierr)
     do isjmp = 1,x1bjump(it)%nsectjumps
         x1bjump(it)%sjmp(isjmp)%njumps = jumpreduce(isjmp)  
     end do
  end if 

  do isjmp = 1,x1bjump(it)%nsectjumps
     if(.not.fill )then
         x1bjump(it)%sjmp(isjmp)%nstart = totjumps
     end if
     totjumps = totjumps + x1bjump(it)%sjmp(isjmp)%njumps

  end do

  if(.not.fill)then
     if(it == 1)then
        totp1bjumps = totjumps
        if(iproc==0 .and. verbose_jumps) write(6,*)totjumps,' proton jumps ' 
     else
        totn1bjumps = totjumps
        if(iproc==0 .and. verbose_jumps)write(6,*)totjumps,' neutron jumps ' 
     end if
  end if

  deallocate(descent%sector)

  nullify(optype)

  if(.not.fill )deallocate( jumpcount, jumpreduce)
  completed1bodyflag = 123
  
!  if(fill .and. it==1 .and. iproc==1)print*,iproc,it,' here we go ',p1b_cop(:)
  
  return
end subroutine master1bodyjumps

!====================================================================
!
!  NOTE: master2bodyjumps SHOULD BE CALLED BEFORE master1bodyjumps
!
subroutine master2bodyjumps(it,fill)

  use system_parameters
  use basis
  use verbosity
  use sectors
  use descendents
  use geneologies
  use chains
  use jumpNbody
  use jump3body
  use nodeinfo
  use flags3body
  use jumplimits
  use flagger
  use timing
  use bmpi_mod
  use butil_mod
  use welder
  use chains
  implicit none
!  include 'binterfaces.inc'
  integer :: ierr
  integer :: aerr

  integer :: it  ! species
  type (groupjump) :: descent
  type (genlink), pointer :: curgen
  integer(8) :: ngeneologies,igen
  integer :: is,fs
  integer :: nbody
  character, pointer :: optype(:)*1
  !   type (handchain), allocatable :: mrchain(:)
  integer :: isjmp
  integer(8):: njumps,nlinks
  integer(8) :: totjumps,totlinks
  integer(8) :: nzero, nzerotot
  logical :: fill
  integer :: iprocs
  integer(8), allocatable :: jumpcount(:),jumpreduce(:)   
                                      ! used to pass jump information in MPI
  integer(8):: minjump,maxjump,startjump,endjump,shiftjump
  real(8) :: descenttime
  integer(8) :: proddescents  ! PRODUCT OF # OF DESCENDENTS; USED TO ESTIMATE TIME FOR CREATING JUMPS
  real(8)    :: rproddescents

!--- VARIABLES FOR DISTRIBUTING WORK ON COUNTING/CREATING JUMPS..
  integer,allocatable :: mysectjumplist(:)
  integer :: nmysectjumps,sjumplistsize, jj
  
!--------- KSM - initialize for -Wuninitialized
! set to values that would cause visible problems if used
   minjump = 100000000
   maxjump = -10000000
!--------------------------------- SET UP DESCENDENTS----------------
!       NOTE FOR FUTURE PARALLELIZATION 
!       CAN DISTRIBUTE -- ONLY GET DESCENDENTS FOR RELEVANT NODES
!

  if ( completed1bodyflag == 123 .and. iproc == 0 .and. .not. fill )then
     print*,' Warning, need to call master2bodyjumps first, then master1bodyjumps. You may get an error '
  end if
 
 !............. SET UP ARRAYS FOR COUNTING JUMPS ACROSS MPI NODES IN FIRST PASS........ 
  if(.not.fill)then
     allocate( jumpcount( x2bjump(it)%nsectjumps ), jumpreduce(x2bjump(it)%nsectjumps) , stat=aerr)
     if(aerr /= 0) call memerror("master2bodyjumps 1")
      jumpcount(:) = 0
      jumpreduce(:) = 0
  end if
  
  nbody = 2
  if(associated(optype))nullify(optype)
  allocate(optype(2*nbody), stat=aerr)
  if(aerr /= 0) call memerror("master2bodyjumps 2")

  optype(1) = 'd'
  optype(2) = 'd'
  optype(3) = 'c'
  optype(4) = 'c'

  descent%nsectors = nsectors(it)
  allocate(descent%sector(nsectors(it)), stat=aerr)
  if(aerr /= 0) call memerror("master2bodyjumps 3")
  call clocker('des','sta')

!$OMP PARALLEL DO SCHEDULE(STATIC)
  do is = 1,nsectors(it)
     call sectordescent(it,is,nbody,descent%sector(is))
  enddo ! is
!$OMP END PARALLEL DO
  call clocker('des','end')

  if(.not.fill) then
     jumptimestart = BMPI_WTIME()  ! MPI clock
  end if
!==========================================================================

  totjumps = 0

  if(fill)then
     if(restrictjumps .and. nproc > 1 )then
        if(it ==1)then
           if(.not.makeppjumps(iproc))return
           minjump = startppjumps(iproc)
           maxjump =  stopppjumps(iproc)
        end if
        if( it ==2)then
           if(.not.makennjumps(iproc))return
           minjump = startnnjumps(iproc)
           maxjump =  stopnnjumps(iproc)
        end if

     else
        minjump = 1
        if(it ==1)then
          maxjump = totp2bjumps
        end if
        if(it ==2)then
         maxjump = totn2bjumps
        end if
     end if

     if ( it == 1 ) then
        if ( .not.allocated( p2b_isd ) ) then
           allocate(p2b_isd( minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master2bodyjumps 10")
        end if
        if ( .not.allocated( p2b_fsd ) ) then
           allocate(p2b_fsd( minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master2bodyjumps 11")
        end if
        if ( .not.allocated( p2b_phase ) ) then
           allocate(p2b_phase(minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master2bodyjumps 12")
        end if
        if(print_jumps)then
              if ( .not.allocated( p2b_isd0 ) ) then
                 allocate(p2b_isd0(minjump:maxjump) , stat=aerr)
                 if(aerr /= 0) call memerror("master2bodyjumps 13")
              end if
              if ( .not.allocated( p2b_fsd0 ) ) then
                 allocate(p2b_fsd0(minjump:maxjump) , stat=aerr)
                 if(aerr /= 0) call memerror("master2bodyjumps 14")
              end if
        endif

       if(threebody)then
          if ( .not.allocated( p2b_dop ) )  then
             allocate(p2b_dop( minjump:maxjump), stat=aerr)
             if(aerr /= 0) call memerror("master2bodyjumps 20")
          end if
          if ( .not.allocated( p2b_cop ) )  then
             allocate(p2b_cop( minjump:maxjump), stat=aerr)
             if(aerr /= 0) call memerror("master2bodyjumps 21")
          end if
       else
           if ( storeXXmesjumps) then
              if ( .not.allocated( p2b_me) ) then
                 allocate(p2b_me ( minjump:maxjump) , stat=aerr)
                 if(aerr /= 0) call memerror("master2bodyjumps 22")
              end if
           else
              if ( .not.allocated( p2b_op ) ) then
                 allocate(p2b_op( minjump:maxjump), stat=aerr)
                 if(aerr /= 0) call memerror("master2bodyjumps 23")
              end if
           end if
       end if
     else
        if ( .not.allocated( n2b_isd ) ) then
           allocate(n2b_isd(minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master2bodyjumps 30")
        end if
        if ( .not.allocated( n2b_fsd ) ) then
           allocate(n2b_fsd( minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master2bodyjumps 31")
        end if
        if ( .not.allocated( n2b_phase)) then
           allocate(n2b_phase(minjump:maxjump), stat=aerr)
           if(aerr /= 0) call memerror("master2bodyjumps 32")
        end if

        if(print_jumps)then
              if ( .not.allocated( n2b_isd0 ) ) then
                 allocate(n2b_isd0(minjump:maxjump) , stat=aerr)
                 if(aerr /= 0) call memerror("master2bodyjumps 33")
              end if
              if ( .not.allocated( n2b_fsd0 ) ) then
                 allocate(n2b_fsd0(minjump:maxjump) , stat=aerr)
                 if(aerr /= 0) call memerror("master2bodyjumps 34")
              end if
        endif

        if(threebody)then
           if ( .not.allocated( n2b_dop ) )  then
              allocate(n2b_dop( minjump:maxjump), stat=aerr)
              if(aerr /= 0) call memerror("master2bodyjumps 50")
           end if
           if ( .not.allocated( n2b_cop ) )  then
              allocate(n2b_cop( minjump:maxjump), stat=aerr)
              if(aerr /= 0) call memerror("master2bodyjumps 51")
           end if
        else
           if ( storeXXmesjumps) then
             if ( .not.allocated( n2b_me) ) then
                allocate(n2b_me ( minjump:maxjump) , stat=aerr)
              if(aerr /= 0) call memerror("master2bodyjumps 52")
             end if
           else
             if ( .not.allocated( n2b_op ) )  then
                allocate(n2b_op( minjump:maxjump), stat=aerr)
              if(aerr /= 0) call memerror("master2bodyjumps 53")
             end if
           end if
        endif
     end if
  endif
!.............................................................................
!.... LOOP OVER SECTOR JUMPS.... DIVIDE AMONG MPI PROCS & OpenMP THREADS.....
!      started in 7.5.0

!........ SET UP ARRAY OF SECTOR JUMPS TO LOOP OVER........

  proddescents = 0
  rproddescents = 0.d0
  
  if(fill .and. restrictjumps .and. compactjumpstorage) call surveysectorjumps(it,2)

  if(.not.allocated(mysectjumplist))then 
	  sjumplistsize=0
	  if( np(1)>1)sjumplistsize=x2bjump(1)%nsectjumps
	  if( np(2)>1)sjumplistsize=bmax(sjumplistsize,x2bjump(2)%nsectjumps)
	  allocate(mysectjumplist( sjumplistsize))
  end if
  mysectjumplist(:) = 0	 
  nmysectjumps      = 0
  if(.not.fill)then   ! counting up the jumps
       if(nproc==1)then
		   nmysectjumps = x2bjump(it)%nsectjumps
		   do isjmp = 1,x2bjump(it)%nsectjumps
			   mysectjumplist(isjmp)=isjmp         ! simplest choice
		   end do
	   else
		   if(nproc > x2bjump(it)%nsectjumps)then  ! few sector jumps than MPI processes--unlikely but must account for it
			   
			   if(iproc < x2bjump(it)%nsectjumps)then
				   nmysectjumps=1
				   mysectjumplist(1)=iproc+1
  		       else
				   nmysectjumps=0
			   end if
		   else
			   
		       if(iproc >0)then ! choose simplest assignment
			      nmysectjumps = x2bjump(it)%nsectjumps/nproc
			      do isjmp = 1,nmysectjumps
					   mysectjumplist(isjmp)=isjmp+(iproc-1)*nmysectjumps
			      end do 
		       else             ! put remainder on root proc
			      nmysectjumps = x2bjump(it)%nsectjumps-(x2bjump(it)%nsectjumps/nproc)*(nproc-1)
			      if(nmysectjumps < 1)then
				     print*,iproc,' whoopsie problem with assigning sector jumps for counting jumps ',nmysectjumps
				      stop
			      end if
			      do isjmp = 1,nmysectjumps
					   mysectjumplist(isjmp)=isjmp+x2bjump(it)%nsectjumps-nmysectjumps
			      end do 
		       end if
	       end if		   
	   end if	  
  else               ! creating the jumps needed for this MPI proc
	  if(restrictjumps .and. nproc > 1)then
   	     do isjmp = 1,x2bjump(it)%nsectjumps
			 if (x2bjump(it)%sjmp(isjmp)%nstart+x2bjump(it)%sjmp(isjmp)%njumps < minjump)cycle
			 if(.not. compactjumpstorage .and. x2bjump(it)%sjmp(isjmp)%nstart > maxjump )cycle  
			 if(compactjumpstorage)then
				 startjump = x2bjump(it)%sjmp(isjmp)%nstart
				 if(it==1)then
					 call fetchjumpshift('P2',startjump)
				 else
					 call fetchjumpshift('N2',startjump)
				 end if
				 if(startjump > maxjump)cycle
			 end if
             nmysectjumps=nmysectjumps+1
             mysectjumplist(nmysectjumps)=isjmp
 	     end do			  
		  
	  else   ! count 'em all
		  nmysectjumps=x2bjump(it)%nsectjumps
  	     do isjmp = 1,x2bjump(it)%nsectjumps
		     mysectjumplist(isjmp)=isjmp         ! simplest choice
	     end do		  
	  end if
  end if
!....... NOW LOOP OVER ASSIGNED SECTOR JUMPS...........

  do jj = 1,nmysectjumps
  	 isjmp = mysectjumplist(jj)

!........ DEFAULT CHOICES......
     shiftjump = 0
	 startjump = 	 x2bjump(it)%sjmp(isjmp)%nstart+1
	 endjump   =     x2bjump(it)%sjmp(isjmp)%nstart+x2bjump(it)%sjmp(isjmp)%njumps
	 if(fill .and. restrictjumps .and. compactjumpstorage)then  ! SET LIMITS ON STORAGE, SHIFT IN STORAGE
		 if(.not.x2bjump(it)%sjmp(isjmp)%containsexons)cycle
		 startjump = x2bjump(it)%sjmp(isjmp)%exstart
		 endjump   = x2bjump(it)%sjmp(isjmp)%exstop
		 shiftjump = x2bjump(it)%sjmp(isjmp)%exshift
		 
	 end if
!	 if(fill .and. it==1 .and. iproc==0)print*,' CHECKING.... ',isjmp,nmysectjumps
!	 if(fill .and. restrictjumps .and. .not. compactjumpstorage)then
!		 if (x2bjump(it)%sjmp(isjmp)%nstart< 0)cycle  ! nstart < 0 signals to skip these jumps
!	 end if

!........... PREPARE TO COUNT JUMPS............
     if(.not.fill ) x2bjump(it)%sjmp(isjmp)%njumps = 0
     if(.not.fill )x2bjump(it)%sjmp(isjmp)%nstart = 0 !totjumps
	 
!.......... CYCLE IF NO CONJUGATE JUMPS...........
     if( x2bjump(it)%csjmp(isjmp)%ncjmps ==0)then
         x2bjump(it)%sjmp(isjmp)%njumps = 0
         x2bjump(it)%sjmp(isjmp)%nstart = 0 !totjumps
         cycle
     endif
	 
     is = x2bjump(it)%isector(isjmp)
     fs = x2bjump(it)%fsector(isjmp)
     proddescents = proddescents + descent%sector(is)%totdescendents* descent%sector(fs)%totdescendents
     rproddescents = rproddescents + dsqrt( real(descent%sector(is)%totdescendents* descent%sector(fs)%totdescendents,8))

     call matchmaker(it,nbody,is,descent%sector(is),nbody,fs,descent%sector(fs), .true., ngeneologies)

     if ( ngeneologies == 0 ) cycle
	 
     njumps =   x2bjump(it)%sjmp(isjmp)%nstart
     curgen =>geneology
	 
     do igen = 1,ngeneologies
        call chaingang(it,nbody*2,curgen,optype)   ! advances curgen
        call weld(it,is,fs,nbody,isjmp,.true., fill , x2bjump(it)%sjmp(isjmp),startjump,endjump,shiftjump,njumps)
     end do  ! igen
     njumps = njumps -x2bjump(it)%sjmp(isjmp)%nstart
	 
!	 if(fill .and. restrictjumps .and. compactjumpstorage)then
!		 njumps = njumps - startjump+1
!	 else		 
!         njumps = njumps -x2bjump(it)%sjmp(isjmp)%nstart
!     end if

     if(.not.fill )then
             x2bjump(it)%sjmp(isjmp)%njumps = njumps
             jumpcount(isjmp) = njumps
     end if

  enddo ! isjmps

!.... END OF LOOP OVER SECTOR JUMPS..........................................  
!.............................................................................
  
  if(.not.fill .and. nproc > 1) then
     jumptime= jumptime+  BMPI_WTIME()-jumptimestart  ! MPI clock
  end if

  if(.not.fill.and. nproc > 1 .and. MPIjumps)then   ! broadcast jumps via MPI
     call BMPI_ALLREDUCE(jumpcount,jumpreduce,x2bjump(it)%nsectjumps, & 
          MPI_SUM,icomm,ierr)
     do isjmp = 1,x2bjump(it)%nsectjumps
         x2bjump(it)%sjmp(isjmp)%njumps = jumpreduce(isjmp)  
     end do
  end if 
!......... ADD UP JUMPS..............
  do isjmp = 1,x2bjump(it)%nsectjumps
     if(.not.fill )then
         x2bjump(it)%sjmp(isjmp)%nstart = totjumps
     end if
     totjumps = totjumps + x2bjump(it)%sjmp(isjmp)%njumps

  end do
  if(.not.fill)then 
     if(it == 1)then
        totp2bjumps = totjumps
        if ( iproc == 0 .and. verbose_jumps ) write(6,*)totjumps,' proton jumps' 
     else
        totn2bjumps = totjumps
        if ( iproc == 0 .and. verbose_jumps) write(6,*)totjumps,' neutron jumps' 
     end if
  end if

  deallocate(descent%sector)
  nullify(optype)
  if(.not.fill )deallocate( jumpcount, jumpreduce)

  return
end subroutine master2bodyjumps


!================================================================
! subroutine set1bsectorjumps
!
!  sets a list of 1body sector jumps to be computed
!
!  IMPORTANT: the following restriction on jumps apply:
!     FOR PROTONS  final parity <= initial parity
!                  final Jz     <= initial Jz
!                  final W       <= initial W
!
!     FOR NEUTRONS final parity <= initial parity  if overall parity = +
!                  final parity >= initial parity  if overall parity = -
!                  final Jz     >= initial Jz
!                 (but no restrictions on W)
!  This has consequences for computing densities
!
!  NB: THIS MIGHT BE CHANGED DEPENDING ON HERMITICITY
!  WHEN CREATING JUMPS FOR 1-BODY DENSITY MATRICES, HERMITICITY TURNED OFF;
!  THIS MEANS NO ORDERING ON INITIAL, FINAL SECTORS
!
subroutine set1bsectorjumps(it,create,hermflag,whermflag,changeq )
   use sectors
   use system_parameters
   use jumpNbody
   use verbosity
   use nodeinfo
   implicit none

   integer :: it  ! species
   logical :: create  ! flag
   logical :: hermflag  !  flag to signal to use general hermiticity
   logical :: whermflag  ! flag to signal use hermiticity on W specifically
   logical :: changeq ! flag for allowing change of q#s
   integer :: nsjmps
   integer :: is,fs
   integer :: ijz,ipar,iw
   integer :: fjz,fpar,fw
   integer :: maxdm  
   integer :: aerr

   nsjmps = 0

   call maxdelJz(it,1,maxdm)
   do is = 1,nsectors(it)
     ijz = xsd(it)%sector(is)%jzX
     ipar= xsd(it)%sector(is)%parX
     iw  = xsd(it)%sector(is)%wX
     do fs = 1,nsectors(it)
   !----------- CHECK IF A ONE-BODY OPERATOR CAN EVEN REACH-----------
        if( abs(ijz - xsd(it)%sector(fs)%jzX) > maxdm)cycle

   ! here is where HERMITICITY comes in
        fjz = xsd(it)%sector(fs)%jzX
        fpar= xsd(it)%sector(fs)%parX
        fw  = xsd(it)%sector(fs)%wX
        if(.not.changeq .and. ( ijz /= fjz .or. ipar /= fpar)) cycle
        if(hermflag)then        
          if(fpar > ipar .and. it == 1)cycle   ! for protons
          if(fpar > ipar .and. iparity == 1 .and. it==2)cycle  ! same order as protons
          if(fpar < ipar .and. iparity ==-1 .and. it==2)cycle  ! opposite order as protons
          if(fpar == ipar)then
              if( fjz > ijz .and. it == 1)cycle     ! protons
              if( fjz < ijz .and. it == 2)cycle  ! opposite order from protons to preserve conjugacy

              if(fjz == ijz)then
                 if(fw > iw .and. whermflag)cycle
              endif
          endif
   ! added in 7.4.6 for fragmented sectors
         if(fpar == ipar .and. fjz==ijz .and. iw==fw .and. is < fs)cycle
        endif
        nsjmps = nsjmps +1
        if(create)then
          x1bjump(it)%isector(nsjmps) = is
          x1bjump(it)%fsector(nsjmps) = fs
          if(is == fs)then
            x1bjump(it)%sjmp(nsjmps)%diag = .true.
          else
            x1bjump(it)%sjmp(nsjmps)%diag = .false.
          endif
        endif
     enddo  ! fs
   enddo  ! is

   if(.not.create)then
      x1bjump(it)%nsectjumps = nsjmps
      if(nsjmps ==0)return
      allocate(x1bjump(it)%isector(nsjmps), x1bjump(it)%fsector(nsjmps) , stat=aerr)
      if(aerr /= 0) call memerror("set1bsectorjumps 1")
      allocate(x1bjump(it)%sjmp(nsjmps), stat=aerr)
      if(aerr /= 0) call memerror("set1bsectorjumps 2")

   endif
   if(create .and. print_sectorjumps .and. iproc==0)then
      print*,nsjmps,' total 1-body sector jumps for ',it
      write(sectorjump_file,*)'1-body jumps ',it
      do is = 1,nsjmps
        write(sectorjump_file,*)x1bjump(it)%isector(is), ' -> ', x1bjump(it)%fsector(is)
      
      enddo  ! is
   endif

   return
end subroutine set1bsectorjumps

!================================================================
! subroutine set2bsectorjumps
!
!  sets a list of 2body sector jumps to be computed; newer version
!
!  IMPORTANT: the following restriction on jumps apply:
!     FOR PROTONS  final parity <= initial parity
!                  final Jz     <= initial Jz
!                  final W       <= initial W
!
!     FOR NEUTRONS final parity <= initial parity  if overall parity = +
!                  final parity >= initial parity  if overall parity = -
!                  final Jz     >= initial Jz
!                 (but no restrictions on W)
!  This has consequences for computing densities
!
!  NB: THIS MIGHT BE CHANGED DEPENDING ON HERMITICITY
!
subroutine set2bsectorjumps(it,create,hermflag,whermflag,changeq,changew,ndm )

   use sectors
   use system_parameters
   use jumpNbody
   use verbosity
   use nodeinfo
   implicit none

   integer :: it  ! speciesq
   logical :: create  ! flag
   logical :: hermflag  !  flag to signal to use general hermiticity
   logical :: whermflag  ! flag to signal use hermiticity on W specifically
   logical :: changeq ! flag for allowing change of q#s
   logical :: changew ! flag for allowing change of W--usually TRUE but set FALSE for J2,T2
   integer :: ndm    ! input into maxdelJz
   integer :: nsjmps
   integer :: is,fs
   integer :: ijz,ipar,iw
   integer :: fjz,fpar,fw
   integer :: maxdm
   integer :: aerr

   nsjmps = 0

   call maxdelJz(it,ndm,maxdm)
   do is = 1,nsectors(it)
     ijz = xsd(it)%sector(is)%jzX
     ipar= xsd(it)%sector(is)%parX
     iw  = xsd(it)%sector(is)%wX
     do fs = 1,nsectors(it)
        !----------- CHECK IF A ONE-BODY OPERATOR CAN EVEN REACH-----------
        if( abs(ijz - xsd(it)%sector(fs)%jzX) > maxdm)cycle

        ! here is where HERMITICITY comes in
        fjz = xsd(it)%sector(fs)%jzX
        fpar= xsd(it)%sector(fs)%parX
        fw  = xsd(it)%sector(fs)%wX
		if(.not.changew .and. iw/=fw)cycle
        if(.not.changeq .and. ( ijz /= fjz .or. ipar /= fpar)) cycle
        if(.not.changeq .and. whermflag .and. iw < fw)cycle
        if(hermflag)then        
          if(fpar > ipar .and. it == 1)cycle   ! for protons
          if(fpar > ipar .and. iparity == 1 .and. it==2)cycle  ! same order as protons
          if(fpar < ipar .and. iparity ==-1 .and. it==2)cycle  ! opposite order as protons
          if(fpar == ipar)then
              if( fjz > ijz .and. it == 1)cycle     ! protons
              if( fjz < ijz .and. it == 2)cycle  ! opposite order from protons to preserve conjugacy

              if(fjz == ijz)then
                 if(fw > iw .and. whermflag)cycle
              endif
          endif
          ! added in 7.4.6 for fragmented sectors
          if(fpar == ipar .and. fjz==ijz .and. iw==fw .and. is < fs)cycle
        endif
        nsjmps = nsjmps +1
        if(create)then
          x2bjump(it)%isector(nsjmps) = is
          x2bjump(it)%fsector(nsjmps) = fs
          if(is == fs)then
            x2bjump(it)%sjmp(nsjmps)%diag = .true.
          else
            x2bjump(it)%sjmp(nsjmps)%diag = .false.
          endif
        endif
     enddo  ! fs
   enddo  ! is

   if(.not.create)then
      x2bjump(it)%nsectjumps = nsjmps
      if(nsjmps ==0)return
      allocate(x2bjump(it)%isector(nsjmps), x2bjump(it)%fsector(nsjmps) , stat=aerr)
      if(aerr /= 0) call memerror("set2bsectorjumps 1")
      allocate(x2bjump(it)%sjmp(nsjmps), stat=aerr)
      if(aerr /= 0) call memerror("set2bsectorjumps 2")
   endif

   if(create .and. print_sectorjumps .and. iproc==0)then
      print*,nsjmps,' total 2-body sector jumps for ',it
      write(sectorjump_file,*)'2-body jumps ',it
      do is = 1,nsjmps
        write(sectorjump_file,*)x2bjump(it)%isector(is), ' -> ', x2bjump(it)%fsector(is)
      enddo  ! is
   endif
   return
end subroutine set2bsectorjumps
!================================================================

!
!  finds conjugate 1-body sector jumps
!
   subroutine masterfindconjumps1b(it,changeq)

   use jumpNbody

   implicit none
   integer it
   logical changeq

   if(changeq)then
      call findconjugatejumps(it,.false.,x1bjump)
      call findconjugatejumps(it,.true.,x1bjump)
   else
      call findconjugatesemidiag(it,.false.,x1bjump)
      call findconjugatesemidiag(it,.true.,x1bjump)
   endif
   return
   end  subroutine masterfindconjumps1b
!====================================================================
!  COMMENT -- WELL, this actually isn't correct
!
   subroutine masterfindconjumps2b(it)

   use jumpNbody

   implicit none
   integer it

   call findconjugatesemidiag(it,.false.,x2bjump)
   call findconjugatesemidiag(it,.true.,x2bjump)

   return
   end  subroutine masterfindconjumps2b
!====================================================================
!
!  subroutine findconjugatejumps
!  for now, this applies only to 1-body jumps
!  in future, modify to apply to couple 1 and 2-body jumps 
!  (e.g., pp 2-body jumps with n 1-body jumps for a ppn 3-body force)
!
   subroutine findconjugatejumps(it,create,xNjumps)

   use jumpdef
   use sectors
   use bsector_mod
   implicit none
   integer :: it
   logical :: create
   type (jumpsect) :: xNjumps(2)

   integer :: itc
   integer :: isjmp,csjmp
   integer :: ixs,fxs,ics,fcs
   integer :: inc,fnc
   logical :: iok,fok
   integer :: nc
   integer :: aerr
!------------- ERROR TRAP VARIABLES --------------
   integer :: dMx,dMy,dparx,dpary,dWx,dWy   
   
   itc = 3 - it

   if(.not.create) then
      allocate(xNjumps(it)%csjmp( xNjumps(it)%nsectjumps), stat=aerr)
      if(aerr /= 0) call memerror("findconjugatejumps 1")
   end if

   do isjmp = 1,xNjumps(it)%nsectjumps    ! loop over sector jumps
     ixs = xNjumps(it)%isector(isjmp)     ! initial sector
     fxs = xNjumps(it)%fsector(isjmp)     ! final sector
     call get_delta_q(it,ixs,fxs,dMx,dparx,dWx)   ! FOR ERROR TRAP
	 
!     if(create)print*,it,isjmp,' : ',ixs,fxs,xNjumps(it)%csjmp(isjmp)%ncjmps
     nc = 0
     do csjmp = 1, xNjumps(itc)%nsectjumps     ! loop over all the sector jumps for conjugate species
       ics = xNjumps(itc)%isector(csjmp)       ! trial initial conjugate sector         
       fcs = xNjumps(itc)%fsector(csjmp)       ! trial initial final sector

       if(ixs == fxs .and. fcs > ics) cycle   ! ENFORCEMENT OF HERMITICITY

!........ NOW CHECK TO SEE IF ics CAN BE A CONJUGATE SECTOR TO ixs 
!         AND IF fcs CAN BE A CONJUGATE SECTOR TO fxs

       iok = .false.   ! FIRST CHECK IF INITIAL SECTOR CONJUGATE
       do inc = 1,xsd(it)%sector(ixs)%ncsectors            
         if(ics == xsd(it)%sector(ixs)%csector(inc))then
           iok = .true.
           exit
         endif
       enddo  ! inc

       fok = .false.    ! CHECK IF FINAL SECTORS CONJUGATE
       do fnc = 1,xsd(it)%sector(fxs)%ncsectors
         if(fcs == xsd(it)%sector(fxs)%csector(fnc))then
           fok = .true.
           exit
         endif
       enddo  ! fnc
       if(iok .and. fok) then
         nc = nc + 1
         if(create)then
           xNjumps(it)%csjmp(isjmp)%cjump( nc ) = csjmp
!.................... TESTING.......................................
           call get_delta_q(itc,ics,fcs,dMy,dpary,dWy)
		   if(dMx/=-dMy .or. dparx/=dpary)then
			   print*,' ??? mismatch in quantum numbers '
			   print*,it,dmx,dmy,dparx,dpary
			   stop
		   end if
         endif
       endif
     enddo  ! csjmp
     if(.not.create)then
       xNjumps(it)%csjmp(isjmp)%ncjmps = nc
       allocate( xNjumps(it)%csjmp(isjmp)%cjump( nc ), stat=aerr)
       if(aerr /= 0) call memerror("findconjugatejumps 2")
     endif
   enddo !isjmp
   return
   end subroutine findconjugatejumps

!==========================================================================
!
!  subroutine findconjugatejumpsALT
!  modified to apply to couple 1 and 2-body jumps 
!  (e.g., pp 2-body jumps with n 1-body jumps for a ppn 3-body force)
!  
!
   subroutine findconjugatejumpsALT(it,create,xNjumps,xCjumps)

   use jumpdef
   use sectors
   implicit none
   integer :: it
   logical :: create
   type (jumpsect) :: xNjumps(2),xCjumps(2)

   integer :: itc
   integer :: isjmp,csjmp
   integer :: ixs,fxs,ics,fcs
   integer :: inc,fnc
   logical :: iok,fok
   integer :: nc
   integer :: aerr
   itc = 3 - it

   if(.not.create) then
      allocate(xNjumps(it)%csjmp( xNjumps(it)%nsectjumps), stat=aerr)
      if(aerr /= 0) call memerror("findconjugatejumpsALT 1")
   end if


   do isjmp = 1,xNjumps(it)%nsectjumps    ! loop over sector jumps
     ixs = xNjumps(it)%isector(isjmp)     ! initial sector
     fxs = xNjumps(it)%fsector(isjmp)     ! final sector
!     if(create)print*,it,isjmp,' : ',ixs,fxs,xNjumps(it)%csjmp(isjmp)%ncjmps
     nc = 0
     do csjmp = 1, xCjumps(itc)%nsectjumps     ! loop over all the sector jumps for conjugate species
       ics = xCjumps(itc)%isector(csjmp)       ! trial initial conjugate sector         
       fcs = xCjumps(itc)%fsector(csjmp)       ! trial initial final sector

       if(ixs == fxs .and. fcs > ics) cycle   ! ENFORCEMENT OF HERMITICITY

!........ NOW CHECK TO SEE IF ics CAN BE A CONJUGATE SECTOR TO ixs 
!         AND IF fcs CAN BE A CONJUGATE SECTOR TO fxs

       iok = .false.   ! FIRST CHECK IF INITIAL SECTOR CONJUGATE

       do inc = 1,xsd(it)%sector(ixs)%ncsectors            
         if(ics == xsd(it)%sector(ixs)%csector(inc))then
           iok = .true.
           exit
         endif
       enddo  ! inc
       fok = .false.    ! CHECK IF FINAL SECTORS CONJUGATE
       do fnc = 1,xsd(it)%sector(fxs)%ncsectors
         if(fcs == xsd(it)%sector(fxs)%csector(fnc))then
           fok = .true.
           exit
         endif
       enddo  ! fnc
       if(iok .and. fok) then
         nc = nc + 1
         if(create)then
           xNjumps(it)%csjmp(isjmp)%cjump( nc ) = csjmp

         endif
       endif
     enddo  ! csjmp
     if(.not.create)then
       xNjumps(it)%csjmp(isjmp)%ncjmps = nc
       allocate( xNjumps(it)%csjmp(isjmp)%cjump( nc ), stat=aerr)
       if(aerr /= 0) call memerror("findconjugatejumpsALT 2")
     endif
   enddo !isjmp
   return
   end subroutine findconjugatejumpsALT

!==========================================================================
!
!  subroutine findconjugatediag
!  for now, this applies only to 2-body (pp and nn) jumps
!  in future, modify to apply to 3-body (e.g. ppp and nnn) jumps 
!
   subroutine findconjugatesemidiag(it,create,xNjumps)

   use jumpdef
   use sectors
   implicit none
   integer :: it
   logical :: create
   type (jumpsect) :: xNjumps(2)

   integer :: itc
   integer :: isjmp,csjmp
   integer :: ixs,fxs,ics,fcs
   integer :: inc,fnc
   logical :: iok,fok
   integer :: nc
   integer :: aerr

   itc = 3 - it

   if(.not.create) then
      allocate(xNjumps(it)%csjmp( xNjumps(it)%nsectjumps), stat=aerr)
      if(aerr /= 0) call memerror("findconjugatesemidiag 1")
   end if

   do isjmp = 1,xNjumps(it)%nsectjumps
     ixs = xNjumps(it)%isector(isjmp)
     fxs = xNjumps(it)%fsector(isjmp)
     nc = 0
!------------- LOOP OVER CONJUGATE SECTORS
     do inc = 1,xsd(it)%sector(ixs)%ncsectors
       ics = xsd(it)%sector(ixs)%csector(inc)
       iok = .false.
       do fnc = 1,xsd(it)%sector(fxs)%ncsectors
         fcs = xsd(it)%sector(fxs)%csector(fnc)
         if(ics == fcs)then
           iok = .true.
           exit
         endif
       enddo  ! fnc
       if(iok)then
         nc = nc + 1
         if(create)then
           xNjumps(it)%csjmp(isjmp)%cjump( nc ) = ics
         endif
       endif
     enddo  ! inc

     if(.not.create)then
       xNjumps(it)%csjmp(isjmp)%ncjmps = nc
       allocate( xNjumps(it)%csjmp(isjmp)%cjump( nc ), stat=aerr)
       if(aerr /= 0) call memerror("findconjugatesemidiag 2")
     endif
   enddo !isjmp
   return
   end subroutine findconjugatesemidiag


!===================================================================
!
!  write information on 1-body jumps needed to model distribution of workload on parallel nodes
!
subroutine write1bjmps4modeling(it,writeitout,nme)

 use sectors
 use jumpNbody
 use verbosity
 implicit none
 integer it,itc

 integer is,fs
 integer isjmp
 integer ic
 integer cs
 integer(kind=8) :: nme        ! equivalent # of matrix elements
 logical writeitout
 integer hfactor  ! to account for hermiticity


! if(.not.print4modelinfo)return
 itc = 3 -it
 nme = 0
!------------- WRITE OUT 1-BODY JUMP INFO ------------

  write(modelinfo_file,*)it,1,x1bjump(it)%nsectjumps

do isjmp = 1,x1bjump(it)%nsectjumps
  is = x1bjump(it)%isector(isjmp)
  fs = x1bjump(it)%fsector(isjmp)
  if(writeitout)then
  write(modelinfo_file,101)is,fs,x1bjump(it)%sjmp(isjmp)%njumps, x1bjump(it)%csjmp(isjmp)%ncjmps

    write(modelinfo_file,101)(x1bjump(it)%csjmp(isjmp)%cjump(ic), ic = 1, x1bjump(it)%csjmp(isjmp)%ncjmps)
   endif
  if(it == 1)then
  do ic = 1,x1bjump(it)%csjmp(isjmp)%ncjmps
     cs = x1bjump(it)%csjmp(isjmp)%cjump(ic)
     if (is == fs .and. x1bjump(it)%sjmp(isjmp)%diag .and. & 
                     x1bjump(itc)%sjmp(cs)%diag)then
         hfactor = 1
     else
         hfactor = 2
     endif
     nme = nme + x1bjump(it)%sjmp(isjmp)%njumps * x1bjump(itc)%sjmp(cs)%njumps*hfactor

  enddo
  end if

101 format(10i7)
enddo  !isjmp

if(it ==1 .and. writeitout)then
  print*,' Hpn creates ',nme,' matrix elements '
endif

return
end subroutine write1bjmps4modeling
!==================================================================================

!
!  write information on 2-body jumps needed to model distribution of workload on parallel nodes
!
subroutine write2bjmps4modeling(it,writeitout,nme)

 use sectors
 use jumpNbody
 use verbosity
 use system_parameters
 implicit none
 integer it
 integer itc
 integer is,fs
 integer isjmp
 integer ic
 integer cs

 integer(kind=8) :: nme,nlocal

 logical writeitout

 nme = 0
 if(np(it) ==0)return

!------------- WRITE OUT 2-BODY JUMP INFO ------------
  itc = 3 -it
  if(writeitout)write(modelinfo_file,*)' species ',it,1,x2bjump(it)%nsectjumps

do isjmp = 1,x2bjump(it)%nsectjumps
  is = x2bjump(it)%isector(isjmp)
  fs = x2bjump(it)%fsector(isjmp)
  if(writeitout)then
	  write(modelinfo_file,*)' conjugate jump list '
  write(modelinfo_file,101)is,fs,x2bjump(it)%sjmp(isjmp)%njumps, x2bjump(it)%csjmp(isjmp)%ncjmps
!---------- THIS PART NEEDS TO BE FIXED
    write(modelinfo_file,101)(x2bjump(it)%csjmp(isjmp)%cjump(ic), ic = 1, x2bjump(it)%csjmp(isjmp)%ncjmps)
101 format(10i7)
  endif
!  write(modelinfo_file,*)x2bjump(it)%csjmp(isjmp)%ncjmps
    nlocal = 0
    do ic = 1,x2bjump(it)%csjmp(isjmp)%ncjmps
       cs = x2bjump(it)%csjmp(isjmp)%cjump(ic)
       nme = nme + x2bjump(it)%sjmp(isjmp)%njumps*xsd(itc)%sector(cs)%nxsd 
       nlocal = nlocal + x2bjump(it)%sjmp(isjmp)%njumps*xsd(itc)%sector(cs)%nxsd 
       if(writeitout)write(modelinfo_file,*)ic,cs, x2bjump(it)%sjmp(isjmp)%njumps,xsd(itc)%sector(cs)%nxsd 
    enddo
   if(writeitout) write(modelinfo_file,*)nlocal,' matrix elements for this sector jump '
enddo  !isjmp

nme = nme*2   ! need to account for Hermiticity back and forth

if(writeitout)then
if(it == 1)then
  print*,' Hpp has ',nme,' matrix elements '
else
  print*,' Hnn has ',nme,' matrix elements '
endif
endif

return
end subroutine write2bjmps4modeling
!=========================================================
!
! Find max change in Jz for an n-body operator
! rather kludgey but no need to make it slick
! not going to fuss about parity, either
!

subroutine maxdelJz(it,nbody,delM)

   use haiku_info
   use spstate
   implicit none

   integer :: it
   integer :: nbody
   integer :: delM

   integer :: list(nbody)  ! list of operators already used
   integer :: ibody

   integer :: Mmax, lastm
   integer :: m,i

   delM =0
   lastm = -1
   do ibody = 1,nbody
      mmax = -1
      do i = 1,nhsps(it)
         m = hspsqn(it,i)%m
         if(m <= mmax) cycle
         !----------- CHECK THIS HAS NOT ALREADY BEEN USED ----------
         if( ibody > 1)then
            if( m == lastm .and. i == list(ibody-1))cycle
         endif
         mmax = m 
         list(ibody) = i
         lastm = m
      enddo
      lastm = mmax
      delM = delM + mmax
   enddo ! ibody

   delM = 2*delM
   return
end subroutine maxdelJz
!================================================================
subroutine setlinkid
   use chains
   implicit none

   currentID = 0
return

end subroutine setlinkid

!==================================================================


 subroutine countchainbase(hchain,ncount)

 use chains
 implicit none
 type (chainbase), pointer :: hchain
 integer, intent(inout) :: ncount
 logical finished
 ncount = ncount + 1
 finished = .false.
 do while(.not.finished)
    ncount = ncount + 1
    if(associated( hchain%next))then
      hchain => hchain%next
    else
       finished = .true.
    endif
 end do

 return
 
 end subroutine countchainbase
!==================================================================

!  routines for creating 3-body jumps
!
!  initiated 9/09 by CWJ @ SDSU
!

!====================================================================
!
!  NOTE: master3bodyjumps SHOULD BE CALLED BEFORE any other
!
subroutine master3bodyjumps(it,fill,changeq)

  use system_parameters
  use basis
  use verbosity
  use sectors
  use descendents
  use geneologies
  use chains
  use jumpNbody
  use jump3body
  use nodeinfo
  use jumplimits
  use flagger
  use bmpi_mod
  use welder
  implicit none

!  include 'binterfaces.inc'
  integer :: ierr

  integer :: it  ! species
  type (groupjump) :: descent
  type (genlink), pointer :: curgen

  integer(8) :: ngeneologies,igen
  integer :: is,fs
  integer :: nbody
  character, pointer :: optype(:)*1
  !   type (handchain), allocatable :: mrchain(:)
  integer :: isjmp
  integer(8):: njumps,nlinks
  integer(8) :: totjumps
  logical :: fill
  logical :: changeq  ! flag to allow changing q#s


  integer(8), allocatable :: jumpcount(:),jumpreduce(:)   
                                      ! used to pass jump information in MPI
  integer(8):: minjump,maxjump,startjump,endjump,shiftjump
  integer :: aerr
!--------------------------------- SET UP DESCENDENTS----------------
!       NOTE FOR FUTURE PARALLELIZATION 
!       CAN DISTRIBUTE -- ONLY GET DESCENDENTS FOR RELEVANT NODES
!

  if(np(it) < 3)return
  minjump = 1
  maxjump = 1
  nbody = 3
  if(associated(optype))nullify(optype)
  allocate(optype(3*nbody), stat=aerr)
  if(aerr /= 0) call memerror("master3bodyjumps 1")

  optype(1) = 'd'
  optype(2) = 'd'
  optype(3) = 'd'
  optype(4) = 'c'
  optype(5) = 'c'
  optype(6) = 'c'
!-------------------- FIGURE OUT WHICH SECTORS WE WANT TO GO TO
  if(.not.fill)then
     call set3bsectorjumps(it,.false.,changeq )
     call set3bsectorjumps(it,.true.,changeq )
  endif
  if(.not.fill)then
      allocate( jumpcount( x3bjump(it)%nsectjumps ), jumpreduce(x3bjump(it)%nsectjumps), stat=aerr )
      if(aerr /= 0) call memerror("master3bodyjumps 2")
      jumpcount(:) = 0
      jumpreduce(:) = 0
  end if
  descent%nsectors = nsectors(it)
  allocate(descent%sector(nsectors(it)), stat=aerr)
  if(aerr /= 0) call memerror("master3bodyjumps 3")

  do is = 1,nsectors(it)
     call sectordescent(it,is,nbody,descent%sector(is))
  enddo ! is
  totjumps = 0
  if(fill)then
     if(restrictjumps .and. nproc > 1 )then
        if(it ==1)then
           if(.not.makepppjumps(iproc))return
           minjump = startpppjumps(iproc)
           maxjump =  stoppppjumps(iproc)
        end if
        if( it ==2)then
           if(.not.makennnjumps(iproc))return
           minjump = startnnnjumps(iproc)
           maxjump =  stopnnnjumps(iproc)
        end if
     else
        minjump = 1
        if(it ==1)then
          maxjump = totp3bjumps
        end if
        if(it ==2)then
         maxjump = totn3bjumps
        end if

     end if


     if ( it == 1 ) then
        if ( .not.allocated( p3b_isd ) ) then
           allocate(p3b_isd( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master3bodyjumps 1")
        end if
        if ( .not.allocated( p3b_fsd ) ) then
           allocate(p3b_fsd( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master3bodyjumps 2")
        end if
        if ( .not.allocated( p3b_op ) )  then
           allocate(p3b_op( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master3bodyjumps 3")
        end if
        if ( .not.allocated( p3b_phase ) )  then
           allocate(p3b_phase( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master3bodyjumps 4")
        end if
     else
        if ( .not.allocated( n3b_isd ) ) then
           allocate(n3b_isd( minjump:maxjump  ), stat=aerr)
           if(aerr /= 0) call memerror("master3bodyjumps 5")
        end if
        if ( .not.allocated( n3b_fsd ) )  then 
           allocate(n3b_fsd( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master3bodyjumps 6")
        end if
        if ( .not.allocated( n3b_op ) )  then
           allocate(n3b_op( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master3bodyjumps 7")
        end if
        if ( .not.allocated( n3b_phase)) then
           allocate(n3b_phase( minjump:maxjump ), stat=aerr)
           if(aerr /= 0) call memerror("master3bodyjumps 8")
        end if
     end if

  endif
  if(fill .and. restrictjumps .and. compactjumpstorage) call surveysectorjumps(it,3)
  do isjmp = 1,x3bjump(it)%nsectjumps
     if( .not.fill .and. nproc > 1 .and. mod(isjmp,nproc) /= iproc .and. MPIjumps)cycle
	 
!........ DEFAULT CHOICES......
	  shiftjump = 0
	  startjump = 	 x3bjump(it)%sjmp(isjmp)%nstart+1
	  endjump   =     x3bjump(it)%sjmp(isjmp)%nstart+x3bjump(it)%sjmp(isjmp)%njumps
	  if(fill .and. restrictjumps .and. compactjumpstorage)then  ! SET LIMITS ON STORAGE, SHIFT IN STORAGE
	 		 if(.not.x3bjump(it)%sjmp(isjmp)%containsexons)cycle
	 		 startjump = x3bjump(it)%sjmp(isjmp)%exstart
	 		 endjump   = x3bjump(it)%sjmp(isjmp)%exstop
	 		 shiftjump = x3bjump(it)%sjmp(isjmp)%exshift
		 
	 end if
     if(fill .and. restrictjumps .and. x3bjump(it)%sjmp(isjmp)%nstart+x3bjump(it)%sjmp(isjmp)%njumps -shiftjump < minjump )cycle  	 
	 if( fill .and. restrictjumps .and. x3bjump(it)%sjmp(isjmp)%nstart-shiftjump > maxjump )cycle  
	 	 
!.......... SKIP IF NO CONJUGATE JUMPS...........
	 if(.not.fill ) x3bjump(it)%sjmp(isjmp)%njumps = 0

!  NO 3-BODY CONJUGATE JUMPS UNTIL 4-BODY IMPLEMENTED
!	 if( x3bjump(it)%csjmp(isjmp)%ncjmps ==0) then
!	          x3bjump(it)%sjmp(isjmp)%njumps = 0
!	          x3bjump(it)%sjmp(isjmp)%nstart = 0 ! totjumps
!	          cycle
!	 endif	 

     is = x3bjump(it)%isector(isjmp)
     fs = x3bjump(it)%fsector(isjmp)
     if(.not.fill)then
        x3bjump(it)%sjmp(isjmp)%nstart = 0
        x3bjump(it)%sjmp(isjmp)%njumps = 0
     end if

     call matchmaker(it,nbody,is,descent%sector(is),nbody,fs,descent%sector(fs), .true., ngeneologies)

     if ( ngeneologies == 0 ) cycle

     if( print_genes ) then
        write(gene_file,*)ngeneologies,' 3-body GENEOLOGIES for ',is,' -> ',fs
     end if
     curgen =>geneology
     njumps =   x3bjump(it)%sjmp(isjmp)%nstart
     do igen = 1,ngeneologies
        call chaingang(it,nbody*2,curgen,optype)
         call weld3b(it,is,fs,nbody,isjmp,.true., fill,x3bjump(it)%sjmp(isjmp),startjump,endjump,shiftjump, njumps)
     end do
     if(verbose_3bodyjumps .and. iproc== 0)write(6,*)isjmp,is,fs,' there are ',njumps , ' 3b jumps ',ngeneologies
     njumps = njumps -x3bjump(it)%sjmp(isjmp)%nstart

     if(.not.fill)then
             x3bjump(it)%sjmp(isjmp)%njumps = njumps
             jumpcount(isjmp) = njumps
     end if

  enddo ! isjmps

  if(.not.fill .and. nproc > 1 .and. MPIjumps)then   ! broadcast jumps via MPI
     call BMPI_ALLREDUCE(jumpcount,jumpreduce,x3bjump(it)%nsectjumps, & 
          MPI_SUM,icomm,ierr)
     call BMPI_BARRIER(icomm,ierr)
     do isjmp = 1,x3bjump(it)%nsectjumps
         x3bjump(it)%sjmp(isjmp)%njumps = jumpreduce(isjmp)  
     end do
  end if 

!......... ADD UP JUMPS..............
  do isjmp = 1,x3bjump(it)%nsectjumps
     if(.not.fill)then
         x3bjump(it)%sjmp(isjmp)%nstart = totjumps
     end if
     totjumps = totjumps + x3bjump(it)%sjmp(isjmp)%njumps

  end do

  if(.not.fill)then 
     if(it == 1)then
        totp3bjumps = totjumps
        if ( iproc == 0 ) write(6,*)totjumps,' proton jumps' 
     else
        totn3bjumps = totjumps
        if ( iproc == 0 ) write(6,*)totjumps,' neutron jumps' 

     end if
  end if

  deallocate(descent%sector)
  nullify(optype)
  if(.not.fill)deallocate( jumpcount, jumpreduce)

  return
end subroutine master3bodyjumps
!========================================================
! subroutine set3bsectorjumps
!
!  sets a list of 3body sector jumps to be computed
!
!  THIS MIGHT BE CHANGED DEPENDING ON HERMITICITY
!
subroutine set3bsectorjumps(it,create,changeq )

use sectors
use jump3body
implicit none

integer it  ! species
logical create  ! flag
logical changeq  ! flag for allowing change of q#s
integer nsjmps
integer is,fs
integer ijz,ipar
integer maxdm
integer :: aerr

nsjmps = 0

call maxdelJz(it,3,maxdm)
do is = 1,nsectors(it)
  ijz = xsd(it)%sector(is)%jzX
  ipar= xsd(it)%sector(is)%parX
  do fs = 1,nsectors(it)
!----------- CHECK IF A ONE-BODY OPERATOR CAN EVEN REACH-----------
!     if( abs(ijz - xsd(it)%sector(fs)%jzX) > maxdm)cycle
! here is where HERMITICITY comes in
     if(fs > is)cycle
!--------- ALLOW ONLY SAME JZ, PARITY, but different W is okay
     if( .not. changeq  .and. & 
        (ijz /= xsd(it)%sector(fs)%jzX .or. ipar /= xsd(it)%sector(fs)%parX)) cycle
     nsjmps = nsjmps +1
     if(create)then
       x3bjump(it)%isector(nsjmps) = is
       x3bjump(it)%fsector(nsjmps) = fs

     endif

  enddo  ! fs
enddo  ! is

if(.not.create)then
   x3bjump(it)%nsectjumps = nsjmps
   if(nsjmps ==0)return
   allocate(x3bjump(it)%isector(nsjmps), x3bjump(it)%fsector(nsjmps), stat=aerr )
   if(aerr /= 0) call memerror("set3bsectorjumps 1")
   allocate(x3bjump(it)%sjmp(nsjmps), stat=aerr)
   if(aerr /= 0) call memerror("set3bsectorjumps 2")
endif



return
end subroutine set3bsectorjumps

!  COMMENT -- WELL, this actually isn't correct
!
   subroutine masterfindconjumps3b(it)

   use jump3body

   implicit none
   integer it

   call findconjugatesemidiag(it,.false.,x3bjump)
   call findconjugatesemidiag(it,.true.,x3bjump)

   return
   end  subroutine masterfindconjumps3b
!====================================================================
!
!  routines to orchestrate finding conjugate 1+2-body sector jumps for PPN/PNN
!  and "trimming" unused sector jumps
!
!  5/2010 CWJ
!
   subroutine masterXXYconjugatejumps(it)

  
   use jumpNbody
   use verbosity
   use nodeinfo
   implicit none
   integer it,itc   ! species labels
   integer nnull    ! # of null sector jumps (no conjugate jumps)
   integer is

   itc = 3 -it
   call findconjugatejumpsALT(it,.false.,x2bjump,x1bjump)
   call findconjugatejumpsALT(it,.true.,x2bjump,x1bjump)
   call findconjugatejumpsALT(itc,.false.,x1bjump,x2bjump)
   call findconjugatejumpsALT(itc,.true.,x1bjump,x2bjump)

!.......... COUNT UP NULL SECTOR JUMPS (no conjugate jumps)..........
  if(iproc ==0 )then
  if(it == 1)then
!     print*,' For PPN '
  else
!     print*,' For PNN '
  endif
  end if

  nnull = 0
  do is = 1,x2bjump(it)%nsectjumps
     if ( x2bjump(it)%csjmp(is)%ncjmps == 0)then
        nnull = nnull +1
        if(print_sectorjumps)then
           write(sectorjump_file,*)it,' null 2-body ',x2bjump(it)%isector(is), ' -> ', x2bjump(it)%fsector(is)
        endif
     endif
  end do ! is
!  if(iproc == 0)print*,nnull,' null 2-body sector jumps (out of ',x2bjump(it)%nsectjumps,')'

  nnull = 0
  do is = 1,x1bjump(itc)%nsectjumps
     if ( x1bjump(itc)%csjmp(is)%ncjmps == 0)then
        nnull = nnull +1
        if(print_sectorjumps)then
           write(sectorjump_file,*)it,' null 1-body ',x1bjump(itc)%isector(is), ' -> ', x1bjump(itc)%fsector(is)
        endif
     endif
  end do ! is
!  if(iproc == 0)print*,nnull,' null 1-body sector jumps (out of ',x1bjump(itc)%nsectjumps,')'

  return
  end  subroutine masterXXYconjugatejumps


!====================================================================
!===================================================================
!
!  write information on 1-body jumps needed to model distribution of workload on parallel nodes
!
subroutine writeXXYjmps4modeling(itx,writeitout,nme)

 use sectors
 use jumpNbody
 use verbosity
 implicit none
 integer itx,ity

 integer is,fs
 integer isjmp
 integer ic
 integer cs
 integer(kind=8) :: nme        ! equivalent # of matrix elements
 logical writeitout


! if(.not.print4modelinfo)return
 ity = 3 -itx
 nme = 0
!------------- WRITE OUT 1-BODY JUMP INFO ------------

  write(modelinfo_file,*)itx,1,x2bjump(itx)%nsectjumps

do isjmp = 1,x2bjump(itx)%nsectjumps
  is = x2bjump(itx)%isector(isjmp)
  fs = x2bjump(itx)%fsector(isjmp)
  if(writeitout)then
  write(modelinfo_file,101)is,fs,x2bjump(itx)%sjmp(isjmp)%njumps, x1bjump(itx)%csjmp(isjmp)%ncjmps

    write(modelinfo_file,101)(x2bjump(itx)%csjmp(isjmp)%cjump(ic), ic = 1, x2bjump(itx)%csjmp(isjmp)%ncjmps)
   endif

  do ic = 1,x2bjump(itx)%csjmp(isjmp)%ncjmps
     cs = x2bjump(itx)%csjmp(isjmp)%cjump(ic)
     nme = nme + x2bjump(itx)%sjmp(isjmp)%njumps * x1bjump(ity)%sjmp(cs)%njumps

  enddo

101 format(10i7)
enddo  !isjmp
nme = nme*2  ! account for hermiticity

if(writeitout)then
if(itx ==1)then
  print*,' Hppn creates ',nme,' matrix elements '
else
  print*,' Hpnn creates ',nme,' matrix elements '
endif
endif

return
end subroutine writeXXYjmps4modeling
!===================================================================

!
!  write information on 2-body jumps needed to model distribution of workload on parallel nodes
!
subroutine writeXXXjmps4modeling(it,writeitout,nme)

 use sectors
 use jump3body
 use verbosity
 use system_parameters
 implicit none
 integer it
 integer itc
 integer is,fs
 integer isjmp
 integer ic
 integer cs

 integer(kind=8) :: nme

 logical writeitout

! if(.not.print4modelinfo)return
 if(np(it) ==0)return

  nme = 0
!------------- WRITE OUT 1-BODY JUMP INFO ------------
  itc = 3 -it
  write(modelinfo_file,*)it,1,x3bjump(it)%nsectjumps

do isjmp = 1,x3bjump(it)%nsectjumps
  is = x3bjump(it)%isector(isjmp)
  fs = x3bjump(it)%fsector(isjmp)
  if(writeitout)then
  write(modelinfo_file,101)is,fs,x3bjump(it)%sjmp(isjmp)%njumps, x3bjump(it)%csjmp(isjmp)%ncjmps
!---------- THIS PART NEEDS TO BE FIXED
    write(modelinfo_file,101)(x3bjump(it)%csjmp(isjmp)%cjump(ic), ic = 1, x3bjump(it)%csjmp(isjmp)%ncjmps)
101 format(10i7)
  endif
    do ic = 1,x3bjump(it)%csjmp(isjmp)%ncjmps
       cs = x3bjump(it)%csjmp(isjmp)%cjump(ic)
       nme = nme + x3bjump(it)%sjmp(isjmp)%njumps*xsd(itc)%sector(cs)%nxsd 
    enddo

enddo  !isjmp
nme = nme*2  ! account for hermiticity
if(writeitout)then
if(it == 1)then
  print*,' Hppp has ',nme,' matrix elements '
else
  print*,' Hnnn has ',nme,' matrix elements '
endif
endif

return
end subroutine writeXXXjmps4modeling
!=================================================================
!===========================================================================
!
!  sometimes it is necessary to reset all the jump information, for example
!  if we have used 3-body (or the skipzeroes options)
!  in order to use the observables J (and possibly T)
!
!  added Dec 2010 by CWJ @ SDSU
!

subroutine reset_jumps4obs

use system_parameters
use flags3body
use jumpNbody
use jump3body
use opbundles
!use distrinfo
use jumplimits
use nodeinfo
use bmpi_mod
use para_main_mod
use operation_stats
!use annexations
implicit none
logical change1bodyqs,change2bodyqs,hermflag,whermflagp,whermflagn
integer ierr

  change1bodyqs = .true.     
  change2bodyqs = .false.     
  hermflag = .true.
  whermflagp = .true.   ! for protons
  whermflagn = .false.  ! for neutrons -- ALWAYS TURN OFF
  if(iproc==0 .and. compactjumpstorage)then
      print*,' Turning off intron deletion in jump storage '
  end if
  compactjumpstorage = .false.   ! don't use introns...??
!  call set_unannexed
  
!----- CONDITIONS FOR RESETTING ----
if(.not.threebody .and. .not.skipzeros)return
 
  if(iproc==0)print*,' RESETTING '
     threebody = .false.   ! reset
     if( allocated(n2b_op))deallocate(n2b_op)
     if( allocated(p2b_op))deallocate(p2b_op)
     if( allocated(n2b_isd))deallocate(n2b_isd)
     if( allocated(n2b_fsd))deallocate(n2b_fsd)
     if( allocated(n2b_phase))deallocate(n2b_phase)
     if( allocated(p2b_isd))deallocate(p2b_isd)
     if( allocated(p2b_fsd))deallocate(p2b_fsd)
     if( allocated(p2b_phase))deallocate(p2b_phase)
     if( allocated(p2b_isd0))deallocate(p2b_isd0)
     if( allocated(n2b_isd0))deallocate(n2b_isd0)
     if( allocated(p2b_fsd0))deallocate(p2b_fsd0)
     if( allocated(n2b_fsd0))deallocate(n2b_fsd0)

     if( allocated(n1b_cop))deallocate(n1b_cop)
     if( allocated(p1b_cop))deallocate(p1b_cop)
     if( allocated(n1b_dop))deallocate(n1b_dop)
     if( allocated(p1b_dop))deallocate(p1b_dop)
     if( allocated(n1b_isd))deallocate(n1b_isd)
     if( allocated(n1b_fsd))deallocate(n1b_fsd)
     if( allocated(n1b_phase))deallocate(n1b_phase)
     if( allocated(p1b_isd))deallocate(p1b_isd)
     if( allocated(p1b_fsd))deallocate(p1b_fsd)
     if( allocated(p1b_phase))deallocate(p1b_phase)
     if( allocated(p1b_isd0))deallocate(p1b_isd0)
     if( allocated(n1b_isd0))deallocate(n1b_isd0)
     if( allocated(p1b_fsd0))deallocate(p1b_fsd0)
     if( allocated(n1b_fsd0))deallocate(n1b_fsd0)

     if( allocated(n3b_op))deallocate(n3b_op)
     if( allocated(p3b_op))deallocate(p3b_op)
     if( allocated(n3b_isd))deallocate(n3b_isd)
     if( allocated(n3b_fsd))deallocate(n3b_fsd)
     if( allocated(n3b_phase))deallocate(n3b_phase)
     if( allocated(p3b_isd))deallocate(p3b_isd)
     if( allocated(p3b_fsd))deallocate(p3b_fsd)
     if( allocated(p3b_phase))deallocate(p3b_phase)


     if(np(1) > 1)then
         call set2bsectorjumps(1,.false.,  hermflag, .true. ,change2bodyqs, .false., 2 )
         call set2bsectorjumps(1,.true.,  hermflag, .true. ,change2bodyqs, .false., 2 )
         call masterfindconjumps2b(1)

        call master2bodyjumps(1,.false.  )

     endif

     if(np(2) > 1)then

         call set2bsectorjumps(2,.false.,  hermflag, .true. ,change2bodyqs, .false., 2 )
         call set2bsectorjumps(2,.true.,  hermflag, .true. ,change2bodyqs, .false., 2 )
       call masterfindconjumps2b(2)

       call master2bodyjumps(2,.false. )
     endif

     if( np(1)*np(2) > 0)then
!................ FIND SECTORS CONNECTED BY 1-BODY JUMPS.........
! NB: logical flags are: fill; hermflag; whermflag 
!
        call set1bsectorjumps(1, .false. , hermflag ,whermflagp , change1bodyqs )
        call set1bsectorjumps(1, .true.  , hermflag, whermflagp , change1bodyqs )

! NOTE that I use "hermiticity" on 1-body jumps but not "W-hermiticity"

        call set1bsectorjumps(2, .false. , hermflag, whermflagn , change1bodyqs )
        call set1bsectorjumps(2, .true.  , hermflag, whermflagn, change1bodyqs )
        call masterfindconjumps1b(1,change1bodyqs)
        call masterfindconjumps1b(2,change1bodyqs)  ! not really needed
        call master1bodyjumps(1,.false.)
        call master1bodyjumps(2,.false. )
!--------------- find conjugate jumps
        call masterfindconjumps1b(1,change1bodyqs)
        call masterfindconjumps1b(2,change1bodyqs)  ! not really needed

     endif  

     call master_op_stat

     if(allocated(opbundle))deallocate(opbundle)
     call distro_opbundles_over_fragments
	 if(iproc==0)print*,' Total opbundles for obs ',nopbundles
     call setMPIjumplimits(.false.)  !NOTE: might need to change logical flag
     if(np(1) > 1)then

        call master2bodyjumps(1,.true.  )

     endif

     if(np(2) > 1)then

      call master2bodyjumps(2,.true. )
     endif

     if( np(1)*np(2) > 0)then
        call master1bodyjumps(1,.true. )
        call master1bodyjumps(2,.true.)

     endif  
	 	 
!	 call set_unannexed

     call BMPI_BARRIER(icomm,ierr)

return

end subroutine reset_jumps4obs
!============================================================================
!
! match initial and final sectors to a sector jump
!
!  CALLED BY:
!   survey_geneologies
!
subroutine findsectorjump(it,nbody,is,fs,isjmp)
	use jumpNbody
	use jump3body
	implicit none
	integer :: it  ! species
	integer :: nbody
	integer :: is,fs  ! initial, final sector
	integer :: isjmp  ! sector jump
	
	integer :: jsjmp
	
	select case (nbody)
	case (1)
    do jsjmp= 1,x1bjump(it)%nsectjumps	   
		if((is == x1bjump(it)%isector(jsjmp) .and.  fs == x1bjump(it)%fsector(jsjmp)))then
			isjmp = jsjmp 
			return
		end if
		if((fs == x1bjump(it)%isector(jsjmp) .and.  is == x1bjump(it)%fsector(jsjmp)))then
			isjmp = -jsjmp 
			return
		end if		
		
	end do
			
	case (2)
    do jsjmp= 1,x1bjump(it)%nsectjumps	   
		if((is == x2bjump(it)%isector(jsjmp) .and.  fs == x2bjump(it)%fsector(jsjmp)))then
			isjmp = jsjmp 
			return 
		end if
		if((fs == x2bjump(it)%isector(jsjmp) .and.  is == x2bjump(it)%fsector(jsjmp)))then
			isjmp = -jsjmp 
			return 
		end if		
	end do
	
	case (3)
    do jsjmp= 1,x3bjump(it)%nsectjumps	   
		if((is == x3bjump(it)%isector(jsjmp) .and.  fs == x3bjump(it)%fsector(jsjmp)))then
			isjmp = jsjmp 
			return 
		end if
		if((fs == x3bjump(it)%isector(jsjmp) .and.  is == x3bjump(it)%fsector(jsjmp)))then
			isjmp = -jsjmp 
			return 
		end if		
	end do
    end select
	print*,' could not find sector jump ',it,nbody,is,fs
	stop
	
end subroutine findsectorjump


!============================================================================

end module jump_mod
