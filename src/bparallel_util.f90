!=================================================================
!
! file BPARALLEL_LIB3.F90
!
! routines to compute MPI distribution for multiple nodes; 
! originally only for one "fragment" (that is, lanczos vectors not broken up)
! all-in-one routine to be broken up in 7.2.13/14...
!
! started 4/2012 by CWJ
! original routines used when Lanczos vectors are NOT broken into fragments
! (supplanting routines handling fragments are found in bparallel_lib2.f90)
!
!===========================================================

module para_util_mod
	use para_bundles_mod
contains
!===============================================================
!
!  goes through opbundles and finds start, stops of each set of jumps
!
!  initiated 9/2012 by CWJ @ SDSU/TRIUMF
!
!  CALLED BY:
!     master_para_distribute (in bparallel_lib1.f90)
!     rebalance_MPI_distro  (in bparallel_lib2.f90)
!  CALLS
!     intron_master
!
!  INPUT:
!    fillredirect -- flag to fill redirection arrays
!
   subroutine setMPIjumplimits(filljumpdirect)
   use nodeinfo
   use jumpNbody
   use jump3body
   use jumplimits
   use flags3body
   use opbundles
   use butil_mod

   implicit none
   logical filljumpdirect

   integer iprocs,procstart,procstop
   integer ibundle
   integer(8) :: sumjumps,maxlocaljumps,localjumps,avgjumps
   integer :: maxproc   ! which MPI process has the most jumps
   integer :: aerr

   if(compactjumpstorage)call intron_master(filljumpdirect,.false.)    ! routine to set up redirection arrays for jumps; do not printout

!---- MODIFIED IN 7.5.4: always loop over all procs. This appears to solve a bug, but it's not clear what the root problem really is

    procstart = 0
	procstop  = nprocs -1

   if( .not. allocated(makep1bjumps))then

      allocate( makep1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 1")
   end if
   if( .not. allocated(maken1bjumps))then
      allocate( maken1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 2")
   end if
   if( .not. allocated(makeppjumps))then
      allocate( makeppjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 3")
   end if
   if( .not. allocated(makennjumps))then
      allocate( makennjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 4")
   end if
   if( .not. allocated(makePPPjumps))then
      allocate( makePPPjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 5")
   end if
   if( .not. allocated(makeNNNjumps))then
      allocate( makeNNNjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 6")
   end if
! KSM - FIX
! Make sure we initialize the whole allocated array because procstop/procstart
! may only be something like 4,4, but we have allocated from 0 to nprocs-1
   makep1bjumps(0:nprocs-1) = .false.
   maken1bjumps(0:nprocs-1) = .false.
   makeppjumps(0:nprocs-1) = .false.
   makennjumps(0:nprocs-1) = .false.
   makePPPjumps(0:nprocs-1) = .false.
   makeNNNjumps(0:nprocs-1) = .false.
! KSM - END FIX
   if(.not.allocated(startp1bjumps))then
      allocate( startp1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 7")
   end if
   if(.not.allocated(stopp1bjumps))then
     allocate( stopp1bjumps(0:nprocs-1), stat=aerr )
     if(aerr /= 0) call memerror("setMPIjumplimits 8")
   end if
   if( .not. allocated(startn1bjumps))then
      allocate( startn1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 9")
   end if
   if( .not. allocated(stopn1bjumps))then
      allocate( stopn1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 10")
   end if
   if( .not. allocated(startppjumps))then
      allocate( startppjumps(0:nprocs-1),  stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 11")
   end if
   if( .not. allocated(stopppjumps))then
      allocate( stopppjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 12")
   end if
    if( .not. allocated(startnnjumps))then
      allocate( startnnjumps(0:nprocs-1),  stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 13")
   end if
   if( .not. allocated(stopnnjumps))then
      allocate( stopnnjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 14")
   end if
   if( .not. allocated(startPPPjumps))then
      allocate( startPPPjumps(0:nprocs-1),  stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 15")
   end if
   if( .not. allocated(stopPPPjumps))then
      allocate( stopPPPjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 16")
   end if
   if( .not. allocated(startNNNjumps))then
      allocate( startNNNjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIjumplimits 17")
   end if
   if( .not. allocated(stopNNNjumps))then
       allocate( stopNNNjumps(0:nprocs-1), stat=aerr )
       if(aerr /= 0) call memerror("setMPIjumplimits 18")
   end if
! KSM - FIX
! initialize complete arrays
   startp1bjumps(0:nprocs-1) = 0
   stopp1bjumps(0:nprocs-1) = 0
   startn1bjumps(0:nprocs-1) = 0
   stopn1bjumps(0:nprocs-1) = 0
   startppjumps(0:nprocs-1) = 0
   stopppjumps(0:nprocs-1) = 0
   startnnjumps(0:nprocs-1) = 0
   stopnnjumps(0:nprocs-1) = 0
   startPPPjumps(0:nprocs-1) = 0
   stopPPPjumps(0:nprocs-1) = 0
   startNNNjumps(0:nprocs-1) = 0
   stopNNNjumps(0:nprocs-1) = 0
! KSM - END FIX
!   maxlocaljumps = 0
!   avgjumps = 0
!   maxproc = 0
   do iprocs = procstart,procstop
!       localjumps = 0
       startp1bjumps(iprocs) = totp1bjumps
       stopp1bjumps(iprocs ) = 0
       startn1bjumps(iprocs) = totn1bjumps
       stopn1bjumps(iprocs ) = 0
       startppjumps(iprocs) = totp2bjumps
       stopppjumps(iprocs ) = 0
       startnnjumps(iprocs) = totn2bjumps
       stopnnjumps(iprocs ) = 0
       startpppjumps(iprocs) = totp3bjumps
       stoppppjumps(iprocs ) = 0
       startnnnjumps(iprocs) = totn3bjumps
       stopnnnjumps(iprocs ) = 0
       makep1bjumps(iprocs) = .false.
       maken1bjumps(iprocs) = .false.
       makePPjumps(iprocs) = .false.
       makeNNjumps(iprocs) = .false.
       makePPPjumps(iprocs) = .false.
       makeNNNjumps(iprocs) = .false.
       if(opbundleend(iprocs) < 1)cycle

! .......... INITIAL SET FOR START/STOP OF JUMPS....
 
       do ibundle = opbundlestart(iprocs),opbundleend(iprocs)
           select case (opbundle(ibundle)%optype )
           case ('PP')
              startppjumps(iprocs) = bmin( startppjumps(iprocs), opbundle(ibundle)%pxstart )      
              stopppjumps(iprocs) = bmax( stopppjumps(iprocs), opbundle(ibundle)%pxend )     
           case ('NN')
              startnnjumps(iprocs) = bmin( startnnjumps(iprocs), opbundle(ibundle)%nxstart )      
              stopnnjumps(iprocs) = bmax( stopnnjumps(iprocs), opbundle(ibundle)%nxend )  

           case ('PPP')
              startpppjumps(iprocs) = bmin( startpppjumps(iprocs), opbundle(ibundle)%pxstart )      
              stoppppjumps(iprocs) = bmax( stoppppjumps(iprocs), opbundle(ibundle)%pxend )  

           case ('NNN')
              startNNNjumps(iprocs) = bmin( startNNNjumps(iprocs), opbundle(ibundle)%nxstart )      
              stopNNNjumps(iprocs) = bmax( stopNNNjumps(iprocs), opbundle(ibundle)%nxend )  

           case ('PN')
              startp1bjumps(iprocs) = bmin( startp1bjumps(iprocs), opbundle(ibundle)%pxstart )      
              stopp1bjumps(iprocs) = bmax( stopp1bjumps(iprocs), opbundle(ibundle)%pxend )  

              startn1bjumps(iprocs) = bmin( startn1bjumps(iprocs), opbundle(ibundle)%nsortstart )      
              stopn1bjumps(iprocs) = bmax( stopn1bjumps(iprocs), opbundle(ibundle)%nsortend ) 

           case ('PPN')
              startppjumps(iprocs) = bmin( startppjumps(iprocs), opbundle(ibundle)%pxstart )      
              stopppjumps(iprocs) = bmax( stopppjumps(iprocs), opbundle(ibundle)%pxend )  

              startn1bjumps(iprocs) = bmin( startn1bjumps(iprocs), opbundle(ibundle)%nsortstart )      
              stopn1bjumps(iprocs) = bmax( stopn1bjumps(iprocs), opbundle(ibundle)%nsortend ) 

           case ('PNN')
              startp1bjumps(iprocs) = bmin( startp1bjumps(iprocs), opbundle(ibundle)%pxstart )      
              stopp1bjumps(iprocs) = bmax( stopp1bjumps(iprocs), opbundle(ibundle)%pxend )  

              startnnjumps(iprocs) = bmin( startnnjumps(iprocs), opbundle(ibundle)%nsortstart )      
              stopnnjumps(iprocs) = bmax( stopnnjumps(iprocs), opbundle(ibundle)%nsortend ) 

           end select
         end do
         makep1bjumps(iprocs) = .false.
         maken1bjumps(iprocs) = .false.
         makePPjumps(iprocs) = .false.
         makeNNjumps(iprocs) = .false.
         makePPPjumps(iprocs) = .false.
         makeNNNjumps(iprocs) = .false.
         if(startppjumps(iprocs) <= stopppjumps(iprocs)) makeppjumps(iprocs) = .true. 
         if(startNNjumps(iprocs) <= stopNNjumps(iprocs)) makeNNjumps(iprocs) = .true. 
         if(startp1bjumps(iprocs) <= stopp1bjumps(iprocs)) makep1bjumps(iprocs) = .true. 
         if(startn1bjumps(iprocs) <= stopn1bjumps(iprocs)) maken1bjumps(iprocs) = .true. 
         if(threebody .and. startPPPjumps(iprocs) <= stoppppjumps(iprocs)) makepppjumps(iprocs) = .true. 
         if(threebody .and. startNNNjumps(iprocs) <= stopNNNjumps(iprocs)) makeNNNjumps(iprocs) = .true. 
! ............MODIFY IF COMPACT STORAGE OF JUMPS........		 
         if(compactjumpstorage)then    ! reset jump storage limits
		  if( makep1bjumps(iprocs)) then
			  startp1bjumps(iprocs) = 1
			  stopp1bjumps(iprocs)  = p1bjumplength(iprocs)
          end if
		  if( maken1bjumps(iprocs)) then
			  startn1bjumps(iprocs) = 1
			  stopn1bjumps(iprocs)  = n1bjumplength(iprocs)
          end if
		  if( makeppjumps(iprocs)) then
			  startPPjumps(iprocs) = 1
			  stopPPjumps(iprocs)  = PPjumplength(iprocs)
		  end if 
		  if( makeNNjumps(iprocs)) then
			  startNNjumps(iprocs) = 1
			  stopNNjumps(iprocs)  = NNjumplength(iprocs)
		  end if 
		  if( makepppjumps(iprocs)) then
			  startPPPjumps(iprocs) = 1
			  stopPPPjumps(iprocs)  = PPPjumplength(iprocs)
		  end if 
		  if( makeNNNjumps(iprocs)) then
			  startNNNjumps(iprocs) = 1
			  stopNNNjumps(iprocs)  = NNNjumplength(iprocs)
          end if

      end if
   end do

   if(nproc == 1 .and. .not.simulateMPI)then    ! set everthing to default
      makep1bjumps(0) = .true.
      maken1bjumps(0) = .true.
      makePPjumps(0) = .true.
      makeNNjumps(0) = .true.

      startp1bjumps(0) = 1
      stopp1bjumps(0)  = totp1bjumps
      startn1bjumps(0) = 1
      stopn1bjumps(0)  = totn1bjumps
      startppjumps(0) = 1
      stopppjumps(0)  = totp2bjumps
      startnnjumps(0) = 1
      stopnnjumps(0)  = totn2bjumps

      if(threebody)then

           makePPPjumps(0) = .true.
           makeNNNjumps(0) = .true.
           startPPPjumps(0) = 1
           stopPPPjumps(0) = totp3bjumps
           startNNNjumps(0) = 1
           stopNNNjumps(0) = totn3bjumps

      else
           makePPPjumps(0) = .false.
           makeNNNjumps(0) = .false.

      end if
      return
   end if
   return
   end subroutine setMPIjumplimits
!===============================================================
!
!  routine to set MPI limits for densities, but only for root node
!
!  CALLED BY
!   density1b_output
!   density1b_from_oldwfn
!
subroutine setMPIlimits4densities
   use nodeinfo
   use jumpNbody
   use jump3body
   use jumplimits
   use flags3body
   use opbundles
   use butil_mod

   implicit none
   integer iprocs,procstart,procstop
   integer ibundle
   logical printout
   integer(8) :: sumjumps,maxlocaljumps,localjumps,avgjumps
   integer :: aerr

!---- always loop over all procs. This appears to solve a bug, but it's not clear what the root problem really is

    procstart = 0
   	procstop  = nprocs -1
   if( .not. allocated(makep1bjumps))then
      allocate( makep1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 1")
   end if
   if( .not. allocated(maken1bjumps))then
      allocate( maken1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 2")
   end if
   if( .not. allocated(makeppjumps))then
      allocate( makeppjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 3")
   end if
   if( .not. allocated(makennjumps))then
      allocate( makennjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 4")
   end if
   if( .not. allocated(makepppjumps))then
      allocate( makePPPjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 5")
   end if
   if( .not. allocated(makeNNNjumps))then
      allocate( makeNNNjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 6")
   end if
   if(.not.allocated(startp1bjumps))then
      allocate( startp1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 10")
  	! initialize complete arrays
  	   startp1bjumps(0:nprocs-1) = 0
   end if
   if(.not.allocated(stopp1bjumps))then
      allocate(  stopp1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 10b")
	    stopp1bjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(startn1bjumps))then
      allocate( startn1bjumps(0:nprocs-1),  stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 11")
      startn1bjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(stopn1bjumps))then
      allocate(  stopn1bjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 11b")
      stopn1bjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(startppjumps))then
      allocate( startppjumps(0:nprocs-1),  stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 12")
      startppjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(stopppjumps))then
      allocate( stopppjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 12b")
      stopppjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(startnnjumps))then
      allocate( startnnjumps(0:nprocs-1),  stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 13")
      startnnjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(stopnnjumps))then
      allocate( stopnnjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 13b")
   stopnnjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(startPPPjumps))then
      allocate( startPPPjumps(0:nprocs-1),  stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 14")
      startPPPjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(stopPPPjumps))then
      allocate(  stopPPPjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 14b")
   stopPPPjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(startNNNjumps))then
      allocate( startNNNjumps(0:nprocs-1), stat=aerr )
      if(aerr /= 0) call memerror("setMPIlimits4densities 15")
      startNNNjumps(0:nprocs-1) = 0
   end if
   if( .not. allocated(stopNNNjumps))then
       allocate( stopNNNjumps(0:nprocs-1), stat=aerr )
       if(aerr /= 0) call memerror("setMPIlimits4densities 15b")
	   stopNNNjumps(0:nprocs-1) = 0
   end if

   makeppjumps(0:nprocs-1) =.false.
   makennjumps(0:nprocs-1) =.false.
   makePPPjumps(0:nprocs-1) =.false.
   makeNNNjumps(0:nprocs-1) =.false.

   startp1bjumps(0:nprocs-1) = 0
   stopp1bjumps(0:nprocs-1) = 0
   startn1bjumps(0:nprocs-1) = 0
   stopn1bjumps(0:nprocs-1) = 0
   startppjumps(0:nprocs-1) = 0
   stopppjumps(0:nprocs-1) = 0
   startnnjumps(0:nprocs-1) = 0
   stopnnjumps(0:nprocs-1) = 0
   startPPPjumps(0:nprocs-1) = 0
   stopPPPjumps(0:nprocs-1) = 0
   startNNNjumps(0:nprocs-1) = 0
   stopNNNjumps(0:nprocs-1) = 0
   
!   maxlocaljumps = 0
!   avgjumps = 0
!   maxproc = 0
   if(iproc==0)then   ! check jumps have been counted up
	   print*, ' total proton 1-body jumps ',totp1bjumps
	   print*, ' total neutron 1-body jumps ',totn1bjumps
   end if
	   
   do iprocs = procstart,procstop
!       localjumps = 0
       startp1bjumps(iprocs) = totp1bjumps
       stopp1bjumps(iprocs ) = 0
       startn1bjumps(iprocs) = totn1bjumps
       stopn1bjumps(iprocs ) = 0
	   
       if(opbundlestart(iprocs) < 1 .or. opbundleend(iprocs) < 1)cycle

! .......... INITIAL SET FOR START/STOP OF JUMPS....
 
       do ibundle = opbundlestart(iprocs),opbundleend(iprocs)
           select case (opbundle(ibundle)%optype )

           case ('P1B')
              startp1bjumps(iprocs) = bmin( startp1bjumps(iprocs), opbundle(ibundle)%pxstart )      
              stopp1bjumps(iprocs) = bmax( stopp1bjumps(iprocs), opbundle(ibundle)%pxend )  
           case ('N1B')
              startn1bjumps(iprocs) = bmin( startn1bjumps(iprocs), opbundle(ibundle)%nsortstart )      
              stopn1bjumps(iprocs) = bmax( stopn1bjumps(iprocs), opbundle(ibundle)%nsortend ) 			  
		   case default
!			  print*, 'OOPS some none onebody opbundle, something wrong !', iprocs, ibundle,opbundle(ibundle)%optype 
			  
		  end select
	      if(startp1bjumps(iprocs) <= stopp1bjumps(iprocs)) makep1bjumps(iprocs) = .true. 
	      if(startn1bjumps(iprocs) <= stopn1bjumps(iprocs)) maken1bjumps(iprocs) = .true. 

! ............MODIFY IF COMPACT STORAGE OF JUMPS........		 
	      if(compactjumpstorage)then    ! reset jump storage limits; TO BE ADDED
			  
		  end if
	  end do   ! ibundle
   end do   ! iprocs

   if(nproc == 1 .and. .not.simulateMPI)then    ! set everthing to default
      makep1bjumps(0) = .true.
      maken1bjumps(0) = .true.

      startp1bjumps(0) = 1
      stopp1bjumps(0)  = totp1bjumps
      startn1bjumps(0) = 1
      stopn1bjumps(0)  = totn1bjumps
  end if
!  if(iproc == 0)then
!	   print*,' MPI LIMITS SET '
!   end if
   return
end subroutine setMPIlimits4densities

!===============================================================
!
!  routine to count up how many jumps STORED on a node
!  INPUT: 
!     iprocs : which (dummy) processor node
!     iunit  : which unit to write to; if 0, don't print
!  OUTPUT
!     localjumps : how many jumps to be stored on node iprocs
!
!  CALLED BY
!      print_jump_distro
!      redistributor_frag
!      rebalance_MPI_distro

   subroutine count_jumps_on_node(iprocs,iunit,localjumps)

   use nodeinfo
   use jumpNbody
   use jump3body
   use jumplimits
   use flags3body
   use opbundles
   
   implicit none
   integer iprocs
   integer :: iunit
   integer(8) localjumps
   integer ibundle


!   if(iproc==0 .and. iprocs==0)print*,' CHECKING LIMITS ON JUMPS PER NODE '
   localjumps = 0
   if(makep1bjumps(iprocs))localjumps = localjumps + stopp1bjumps(iprocs)-startp1bjumps(iprocs)+1
   if(maken1bjumps(iprocs))localjumps = localjumps + stopn1bjumps(iprocs)-startn1bjumps(iprocs)+1
   if(makeppjumps(iprocs))localjumps = localjumps + stopppjumps(iprocs)-startPPjumps(iprocs)+1
   if(makennjumps(iprocs))localjumps = localjumps + stopnnjumps(iprocs)-startNNjumps(iprocs)+1
   if(makePPPjumps(iprocs))localjumps = localjumps + stoppppjumps(iprocs)-startpppjumps(iprocs)+1
   if(makeNNNjumps(iprocs))localjumps = localjumps + stopNNNjumps(iprocs)-startNNNjumps(iprocs)+1

!--------- WRITE OUT DATA ----------
   if (iunit > 0 .and. iproc ==0)then
      if( makep1bjumps(iprocs))write(iunit,*)' P1B ',stopp1bjumps(iprocs),startp1bjumps(iprocs)
   	  if( maken1bjumps(iprocs))write(iunit,*)' N1B ',stopn1bjumps(iprocs),startn1bjumps(iprocs)
   	  if( makeppjumps(iprocs))write(iunit,*)' PP ',stopppjumps(iprocs)-startppjumps(iprocs)+1
   	  if( makennjumps(iprocs))write(iunit,*)' NN ',stopnnjumps(iprocs)-startnnjumps(iprocs)+1
   	  if( makepppjumps(iprocs))write(iunit,*)' PPP ',stoppppjumps(iprocs)-startppjumps(iprocs)+1
      if( makennnjumps(iprocs))write(iunit,*)' NNN ',stopnnnjumps(iprocs)-startnnjumps(iprocs)+1
end if
   return

   end subroutine count_jumps_on_node

!===============================================================
!
!  writes information on operations
!
! CALLED BY master_para_distribute
!
!  NEEDS TO BE REVISED
!
subroutine print_op_stats
    use operation_stats
	use nodeinfo
	use io
	use fragments
	
	implicit none
	logical :: print_frag_stats = .true.  ! flag to print out information on frag-to-frag operations
	integer :: ifrag,jfrag
	
	if(iproc/=0)return
	
	write(6,*)' '
	write(6,*)' Total # of operations = ',allops
	write(6,*)' '
	write(logfile,*)' Total # of operations = ',allops
	
!............. THE FOLLOWING IS OPTIONAL......

    if(print_frag_stats)then
		write(logfile,*)' '
		write(logfile,*)' Information on distribution of work between fragments '
		do ifrag = 1,nfragments
			write(logfile,*)' Fragment ',ifrag,' -> :'
				write(logfile,'(5(1i5,1i18))')(jfrag,int(opfragstat(ifrag,jfrag)%optot,kind=8),jfrag=ifrag,nfragments)
         end do		
		
	end if	
	
	return
end subroutine print_op_stats
!===============================================================
!
! COMPUTES THE NUMBER OF OPERATIONS AS WELL THE JUMPS OVER WHICH THEY ARE DISTRIBUTED
!
!  ADDED 7.7.9
!   replaces print_ops_distro, print_jump_distro, print_opbundle_distro
!   NEW OUTPUT:  distrodata.bigstick 
!       REPLACES jumpdistro.bigstick,procjumpdistro.bigstick, opsdistro.bigstick,  
!
!  CALLED BY: 
!     master_para_distribute
!
subroutine print_all_distro
   use nodeinfo
   use opbundles
   use operation_stats
   use butil_mod
   use program_info
   use fragments
   use jumplimits,only:maxjumpmemory
   use io,only:logfile
   
   implicit none

   integer :: iob,jnode
   integer(8) :: localops,maxlocalops,localjumps
   integer(8) :: startz,stopz,laststart
   integer :: alldistrofile=97
   real    :: localwt
   real(8) :: maxwtops, totwtops,wtlocalops
   character(10) :: date,time,zone
   integer :: datetimeval(8)
   integer(8) :: sumjumps,maxlocaljumps,avgjumps
   integer :: maxproc
   integer :: bytesPerJump = 21

   if(iproc /= 0)return
   print*, ' Printing out operations distro'

   open(unit=alldistrofile,file='distrodata.bigstick',status='unknown')
   write(alldistrofile,*)' BIGSTICK Version ',version,lastmodified
   call date_and_time(date,time,zone,datetimeval)
   write(alldistrofile,*)' Run date: ',date(1:4),'-',date(5:6),'-',date(7:8)
   write(alldistrofile,*)'# bundles:    proton    x    neutron       total jumps'
!   write(alldistrofile,*)' Distribution of operations across MPI processes '
   maxlocalops = 0
   maxwtops = 0.0
   totwtops = 0.0 
   avgjumps = 0
   maxlocaljumps = 0
   maxproc = -1
   do jnode = 0,nprocs-1	 
      localops = 0
      if(opbundleend(jnode) < 1)cycle
      wtlocalops = 0.0
      do iob = opbundlestart(jnode),opbundleend(jnode)
         if(usewts)then
             select case (opbundle(iob)%optype)
                case('SPE')
                   localwt = opwtSPE
                case('PP ')
                    localwt = opwtPP
                   if (opbundle(iob)%hchar == 'b') localwt = opwtPPb
                case('NN ')
                    localwt = opwtNN
                case('PN ')
                    localwt = opwtPN
                   if (opbundle(iob)%hchar == 'b') localwt = opwtPNb
                case('PPP')
                    localwt = opwtPPP
                case('NNN')
                    localwt = opwtNNN
                case('PPN')
                    localwt = opwtPPN
                case('PNN')
                    localwt = opwtPNN
                case default
                    if(iproc == 0)print*,' OOPS wrong optype ',draft_opbundle(iob)%optype
                    stop
             end select
         else
             localwt = 1.d0
         end if
         select case ( opbundle(iob)%optype )
            case('PP ', 'PPP','PPN', 'SPE')
               startz = opbundle(iob)%pxstart
               stopz =  opbundle(iob)%pxend

            case('NN ','PN ','NNN','PNN')
               startz = opbundle(iob)%nxstart
               stopz =  opbundle(iob)%nxend

            case default
               startz = 0 ! prevent compiler warnings
               stopz = 0
         end select
         localops = localops + int(stopz-startz+1,8)*int(opbundle(iob)%min_nop,8)
         wtlocalops = wtlocalops+localwt*real(stopz-startz+1,8)*real(opbundle(iob)%min_nop,8)
         totwtops = totwtops + localwt*real(stopz-startz+1,8)*real(opbundle(iob)%min_nop,8)
      end do
      maxlocalops = bmax(maxlocalops, localops)
      maxwtops = bmax(maxwtops,wtlocalops)
      write(alldistrofile,*)' processor ',jnode !,' total operations ',localops,', time (ns) for ops ',wtlocalops
!      write(alldistrofile,*)' op bundles ',opbundlestart(jnode), opbundleend(jnode)
 
write(alldistrofile,*)' initial, final fragments : ',nodal(jnode)%ifragment, nodal(jnode)%ffragment

      do iob = opbundlestart(jnode),opbundleend(jnode)

         select case ( opbundle(iob)%optype )
             case('SPE')
                localjumps = 0
             case('PP','PPP')
                localjumps = int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)
             case('NN','NNN')
                localjumps = int( opbundle(iob)%nsortend -opbundle(iob)%nsortstart+1,8)
             case('PN','PPN')
                localjumps = int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)    & 
               + int( opbundle(iob)%nsortend -opbundle(iob)%nsortstart+1,8)

         end select
!         if(opbundle(iob)%optype=='PN')then
!             write(alldistrofile,'(1x,a3,1x,a,i14,a12,i12,a12,i12,a12)')opbundle(iob)%optype,opbundle(iob)%hchar, & 
!             int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)* & 
!             int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8), ' operations,', & 
!             int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8), ' P1b jumps, ',&		 
!             int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8), ' N1b jumps  '		
!			    write(distrofile,*)opbundle(iob)%pxend ,opbundle(iob)%pxstart, opbundle(iob)%nsortend ,opbundle(iob)%nsortstart
!         else
!             write(alldistrofile,'(1x,a3,1x,a,i14,a12,i12,a7)')opbundle(iob)%optype,opbundle(iob)%hchar, & 
!             int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)* & 
!             int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8), ' operations,',localjumps,' jumps '
!         end if
		 
		 write(alldistrofile,999)iob,opbundle(iob)%optype,opbundle(iob)%hchar, & 
         int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8), & 
         int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8),  & 
!         int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)* & 
!         int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8), &
		  localjumps
      end do
	  
!--------------- JUMPS STORAGE --------------------------------	  
      sumjumps = 0
      call count_jumps_on_node(jnode,0,sumjumps)
      write(alldistrofile,*)' Node ',jnode,', jump storage = ', &
         1.0e-6*bytesPerJump*sumjumps,' Mb '!, localjumps*1.0e-6*bytesPerJump
      if(sumjumps > maxlocaljumps)maxproc = jnode
      maxlocaljumps = bmax(maxlocaljumps,sumjumps)
      avgjumps = avgjumps + sumjumps	
999 format(i8,2x,a3,a1,3i15)	  


   end do   ! jnode

   totwtops = totwtops/float(nprocs)
   print*,' '
   print*,' Estimated efficiency of MPI distribution is approx ', & 
              totwtops/maxwtops
   print*,' (Detailed information on distribution across MPI in distrodata.bigstick)'
   print*,' '
   maxjumpmemory = 1.0e-6*bytesPerJump*maxlocaljumps
	
   avgjumps = int(avgjumps/real(nprocs,kind=8),8)
   write(6,*)' max local storage of jumps is ', &
       1.0e-6*bytesPerJump*maxlocaljumps,' Mb, on process ',maxproc
   write(6,*)' avg storage of jumps is ',1.0e-6*bytesPerJump*avgjumps

   write(logfile,*)' max local storage of jumps is ', &
       1.0e-6*bytesPerJump*maxlocaljumps,' Mb, on process ',maxproc
   write(logfile,*)' avg storage of jumps is ',1.0e-6*bytesPerJump*avgjumps
   
   close(alldistrofile)
   return
 end subroutine print_all_distro 

!===============================================================
!
! COMPUTES THE NUMBER OF OPERATIONS AS WELL THE JUMPS OVER WHICH THEY ARE DISTRIBUTED
!
!  CALLED BY: 
!     master_para_distribute
!
subroutine print_ops_distro
   use nodeinfo
   use opbundles
!   use tribution
   use operation_stats
   use butil_mod
   implicit none

   integer :: iob,jnode
   integer(8) :: localops,maxlocalops,localjumps
   integer(8) :: startz,stopz,laststart
   integer :: distrofile=39
   real    :: localwt
   real(8) :: maxwtops, totwtops,wtlocalops

   if(iproc /= 0)return
   print*, ' Printing out operations distro'

   open(unit=distrofile,file='opsdistro.bigstick',status='unknown')
   write(distrofile,*)' Distribution of operations across MPI processes '
   maxlocalops = 0
   maxwtops = 0.0
   totwtops = 0.0
   do jnode = 0,nprocs-1
      localops = 0
      if(opbundleend(jnode) < 1)cycle
      wtlocalops = 0.0
      do iob = opbundlestart(jnode),opbundleend(jnode)
         if(usewts)then
             select case (opbundle(iob)%optype)
                case('SPE')
                   localwt = opwtSPE
                case('PP ')
                    localwt = opwtPP
                   if (opbundle(iob)%hchar == 'b') localwt = opwtPPb
                case('NN ')
                    localwt = opwtNN
                case('PN ')
                    localwt = opwtPN
                   if (opbundle(iob)%hchar == 'b') localwt = opwtPNb
                case('PPP')
                    localwt = opwtPPP
                case('NNN')
                    localwt = opwtNNN
                case('PPN')
                    localwt = opwtPPN
                case('PNN')
                    localwt = opwtPNN
                case default
                    if(iproc == 0)print*,' OOPS wrong optype ',draft_opbundle(iob)%optype
                    stop
             end select
         else
             localwt = 1.d0
         end if
         select case ( opbundle(iob)%optype )
            case('PP ', 'PPP','PPN', 'SPE')
               startz = opbundle(iob)%pxstart
               stopz =  opbundle(iob)%pxend

            case('NN ','PN ','NNN','PNN')
               startz = opbundle(iob)%nxstart
               stopz =  opbundle(iob)%nxend

            case default
               startz = 0 ! prevent compiler warnings
               stopz = 0
         end select
         localops = localops + int(stopz-startz+1,8)*int(opbundle(iob)%min_nop,8)
         wtlocalops = wtlocalops+localwt*real(stopz-startz+1,8)*real(opbundle(iob)%min_nop,8)
         totwtops = totwtops + localwt*real(stopz-startz+1,8)*real(opbundle(iob)%min_nop,8)
      end do
      maxlocalops = bmax(maxlocalops, localops)
      maxwtops = bmax(maxwtops,wtlocalops)
      write(distrofile,*)' processor ',jnode,' total operations ',localops,', time (ns) for ops ',wtlocalops
      write(distrofile,*)' op bundles ',opbundlestart(jnode), opbundleend(jnode)

      do iob = opbundlestart(jnode),opbundleend(jnode)

         select case ( opbundle(iob)%optype )
             case('SPE')
                localjumps = 0
             case('PP','PPP')
                localjumps = int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)
             case('NN','NNN')
                localjumps = int( opbundle(iob)%nsortend -opbundle(iob)%nsortstart+1,8)
             case('PN','PPN')
                localjumps = int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)    & 
               + int( opbundle(iob)%nsortend -opbundle(iob)%nsortstart+1,8)

         end select
         if(opbundle(iob)%optype=='PN')then
             write(distrofile,'(1x,a3,1x,a,i14,a12,i12,a12,i12,a12)')opbundle(iob)%optype,opbundle(iob)%hchar, & 
             int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)* & 
             int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8), ' operations,', & 
             int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8), ' P1b jumps, ',&		 
             int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8), ' N1b jumps  '		
!			    write(distrofile,*)opbundle(iob)%pxend ,opbundle(iob)%pxstart, opbundle(iob)%nsortend ,opbundle(iob)%nsortstart
         else
             write(distrofile,'(1x,a3,1x,a,i14,a12,i12,a7)')opbundle(iob)%optype,opbundle(iob)%hchar, & 
             int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)* & 
             int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8), ' operations,',localjumps,' jumps '
         end if
      end do

      write(distrofile,*)' - - - - - - - - - - - - - - - - - '
   end do   ! jnode

   totwtops = totwtops/float(nprocs)
   print*,' '
   print*,' Estimated efficiency of MPI distribution is approx ', & 
              totwtops/maxwtops
   print*,' '
   close(distrofile)
   return
 end subroutine print_ops_distro 
   
!===============================================================
!
! routine to print information on storage of jumps over MPI processes
!
!  INPUT : iunit
!     normally iunit = 64
!     but if iunit = 62, then print out to an alternate file (for debugging purposes only)
!
!  CALLED BY: master_para_distribute
!
subroutine print_jump_distro(iunit)
   use nodeinfo
   use jumplimits
   use butil_mod
	
	implicit none
	integer :: iunit      ! output unit
   integer(8) :: sumjumps,maxlocaljumps,localjumps,avgjumps
	integer :: iprocs,maxproc
   integer :: bytesPerJump = 21
	
	if(iproc/=0 .or. nprocs == 1)return
	print*,' Printing out jump storage distro '
	if(iunit==64)then
        open(unit=64, file='jumpdistro.bigstick',status='unknown')
   end if
	if(iunit==62)then    ! used only for debugging
        open(unit=62, file='jumpdistro.alt.bigstick',status='unknown')
   end if
	write(iunit,*)' Distribution of jumps stored per MPI process '
	write(iunit,*)' (NOTE: More jumps may be stored than used due to sorting limits )'
	avgjumps = 0
	maxlocaljumps = 0
	do iprocs = 0,nprocs -1
      sumjumps = 0
      call count_jumps_on_node(iprocs,iunit,sumjumps)
      write(iunit,*)' Node ',iprocs,', jump storage = ', &
         1.0e-6*bytesPerJump*sumjumps,' Mb '!, localjumps*1.0e-6*bytesPerJump
      write(iunit,*)' -- -- -- -- '
      write(iunit,*)' '
      if(sumjumps > maxlocaljumps)maxproc = iprocs
      maxlocaljumps = bmax(maxlocaljumps,sumjumps)
      avgjumps = avgjumps + sumjumps		
	end do
   ! KSM:  Some array element sizes have changed.  Need to update
   maxjumpmemory = 1.0e-6*bytesPerJump*maxlocaljumps
	
   avgjumps = int(avgjumps/real(nprocs,kind=8),8)
   write(6,*)' max local storage of jumps is ', &
       1.0e-6*bytesPerJump*maxlocaljumps,' Mb, on process ',maxproc
   write(6,*)' avg storage of jumps is ',1.0e-6*bytesPerJump*avgjumps

   write(iunit,*)' max local storage of jumps is ', &
       1.0e-6*bytesPerJump*maxlocaljumps,' Mb, on process ',maxproc
   write(iunit,*)' avg storage of jumps is ',1.0e-6*bytesPerJump*avgjumps
	close(iunit)
	call print_opbundle_distro
	return
end subroutine print_jump_distro
!===============================================================
!
!
! CALLED BY:
!    redistributor_frag
!
subroutine get_jumpcount(iob,startX,stopX,startY,stopY,reset,jumpcount)
   use opbundles
   use butil_mod
   implicit none
   integer iob
   integer(8) :: startX,stopX,startY,stopY,jumpcount
   logical :: reset

   if(reset)then
      stopX = 0
      stopY = 0
      startX = opbundle(iob)%pxstart
      startY = opbundle(iob)%nsortstart

   end if
   select case (opbundle(iob)%optype )
         case( 'PP', 'PPP')
             startX = bmin(startX,opbundle(iob)%pxstart)
             stopX  = bmax(stopX,opbundle(iob)%pxend)
             startY = 1
             stopY  = 0

         case( 'NN', 'NNN')
             startY = bmin(startY,opbundle(iob)%nxstart)
             stopY = bmax(stopY,opbundle(iob)%nxend)
             startX = 1
             stopX  = 0

         case( 'PN', 'PPN', 'PNN')
             startX = bmin(startX,opbundle(iob)%pxstart)
             stopX  = bmax(stopX,opbundle(iob)%pxend)
             startY = bmin(startY,opbundle(iob)%nsortstart)
             stopY = bmax(stopY,opbundle(iob)%nsortend)

         case('SPE')
             startY = 1
             stopY  = 0
             startX = 1
             stopX  = 0
   end select

   jumpcount = stopX +1-startX + stopY+1-startY
   return
end subroutine get_jumpcount

!===============================================================
!
!  a subroutine to check we didn't lose any operations in 
!  splitting up opbundles
!
!  CALLED BY:
!      split_draft_opbundles_g
!      split_draft_opbundles_frag  (needs to be revisited for this)
!
subroutine check_ops(optarget,nob_draft)
   use opbundles
   use nodeinfo
   implicit none
   character(3) :: optarget
   integer nob_draft
   integer ibundle
   integer(8) :: nops, nops_draft

   nops_draft = 0

   do ibundle = 1,nob_draft
      if( draft_opbundle(ibundle)%optype == optarget)then

          nops_draft = nops_draft + & 
            int(1-draft_opbundle(ibundle)%pxstart+ draft_opbundle(ibundle)%pxend,8) & 
          *int( 1-draft_opbundle(ibundle)%nxstart+ draft_opbundle(ibundle)%nxend,8)
      end if

   end do
   nops = 0

   do ibundle = 1,nopbundles
      if( opbundle(ibundle)%optype == optarget)then

          nops = nops + & 
            int(1-opbundle(ibundle)%pxstart+ opbundle(ibundle)%pxend,8) & 
          *int( 1-opbundle(ibundle)%nxstart+ opbundle(ibundle)%nxend,8)
      end if

   end do
   if(nops /= nops_draft .and. iproc==0)then
      print*,' did not match # of operations ',optarget
      print*,nops_draft, nops
   do ibundle = 1,nob_draft    
	   if(draft_opbundle(ibundle)%optype==optarget)then
           print*,ibundle,' draft ',draft_opbundle(ibundle)%pxstart,draft_opbundle(ibundle)%pxend,& 
		                 draft_opbundle(ibundle)%nxstart,draft_opbundle(ibundle)%nxend
      end if
    end do

   do ibundle = 1,nopbundles    
      if(opbundle(ibundle)%optype==optarget)then
           print*,ibundle,opbundle(ibundle)%pxstart,opbundle(ibundle)%pxend,opbundle(ibundle)%nxstart,opbundle(ibundle)%nxend
      end if
    end do

      stop
   end if
   return

end subroutine check_ops
!========================================================
!
!  Routines for storing Lanczos vectors in core memory
!
!  Break up each vectors into "pieces"; 
!  on any given core store the same piece of every vector
!  The pieces are of length Lpiece << dimbasis
!  (the last piece will be shorter)
!
!  I.e., if a Lanczos vector looks like this:
!
! vec = (v_1 v_2 v_3 ... v_N)
!
!     then we break up each Lanczos vector
!       (v_1 v_2 ... v_Lpiece |  v_Lpiece+1  v_Lpiece+2 ... V_2*Lpiece |
!         V_2*Lpiece+1 V_2*Lpiece+2 ... V_3*Lpiece | V_3*Lpiece+1...)
!
! On compute node 1 we store the first "piece"
!     (v_1 v_2 ... v_Lpiece )  for ALL lanczos vectors
! One compute node 2 we store the second "piece"
!    (v_Lpiece+1  v_Lpiece+2 ... V_2*Lpiece )  for ALL Lanczos vectors
! and so on. 
! 
! If the first Lanczos vector is stored as (v1_1, v1_2, v1_3...)
! and the second Lanczos vector is stored as (v2_1, v2_2, v2_3...),
! then on each node we store the pieces sequentially (in vecstore). 
! Specifically on node 1 we store the first piece for all the lanczos vectors
! so that vecstore = (v1_1 v1_2 v1_3.. v1_Lpiece, v2_1 v2_2, ... v2_Lpiece...
!   v3_1 v3_2 ... v3_Lpiece, v4_1...) 
! while on node 2 we store the second piece of all the lanczos vectors
! vecstore = (v1_Lpiece+1 v1_Lpiece+2...v1_2*Lpiece, v2_Lpiece+1 ...etc)
! and so on. 
!
! This aids in reorthogonalization by minimizing communication:
! The basic step is the dot product, say between vi and vj
!           dot(i,j) = sum_k vi_k vj_k 
! In a naive implementation one would have to communicate entire lanczos vectors,
! for example, we might store v1 on node 1, v2 on node 2, and so on.
! But when we reorthogonalize vector vi we would have to send v1 to node i
! then send v2 to node i, etc.. 
!
! Instead we break up the dot product. On node 1 we carry out
!     dotpiece1(i,j) = sum_k=1,Lpiece vi_k vj_k
! while on node 2 we carry out
!     dotpiece2(i,j) = sum_k=Lpiece+2, 2*Lpiece vi_k vj_k
! and so on.  This is done very quickly. The dotpieces are 
! then transmited to node 0, where
!       dot(i,j) = sum_m dotpiecem(i,j)
! This combined dot product is then broadcast to all the nodes, that carry out
! the reorthogonalization, e.g.
!     vi_k = vi_k - sum_j dot(i,j) vj_k
! but with k =1, Lpiece on node 1, k = Lpiece+1, 2*Lpiece on node 2, and so on
! The dotpieces and dot products are very small packets of data. 
! After this one must normalize the vector in a similar fashion. 
! The final vector must be broadcast to start the next Lanczos iteration,
! but this is only one broadcast of a single large vector instead of multiple 
! transmissions.
!  -- CWJ Sept/Nov 2012
!     (Note: Oct 2012 confirmed with P. Maris that MFD uses similar technology)


!  IMPORTANT!  NEED TO REWRITE SO TAKING DOT PRODUCT AGAINST EXISTING VECTOR!!

!=======================================================
!
!  COMPUTE DISTRIBUTION OF LANCZOS VECTORS ACROSS NODES
!  as broken up into "pieces"
!
!  CALLED BY: lanczos_master in file blanczoslib1.f90

  subroutine distribute_lanczos_pieces
  use nodeinfo
  use basis
  use localvectors
  use flagger
  use precisions
  use lanczos_info
  use bmpi_mod
  use butil_mod
  implicit none

  integer(4) :: ierr
  real(kind=8)  :: maxmem
  integer (kind=8):: vstart,vend,vlength
  integer (kind=8) :: maxpiecedim
  integer(kind=8) :: leftover
  integer tag
  integer :: stat(MPI_STATUS_SIZE)
  integer :: iprocs
  integer :: aerr

  if(.not.storelanczosincoreMPI .or. nproc ==1)then
     Lstart = 1
     return
  end if
  if(asklancmemstore)then
        if(iproc ==0)then
           print*,' '
           print*,' Enter max memory/node to store Lanczos vectors (in Gb) '
           read*,maxmem
           maxpiecedim = int(maxmem/(lanc_prec)*1.e9,kind=8)
        end if
        call BMPI_BCAST(maxpiecedim,1,0,icomm,ierr) 
     else
        maxpiecedim = maxpiecedimdefault
        if(iproc ==0)then
           print*,' '
           print*,' Assuming max memory per node to store Lanczos vectors ',maxpiecedim*lanc_prec*1.e-9,' Gb '
           print*,' '
        end if
  end if

  Lpiece= dimbasis/nproc    ! default size
  if(Lpiece*nprocs < dimbasis)Lpiece = Lpiece +1
!............. CHECK THEY WILL FIT..............
  if (Lpiece*niter > maxpiecedim)then        ! exceeds memory
     if(iproc ==0)then
          print*,' Cannot store all lanczos vectors in memory '
          print*,' Would require ',Lpiece*niter*lanc_prec*1.e-9,' Gb'
          print*,' One option is to reduce # of iterations '
     end if
     call BMPI_ABORT(icomm,101,ierr)
  end if

  if( Lpiece < minpiece .and. minpiece*niter <= maxpiecedim)then  ! rescale
!      Lpiece = minpiece                                           ! this allows for efficiency
  endif

!................... ALLOCATE.........................

  vstart = Lpiece*iproc        ! starting point for the vector to be stored
  vlength = bmin(Lpiece,dimbasis-vstart)
  vend   = vstart+vlength

  Ldim = dimbasis/nproc
  buff_max = Ldim + 1
  leftover = dimbasis - nproc*Ldim
  if(leftover == 0)then
     Lstart = iproc*Ldim + 1
  else
     if(iproc < leftover)then
        Ldim = Ldim + 1
        Lstart = iproc*Ldim + 1
     else
        Lstart = leftover*(Ldim + 1) + (iproc - leftover)*Ldim + 1
     end if
  end if

  allocate(Lvec(Ldim,niter), stat=aerr)
  if(aerr /= 0) call memerror("subroutine distribute_lanczos_pieces 1")

!---- COLLECT Lstart INFO FROM EACH NODE ON ROOT NODE----

  
  if(iproc /= 0)then
       tag = iproc
       call BMPI_SEND(Lstart,1,0,tag,icomm,ierr) 
  end if

  if(iproc==0)then
      allocate(Lstart_r(0:nproc) , stat=aerr)
      if(aerr /= 0) call memerror("subroutine distribute_lanczos_pieces 2")
      Lstart_r(0) = 1
      Lstart_r(nproc) = dimbasis+1
      do iprocs = 1,nproc-1
          tag =iprocs
          call BMPI_Recv(Lstart_r(iprocs),1,iprocs,tag,icomm,stat,  ierr) 
      end do

  end if

  npiece = int(dimbasis / Lpiece, 4)
  if(dimbasis > npiece * Lpiece)npiece = npiece + 1
  if(npiece > nprocs)then
     if(iproc ==0)then
        print*,' Error in distributing lanczos vectors '
        print*,npiece, nproc
        call BMPI_ABORT(icomm,101,ierr)
     end if
  end if
  if(iproc ==0)then
     print*,' '
     print*,' Storage of Lanczos vectors distributed up across ',npiece,' nodes '
     print*,' Memory per node = ',Lpiece*niter*lanc_prec*1.e-9,' Gb '
     print*,' '
  end if
  return

  end subroutine distribute_lanczos_pieces
!=======================================================
!
!  MODEL DISTRIBUTION OF LANCZOS VECTORS ACROSS NODES
!  as broken up into "pieces"
!

  subroutine model_distribute_lanczos_pieces
  use nodeinfo
  use basis
  use localvectors
  use flagger
  use precisions
  use bmpi_mod
  implicit none

  integer(4) :: ierr
  real  :: maxmem
  integer (kind=8):: vstart,vend,vlength
  integer (kind=8) :: maxpiecedim
  integer :: niter

  if(iproc ==0)then
    print*,' '
    print*,' Enter max number of Lanczos iterations '
    read*,niter
    print*,' '
  end if
  if(asklancmemstore)then
        if(iproc ==0)then
           print*,' '
           print*,' Enter max memory/node to store Lanczos vectors (in Gb) '
           read*,maxmem
           maxpiecedim = int(maxmem/(lanc_prec)*1.e9,kind=8)
        end if
        call BMPI_Bcast(maxpiecedim,1,0,icomm,ierr) 
     else
        maxpiecedim = maxpiecedimdefault
        if(iproc ==0)then
           print*,' '
           print*,' Assuming max memory per node to store Lanczos vectors ',maxpiecedim*lanc_prec*1.e-9,' Gb '
           print*,' '
        end if
  end if

  Lpiece= dimbasis/nprocs    ! default size
  if(Lpiece*nprocs < dimbasis)Lpiece = Lpiece +1

!............. CHECK THEY WILL FIT..............
  if (Lpiece*niter > maxpiecedim)then        ! exceeds memory
     if(iproc ==0)then
          print*,' Cannot store all lanczos vectors in memory '
          print*,' Would require ',Lpiece*niter*lanc_prec*1.e-9,' Gb'
          print*,' One option is to reduce # of iterations '
          print*,Lpiece
     end if
     return
  end if

  if( Lpiece < minpiece .and. minpiece*niter <= maxpiecedim)then  ! rescale
      Lpiece = minpiece                                           ! this allows for efficiency
  endif

!................... ALLOCATE.........................

  npiece = int(dimbasis / Lpiece, 4)
  if(dimbasis > npiece * Lpiece)npiece = npiece + 1
  if(npiece > nprocs)then
     if(iproc ==0)then
        print*,' Error in distributing lanczos vectors '
        print*,npiece, nprocs
        print*,dimbasis,npiece*Lpiece
        call BMPI_ABORT(icomm,101,ierr)
     end if
  end if
  if(iproc ==0)then
     print*,' Storage of Lanczos vectors distributed up across ',npiece,' nodes '
     print*,' Memory per node = ',Lpiece*niter*lanc_prec*1.e-9,' Gb '
     print*,' '
  end if
  return

  end subroutine model_distribute_lanczos_pieces
!===============================================================
end module para_util_mod
