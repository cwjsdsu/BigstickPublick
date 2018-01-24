!
!
!  BINTRON.f90
!
!  *  routines used to study discontinuous storage of jumps 
!  which may contribute to large MPI memory burdens
!  *  initiated in 7.4.0 by CWJ Jan 2015
!
!  - - - - - - - - - - - - - - - - - -
!  steps:
!  * for a given MPI process, set up start and stop for each set of jumps used
!  * look for initial, final, and overlaps and compute discontinuities ("introns")
!  * calculate size of discontinuities/"introns" relative to total size
! 
!  Remedy: introduce shifts in storage of jumps, with corresponding shifts in opbundles
!
!  The ultimate products are arrays found in module jumplimits (in bmodules_parallel.f90)
!  P1Bjumpcut(:), P1Bendcut(:),P1Bjumpshift(:), etc (for N1B, PP, NN, PPP, NNN)
!  for a given MPI process, these index, respectively, the start and end of a set of jumps
!  and the shift needed to have all "required" jumps contiguous, starting from 1.
!  These are ultimately used in routine fetchjumpshift
!
!========================================================================= 
!
!   master routine to search for discontinuities in jump storage = "introns"
!   will always fill arrays PPjumplength(iproc) etc
!   option to fill arrays PPjumpcut(), PPjumpshift(),etc which are unique to each MPI process
!
! INPUT:
!    fillup--logical flag to (ultimately) fill up Pjumpcut/Xjumpshift arrays
!           should NOT be called if modeling
!    printout -- logical flag, writes to introns.bigstick
!   IMPORTANT: Should not do both fillup and printout simultaneously
!  
!  CALLED BY: 
!    setMPIjumplimits 
!    distro_opbundles_over_fragments
!
!  SUBROUTINES CALLED:
!    measure_introns

subroutine intron_master(fillup,printout)
   use flags3body
   use jumplimits
   use nodeinfo
   use menu_choices
   use bmpi_mod
   implicit none
   logical fillup,printout 
   integer :: procstart,procstop
   integer iprocs
   integer aerr,ierr

   call clocker('int','sta')
   if(iproc==0 .and. printout)print*,' Printing intron deletions to file introns.bigstick '
!   if(iproc==0 .and. fillup)print*,' Final filling of intron shifts '
   if(menu_char=='m' .and. fillup)then
	   
	   if(iproc==0)then
		   print*,' bad call, when modeling set fillup = .false. in call to intron_master '
	   endif
	   fillup= .false.
   end if
   if(fillup .and. printout)then
	   if(iproc==0)print*,' WARNING misuse of intron master, cannot fill up and print simultaneously '
	   call BMPI_Abort(icomm,101,ierr)
   end if

   if(.not.allocated(P1Bjumplength))then
       allocate(P1Bjumplength(0:nprocs-1),stat=aerr)
	   if(aerr /= 0) call memerror("intron_master 1")
       allocate(N1Bjumplength(0:nprocs-1),stat=aerr)
	   if(aerr /= 0) call memerror("intron_master 2")
       allocate(PPjumplength(0:nprocs-1),stat=aerr)
	   if(aerr /= 0) call memerror("intron_master 3")
       allocate(NNjumplength(0:nprocs-1),stat=aerr)
	   if(aerr /= 0) call memerror("intron_master 4")
	   if (threebody)then
	       allocate(PPPjumplength(0:nprocs-1),stat=aerr)
		   if(aerr /= 0) call memerror("intron_master 5")
	       allocate(NNNjumplength(0:nprocs-1),stat=aerr)
		   if(aerr /= 0) call memerror("intron_master 6")		   
	   end if
   end if
   if(iproc == 0 .and. printout)open(unit=63,file='introns.bigstick',status="unknown")
   if(.not.fillup)then
	   procstart = 0
	   procstop  = nprocs-1
   else
	   procstart = iproc
	   procstop = iproc
   end if
   do iprocs = procstart,procstop
      call measure_introns(iprocs,'PN ',1,fillup,printout)
      call measure_introns(iprocs,'PN ',2,fillup,printout)

      call measure_introns(iprocs,'NN ',2,fillup,printout)
      call measure_introns(iprocs,'PP ',1,fillup,printout)
      call measure_introns(iprocs,'PPP',1,fillup,printout)
      call measure_introns(iprocs,'NNN',2,fillup,printout)

      call measure_introns(iprocs,'PPN',1,fillup,printout)
      call measure_introns(iprocs,'PPN',2,fillup,printout)
      call measure_introns(iprocs,'PNN',1,fillup,printout)
      call measure_introns(iprocs,'PNN',2,fillup,printout)
   end do
   if(iproc ==0 .and. printout)close(63)
   call clocker('int','end')
   
   return

end subroutine intron_master
!========================================================================= 
!
!   routine to search for discontinuities in jump storage = "introns"
!
! INPUT:
!    iprocs -- MPI process (actual or modeled)
!    optype -- which kind of operation (PP etc)
!    it : species 
!    fillup--logical flag to (ultimately) fill up Pjumpcut/Xjumpshift arrays
!           should NOT be called if modeling
!    printout -- writes to introns.bigstick
!
!  CALLED BY: intron_master
!
!  SUBROUTINES CALLED:
!    survey_jumpstorage
!    map_introns

subroutine measure_introns(iprocs,optype,it,fillup,printout)
   use menu_choices
   use jumplimits
   use jumpNbody
   use jump3body
   use nodeinfo
   implicit none 

   interface
      subroutine survey_jumpstorage(iprocs,optype,it,njumpsets,jstarts,jstops)  !INTERFACE
      implicit none
!........INPUT ..............
      integer :: iprocs  ! dummy label for "current" MPI process (may be used for modeling)
      character(3) :: optype   ! label controls type of operation
      integer it
!....... OUTPUT...............
      integer :: njumpsets  ! # of sets of jump storage
      integer(8), pointer :: jstarts(:), jstops(:)  ! start and stop of jump storage
      end subroutine survey_jumpstorage       ! INTERFACE

      subroutine map_introns(njumpsets,jstarts,jstops,nexons,exstarts,exstops) ! INTERFACE
      implicit none
!........INPUT ..............
      integer :: njumpsets  ! # of sets of jump storage
      integer(8), pointer :: jstarts(:), jstops(:)  ! start and stop of jump storage
!....... OUTPUT...............
      integer ::  nexons 
      integer(8), pointer :: exstarts(:), exstops(:)
      end subroutine map_introns  ! INTERFACE

   end interface

!................INPUT................
   integer iprocs  ! which processor
   character(3) :: optype   ! what time of operation
   integer    :: it   ! what species for mixed-species operations
   logical    :: fillup    !  whether or not to fill out all arrays
   logical    :: printout   ! prints out introns
!............... INTERNAL.....
   integer(8),pointer :: jstarts(:), jstops(:),exstarts(:),exstops(:)
   integer :: njumpsets, nexons
   integer(8) :: storagesize,intronsize,exonsize
   integer(8), pointer :: cutstart(:) => NULL()
   integer(8), pointer :: cutstop(:) => NULL()
   integer(8), pointer :: jumpshift(:) => NULL()
   integer :: i 
   logical :: fill
   integer :: aerr
   logical :: printdetails = .false.  ! used only for debugging

   call survey_jumpstorage(iprocs,optype,it,njumpsets,jstarts,jstops)
   if(njumpsets==0)return
   if(associated(exstarts))nullify(exstarts)
   allocate(exstarts(njumpsets), stat=aerr)
   if(aerr /= 0) call memerror("measure introns 1")
   if(associated(exstops))nullify(exstops)
   allocate( exstops(njumpsets) , stat=aerr)
   if(aerr /= 0) call memerror("measure introns 2")
   
   call map_introns(njumpsets,jstarts,jstops,nexons,exstarts,exstops)
!..........
   fill = (nproc > 1 .and. menu_char/='m' .and. fillup)
   
   if(fill)then
      select case (optype)

         case('PN ')
            if(it==1)then
				if(allocated(P1Bjumpcut))deallocate(P1Bjumpcut)
				if(allocated(P1Bendcut))deallocate(P1Bendcut)
				if(allocated(P1Bjumpshift))deallocate(P1Bjumpshift)

               allocate(P1Bjumpcut(nexons+1),P1Bjumpshift(nexons),P1bendcut(nexons+1),stat=aerr )
			   if(aerr /= 0) call memerror("measure introns 3")
			   cutstop => P1Bendcut
               cutstart => P1Bjumpcut
               jumpshift=> P1Bjumpshift				   
            else
				if(allocated(N1Bjumpcut))deallocate(N1Bjumpcut)
				if(allocated(N1Bendcut))deallocate(N1Bendcut)
				if(allocated(N1Bjumpshift))deallocate(N1Bjumpshift)
               allocate(N1Bjumpcut(nexons+1),N1Bjumpshift(nexons),N1bendcut(nexons+1),stat=aerr )
			   if(aerr /= 0) call memerror("measure introns 4")
			   cutstop => N1Bendcut
               cutstart => N1Bjumpcut
               jumpshift=> N1Bjumpshift

            end if

         case('PPN')
            if(it==1)then
				if(allocated(PPjumpcut))deallocate(PPjumpcut)
				if(allocated(PPjumpshift))deallocate(PPjumpshift)
				if(allocated(PPendcut))deallocate(PPendcut)
				
               allocate(PPjumpcut(nexons+1),PPjumpshift(nexons),PPendcut(nexons+1),stat=aerr )
			   if(aerr /= 0) call memerror("measure introns 5")
			   
               cutstart => PPjumpcut
               jumpshift=> PPjumpshift
     		   cutstop => PPendcut
			   

            else
				if(allocated(N1Bjumpcut))deallocate(N1Bjumpcut)
				if(allocated(N1Bjumpshift))deallocate(N1Bjumpshift)
				if(allocated(N1Bendcut))deallocate(N1Bendcut)
				
               allocate(N1Bjumpcut(nexons+1),N1Bjumpshift(nexons),N1bendcut(nexons+1),stat=aerr )
			   if(aerr /= 0) call memerror("measure introns 6")
			   cutstop => N1Bendcut			   
               cutstart => N1Bjumpcut
               jumpshift=> N1Bjumpshift

            end if
         case('PNN')
            if(it==1)then
				if(allocated(P1Bjumpcut))deallocate(P1Bjumpcut)
				if(allocated(P1Bjumpshift))deallocate(P1Bjumpshift)
				if(allocated(P1Bendcut))deallocate(P1Bendcut)
				
               allocate(P1Bjumpcut(nexons+1),P1Bjumpshift(nexons),P1bendcut(nexons+1),stat=aerr )
			   if(aerr /= 0) call memerror("measure introns 7")
			   
               cutstart => P1Bjumpcut
               jumpshift=> P1Bjumpshift
			   cutstop => P1Bendcut

            else
			    if(allocated(NNjumpcut))deallocate(NNjumpcut)
			    if(allocated(NNjumpshift))deallocate(NNjumpshift)
				if(allocated(NNendcut))deallocate(NNendcut)
				
               allocate(NNjumpcut(nexons+1),NNjumpshift(nexons),NNendcut(nexons+1),stat=aerr )

			   if(aerr /= 0) call memerror("measure introns 8")
			   cutstop => NNendcut  
               cutstart => NNjumpcut
               jumpshift=> NNjumpshift

            end if

         case('PP ')
   		    if(allocated(PPjumpcut))deallocate(PPjumpcut)
		    if(allocated(PPjumpshift))deallocate(PPjumpshift)
			if(allocated(PPendcut))deallocate(PPendcut)
			
            allocate(PPjumpcut(nexons+1),PPjumpshift(nexons),PPendcut(nexons+1),stat=aerr )
		    if(aerr /= 0) call memerror("measure introns 9")
            cutstart => PPjumpcut
            jumpshift=> PPjumpshift
		    cutstop => PPendcut
			

         case('NN ')
		    if(allocated(NNjumpcut))deallocate(NNjumpcut)
		    if(allocated(NNjumpshift))deallocate(NNjumpshift)
			if(allocated(NNendcut))deallocate(NNendcut)
            allocate(NNjumpcut(nexons+1),NNjumpshift(nexons),NNendcut(nexons+1),stat=aerr )
		    if(aerr /= 0) call memerror("measure introns 10")
 		    cutstop => NNendcut  
            cutstart => NNjumpcut
            jumpshift=> NNjumpshift

         case('PPP')
   		    if(allocated(PPPjumpcut))deallocate(PPPjumpcut)
		    if(allocated(PPPjumpshift))deallocate(PPPjumpshift)
			if(allocated(PPPendcut))deallocate(PPPendcut)
			
            allocate(PPPjumpcut(nexons+1),PPPjumpshift(nexons),PPPendcut(nexons+1),stat=aerr )

		    if(aerr /= 0) call memerror("measure introns 11")
            cutstart => PPPjumpcut
            jumpshift=> PPPjumpshift
		    cutstop => PPPendcut
			

         case('NNN')
   	        if(allocated(NNNjumpcut))deallocate(NNNjumpcut)
	        if(allocated(NNNjumpshift))deallocate(NNNjumpshift)
			if(allocated(NNNendcut))deallocate(NNNendcut)
			
            allocate(NNNjumpcut(nexons+1),NNNjumpshift(nexons),NNNendcut(nexons+1),stat=aerr )
		    if(aerr /= 0) call memerror("measure introns 12")
 		    cutstop => NNNendcut  	
            cutstart => NNNjumpcut
            jumpshift=> NNNjumpshift

       end select
   end if
!......... NOW ANALYZE FOR FRACTION OF INTRONS.......
!           AND SET UP REDIRECTION ARRAYS
!  for jump type X, if native (original) jump label i,  for some index n Xjumpcut(n)<= i < Xjumpcut(n+1)
!  new jump label iprime = i - Xjumpshift(n) 
!  new indices go from 1 to Xjumplength(iproc)

   storagesize = jstops(njumpsets)-jstarts(1)+1
   intronsize = 0
   exonsize   = exstops(1) - exstarts(1)+1
   if(fill)cutstart(1) = exstarts(1)
   if(fill)cutstop(1)=   exstops(1)
   if(fill)jumpshift(1) = exstarts(1)-1

   if(nexons>1)then
       do i = 2,nexons
          if(exstarts(i) - exstops(i-1) < 2)then
              print*,' some problems with exons/introns '
			  print*,i,nexons
			  print*,exstarts
			  print*,exstops
              stop
          end if
          intronsize = intronsize + exstarts(i)-exstops(i-1)-1
          if(fill)cutstart(i) = exstarts(i)
		  if(fill)cutstop(i)  =exstops(i)
          if(fill)jumpshift(i) = exstarts(i)-exonsize-1
          exonsize   = exonsize   + exstops(i) -exstarts(i) +1
		  
       end do
   end if
   if(iproc==0 .and. printout)then
	     write(63,6363)iprocs,optype, it, real(intronsize)*13.e-6,storagesize*13.e-6
6363 format(' On MPI proc ',i7,a4,i2,f12.1,' Mb deleted out of ',f10.1)		 
6364 format(' On MPI proc ',i7,a4,i2,i12,' introns deleted out of ',i12)		 

	 end if
!....................................................
   select case (optype)
         case('PN ')
            if(it==1)then
               p1bjumplength(iprocs) = exonsize
               if(fill) cutstart(nexons+1) = totp1bjumps+1
               if(fill) cutstop(nexons+1) = totp1bjumps+1
            else
               n1bjumplength(iprocs) = exonsize
               if(fill)cutstart(nexons+1) = totn1bjumps+1
               if(fill)cutstop(nexons+1) = totn1bjumps+1
            end if
         case('PNN')
            if(it==1)then
               p1bjumplength(iprocs) = exonsize
               if(fill)cutstart(nexons+1) = totp1bjumps+1
               if(fill)cutstop(nexons+1) = totp1bjumps+1
            else
               nnjumplength(iprocs) = exonsize
               if(fill)cutstart(nexons+1) = totn2bjumps+1
               if(fill)cutstop(nexons+1) = totn2bjumps+1
            end if
         case('PPN')
            if(it==1)then
               ppjumplength(iprocs) = exonsize
               if(fill)cutstart(nexons+1) = totp2bjumps+1
               if(fill)cutstop(nexons+1) = totp2bjumps+1
            else
               n1bjumplength(iprocs) = exonsize
               if(fill)cutstart(nexons+1) = totn1bjumps+1
               if(fill)cutstop(nexons+1) = totn1bjumps+1
            end if

         case('PP')
            ppjumplength(iprocs) = exonsize
            if(fill)cutstart(nexons+1) = totp2bjumps+1
            if(fill)cutstop(nexons+1) = totp2bjumps+1
         case('NN')
            nnjumplength(iprocs) = exonsize
            if(fill)cutstart(nexons+1) = totn2bjumps+1
            if(fill)cutstop(nexons+1) = totn2bjumps+1
         case('PPP')
            pppjumplength(iprocs) = exonsize
            if(fill)cutstart(nexons+1) = totp3bjumps+1
            if(fill)cutstop(nexons+1) = totp3bjumps+1
         case('NNN')
            nnnjumplength(iprocs) = exonsize
            if(fill)cutstart(nexons+1) = totn3bjumps+1
            if(fill)cutstop(nexons+1) = totn3bjumps+1
   end select

   nullify(jstarts,jstops,exstarts,exstops)
   return

end subroutine measure_introns

!========================================================================= 
!
! for a given MPI process and optype, set up and fill arrays for jump start, stops
!
! INPUT: 
!    iprocs = MPI process (may be dummy for modeling)
!    optype = type of operation, i.e. PN, PP, etc.
!    it     = species for mixed operator, e.g. P of PN, N of PPN, etc.
!
! OUTPUT:
!    njumpsets:  # of sets of jumps of type "optype"
!    jstarts(:),jstops(:): arrays of starts and stops, i = 1,njumpsets
!
! CALLED BY: measure_introns
!
subroutine survey_jumpstorage(iprocs,optype,it,njumpsets,jstarts,jstops)

   use opbundles
   use nodeinfo
   use jumplimits
   use flags3body
   implicit none
!........INPUT ..............
   integer :: iprocs  ! dummy label for "current" MPI process (may be used for modeling)
   character(3) :: optype   ! label controls type of operation
   integer :: it
!....... OUTPUT...............
   integer :: njumpsets  ! # of sets of jump storage
   integer(8), pointer :: jstarts(:), jstops(:)  ! start and stop of jump storage
!........ INTERNAL.................
   integer :: ibundle
   integer :: nset
   integer :: aerr
   
!...... COUNT UP
   njumpsets = 0
   do ibundle = opbundlestart(iprocs),opbundleend(iprocs)
       if(opbundle(ibundle)%optype /= optype)cycle
       njumpsets = njumpsets+1
   end do ! ibundle
   if(njumpsets==0)return
!....... ALLOCATE......
   if(associated(jstarts))nullify(jstarts)
   if(associated(jstops))nullify(jstops)
   allocate( jstarts(njumpsets), stat=aerr)
   if(aerr /= 0) call memerror("survey jump storage 1")
   allocate( jstops(njumpsets), stat=aerr)
   if(aerr /= 0) call memerror("survey jump storage 2")
   nset = 0
   do ibundle = opbundlestart(iprocs),opbundleend(iprocs)
       if(opbundle(ibundle)%optype /= optype)cycle
       nset = nset+1
       select case(optype)
          case('PP','PPP')
             jstarts(nset) = opbundle(ibundle)%pxstart
             jstops(nset) = opbundle(ibundle)%pxend

          case('NN','NNN')
             jstarts(nset) = opbundle(ibundle)%nxstart
             jstops(nset) = opbundle(ibundle)%nxend
          case('PN','PPN','PNN')
             if(it == 1)then
                jstarts(nset) = opbundle(ibundle)%pxstart
                jstops(nset) = opbundle(ibundle)%pxend
             else

               jstarts(nset) = opbundle(ibundle)%nsortstart
               jstops(nset) = opbundle(ibundle)%nsortend
             end if

          case default 
             print*,iproc,' Wrong kind of optype ',optype
             stop
       end select

   end do ! ibundle  

   return

end subroutine survey_jumpstorage
!========================================================================= 
!
! for a given set of jump start and stop arrays, create a master map
! (actually maps "exons" stored jumps, but same idea)
!
! INPUT:
!   njumpsets = original number of sets of jumps of a given type on a given MPI process
!   jstarts(:), jstops(:) : arrays of starts and stops, i = 1,njumpsets
!         (computed in routine survey_jumpstorage)
! OUTPUT:
!   nexons: = final number of continous sets of jump storage
!   exstarts(i), exstops(i): where "exons" (actually storage) start and stop
!
! CALLED BY: measure_introns
!
subroutine map_introns(njumpsets,jstarts,jstops,nexons,exstarts,exstops)
   use butil_mod
   implicit none
!........INPUT ..............

   integer :: njumpsets  ! # of sets of jump storage
   integer(8), pointer :: jstarts(:), jstops(:)  ! start and stop of jump storage
!....... OUTPUT...............
   integer ::  nexons 
   integer(8), pointer :: exstarts(:), exstops(:)

!........INTERNAL.......
   integer :: i,j,istart
   integer(8) :: tmpstart,tmpstop

!....... SORT BY START.....

   do i = 1,njumpsets-1
      tmpstart = jstarts(i)
      istart = i
      do j = i+1,njumpsets
         if( jstarts(j) < tmpstart)then
            istart = j
            tmpstart = jstarts(j)
         end if
      end do  ! j
      if(istart/=i)then ! swap
         jstarts(istart) = jstarts(i)
         jstarts(i) = tmpstart
         tmpstop = jstops(istart)
         jstops(istart) = jstops(i)
         jstops(i) = tmpstop
      end if
   end do  !i

!.......... NOW CHECK FOR OVERLAPS.......
   exstarts(1) = jstarts(1)
   exstops(1)  = jstops(1)
   nexons = 1

   do i = 2,njumpsets
      if(jstarts(i) > exstops(nexons)+1)then ! start new exons
         nexons = nexons+1
         exstarts(nexons)= jstarts(i)
         exstops(nexons) = jstops(i)

      else
         exstops(nexons) = bmax(exstops(nexons),jstops(i) )
      end if
   end do

!...........
end subroutine map_introns
!===========================================================
!
! routine to compute shift in jumps  added 7.4.2
! these are specific to a given MPI process and will not be the same across processes
!
! INPUT:
!   jumptype = 2-char code for type of jump arrays
!   jumpindx    :  jump label (changed upon return)
!
!  for jump type X, if native (original) jump label i,  for some index n Xjumpcut(n)<= i < Xjumpcut(n+1)
!  new jump label iprime = i - Xjumpshift(n) 
!
! CALLED BY:
!

subroutine fetchjumpshift(jumptype,jumpindx)
   use jumplimits
   use nodeinfo
   
   implicit none
   character(2) :: jumptype
   integer(8) :: jumpindx
   integer n
   integer(8), pointer :: Xjumpcut(:), Xjumpshift(:),Xendcut(:)

   if(.not.compactjumpstorage)return

   select case(jumptype)
     case('P1')
        Xjumpcut => P1Bjumpcut
        Xendcut  => P1Bendcut
        Xjumpshift => P1Bjumpshift
     case('N1')
        Xjumpcut => N1Bjumpcut
        Xendcut  => N1Bendcut
        Xjumpshift => N1Bjumpshift
     case('P2')
        Xjumpcut => PPjumpcut
        Xendcut  => PPendcut
        Xjumpshift => PPjumpshift
     case('N2')
        Xjumpcut => NNjumpcut
        Xendcut  => NNendcut
        Xjumpshift => NNjumpshift
     case('P3')
        Xjumpcut => PPPjumpcut
        Xendcut  => PPPendcut
        Xjumpshift => PPPjumpshift
     case('N3')
        Xjumpcut => NNNjumpcut
        Xendcut  => NNNendcut
        Xjumpshift => NNNjumpshift
     case default
        print *, "bad case" ! helps uninitialized var detection
        stop 1
   end select    
  
!........ find index

   do n = 1,100000
      if(jumpindx < Xjumpcut(n))exit
   end do
   n=n-1
   if(n < 1)then
	   jumpindx=-1
	   return
   end if
   if(jumpindx > Xendcut(n))then
	   jumpindx = -1
   else
       jumpindx =jumpindx - Xjumpshift(n)
   end if
   return
end subroutine fetchjumpshift
!===========================================================
!
! routine to update opbundles 
! added in 7.5.2
!
!  CALLED BY MAIN ROUTINE AFTER JUMPMASTER
!

subroutine updateopbundles
   use jumplimits
   use nodeinfo
   use opbundles
   use bmpi_mod
   implicit none
   interface 
       subroutine fetchjumpshift(jumptype,jumpindx)        ! INTERFACE
       use jumplimits
       implicit none
       character(2) :: jumptype
       integer(8) :: jumpindx
       end subroutine fetchjumpshift                       ! INTERFACE
   end interface
   integer :: procstart,procstop,iprocs
   integer :: ibundle
   integer(8)::startjump,stopjump
   integer :: ierr
   
   if(.not.compactjumpstorage)return
!........ OPTION TO SIMULATE MPI ON A SINGLE "CORE".....
   if((simulateMPI .and. nproc ==1))then
      procstart = 0
      procstop  = nprocs -1
   else
      procstart = iproc
      procstop  = iproc
   end if   
!   print*,iproc,' Updating opbundles...'
   do iprocs = procstart,procstop
      do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
 
          select case( opbundle(ibundle)%optype )
          
          case ('PN')
             startjump = opbundle(ibundle)%pxstart
             stopjump  = opbundle(ibundle)%pxend
             call fetchjumpshift('P1',startjump)
             call fetchjumpshift('P1',stopjump)
! print*,iproc,ibundle,opbundle(ibundle)%pxstart,startjump
             opbundle(ibundle)%pxstart = startjump
             opbundle(ibundle)%pxend   =stopjump
             startjump = opbundle(ibundle)%nxstart
             stopjump  = opbundle(ibundle)%nxend
             call fetchjumpshift('N1',startjump)
             call fetchjumpshift('N1',stopjump)
             opbundle(ibundle)%nxstart = startjump
             opbundle(ibundle)%nxend   =stopjump     
             startjump = opbundle(ibundle)%nsortstart
             stopjump  = opbundle(ibundle)%nsortend
             call fetchjumpshift('N1',startjump)
             call fetchjumpshift('N1',stopjump)
             opbundle(ibundle)%nsortstart = startjump
             opbundle(ibundle)%nsortend   =stopjump            
             
          case ('PP')
             startjump = opbundle(ibundle)%pxstart
             stopjump  = opbundle(ibundle)%pxend
             call fetchjumpshift('P2',startjump)
             call fetchjumpshift('P2',stopjump)
             opbundle(ibundle)%pxstart = startjump
             opbundle(ibundle)%pxend   =stopjump

          case ('NN')
             startjump = opbundle(ibundle)%nxstart
             stopjump  = opbundle(ibundle)%nxend
             call fetchjumpshift('N2',startjump)
             call fetchjumpshift('N2',stopjump)
             opbundle(ibundle)%nxstart = startjump
             opbundle(ibundle)%nxend   =stopjump
		  
          case ('PPN')
             startjump = opbundle(ibundle)%pxstart
             stopjump  = opbundle(ibundle)%pxend
             call fetchjumpshift('P2',startjump)
             call fetchjumpshift('P2',stopjump)
! print*,iproc,ibundle,opbundle(ibundle)%pxstart,startjump
             opbundle(ibundle)%pxstart = startjump
             opbundle(ibundle)%pxend   =stopjump
             startjump = opbundle(ibundle)%nxstart
             stopjump  = opbundle(ibundle)%nxend
             call fetchjumpshift('N1',startjump)
             call fetchjumpshift('N1',stopjump)
             opbundle(ibundle)%nxstart = startjump
             opbundle(ibundle)%nxend   =stopjump     
             startjump = opbundle(ibundle)%nsortstart
             stopjump  = opbundle(ibundle)%nsortend
             call fetchjumpshift('N1',startjump)
             call fetchjumpshift('N1',stopjump)
             opbundle(ibundle)%nsortstart = startjump
             opbundle(ibundle)%nsortend   =stopjump   
			 
           case ('PNN')
                startjump = opbundle(ibundle)%pxstart
                stopjump  = opbundle(ibundle)%pxend
                call fetchjumpshift('P1',startjump)
                call fetchjumpshift('P1',stopjump)
   ! print*,iproc,ibundle,opbundle(ibundle)%pxstart,startjump
                opbundle(ibundle)%pxstart = startjump
                opbundle(ibundle)%pxend   =stopjump
                startjump = opbundle(ibundle)%nxstart
                stopjump  = opbundle(ibundle)%nxend
                call fetchjumpshift('N2',startjump)
                call fetchjumpshift('N2',stopjump)
                opbundle(ibundle)%nxstart = startjump
                opbundle(ibundle)%nxend   =stopjump     
                startjump = opbundle(ibundle)%nsortstart
                stopjump  = opbundle(ibundle)%nsortend
                call fetchjumpshift('N2',startjump)
                call fetchjumpshift('N2',stopjump)
                opbundle(ibundle)%nsortstart = startjump
                opbundle(ibundle)%nsortend   =stopjump   
				
	            case ('PPP')
	               startjump = opbundle(ibundle)%pxstart
	               stopjump  = opbundle(ibundle)%pxend
	               call fetchjumpshift('P3',startjump)
	               call fetchjumpshift('P3',stopjump)
	               opbundle(ibundle)%pxstart = startjump
	               opbundle(ibundle)%pxend   =stopjump

	            case ('NNN')
	               startjump = opbundle(ibundle)%nxstart
	               stopjump  = opbundle(ibundle)%nxend
	               call fetchjumpshift('N3',startjump)
	               call fetchjumpshift('N3',stopjump)
	               opbundle(ibundle)%nxstart = startjump
	               opbundle(ibundle)%nxend   =stopjump
	       end select
		  
      end do ! ibundle
   end do ! iprocs
   return
end subroutine updateopbundles
!===========================================================
! added in 7.5.3
!
! checks sectorjumps to see if "exons" created within a sectorjump
! are contiguous
!
! INPUT:
!   it = species (proton/neutron)
!   nbody = 1, 2, or 3-body
!
!  CALLED BY:
!   master1bodyjumps
!   master2bodyjumps
!   master3bodyjumps

subroutine surveysectorjumps(it,nbody)
	use jumpNbody
	use jump3body
	use jumplimits
	use nodeinfo
   use butil_mod
	implicit none
	integer :: it  ! species (proton/neutron)
	integer :: nbody   ! rank of jump operator
	
	type (jumpsect), pointer :: xNjump
	integer :: isjmp
   integer :: nex
   integer(8), pointer :: Xjumpcut(:), Xjumpshift(:),Xendcut(:)
	
	logical :: containsexons, contiguous
	integer(8) :: jumpstart,jumpend
	integer :: nexons

	select case (nbody)
	   case (1)
	       xNjump => x1bjump(it)

	       if(it==1)then
			   if(.not.makep1bjumps(iproc))return
		       Xjumpcut => P1Bjumpcut
		       Xendcut  => P1Bendcut
		       Xjumpshift => P1Bjumpshift
			   nexons = size(P1bjumpshift)
		   else
			   if(.not.maken1bjumps(iproc))return

		       Xjumpcut => N1Bjumpcut
		       Xendcut  => N1Bendcut
		       Xjumpshift => N1Bjumpshift
			   nexons = size(N1bjumpshift)
			   
	       end if
	
	   case (2)
           xNjump => x2bjump(it)
	       if(it==1)then
			   if(.not.makeppjumps(iproc))return
			   
		       Xjumpcut => PPjumpcut
		       Xendcut  => PPendcut
		       Xjumpshift => PPjumpshift
			   nexons = size(PPjumpshift)
			   
		   else
			   if(.not.makeNNjumps(iproc))return
			   
		       Xjumpcut => NNjumpcut
		       Xendcut  => NNendcut
		       Xjumpshift => NNjumpshift
			   nexons = size(NNjumpshift)
			   
	       end if
		   
	   case(3)
           xNjump => x3bjump(it)
	       if(it==1)then
			   if(.not.makepppjumps(iproc))return
			   
		       Xjumpcut => PPPjumpcut
		       Xendcut  => PPPendcut
		       Xjumpshift => PPPjumpshift
			   nexons = size(PPPjumpshift)
			   
		   else
			   if(.not.makennnjumps(iproc))return
			   
		       Xjumpcut => NNNjumpcut
		       Xendcut  => NNNendcut
		       Xjumpshift => NNNjumpshift
			   nexons = size(NNNjumpshift)
			   
	       end if		   
      case default
        print *, "stop in subroutine surveysectorjumps" ! helps uninitialized var detection
        stop 1
    end select
	do isjmp = 1,xNjump%nsectjumps
		containsexons=.false.
		contiguous   =.true.
		
		jumpstart = xNjump%sjmp(isjmp)%nstart +1
		jumpend   = xNjump%sjmp(isjmp)%nstart +xNjump%sjmp(isjmp)%njumps
				
!...... CHECK IF THERE ARE EXONS, AND HOW MANY, BETWEEN THESE TWO LIMITS.........	
!       3 minimal cases:  jumpstart  <= Xjumpscut <= jumpend
!                 jumpstart <= Xendcut <= jumpend
!                  Xjumpcut <= jumpstart  <= Xendcut
!
        do 	nex = 1,nexons
			
		    if( (jumpstart <= Xjumpcut(nex) .and. Xjumpcut(nex) <= jumpend ).or. & 
			   (jumpstart <= Xendcut(nex) .and. Xendcut(nex) <= jumpend ) .or. &
			   (jumpstart >= Xjumpcut(nex) .and. Xendcut(nex) >= jumpend ))then
			   containsexons=.true.
			   exit
		   end if 
		end do  ! nex
		
		if(containsexons)then  !check if another exon exists
			if(Xendcut(nex)<jumpend .and. Xjumpcut(nex+1) <=jumpend)then
				contiguous=.true.
			    print*,iproc,' WHOA non-contigous introns',nbody,it,isjmp, jumpstart,jumpend ,nex,nexons
				print*,Xjumpcut(:)
				print*,Xendcut(:)
 			end if
			xNjump%sjmp(isjmp)%containsexons=containsexons
			xNjump%sjmp(isjmp)%exstart      =bmax(jumpstart, Xjumpcut(nex))
			xNjump%sjmp(isjmp)%exstop       =bmin(jumpend, Xendcut(nex))
			xNjump%sjmp(isjmp)%exshift        = Xjumpshift(nex)
		else
			xNjump%sjmp(isjmp)%containsexons=containsexons
			xNjump%sjmp(isjmp)%exstart      = jumpstart
			xNjump%sjmp(isjmp)%exstop       = jumpend 
			xNjump%sjmp(isjmp)%exshift        = 0		
			
		end if
		
	end do   ! isjmp
	
	return
	
end subroutine surveysectorjumps
!===========================================================
