!
!  BXYBMELIB.f90
!
! library of routines to reduce storage of PN/PPN/PNN matrix elements
! uses "geneologies" to determine the matrix elements of required operators
!
! initiated in V7.3.6 Nov 2014 by CWJ @ SDSU
!
!  SOME NOTES ON STORAGE
!
!  Each operator X^+X or Y^+Y (this includes X^+X+ XX etc) have allowed dM,dpar,dW
!  For a given opbundle (which connects initial and final sectors with fixed quantum #s)
!  dMx,dMy, etc are fixed. 
!  We can in principle restrict further because depending on the initial and final 
!  sectors, only some operators may be available.  For "middle" sectors however
!  this is unlikely to be a big restriction so it is not clear how much one would gain.
!  
!  For dM/dpar, this is straightforward. The problem comes with dW. 
!  With different dWx, dWy, the indexing might become more complicated. 
!  (This is what defeated us on previous set ups.) The issue is that different dWy
!  and dWx are associated, making blocking difficult. 
!
!  What one can set up, however, is a simple opbundle-dependent transformation:
!
!  matrix element = xyme(indx)
!  where indx = xoplabel + a*yoplabel + b
!
!  where a and b depend on the opbundle. While this adds additional operations, 
!  it reduces lookup in arrays so ultimately should be faster.
!
!  The first step in testing this is to set up an indexing scheme and see how 
!  many matrix elements are required.  Need to make sure we only construct the
!  required operators.  We probably need to poll the geneologies for this 
!  information anyway.
!
!  Notes Jan 2015 (7.4.0):
!  With regards to operator construction: we first construct "grouples" which is 
!  operators labled by "groups" (which in BIGSTICK means creation/annhilation operators
!  labled by m,par,w, that is, things that can change Abelian quantum numbers)
!  While we can count up how many operators are in a grouple, this will lead to overcounting
!  because different geneologies can have the same or overlapping grouples.
!  Hence we need to eliminate duplicate grouples.
!
!  ADDITIONAL NOTES 1/12/2015:  May need to correlate the proton and neutron grouples;
!  We appear to be getting far too many matrix elements.   The solution is probably to use
!  ntuple structures. Still, it's not going to be easy.
!
!  NOTES IN VERSION 7.5.5 (Aug 2015): Revisiting this issue. Create more data structures.
!  For a given MPI node, loop over all opbundles to 
!    -- find all allowed combinations of dM,dpar,dW
!    -- within those, find all grouples; the problem may be that not all grouples can
!       connect to their conjugates, if we have very different starting Ws. Presumably in 
!       large calculations this is not a problem, but still have to check it out
!
! NOTES 7.5.7 (Oct 2015): need to either index with initial W (most likely) or initial M
!
!==================================================================================
!
! master routine
!
!

module bxybme_mod
contains

subroutine new_setup4XYmes


end subroutine new_setup4XYmes

!==================================================================================
!
! count up XY matrix elements for storage
!
! iprocs = process; can be physical process or can loop over processes for modeling
!
!
!  CALLED BY:
!
!  SUBROUTINES CALLED:
!    survey_geneologies
!    conjugate_grouple_W
!
subroutine new_countupXYmes(iprocs)

   use ntuple_info
   use butil_mod
   implicit none
!!    interface
!! 	   subroutine count_create_pairs_in_grouple(it,create,tmpwgtuple,npairs) ! INTERFACE
!! 	      use nodeinfo
!! 	      use spstate
!! 	      use ntuple_info
!! 
!! 	      implicit none
!! 
!! 	      integer it  !  species
!! 	      logical create
!! 	      type (wgtuple),target  :: tmpwgtuple	   
!! 	      integer :: npairs
!!    	   end subroutine count_create_pairs_in_grouple         ! INTERFACE
!!    end interface
   integer :: iprocs      ! "local" MPI process

   integer :: Mx,My, parx,pary, iWx,iWy,Wx,Wy   ! quantum numbers
   integer :: npairX, npairY,npair
   integer(8) :: nXYmes    ! local variable for # of proton-neutron matrix elements
   integer(8) :: maxlocal
   integer(8) :: partialmes
   
   nXYmes = 0
!...... FIND ALL THE POSSIBLE "GROUPLES" OF QUANTUM NUMBERS.........

   call survey_geneologies(1, 1,iprocs)  ! find proton one-body operators
   call survey_geneologies(2, 1,iprocs)  ! find neutron one-body operators
   call findMinMaxWdops(1,1)
   call findMinMaxWdops(2,1)
   call conjugate_grouple_W(1)
   maxlocal = 0
!...... LOOP OVER QUANTUM NUMBERS.....
   do Mx = XXgrouples(1)%jzmin, XXgrouples(1)%jzmax
      My = -Mx
      if( XXgrouples(2)%jzmin > My .or. XXgrouples(2)%jzmax < My)cycle
!      print*, MX,XXgrouples(1)%Mx(Mx)%parmin,XXgrouples(1)%Mx(Mx)%parmax
      do parx = XXgrouples(1)%Mx(Mx)%parmin,XXgrouples(1)%Mx(Mx)%parmax
         pary = parx
         if( XXgrouples(2)%Mx(my)%parmin > pary .or. XXgrouples(2)%Mx(my)%parmax < pary)cycle
         do iWx = 1,XXgrouples(1)%Mx(Mx)%par(parx)%ndWx
            Wx = XXgrouples(1)%Mx(Mx)%par(parx)%dWlist(iWx)
			npairX=0				
            call count_create_pairs_in_grouple(1,.false.,XXgrouples(1)%Mx(mx)%par(parx)%Wx(Wx),npairX)
!            print*,Mx,XXgrouples(1)%Mx(mx)%par(parx)%Wx(Wx)%ngops,XXgrouples(1)%Mx(mx)%par(parx)%Wx(Wx)%ngopsmax,npairX
!.............. FOR NOW ASSUME ALL CONJUGATE dWy ARE ALLOWED.....
              partialmes = 0
             do Wy = XXgrouples(1)%Mx(Mx)%par(parx)%Wx(Wx)%minWconj, & 
				 XXgrouples(1)%Mx(Mx)%par(parx)%Wx(Wx)%maxWconj
                 call count_create_pairs_in_grouple(1,.false.,XXgrouples(2)%Mx(my)%par(pary)%Wx(Wy),npairY)
!				 print*,npairX,npairY
                nXYmes = nXYmes + npairX*npairY
				partialmes=partialmes+npairX*npairY
!				if(npairy*npairx>0)print*,mx,parx,wx,wy,npairx*npairy
!                print*,Mx,Wx,npairX*npairY,nXYmes
                maxlocal = bmax(maxlocal,npairX*npairY)

             end do
! 			print*,Mx,parx,Wx,XXgrouples(1)%Mx(mx)%par(parx)%Wx(Wx)%minWdops,XXgrouples(1)%Mx(mx)%par(parx)%Wx(Wx)%maxWdops,partialmes
			 
!            do iWy = 1,XXgrouples(2)%Mx(My)%par(pary)%ndWx
!                Wy = XXgrouples(2)%Mx(My)%par(pary)%dWlist(iWy)
!                call count_create_pairs_in_grouple(1,.false.,XXgrouples(2)%Mx(my)%par(pary)%Wx(Wy),npairY)
!                nXYmes = nXYmes + npairX*npairY
!                print*,Mx,parx,npairX,npairY,nXYmes
!            end do  ! Wy
         end do ! Wx
      end do  ! parx
   end do ! Mx
   print*,' I expect ',nXYmes,' PN matrix elements '
   print*,' Max local = ',maxlocal
   return
end subroutine new_countupXYmes
!
!==================================================================================
!
!  a master routine to loop over geneologies and extract the quantum number of 
!  operators
!
!  CALLED BY:
!
! SUBROUTINES CALLED:
!   sector_descent
!   matchmaker
!   interpreter_of_geneologies
!   setup_grouples
!   get_delta_q
!   findsectorjump
!

subroutine survey_geneologies(it, nbody,iprocs)
   use spstate
   use nodeinfo
   use sectors
   use descendents
   use geneologies  
!   use jumpNbody 
   use spstate
   use ntuple_info
   use opbundles
   use jump_mod
   use butil_mod
   use bsector_mod
   use marriage
   implicit none

   integer :: it ! = species
   integer :: nbody ! = rank of operator
   integer :: iprocs  ! which MPI processor (either assigned or modeled)
   type (genlink), pointer :: curgen
   type (groupjump) :: descent
   integer(8) :: ngeneologies
   integer :: igen, ngeneologies4
   integer is,fs
   character, pointer :: optype(:)*1
   integer isjmp
   integer :: ibundle
   integer dw,dpar,dm,idW
   integer(8) :: ngops
   integer(8) :: ngopsall,ngopsmaxall
   integer :: aerr
   type (gtuple),pointer :: tmpgtuple
!   type (jumpsect), pointer :: xNbjump
   optype => NULL()

   select case (nbody)
    case(1)
      tmpgtuple => XXgrouples(it)
!      xNbjump => x1bjump(it)
    case(2)
      tmpgtuple => XXXXgrouples(it)
!      xNbjump => x2bjump(it)
    case default
      stop 2
   end select

   if(associated(optype))nullify(optype)
   if(.not.associated(optype)) then
      allocate(optype(2*nbody), stat=aerr)
      if(aerr /= 0) call memerror("survey_geneologies 1")
   end if

   select case(nbody)
      case (1)
      optype(1) = 'd'
      optype(2) = 'c'
   end select
!-------------------- FIGURE OUT WHICH SECTORS WE WANT TO GO TO
! NB: taken outside subroutine May 2010 by CWJ
!     in order to expedite 3-body interactions

   descent%nsectors = nsectors(it)
   allocate(descent%sector(nsectors(it)), stat=aerr)
   if(aerr /= 0) call memerror("survey_geneologies 2")
   call clocker('des','sta')
!$OMP PARALLEL DO SCHEDULE(STATIC)
   do is = 1,nsectors(it)
!--- CHECK WHETHER SECTORS BELONG TO THIS PROCESSOR --

      call sectordescent(it,is,nbody,descent%sector(is))
   enddo ! is
!$OMP END PARALLEL DO

   call setup_grouples(it,nbody,0)
   if(iproc==0)print*,' Grouples set up! '
   ngopsall = 0
   ngopsmaxall = 0
!--------- GO THROUGH AND STORE GROUPLES, CHECKING FOR DUPLICATES------
!   do isjmp = 1,x1bjump(it)%nsectjumps	   
!     is = x1bjump(it)%isector(isjmp)
!     fs = x1bjump(it)%fsector(isjmp)
   do ibundle=opbundlestart(iprocs),opbundleend(iprocs)
  	 select case (it)
  	   case (1)
  	      if(nbody==1 .and. opbundle(ibundle)%optype/='PN' .and.opbundle(ibundle)%optype/='PNN') cycle
  	      if(nbody==2 .and. opbundle(ibundle)%optype/='PPN' ) cycle
          is = opbundle(ibundle)%isector
          fs = opbundle(ibundle)%fsector	
  	   case(2)
          if(nbody==1 .and. opbundle(ibundle)%optype/='PN' .and.opbundle(ibundle)%optype/='PPN') cycle
          if(nbody==2 .and. opbundle(ibundle)%optype/='PNN' ) cycle
          is = opbundle(ibundle)%insector
          fs = opbundle(ibundle)%fnsector
   	 end select
     call get_delta_q(it,is,fs,dm,dpar,dW)
     call findsectorjump(it,nbody,is,fs,isjmp)
!--- CHECK WHETHER SECTORS BELONG TO THIS PROCESSOR --
     call matchmaker(it,nbody,is,descent%sector(is),nbody,fs,descent%sector(fs),.false.,  ngeneologies)
     if(ngeneologies == 0)cycle
     ngopsmaxall = ngopsmaxall+ngeneologies
     curgen =>geneology

     ngeneologies4 = bint8to4(ngeneologies, 301)  ! test size of ngeneologies
     do igen = 1, ngeneologies4
        call interpreter_of_geneologies(it,dm,dpar,dw,nbody,igen,curgen,isjmp)
        curgen => curgen%next
     end do  ! igen
!	 call Wsort_grouples(it,nbody  , tmpgtuple%Mx(dm)%par(dpar)%Wx(dw))
   end do  ! isjmp
!...................TALLY UP ALL GROUP OPERATORS................
   ngops = 0  
   do dM = tmpgtuple%jzmin,  tmpgtuple%jzmax
      do dpar = tmpgtuple%Mx(dm)%parmin,tmpgtuple%Mx(dm)%parmax
         do idW = 1,tmpgtuple%Mx(dm)%par(dpar)%ndWx
            dW = tmpgtuple%Mx(dm)%par(dpar)%dWlist(idW)
!			print*,dm, tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%ngops
			ngops = ngops + tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%ngops
		 end do
	  end do
   end do
   if(iproc==0)print*,' total # of geneologies/group operators = ',ngops,' out of ',ngopsmaxall
   return
end subroutine survey_geneologies
!==================================================================================
! ROUTINE findsectorjump moved to bjumplib_master.f90

!
!  set up "grouples" which are listings of group creation/destruction operators
!  by looping over sector jumps, constructs data for changes in M, parity, and W
!  and allocates memory for
!
!  CALLED BY
!    survey_geneologies
!
subroutine setup_grouples(it,nbody,iprocs)
   use sectors
   use nodeinfo
   use ntuple_info
!   use jumpNbody
   use opbundles
   use bsector_mod
   implicit none
   integer :: it  ! species
   integer :: nbody ! = 2 or 3 
   integer :: iprocs ! to be used to restrict on a particular MPI process

   integer :: isjmp  ! sector jump
   integer :: is,fs  ! initial,final sectors
   integer :: dM,dpar,dW

   integer :: mmin,mmax
   integer :: parmin,parmax
   integer :: m
   integer :: par
   integer :: wgmin,wgmax,W,nw
   logical, allocatable :: wglist(:)

   integer ibundle
   type (gtuple),pointer :: tmpgtuple
!   type (jumpsect), pointer :: xNbjump
   integer :: aerr

   select case (nbody)
    case(1)
      tmpgtuple => XXgrouples(it)
!      xNbjump => x1bjump(it)
    case(2)
      tmpgtuple => XXXXgrouples(it)
!      xNbjump => x2bjump(it)
    case default
      stop 2
   end select
   mmin = 10000
   mmax = -10000
!   do isjmp = 1,xNbjump%nsectjumps
!     is = xNbjump%isector(isjmp)
!     fs = xNbjump%fsector(isjmp)
!     call get_delta_q(it,is,fs,dm,dpar,dW)
!     mmin = MIN(mmin,dm)
!     mmax = MAX(mmax,dm)
!   end do
   do ibundle=opbundlestart(iprocs),opbundleend(iprocs)
	   select case (it)
	   case (1)
        is = opbundle(ibundle)%isector
        fs = opbundle(ibundle)%fsector	
	   case(2)
        is = opbundle(ibundle)%insector
        fs = opbundle(ibundle)%fnsector
 	 end select
     call get_delta_q(it,is,fs,dm,dpar,dW)
     mmin = MIN(mmin,dm)
     mmax = MAX(mmax,dm)
   end do   
   
   tmpgtuple%jzmin = mmin
   tmpgtuple%jzmax = mmax

   allocate( tmpgtuple%Mx(mmin:mmax) , stat=aerr)
   
   do M = mmin,mmax
      parmin =  2
      parmax = -2
!      do isjmp = 1,xNbjump%nsectjumps
!         is = xNbjump%isector(isjmp)
!         fs = xNbjump%fsector(isjmp)
      do ibundle=opbundlestart(iprocs),opbundleend(iprocs)
         select case (it)
         case (1)
            is = opbundle(ibundle)%isector
            fs = opbundle(ibundle)%fsector	
         case(2)
            is = opbundle(ibundle)%insector
            fs = opbundle(ibundle)%fnsector
         end select
         call get_delta_q(it,is,fs,dm,dpar,dW)
         if( dm /= m)cycle
         parmin = MIN(dpar,parmin)
         parmax = MAX(dpar,parmax)
      end do
      tmpgtuple%Mx(m)%parmin = parmin
      tmpgtuple%Mx(m)%parmax = parmax
      do par=parmin,parmax      
         Wgmin = 10000
         Wgmax = 0
!         do isjmp = 1,xNbjump%nsectjumps
!            is = xNbjump%isector(isjmp)
!            fs = xNbjump%fsector(isjmp)
         do ibundle=opbundlestart(iprocs),opbundleend(iprocs)
            select case (it)
            case (1)
               is = opbundle(ibundle)%isector
               fs = opbundle(ibundle)%fsector	
            case(2)
               is = opbundle(ibundle)%insector
               fs = opbundle(ibundle)%fnsector
            end select
            call get_delta_q(it,is,fs,dm,dpar,dW)
            if( dm /= m)cycle
            if(par/=dpar) cycle
            Wgmin = MIN(wgmin,dw)
            wgmax = MAX(wgmax,dw)
         end do
         tmpgtuple%Mx(m)%par(par)%maxdWx = wgmax
         tmpgtuple%Mx(m)%par(par)%mindWx = wgmin
         allocate( wglist(wgmin:wgmax) , stat=aerr)
         if(aerr /= 0) then
            call memerror("setup_grouples 1")
            stop 5
         end if
         wglist(:) = .false.
!         do isjmp = 1,xNbjump%nsectjumps
!            is = xNbjump%isector(isjmp)
!            fs = xNbjump%fsector(isjmp)
		 do ibundle=opbundlestart(iprocs),opbundleend(iprocs)
		 	select case (it)
		 	   case (1)
		         is = opbundle(ibundle)%isector
		         fs = opbundle(ibundle)%fsector	
		 	   case(2)
		         is = opbundle(ibundle)%insector
		         fs = opbundle(ibundle)%fnsector
		  	end select
            call get_delta_q(it,is,fs,dm,dpar,dW)
            if( dm /= m)cycle
            if(par/=dpar) cycle
            wglist(dw) = .true.
         end do
         Nw = 0
         do w = wgmin,wgmax
            if(wglist(w))nw = nw+1
         end do
         tmpgtuple%Mx(M)%par(par)%ndWx = nW
         allocate( tmpgtuple%Mx(M)%par(par)%dWlist(nw), stat=aerr)
         if(aerr /= 0) call memerror("setup_grouples 2")
         allocate( tmpgtuple%Mx(M)%par(par)%Wx(wgmin:wgmax), stat=aerr)
         if(aerr /= 0) call memerror("setup_grouples 3")
          Nw = 0
         do w = wgmin,wgmax
            tmpgtuple%Mx(M)%par(par)%Wx(W)%ngops = 0
            if(wglist(w))then
               nw = nw+1
               tmpgtuple%Mx(m)%par(par)%dWlist(nw) = w
            end if
         end do   

         deallocate(wglist)
      end do

   end do

end subroutine setup_grouples
!==================================================================================
!
! this subroutine goes through the proton quantum numbers and finds the allowed neutron conjugate dW
! ALSO: finds maximum sum of Wdops allowed
!
!
subroutine conjugate_grouple_W(nbody)

   use sectors
   use nodeinfo
   use ntuple_info
   use jumpNbody
   use bsector_mod
   implicit none

   integer :: nbody
   integer :: isjmp  ! sector jump
   integer :: signjmp  ! whether sector jumps is pos or neg
   integer :: is,fs  ! initial,final sectors
   integer :: dM,dpar,idW,dw
   integer :: ncjmps,cjmp,ic,cW
   integer :: dcm, dcpar, dcW
   integer :: aerr

   type (gtuple),pointer :: ptmpgtuple,ntmpgtuple
   type (jumpsect), pointer :: pNbjump,nNbjump
   type(sectorjump4grouples),pointer :: cursjmp,curNsj
   type (grouple), pointer :: curgroople, curNgr
   integer igp,ins,ign,iii
   integer pWdops,nWdops,Wdopsum
   logical :: foundjump

   select case (nbody)
    case(1)
      ptmpgtuple => XXgrouples(1)
      pNbjump => x1bjump(1)
      ntmpgtuple => XXgrouples(2)
      nNbjump => x1bjump(2)
    case(2)
      ptmpgtuple => XXXXgrouples(1)
      pNbjump => x2bjump(1)
      ntmpgtuple => XXXXgrouples(2)
      nNbjump => x2bjump(2)
    case default
      print *, "conjugate_grouple_W: bad nbody"
      stop 1
   end select

   do dM = ptmpgtuple%jzmin,  ptmpgtuple%jzmax
      do dpar = ptmpgtuple%Mx(dm)%parmin,ptmpgtuple%Mx(dm)%parmax
         do idW = 1,ptmpgtuple%Mx(dm)%par(dpar)%ndWx
            dW = ptmpgtuple%Mx(dm)%par(dpar)%dWlist(idW)
			ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%minWconj=10000
			ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%maxWconj=-10000
			ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%maxWdopsum = 0
            curgroople => ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%groop
		    pWdops = curgroople%Wdops	 
			
			do igp = 1, ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%ngops
				
               cursjmp=> curgroople%sj
			   do ins = 1,curgroople%nsjmps
                  isjmp = cursjmp%sjmp    
				  if( isjmp < 0)then
					  signjmp = -1
					  isjmp = -isjmp
				  else
					  signjmp = 1
				  end if
!---------------- OBTAIN CONJUGATE JUMPS
                  ncjmps = pNbjump%csjmp(isjmp)%ncjmps
                  do ic = 1,ncjmps
                     cjmp = pNbjump%csjmp(isjmp)%cjump(ic)
                     is = nNbjump%isector(cjmp)
                     fs = nNbjump%fsector(cjmp)
                     call get_delta_q(2,is,fs,dcm,dcpar,dcW)
					 if(dcm/=-signjmp*dm .or. dcpar/=dpar)then   ! ERROR TRAP
						 print*,nbody,' wrong conjugate ',dm,dcm,dpar,dcpar
						 print*,signjmp*isjmp,cjmp
						 stop
					 end if
				!	 print*,' Okay ',dm,dcm,isjmp,cjmp
                     ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%minWconj= & 
			             MIN(dcW,ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%minWconj)
                     ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%maxWconj= & 
   			             MAX(dcW,ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%maxWconj)
!------------------ SEARCH THROUGH NEUTRON GROUPLES FOR CONJUGATE Wdop
!                  Wdops is the W value on the destruction operator(s).
!                  In order to further restrict storage of matrix elements,
!                  find the sum of Wdop for protons and neutrons

                     curNgr => ntmpgtuple%Mx(dcm)%par(dcpar)%Wx(dcW)%groop
					 foundjump=.false.
					 do ign = 1,	 ntmpgtuple%Mx(dcm)%par(dcpar)%Wx(dcW)%ngops  ! loop over neutron grouples
						 curNsj => curNgr%sj
						 do iii = 1,curNgr%nsjmps   ! for each neutron groople,loop through associated sector jumps
!							 print*,curNsj%sjmp
							 if(curNsj%sjmp==cjmp)then
								 foundjump=.true.
								 nWdops = curNgr%Wdops
								 exit
							 end if
							 curNsj => curNsj%next
						 end do
						 if(foundjump)exit
						 curNgr=>curNgr%next						 
					 end do ! ign
					 if(.not.foundjump)then
						 print*,' could not find sector jumps '
						 print*,dm,dcm
						 print*,isjmp
						 print*, ntmpgtuple%Mx(dcm)%par(dcpar)%Wx(dcW)%ngops 
						 stop
					 end if
				     ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%maxWdopsum =  & 
					 MAX( ptmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%maxWdopsum, pWdops+nWdops)
						 
                   end do
				   cursjmp=>cursjmp%next
				end do
			   curgroople=>curgroople%next
		    end do
         end do ! idW
      end do  ! dpar
   end do  ! dm

   return
end subroutine conjugate_grouple_W

!==================================================================================
!
!  this subroutine looks at a "geneology" and extracts information 
!  on the operators needed
!
!  from a given n-body geneology, we extract the "group" (i.e. quantum numbers M,par,W)
!  of creation and destruction operators
!
!  ADDED IN 7.4.0: Check that grouples are unique
!  ADDED IN 7.5.9: rewrote so that we have linked lists of grouples, so that we don't have to 
!                  count ahead
!
!  NOTES: In a given geneology (e.g., curgen), curgen%gen(i) lists the blockhop for the 
!         operation for the ith operator; if < 0 then left blockhop, if > 0 then right blockhop
!         lists first the annihilation blockhops (on the initial sector)
!         then the creation blockhops (to get to the final sector)
!
!         To get the groups for the destruction operators, 
!         = bunny(it* sgn( curgen%gen(1:nbody)  )%dblockhop( abs(curgen%gen(1:nbody)))%group
!                and for creation operators
!         = bunny(it* sgn( curgen%gen(nbody+1:2*nbody)  )%cblockhop( abs(curgen%gen(1+nbody:2*nbody)))%group
!         where the derived type bunny is found in module bugs. 
!         These group labels are all positive; to distinguish between left-handed (jz < 0) and 
!         right-hand (jz >= 0) groups, add in a sign
!
!  CALLED BY
!     survey_geneologies
!
!  SUBROUTINES CALLED
!
subroutine interpreter_of_geneologies(it,m,par,w,nbody,pastgen,curgen,isjmp)
   use nodeinfo
   use sectors
   use descendents
   use geneologies   
   use hoppy
   use blocks
   use spstate
   use ntuple_info
   implicit none
   integer :: it ! = species
   integer :: m,par,w !  quantum numbers
   integer :: nbody ! = rank of operator
   integer :: isjmp  ! sector jump
   integer :: pastgen  ! which geneologies NOT NEEDED????
   type (genlink), pointer :: curgen
   type (grouple), pointer :: curgroople
   type(sectorjump4grouples),pointer :: cursjmp
   integer :: igroop
   integer :: cgrops(nbody),dgrops(nbody)   ! temporary arrays for storing creation and destruction operators
   integer ib,jb   ! counting up creation/destruction operators
   integer handsignd,handsignc
   type (wgtuple),pointer :: tmpwgtuple
   logical duplicate     ! flag for duplicate grouples
   integer :: Wdops
   integer :: is

   select case (nbody)
    case(1)
      tmpwgtuple => XXgrouples(it)%Mx(m)%par(par)%Wx(w)
    case(2)
      tmpwgtuple => XXXXgrouples(it)%Mx(m)%par(par)%Wx(w)
    case default
      print *, "interpreter_of_geneologies: bad nbody"
      stop 1
   end select
   Wdops = 0
   do ib = 1,nbody   !loop over ngen destruction operators
    handsignd = 1
    if(curgen%gen(ib) < 0)handsignd= -1
    dgrops(ib) = handsignd*bunny(it*handsignd)%dblockhop( abs(curgen%gen(ib)))%group
! ....... SET W FOR THIS/THESE OPERATORS.......
     Wdops = Wdops +group(handsignd*it)%W(handsignd* dgrops(ib)) 
   end do
   do ib = nbody+1,2*nbody  ! loop over creation operators
    handsignc = 1
    if(curgen%gen(ib) < 0)handsignc= -1
    cgrops(ib-nbody) = handsignc*bunny(it*handsignc)%cblockhop( abs(curgen%gen(ib)))%group
   end do
 
!------ SEARCH THROUGH LINKED LIST TO SEE IF GROUPLE ALREADY EXISTS ------------
   duplicate = .false.
   curgroople => tmpwgtuple%groop    ! start at beginning of linked list and search through it
   if(tmpwgtuple%ngops>0)then
      do igroop = 1, tmpwgtuple%ngops 
		 if(igroop > 1)curgroople => curgroople%next	  
	     duplicate = .true.
	     do ib = 1,nbody
	           if ( curgroople%dgrops(ib)/= dgrops(ib) & 
	           .or. curgroople%cgrops(ib)/= cgrops(ib))duplicate=.false.
	     end do ! igen
	     if(duplicate)then
	             exit
	     end if
		 duplicate = .false.
      end do ! igroop
   end if
   if(.not.duplicate)then   ! make new
	   tmpwgtuple%ngops = tmpwgtuple%ngops + 1       
	   if(tmpwgtuple%ngops /=1)then
		   if(associated(curgroople%next))nullify(curgroople%next)
  	       allocate(curgroople%next)
		   nullify(curgroople%next%next)
	       curgroople => curgroople%next
	   end if
	   allocate(curgroople%cgrops(nbody),curgroople%dgrops(nbody))
	   do ib = 1,nbody
		   curgroople%dgrops(ib)=dgrops(ib)
		   curgroople%cgrops(ib)=cgrops(ib)
		   curgroople%Wdops = Wdops
	   end do
	   curgroople%nsjmps  = 1
	   curgroople%sj%sjmp = isjmp
	   
   else     ! UPDATE INFO ON SECTOR JUMPS
	   cursjmp => curgroople%sj
	   do is = 1,curgroople%nsjmps-1
		   cursjmp => cursjmp%next
	   end do ! is
	   allocate(cursjmp%next)
	   cursjmp => cursjmp%next
	   cursjmp%sjmp = isjmp
	   curgroople%nsjmps = curgroople%nsjmps+1
   end if

   return
end subroutine interpreter_of_geneologies
!==================================================================================
!
!  subroutine to count up pair operators (X^+X) in a grouple
!  CALLED BY
!  
!  SUBROUTINES CALLED:
!   Wsort_grouples
subroutine count_create_pairs_in_grouple(it,create,tmpwgtuple,npairs)
   use nodeinfo
   use spstate
   use ntuple_info

   implicit none

   integer it  !  species
   logical create
   type (wgtuple),target  :: tmpwgtuple	   
   integer :: npairs

   type (grouple), pointer :: curgroople

   integer n
   integer handsignd,handsignc
   integer :: dgroup,cgroup,ndgroup,ncgroup
   npairs = 0
   curgroople => tmpwgtuple%groop
   do n = 1, tmpwgtuple%ngops
      cgroup = curgroople%cgrops(1)
      dgroup = curgroople%dgrops(1)
      if(cgroup < 0)then
         handsignc = -1
      else
         handsignc = 1
      end if
      if(dgroup < 0)then
         handsignd = -1
      else
         handsignd = 1
      end if
      if(.not.create)then  ! just count up
       
         ncgroup = group(it*handsignc)%fin(abs(cgroup)) - group(it*handsignc)%start(abs(cgroup))+1
         ndgroup = group(it*handsignd)%fin(abs(dgroup)) - group(it*handsignd)%start(abs(dgroup))+1
         if(ncgroup < 1 .or. ndgroup < 1)then
           print*,' no elements in group ',ncgroup,ndgroup
           stop
         end if
         npairs = npairs+ncgroup*ndgroup
         if(it==1)then 
!              print*,n,' c d , pairs ',cgroup,ncgroup,dgroup,ndgroup
         end if
      else  ! create the pairs, unleash the kraken, etc etc etc
 
      end if
	  curgroople=>curgroople%next
   end do
   return
end subroutine count_create_pairs_in_grouple
!==================================================================================
!
! goes through and finds min, max W on the destruction operator for grouples;
! this will be necessary for more restricted storage
! 
! CALLED BY:
!

subroutine findMinMaxWdops(it,nbody)
	use nodeinfo
	use ntuple_info
	implicit none
	integer :: it   ! species
	integer :: nbody ! particle rank of operator (1, 2 body...)
!...............	

   integer :: isjmp  ! sector jump
   integer :: is,fs  ! initial,final sectors
   integer :: dM,dpar,idW,dw
   integer :: ncjmps,cjmp,ic,cW
   integer :: dcm, dcpar, dcW
   integer :: aerr

   type (gtuple),pointer :: tmpgtuple
   type(sectorjump4grouples),pointer :: cursjmp
   type (grouple), pointer :: curgroople
   integer igp,ins

   select case (nbody)
    case(1)
      tmpgtuple => XXgrouples(it)
    case(2)
      tmpgtuple => XXXXgrouples(it)
    case default
      print *, "findMinMaxWdops case"
      stop 1
   end select

   do dM = tmpgtuple%jzmin,  tmpgtuple%jzmax
      do dpar = tmpgtuple%Mx(dm)%parmin,tmpgtuple%Mx(dm)%parmax
         do idW = 1,tmpgtuple%Mx(dm)%par(dpar)%ndWx
            dW = tmpgtuple%Mx(dm)%par(dpar)%dWlist(idW)
			tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%minWdops=10000
			tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%maxWdops=-10000
            curgroople => tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%groop
			do igp = 1,tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%ngops

                tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%minWdops= & 
			             MIN(curgroople%Wdops,tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%minWdops)
                tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%maxWdops= & 
   			             MAX(curgroople%Wdops,tmpgtuple%Mx(dm)%par(dpar)%Wx(dW)%maxWdops)
			   curgroople=>curgroople%next
		    end do
         end do ! idW
      end do  ! dpar
   end do  ! dm	
	
    return
end subroutine findMinMaxWdops

!==================================================================================
!
! sorts grouples on their creation W; simple bubble sort
!
!  CALLED BY: 
!     count_create_pairs_in_grouple
!
subroutine Wsort_grouples(it,nbody,tmpwgtuple)
   use sectors
   use nodeinfo
   use ntuple_info
   use spstate
   implicit none

   integer :: it   ! species
   integer :: nbody
   type (wgtuple) :: tmpwgtuple

   integer :: isjmp  ! sector jump
   integer :: is,fs  ! initial,final sectors
   integer :: dM,dpar,idW,dw
   integer :: ncjmps,cjmp,ic,cW
   integer :: dcm, dcpar, dcW

   integer, allocatable :: wlist(:)

   integer :: i,igop,igroup,isgn,jgop
   integer wnow,iplace
   integer :: itmp
   integer :: aerr



!   write(6,'(15i4)')wlist
!......... THE FINAL STEP IS TO FIGURE OUT A MIN AND MAX

   return
end subroutine Wsort_grouples


!==================================================================================
end module bxybme_mod
