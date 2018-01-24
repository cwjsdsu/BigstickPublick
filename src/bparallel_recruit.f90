!
!  New MPI distribution routines -- started in 7.6.2 (Jan 2016) -- CWJ @ SDSU
!
!  Old version of routines attempted to divide up work.  In some cases, however,
!  some opbundles required much more jump storage than others. 
!  (This is understood: it is because of "narrow" opbundles.  Every opbundle 
!   loops over protons and neutrons in some way, either jumps or states, and 
!   if the dimensions of one loop is particularly small, it is a narrow opbundle.
!   For example, suppose for NN operations, the # of conjugate proton Slater 
!   determinants to loop over is only 1 or 2.  This is not a particularly 
!   efficient use of factorization, but a case we sometimes have to deal with.)
!  The old algorithm divided up processes evenly at first, and then attempted to 
!  readjust. In the case of very narrow opbundles, this algorithm fails. In
!  particular, the user could ask for too few MPI processes.
!
!  Instead we make a recruitment model. This has the advantage of determining 
!  up front the minimum number of MPI processes required.
!
!  The steps are:
!  -- Determine number of fragments (done prior to entering these routines)
!  -- Count up operations and jump storage for each fragment-to-fragment (f2f)
!     matvec block
!  ->  NOTE: # of operations for f2f is computed in routine getfragmentstats in bopstat.f90 <-
!  -- From jump storage limits, determine minimum # of MPI procs required for 
!     each f2f block.
!  -- Ask user for actual # of MPI procs
!  -- Recruit MPI procs to f2fs to best load balance
!

module recruiter
	
	implicit none
	
!.......... DATA USED TO SET UP DISTRIBUTION.......
    

    	
!........... NOW SUBROUTINES.......................
	
   contains
   
!============================================================================   

!subroutine masterrecruitMPIprocs

!end subroutine masterrecruitMPIprocs

!============================================================================
!
!  count up the # the amount of jump storage 
!  in an f2f (fragment-to-fragment) matvec block
!
!  This is done by figuring out what initial and final sectors, and 
!  thus what sector jumps, belong to an f2f block
!  and then look up the # of jumps in each sector block
!
!  
subroutine f2fjumpstats
	use fragments
	use flags3body
	use jumpdef
	use jumpNbody
	use jump3body
	use nodeinfo
	use io
	
	implicit none
	integer :: ifrag,ffrag   ! initial and final fragments
	logical, allocatable :: thissectorjump(:)    ! logical array of sector jumps to count up
	integer(4) :: maxsectorjumps
	integer(8) :: isstart,isstop,fsstart,fsstop ! sector start/stop
	integer(8) :: is,fs
	integer(8) :: isjmp    ! sector jump label
	integer(8), pointer :: mynjumps
	integer it      ! species
	integer nbody   ! whether 1,2, or 3-body jumps
	type (jumpsect), pointer :: xNbjump
		integer :: aerr
		
	integer(8) :: totjumpstore,sumjumpstore    ! # of bytes for each jump; allows for different amount of memory...
	
	logical :: localverbose = .true.
	
!	if(localverbose .and. iproc == 0) print*,' In f2fjumpstats'
!print*,onebodyonly,' one body only '

    if(.not.allocated(f2fjumpstat))then
			allocate(f2fjumpstat(nfragments,nfragments),stat=aerr)
			if(aerr/=0)       call memerror("f2fjumpstats ")

				
	end if
			
	f2fjumpstat(:,:)%nX1bjump(1)=0
	f2fjumpstat(:,:)%nXX2bjump(1)=0
	f2fjumpstat(:,:)%nXXX3bjump(1)=0
	f2fjumpstat(:,:)%nX1bjump(2)=0
	f2fjumpstat(:,:)%nXX2bjump(2)=0
	f2fjumpstat(:,:)%nXXX3bjump(2)=0	
	
!.............. SET UP TEMPORARY LOGICAL ARRAY OF INCLUDED SECTOR JUMPS...............	

    maxsectorjumps = max(x1bjump(1)%nsectjumps,x1bjump(2)%nsectjumps)
    maxsectorjumps = max(maxsectorjumps,x2bjump(1)%nsectjumps )
    maxsectorjumps = max(maxsectorjumps,x2bjump(2)%nsectjumps )
	if(threebody)then
	    maxsectorjumps = max(maxsectorjumps,x3bjump(1)%nsectjumps )
	    maxsectorjumps = max(maxsectorjumps,x3bjump(2)%nsectjumps )		
	end if

    allocate(thissectorjump(maxsectorjumps))
	sumjumpstore = 0

!.................. LOOP OVER f2f MATVEC BLOCKS.....................................
	do ifrag = 1,nfragments
		do ffrag = ifrag,nfragments
			
			do Nbody = 1,3
				if(Nbody==2  .and. onebodyonly) cycle
				if(Nbody==3 .and. .not. threebody)cycle
				do it = 1,2
!................. LOOP OVER INITIAL PROTON SECTORS BELONGING TO INITIAL FRAGMENT

                   thissectorjump(:)=.false.
			 	   if(it==1)then
					   isstart = fragmentlist(ifrag)%ssectorstart
					   isstop  = fragmentlist(ifrag)%ssectorend
					   fsstart = fragmentlist(ffrag)%ssectorstart
					   fsstop  = fragmentlist(ffrag)%ssectorend
				   else
					   isstart = fragmentlist(ifrag)%csectorstart
					   isstop  = fragmentlist(ifrag)%csectorend
					   fsstart = fragmentlist(ffrag)%csectorstart
					   fsstop  = fragmentlist(ffrag)%csectorend					   
				   end if
				   
				   select case (Nbody)
					   
				      case (1) 
				         xNbjump => x1bjump(it)
						 f2fjumpstat(ifrag,ffrag)%nX1bjump(it) = 0
						 mynjumps => f2fjumpstat(ifrag,ffrag)%nX1bjump(it)
   				      case (2) 
   				         xNbjump => x2bjump(it)	
						 f2fjumpstat(ifrag,ffrag)%nXX2bjump(it) = 0
						 mynjumps => f2fjumpstat(ifrag,ffrag)%nXX2bjump(it)		
   				      case (3) 
   				         xNbjump => x3bjump(it)
						 f2fjumpstat(ifrag,ffrag)%nXXX3bjump(it) = 0
						 mynjumps => f2fjumpstat(ifrag,ffrag)%nXXX3bjump(it)
				   end select
			
                   do is = isstart,isstop
			
!................. LOOP OVER FINAL SECTORS BELONGING TO FINAL FRAGMENT

                      do fs = fsstart,fsstop
!.................. find sector jumps.........................................
                         do isjmp = 1,xNbjump%nsectjumps
						    if( (is== xNbjump%isector(isjmp) .and. fs == xNbjump%fsector(isjmp)) .or.  & 
						       (is== xNbjump%fsector(isjmp) .and. fs == xNbjump%isector(isjmp)) )then
						          thissectorjump(isjmp)=.true.
					        end if
	
                         end do ! isjmp
                      end do  ! fs
                  end do  ! is
!..................... NOW COUNT UP JUMPS FOR ASSOCIATED SECTOR JUMPS............................	
                  do isjmp = 1,xNbjump%nsectjumps
			          if(thissectorjump(isjmp))then
						  mynjumps = mynjumps + xNbjump%sjmp(isjmp)%njumps
					  end if
				  end do  !isjump
		     end do ! it
		  end do ! Nbody
!..................... COMPUTE TOTAL JUMP STORAGE......................................
          totjumpstore = 0
		  totjumpstore = bytesper1Bjump*(f2fjumpstat(ifrag,ffrag)%nX1bjump(1)+f2fjumpstat(ifrag,ffrag)%nX1bjump(2) )
		  totjumpstore = totjumpstore+ bytesper2Bjump*(f2fjumpstat(ifrag,ffrag)%nXX2bjump(1)+f2fjumpstat(ifrag,ffrag)%nXX2bjump(2) )
		  if(threebody)then
			  totjumpstore = totjumpstore+ bytesper3Bjump*(f2fjumpstat(ifrag,ffrag)%nXXX3bjump(1)+f2fjumpstat(ifrag,ffrag)%nXXX3bjump(2) )
		  end if
          f2fjumpstat(ifrag,ffrag)%totjumpstorage = totjumpstore
          f2fjumpstat(ffrag,ifrag)%totjumpstorage = f2fjumpstat(ifrag,ffrag)%totjumpstorage
		  sumjumpstore = sumjumpstore + totjumpstore
		  if(ifrag /= ffrag) sumjumpstore = sumjumpstore + totjumpstore  ! take care of off-diagonal
!......................................................................................		  
			
		end do ! ffrag
		
	end do ! ifrag
    if(iproc==0)write(logfile,*)' Total storage of jumps = ',	sumjumpstore*1.0e-9,' Gb '
!----- (OPTIONAL ) printout

    if(iproc==0 .and. localverbose)then
		print*,' Initial analysis of jump storage: '
		print*,' Total storage of jumps = ',sumjumpstore*1.0e-9,' Gb '
	end if
!...........	
	
	deallocate(thissectorjump)
	return
end subroutine f2fjumpstats
!============================================================================
!
! routine to determine the minimal # of MPI processes/nodes required
! for each "f2f" (fragment-to-fragment) set of operations;
! uses information f2fjumpstat(i,f)%totjumpstorage to determine this
!  as well as "maxjumpmemory"
!
subroutine minprocs4f2f
	use fragments
	use nodeinfo
	use flagger
	use jumplimits
	use io
	implicit none
	
	integer(4) :: jprocs
	integer(4) :: ifrag,ffrag
	integer(4) :: nproctmp,nprocsum
	
	if(iproc==0)then
		print*,' Maximum memory for storage of jumps = ',maxjumpmemory_default  ! NOTE THIS MIGHT CHANGE
	end if
	
	nprocsum = 0
	if(iproc==0)write(logfile,*)' Information on jump storage between fragments  (Gb)'
	do ifrag = 1,nfragments
		if(iproc==0)write(logfile,*)' Fragment ',ifrag,' -> : '
		write(logfile,'(5(1i5,1f12.3))')(ffrag,(f2fjumpstat(ifrag,ffrag)%totjumpstorage*1.e-9),ffrag=ifrag,nfragments)
		do ffrag = 1,nfragments
			nproctmp = int(f2fjumpstat(ifrag,ffrag)%totjumpstorage*4e-9/maxjumpmemory_default)+1
			f2fjumpstat(ifrag,ffrag)%minprocs = nproctmp
			nprocsum = nprocsum+nproctmp
			if(iproc==0)print*,' fragment-to-fragment ',ifrag,ffrag,' minimum of ',nproctmp,' MPI procs ',& 
			               f2fjumpstat(ifrag,ffrag)%totjumpstorage
		end do  ! ffrag
	end do  ! ifrag
	if(iproc==0)print*,' Need a minimum of ',nprocsum,' MPI processes '
	
	return
end subroutine minprocs4f2f	

!============================================================================
!
! DRAFT (7.6.8) routine to allocate MPI procs to frag-to-frag blocks;
! used initially for applying one-body
! Distribution based upon jump storage, not work efficiency
!
! AFTER this routine is called, call routine "setnodaltribution"
!
! CALLED BY:
!   master_para_distribute_onebody
!

subroutine setnodal4applyonebody_DRAFT
	
use nodeinfo
use fragments
use io
use onebodypot, only: meanie

implicit none
	
integer :: ifrag, ffrag

real(8)    :: sumjumpstore,avgjumpstore
integer  :: sumprocs,tempprocs,jprocs
integer  :: procstart,procend


if(.not.allocated(opfragstat))then
	allocate(opfragstat(nfragments,nfragments))
end if

!......... CHECK IF ALREADY DISTRIBUTED..........
!
!  IF ALREADY DISTRIBUTED THEN EXTRACT DISTRIBUTION FROM nodal
!
if(.not.auto_readin .and. .not. meanie)then
	if(.not.allocated(nodal))then
		print*,' WARNING! Some error, nodal has not been allocated '
		print*,' (I am here in setnodal4applyonebody_DRAFT)'
		stop
	end if
	
	if(nfragments == 1)then
		opfragstat(1,1)%nnodes=nproc
		return
	end if
	do ifrag = 1,nfragments
		do ffrag = 1,nfragments
			tempprocs = 0
			procstart = nprocs -1
			procend   = -1
			do jprocs = 0,nprocs-1
				if(nodal(jprocs)%ifragment==ifrag .and. nodal(jprocs)%ffragment==ffrag)then
					tempprocs = tempprocs+1
					procstart = min(procstart,jprocs)
					procend = max(procend,jprocs)
					
				end if
			end do
			opfragstat(ifrag,ffrag)%nnodes = tempprocs
!............... CHECK CONTIGUITY.............................
            do jprocs = procstart,procend
				if(nodal(jprocs)%ifragment/=ifrag .or. nodal(jprocs)%ffragment/=ffrag)then
					if(iproc==0)then
						print*,' WHOOPS ! Non-contiguous MPI proc for fragments  '
						print*,' for fragments ',ifrag,ffrag,' start start/stop ',ifrag,ffrag
						print*,' but on proc ',jprocs,' it is ', nodal(jprocs)%ifragment,nodal(jprocs)%ffragment
					end if
				end if
            end do

!........................................................			
			
		end do ! ffrag
	
	end do ! ifrag
	return       ! exit
end if

!............ OTHERWISE DISTRIBUTE NODES................

sumjumpstore = 0
do ifrag = 1,nfragments
	do ffrag = 1,nfragments
		sumjumpstore = sumjumpstore +  f2fjumpstat(ifrag,ffrag)%totjumpstorage
	end do ! ffrag
end do ! ifrag

avgjumpstore = sumjumpstore/float(nprocs)
if(iproc == 0)then
	print*,' total jumpstorage = ',sumjumpstore,' over ',nprocs,' processors '
	print*,' Avergage jumpstore per proc = ',avgjumpstore
end if

1001 continue

sumprocs = 0

do ifrag = 1,nfragments
	do ffrag = ifrag,nfragments
		if(f2fjumpstat(ifrag,ffrag)%totjumpstorage==0)then
			opfragstat(ifrag,ffrag)%nnodes=0
			opfragstat(ffrag,ifrag)%nnodes=0
			cycle
		end if
		tempprocs = nint(f2fjumpstat(ifrag,ffrag)%totjumpstorage/ avgjumpstore)
		if(tempprocs == 0)tempprocs = 1
		opfragstat(ifrag,ffrag)%nnodes = tempprocs
		sumprocs = sumprocs + tempprocs
		if(ifrag/=ffrag)then
     		opfragstat(ffrag,ifrag)%nnodes = tempprocs
			sumprocs = sumprocs + tempprocs
		end if
		
	end do
	
end do 
if(sumprocs /= nprocs) then ! MUST REBALANCE
	if(iproc == 0)then
		print*,' (crude) Rebalancing distribution ',nprocs, sumprocs
	end if
	
	avgjumpstore = avgjumpstore*((float(nprocs)+0.3*(sumprocs-nprocs))/float(nprocs))
	
   goto 1001	
	
end if


return

end subroutine setnodal4applyonebody_DRAFT

!============================================================================
!
! DRAFT (7.6.8) routine to allocate opbundles to MPI procs
! assuming initial and final fragments have already been assigned
! used initially for applying one-body
! Distribution based upon jump storage, not work efficiency
!


!============================================================================

end module recruiter