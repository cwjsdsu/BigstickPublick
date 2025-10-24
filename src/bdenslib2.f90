!==================================
!
! file BDENSLIB2.f90
!
!  SUPPORT ROUTINES for one-body densities and applying one-body operators
!  count up operations, set up opbundles, and distribute work (if in MPI)
!
!
!==============================================================


!
! subroutine master_op_stat_den
!
! master subroutine for tracking statistics 
!        on operations between sectors
!
! started 8/2011 by CWJ @ SDSU
!
!  CALLED BY
!    density1b_output
!    density1b_from_oldwfn
!
!  SUBROUTINES CALLED:
!    count_Xop_stats
!    count_total_ops
!    subsectoropstats -- slated for obsolescence
!    getfragmentstats

! May 2020: THIS ROUTINE SEEMS TO NOT FILL any ops
subroutine master_op_stat_den

  use operation_stats
  use sectors
!  use flags3body
  implicit none

  interface
     subroutine count_Xop_stats(itx, is,fs,nops)   
     implicit none
     integer    :: itx
     integer    :: is,fs
     integer(8) :: nops
     end subroutine count_Xop_stats

  end interface

  integer is,fs
  integer(8) nops
  integer :: aerr

  if(.not.allocated(opstat)) then
     allocate ( opstat( nsectors(1), nsectors(1) ) , stat=aerr )
     if(aerr /= 0) call memerror("master_op_stat_den")
  end if
  opstat(:,:)%nopP = 0
  opstat(:,:)%nopN = 0
  opstat(:,:)%nopPP = 0
  opstat(:,:)%nopNN = 0
  opstat(:,:)%nopPPP = 0
  opstat(:,:)%nopNNN = 0
  opstat(:,:)%nopPN = 0
  opstat(:,:)%nopPPN = 0
  opstat(:,:)%nopPNN = 0  
  opstat(:,:)%nopSPE = 0

  do is = 1,nsectors(1)
     do fs = 1,is    ! this ordering is mandatory in bigstick

        call count_Xop_stats(1,is,fs,nops)
        opstat(is,fs)%nopP = nops
        opstat(fs,is)%nopP = nops

        call count_Xop_stats(2,is,fs,nops)
        opstat(is,fs)%nopN =  nops
        opstat(fs,is)%nopN =  nops
		
     end do
  end do
  call count_total_ops
!  call subsectoropstats  ! MAY GO AWAY
!  call getfragmentstats1body


  return

end subroutine master_op_stat_den

!===============================================================
!
! subroutine count_Xop_stats
!
! counts operations between two sectors
!
! INPUT:
!   is = initial (proton) sector 
!   fs = final (proton) sector
!  itx = species (P or N )
!
! OUTPUT:

!   nops  = # of noperations
!
subroutine count_Xop_stats(itx, is,fs,nops)

  use sectors
  use jumpNbody
  use system_parameters

  implicit none
  integer    :: itx
  integer    :: is,fs
  integer(8) :: nops

  integer    :: xsj  ! sector jump
  integer    :: ity
  integer(8) :: ncsd    ! # of conjugate SDs
  integer    :: csj    ! conjugate sector(jumps)
  integer    :: ic
  integer    :: cs
  integer(8) :: njumps

  ity = 3 -itx

  nops = 0
  if(np(itx) < 2)return

!  if(itx == 2 .and. is /= fs)return   ! NN can't change proton sector

!..... LOOP OVER SECTOR JUMPS......... 

  do xsj = 1,x2bjump(itx)%nsectjumps
     if ( itx == 1)then   ! p -sector jumps
         if(    is == x1bjump(itx)%isector(xsj) .and.  & 
                fs == x1bjump(itx)%fsector(xsj) ) then
            njumps = x1bjump(itx)%sjmp(xsj)%njumps
            ncsd = 0

!.............COUNT UP CONJUGATE SPECTATOR SDs

            do ic = 1,x1bjump(itx)%csjmp(xsj)%ncjmps  ! loop over conjugate sectors
                 cs = x1bjump(itx)%csjmp(xsj)%cjump(ic)
                 ncsd = ncsd + xsd(ity)%sector(cs)%nxsd
            end do
            nops = nops + njumps*ncsd
         endif

     else   ! N sector-jumps
         do ic = 1,x1bjump(itx)%csjmp(xsj)%ncjmps
            cs = x1bjump(itx)%csjmp(xsj)%cjump(ic)
            njumps = x1bjump(itx)%sjmp(xsj)%njumps

            if(cs == is .and. cs == fs) then
              ncsd =xsd(ity)%sector(cs)%nxsd
              nops = nops + njumps*ncsd
            end if

         end do

     end if

  end do ! xsj

  return

end subroutine count_Xop_stats
!=================================================================
!
!  master routine for parallel distribution for work for one-body density matrices,
!  apply one-body operators; added in 7.6.8 (July 2016)
!
!  CALLED BY:
!     density1b_output
!     applicator1b
!
!  SUBROUTINES CALLED:
!     master_op_stat_den
!     density_bundles_setup
!     f2fjumpstats
!     setnodal4applyonebody_DRAFT  :  assign NUMBER of MPI procs/nodes for each set 
!                                   of initial/final fragments  (in variable opfragstat)
!     distro_opbundles_over_fragments_onebody  : set up initial/final starting opbundles for each MPI process
!     setnodaltribution           :  using informatin from opfragstat, assign initial/final fragments
!                                    to each MPI proc/node
!     setMPIlimits4densities
!
subroutine master_para_distribute_onebody
	
	use jumpNbody
!	use recruiter
	use para_main_mod
	use para_util_mod
	
	use fragments
	use nodeinfo
	use io
	use opbundles
	implicit none
	logical :: draft
	
	onebodyonly = .true.
	if(.not.allocated(opbundlestart))allocate(opbundlestart(0:nprocs-1))
	if(.not.allocated(opbundleend))allocate(opbundleend(0:nprocs-1))
	
    call master_op_stat_den()    	

	call density_bundles_setup(.false.)		
    if(nprocs > nopbundles)then
 	   draft = .true.
       call density_bundles_setup(.true.)		

    else
       draft = .false.
    end if
!	call f2fjumpstats	
    call getfragmentstats1body(draft)

   	call draftnodal4applyonebody	
	if(draft)call split_density_opbundles  ! ADDED  7.9.5
	
	call distro_opbundles_over_fragments_onebody_ALT		
	
	if(auto_readin  .or. modeldensities)then    ! reading in an established wavefunction
		
    	call  setnodaltribution_alt
	else
		call checknodaltribution
	end if	
	
    call setMPIlimits4densities()   
	
	return
end subroutine master_para_distribute_onebody
!========================================================================
!
!  routine to handle actual distribution of (one-body) opbundles for parallel work
!  added in 7.6.8 (July 2016)
!
!  overall goal is to divide up the jumps as evenly as possible 
!
!  assign start/stop of opbundles on each MPI process /node
!  specifically, using information in opfragstat, fill arrays opbundlestart/opbundleend
!
subroutine distro_opbundles_over_fragments_onebody
	
    use fragments
!	use recruiter	
    use para_util_mod

	use opbundles
	use nodeinfo
	use bmpi_mod
	use io
    implicit none
	
	integer :: ifrag,ffrag
	integer :: myprocs,localproc,testproc,sumprocs
	integer :: ibundle, startopbundle,stopopbundle
	integer(8) :: localjumps,avgjumps,totjumps,jumpssofar
	integer :: ierr
	logical :: fixedprocs  ! = .true. if MPI processes/nodes already assigned

    if(nprocs==1)then
		if(.not.allocated(opbundlestart))print*,' WARNING opbundle start /stop not set '
		opbundlestart(0)=1
		opbundleend(0)=nopbundles
		return
	end if
	
	fixedprocs = .false.
	
!	if(iproc==0)print*,' nodal ',allocated(nodal),' allocated?  also: ',auto_readin
	if(.not.auto_readin .and. .not. modeldensities)then
		fixedprocs = .true.
		if(.not.allocated(nodal))then
			print*,' WARNING!  nodal is not allocated (1)'
			print*,' but it should already be set '
			stop
		end if
		
	end if
	localproc = 0
	
	opbundlestart(:)=0
	opbundleend(:)  =0
	sumprocs = 0
	do ifrag = 1,nfragments
		do ffrag = 1,nfragments
!................ NOW FOR THESE FRAGMENTS DISTRIBUTE THE ONE-BODY JUMPS EVENLY......
!                 FIND START, STOP OF OPBUNDLES
            if(fixedprocs)localproc = sumprocs 
            myprocs = opfragstat(ifrag,ffrag)%nnodes
!			if(fixedprocs .and. iproc==0)print*,' fragments ',ifrag,ffrag,' go from proc ',sumprocs,sumprocs + myprocs-1
			if(myprocs==0)cycle
			startopbundle = nopbundles+1
			stopopbundle = -1
			totjumps = 0
			do ibundle = 1,nopbundles
				if(opbundle(ibundle)%ifragment==ifrag .and. opbundle(ibundle)%ffragment==ffrag)then
					startopbundle = min(startopbundle,ibundle)
					stopopbundle = max(stopopbundle,ibundle)
					totjumps = totjumps + opbundle(ibundle)%njumps
				end if
			end do
!................. CHECK THERE ARE ACTUALLY FRAGMENTS ASSIGNED....			
            if(stopopbundle == -1)then
				if(fixedprocs)sumprocs = sumprocs + opfragstat(ifrag,ffrag)%nnodes
				if(fixedprocs .and. localproc > sumprocs)then
					if(iproc==0)print*,' somehow missed my mark too many procs ',localproc,sumprocs
				end if
			    cycle
			end if
!............. CHECK CONTIGUOUS.................
			
            do ibundle=startopbundle,stopopbundle
				if(opbundle(ibundle)%ifragment/=ifrag .and. opbundle(ibundle)%ffragment/=ffrag)then
					if(iproc == 0)then
 					  print*,' non contiguous opbundles assigned to fragments '
					   print*,' expect ',ifrag,ffrag,' between ',startopbundle,stopopbundle
					   print*,' but bundle ',ibundle,' has fragments ',opbundle(ibundle)%ifragment,opbundle(ibundle)%ffragment
				   end if
#ifdef _MPI
				   call BMPI_ABORT(MPI_COMM_WORLD,101,ierr)
#endif
				   stop
				end if
			end do			
			avgjumps = totjumps/myprocs
						
			localjumps = 0
			opbundlestart(localproc)=startopbundle
			if(localproc > 0 .and. localproc < nprocs-1)opbundleend(localproc-1)=startopbundle-1
			testproc = 1
			jumpssofar = 0
			do ibundle = startopbundle,stopopbundle
				if(opbundle(ibundle)%ifragment/=ifrag .or. opbundle(ibundle)%ffragment/=ffrag)then
					if(iproc==0)then
					   print*,' error in assigning bundles '
					   print*,ifrag,ffrag
					   print*,ibundle
				    endif
				end if

				localjumps = localjumps + opbundle(ibundle)%njumps
				jumpssofar = jumpssofar + opbundle(ibundle)%njumps
				
				if( localjumps > avgjumps)then
					if(testproc < myprocs-1)then   ! reset
					   opbundleend(localproc) = ibundle
					   localproc = localproc+1					
					   testproc = testproc+1					
					   if(ibundle < stopopbundle)then
					  	  opbundlestart(localproc)=ibundle+1
					   end if
					   localjumps = 0
!.... AFTER THIS POINT, READJUST THE DISTRIBUTION....		

                       avgjumps = (totjumps-jumpssofar)/(myprocs-testproc)
				   else
					   avgjumps=localjumps
				   end if
			   end if
			end do
			opbundleend(localproc)=stopopbundle
!............. CHECK IF NOT ALL PROCS USED...........................			

			localproc = localproc+1
			if(fixedprocs)sumprocs = sumprocs + opfragstat(ifrag,ffrag)%nnodes
			if(fixedprocs .and. localproc > sumprocs)then
				if(iproc==0)print*,' somehow missed my mark too many procs ',localproc,sumprocs
			end if
				
			
		end do   ! ffrag		
	end do  ! ifrag

	return
end subroutine distro_opbundles_over_fragments_onebody

!========================================================================
!
!  ADDED in 7.9.5 (April 2020)
!  the problem is, on large machine only a few processes are computing 1-b densities
!  To distribute better, split.
!  Here, we don't try for a optimum distribution, merely a better one.
!
!  Basic algorithm: 
!   -- find avg # of operations
!   -- if an opbundle has > 1.5 x avg, then split
!
subroutine split_density_opbundles
    use fragments
!	use recruiter	
    use para_util_mod

	use opbundles
	use nodeinfo
	use bmpi_mod
	use io
    implicit none
	
	integer :: ifrag,ffrag
	integer :: myprocs,sumprocs
	integer :: ibundle,jbundle, startopbundle,stopopbundle
	integer :: k
	integer(8) :: avgjumps,totjumps,splitjumps,localjumps,jumpbunch
	
	integer(4), allocatable :: split(:)	
	real(4) :: splitratio = 0.6   ! factor to determine if splitting an opbundle
	integer :: maxsplit
	integer :: ierr
	integer :: nopbundles_new
	
	allocate(split(nopbundles))
	
	if(iproc==0)then
		print*,'  '
		print*,' * * * * * '
		print*,' Splitting up opbundles for densities '
		print*,' '
	end if
	split = 0
	sumprocs = 0
	do ifrag = 1,nfragments
		do ffrag = 1,nfragments
!................ NOW FOR THESE FRAGMENTS DISTRIBUTE THE ONE-BODY JUMPS EVENLY......
!                 FIND START, STOP OF OPBUNDLES
 !           if(fixedprocs)localproc = sumprocs 
            myprocs = opfragstat(ifrag,ffrag)%nnodes
!			print*,' nodes ',ifrag,ffrag,myprocs
!			if(fixedprocs .and. iproc==0)print*,' fragments ',ifrag,ffrag,' go from proc ',sumprocs,sumprocs + myprocs-1
			if(myprocs==0)cycle
			startopbundle = nopbundles+1
			stopopbundle = -1
			totjumps = 0
			do ibundle = 1,nopbundles
				if(opbundle(ibundle)%ifragment==ifrag .and. draft_opbundle(ibundle)%ffragment==ffrag)then
					startopbundle = min(startopbundle,ibundle)
					stopopbundle = max(stopopbundle,ibundle)
					
					totjumps = totjumps + draft_opbundle(ibundle)%nops
				end if
			end do

!............. CHECK CONTIGUOUS.................
			
            do ibundle=startopbundle,stopopbundle
				if(draft_opbundle(ibundle)%ifragment/=ifrag .and. draft_opbundle(ibundle)%ffragment/=ffrag)then
					if(iproc == 0)then
 					  print*,' non contiguous opbundles assigned to fragments '
					   print*,' expect ',ifrag,ffrag,' between ',startopbundle,stopopbundle
					   print*,' but bundle ',ibundle,' has fragments ',draft_opbundle(ibundle)%ifragment,draft_opbundle(ibundle)%ffragment
				   end if
#ifdef _MPI
				   call BMPI_ABORT(MPI_COMM_WORLD,101,ierr)
#endif
				   stop
				end if
			end do			
			avgjumps = totjumps/(myprocs+1)
!			print*,ifrag,ffrag,' avg jumps ',avgjumps
			splitjumps = int(avgjumps*splitratio,kind=8)
			
            do ibundle=startopbundle,stopopbundle
				if(draft_opbundle(ibundle)%nops > splitjumps .and.  &
				   draft_opbundle(ibundle)%njumps >= draft_opbundle(ibundle)%nops/splitjumps )then
					split(ibundle)= draft_opbundle(ibundle)%nops / splitjumps
					
				end if
				
			end do 
			
			
		end do ! ffrag
	end do ! ifrag
!............. ALLOCATE ACTUAL OPBUNDLES
	nopbundles_new = nopbundles
	do ibundle= 1,nopbundles
		nopbundles_new = nopbundles_new+split(ibundle)
	end do
	if(iproc==0 .and. nopbundles_new > nopbundles)then
		write(6,*)' Increasing # of density opbundles to ',nopbundles_new
		write(6,*)' * * * * * * * '
		print*,' '
	end if
    if(nopbundles_new > nopbundles)then
		if(allocated(opbundle))deallocate(opbundle)
		allocate(opbundle(nopbundles_new))
	end if
!............ SPLIT	AND COPY OVER
    jbundle = 0
	do ibundle = 1,nopbundles
!		print*,ibundle,' bundle ops ',draft_opbundle(ibundle)%nops,split(ibundle),draft_opbundle(ibundle)%nops/(1+split(ibundle))
		jbundle = jbundle +1
		call copy_draft2bundle(ibundle,jbundle)
		
		if(split(ibundle)>0)then
!............... SET UP SPLIT POINTS.........			
			localjumps = draft_opbundle(ibundle)%njumps
			jumpbunch = localjumps/(split(ibundle)+1)
			do k = 1,split(ibundle)
				jbundle = jbundle + 1
				call copy_draft2bundle(ibundle,jbundle)
				
				if(draft_opbundle(ibundle)%optype=='P1B')then
					opbundle(jbundle-1)%pxend = opbundle(jbundle-1)%pxstart + jumpbunch-1
					opbundle(jbundle-1)%njumps = jumpbunch
					opbundle(jbundle-1)%nops = jumpbunch*(opbundle(jbundle-1)%nxend -opbundle(jbundle-1)%nxstart+1 )
					opbundle(jbundle)%pxstart = opbundle(jbundle-1)%pxend+1
					
					if(k==split(ibundle))then  ! finish data on last new ibundle
						opbundle(jbundle)%njumps =opbundle(jbundle)%pxend -  opbundle(jbundle)%pxstart +1
						if(opbundle(jbundle)%njumps < 0)then
							print*,' Somehow got negative (proton) jumps '
							stop
						end if						
						opbundle(jbundle)%nops = opbundle(jbundle)%njumps*(opbundle(jbundle)%nxend -opbundle(jbundle)%nxstart+1 )	
										
					end if
!					print*,' P1B bundle ',jbundle-1, opbundle(jbundle-1)%nops , ibundle, draft_opbundle(ibundle)%nops 

				else
					if(draft_opbundle(ibundle)%optype/='N1B')then    ! ERROR TRAP
						print*,' wrong optype ',draft_opbundle(ibundle)%optype
						stop
					end if
!					print*,' splitting N1B bundle ',ibundle
					opbundle(jbundle-1)%nxend = opbundle(jbundle-1)%nxstart + jumpbunch-1
					opbundle(jbundle-1)%njumps = jumpbunch
					opbundle(jbundle-1)%nops = jumpbunch*(opbundle(jbundle-1)%pxend -opbundle(jbundle-1)%pxstart+1 )
					opbundle(jbundle)%nxstart = opbundle(jbundle-1)%nxend+1
					
					if(k==split(ibundle))then  ! finish data on last new ibundle
						opbundle(jbundle)%njumps =opbundle(jbundle)%nxend -  opbundle(jbundle)%nxstart +1
						if(opbundle(jbundle)%njumps < 0)then
							print*,' Somehow got negative (neutron) jumps '
							stop
						end if						
						opbundle(jbundle)%nops = opbundle(jbundle)%njumps*(opbundle(jbundle)%pxend -opbundle(jbundle)%pxstart+1 )					
					end if					
					
				end if
							
			
			end do
			
		end if
		
	end do ! ibundle
	nopbundles=nopbundles_new
	
	return	
	
end subroutine split_density_opbundles
!========================================================================
!
! copy over the basic information
!
subroutine copy_draft2bundle(ibundle,jbundle)
	use opbundles
	implicit none
	integer,intent(IN) :: ibundle,jbundle
	opbundle(jbundle)%isector = draft_opbundle(ibundle)%isector
	opbundle(jbundle)%fsector = draft_opbundle(ibundle)%fsector
	opbundle(jbundle)%insector = draft_opbundle(ibundle)%insector
	opbundle(jbundle)%fnsector = draft_opbundle(ibundle)%fnsector
	opbundle(jbundle)%ifragment = draft_opbundle(ibundle)%ifragment
	opbundle(jbundle)%ffragment = draft_opbundle(ibundle)%ffragment
	
	opbundle(jbundle)%optype = draft_opbundle(ibundle)%optype
	opbundle(jbundle)%hchar = draft_opbundle(ibundle)%hchar
	
	opbundle(jbundle)%nsortstart = draft_opbundle(ibundle)%nsortstart
	opbundle(jbundle)%nsortend = draft_opbundle(ibundle)%nsortend
	
	opbundle(jbundle)%pxstart = draft_opbundle(ibundle)%pxstart
	opbundle(jbundle)%pxend = draft_opbundle(ibundle)%pxend
	opbundle(jbundle)%nxstart = draft_opbundle(ibundle)%nxstart
	opbundle(jbundle)%nxend = draft_opbundle(ibundle)%nxend
	
	opbundle(jbundle)%cstride = draft_opbundle(ibundle)%cstride
	opbundle(jbundle)%nops = draft_opbundle(ibundle)%nops
	opbundle(jbundle)%min_nop = draft_opbundle(ibundle)%min_nop
	opbundle(jbundle)%nsetops = draft_opbundle(ibundle)%nsetops
	opbundle(jbundle)%njumps = draft_opbundle(ibundle)%njumps
	
	
	return
	
end subroutine copy_draft2bundle

!========================================================================
!
!  alternate routine to handle actual distribution of (one-body) opbundles for parallel work
!  added in 7.9.5 (April 2019)
!
!  overall goal is to divide up the jumps as evenly as possible 
!
!  assign start/stop of opbundles on each MPI process /node
!  specifically, using information in opfragstat, fill arrays opbundlestart/opbundleend
!
subroutine distro_opbundles_over_fragments_onebody_ALT
	
    use fragments
!	use recruiter	
    use para_util_mod
	use opbundles
	use nodeinfo
	use bmpi_mod
	use io
    implicit none
	
	integer :: ifrag,ffrag
	integer :: myprocs,localproc,testproc,sumprocs
	integer :: ibundle, startopbundle,stopopbundle
	integer(8) :: localjumps,avgjumps,totjumps,jumpssofar
	integer :: ierr
	logical :: fixedprocs  ! = .true. if MPI processes/nodes already assigned


    if(nprocs==1)then
		if(.not.allocated(opbundlestart))print*,' WARNING opbundle start /stop not set '
		opbundlestart(0)=1
		opbundleend(0)=nopbundles
		return
	end if
	
	fixedprocs = .false.
	
!	if(iproc==0)print*,' nodal ',allocated(nodal),' allocated?  also: ',auto_readin
	if(.not.auto_readin .and. .not. modeldensities)then
		fixedprocs = .true.
		if(.not.allocated(nodal))then
			print*,' WARNING!  nodal is not allocated (2) '
			print*,' but it should already be set '
			stop
		end if
		
	end if
		
	localproc = 0
	
	opbundlestart(:)=0
	opbundleend(:)  =0
	sumprocs = 0
	do ifrag = 1,nfragments
		do ffrag = 1,nfragments
!................ NOW FOR THESE FRAGMENTS DISTRIBUTE THE ONE-BODY JUMPS EVENLY......
!                 FIND START, STOP OF OPBUNDLES
            if(fixedprocs)localproc = sumprocs 
            myprocs = opfragstat(ifrag,ffrag)%nnodes
!			if(fixedprocs .and. iproc==0)print*,' fragments ',ifrag,ffrag,' go from proc ',sumprocs,sumprocs + myprocs-1
			if(myprocs==0)cycle
			startopbundle = nopbundles+1
			stopopbundle = -1
			totjumps = 0
			do ibundle = 1,nopbundles
				if(opbundle(ibundle)%ifragment==ifrag .and. opbundle(ibundle)%ffragment==ffrag)then
					startopbundle = min(startopbundle,ibundle)
					stopopbundle = max(stopopbundle,ibundle)
					
					totjumps = totjumps + opbundle(ibundle)%nops
				end if
			end do
!................. CHECK THERE ARE ACTUALLY FRAGMENTS ASSIGNED....			
            if(stopopbundle == -1)then
				if(fixedprocs)sumprocs = sumprocs + opfragstat(ifrag,ffrag)%nnodes
				if(fixedprocs .and. localproc > sumprocs)then
					if(iproc==0)print*,' somehow missed my mark too many procs ',localproc,sumprocs
				end if
			    cycle
			end if
!............. CHECK CONTIGUOUS.................
			
            do ibundle=startopbundle,stopopbundle
				if(opbundle(ibundle)%ifragment/=ifrag .and. opbundle(ibundle)%ffragment/=ffrag)then
					if(iproc == 0)then
 					  print*,' non contiguous opbundles assigned to fragments '
					   print*,' expect ',ifrag,ffrag,' between ',startopbundle,stopopbundle
					   print*,' but bundle ',ibundle,' has fragments ',opbundle(ibundle)%ifragment,opbundle(ibundle)%ffragment
				   end if
#ifdef _MPI
				   call BMPI_ABORT(MPI_COMM_WORLD,101,ierr)
#endif
				   stop
				end if
			end do			
			avgjumps = totjumps/(myprocs+1)
			
!			print*,' avg ops ',avgjumps
			
			localjumps = 0
			opbundlestart(localproc)=startopbundle
			if(localproc > 0 .and. localproc < nprocs-1)opbundleend(localproc-1)=startopbundle-1
			testproc = 1
			jumpssofar = 0
			do ibundle = startopbundle,stopopbundle
				if(opbundle(ibundle)%ifragment/=ifrag .or. opbundle(ibundle)%ffragment/=ffrag)then
					if(iproc==0)then
					   print*,' error in assigning bundles '
					   print*,ifrag,ffrag
					   print*,ibundle
				    endif
				end if

!				localjumps = localjumps + opbundle(ibundle)%nops
				jumpssofar = jumpssofar + opbundle(ibundle)%nops
!				print*,ibundle,opbundle(ibundle)%nops
                if(jumpssofar > avgjumps*testproc)then			
!					print*,ibundle,jumpssofar, ' resetting  THIS NEEDS REVISION '
!					if(iproc==0)print*,' distro den1b ',jumpssofar/testproc
!				if( localjumps > avgjumps)then
					if(testproc < myprocs)then   ! reset
					   opbundleend(localproc) = ibundle
					   localproc = localproc+1					
					   testproc = testproc+1					
					   if(ibundle < stopopbundle)then
					  	  opbundlestart(localproc)=ibundle+1
					   end if
!					   localjumps = 0
!.... AFTER THIS POINT, READJUST THE DISTRIBUTION....		

 !                      avgjumps = (totjumps-jumpssofar)/(myprocs-testproc)
				   else
!					   avgjumps=localjumps
				   end if
			   end if
			end do
			opbundleend(localproc)=stopopbundle
!............. CHECK IF NOT ALL PROCS USED...........................			

			localproc = localproc+1
			if(fixedprocs)sumprocs = sumprocs + opfragstat(ifrag,ffrag)%nnodes
			if(fixedprocs .and. localproc > sumprocs)then
				if(iproc==0)print*,' somehow missed my mark too many procs ',localproc,sumprocs
			end if
				
			
		end do   ! ffrag		
	end do  ! ifrag

	return
end subroutine distro_opbundles_over_fragments_onebody_ALT

!========================================================================
!
! routines to create "opbundles" for 1-body density matrices for parallelization
!
!========================================================================
!
! NOTES for 7.3.8  Dec 2014
!  (1) Look at assigned initial, final fragments
!  (2) within each ifrag,ffrag, identify relevant 1-body sector jumps; may need to look at F/R
!  (3) create opbundles for each sector jump; may need to look at F/R
!  (4) assume density matrix goes so fast you don't need to divide up work
!
!  in set1bsectorjumps, hermiticity is not used;
!  it is used, however, in old routine density1b_sector_pjumps and density1b_sector_njumps
!  when skiph=x1bjump(it)%sjmp(spjmp)%diag  is false which happens when in different sectors;
!  I worry this is double counting but it appears not to be
!  it appears opchar = 'h' is not used but need to confirm
!
!  Assumes fragments have already been distributed and assigned
!
!  CALLED BY:
!    applicator1b
!    density1b_output 
!    density1b_from_oldwfn
!  
!
   subroutine density_bundles_setup(draft)

   use opbundles
   use fragments
   use verbosity
   use nodeinfo
   use operation_stats
   use io
   implicit none
   integer ifrag,ffrag
   integer nob,nob_draft
   integer(8) :: temptot
   logical draft
   integer ibundle
   integer :: aerr

   nob = 0
   
!   if(nproc > nopbundles)then
!	   draft = .true.
!   else
!      draft = .false.
!   end if
   do ifrag = 1,nfragments
       do ffrag = 1,nfragments
              call count_create_draft_opbundles_den(draft,ifrag,ffrag,.false.,'P1B','f',nob)
              call count_create_draft_opbundles_den(draft,ifrag,ffrag,.false.,'P1B','h',nob)
              call count_create_draft_opbundles_den(draft,ifrag,ffrag,.false.,'N1B','f',nob)
              call count_create_draft_opbundles_den(draft,ifrag,ffrag,.false.,'N1B','h',nob)
	       if(ifrag <= ffrag)then
              call count_create_draft_opbundles_den(draft,ffrag,ifrag,.false.,'P1B','b',nob)
              call count_create_draft_opbundles_den(draft,ffrag,ifrag,.false.,'N1B','b',nob)
           end if
       end do

   end do
   if(iproc==0)then
       print*,' '
       print*,nob,' opbundles needed for density matrices '
       print*,' '
	   if(nprocs > nob*.25 .and. nprocs > 1)then
		   print*,' * * * WARNING * * * '
		   print*,' May have too many MPI ranks to distribute opbundles; could crash '
		   print*,' You can compute wfn only and then use option dx/dxp to generate densities '
		   print*,' * * * WARNING * * * '
		   print*,' '
		   write(logfile,*)' WARNING '
		   write(logfile,*)' May have too many MPI ranks to distribute opbundles '
		   write(logfile,*)' You can compute wfn only and then use option dx/dxp to generate densities '
		   write(logfile,*)' WARNING '
       end if
	   
   end if
   if(draft)then
       if(nob > nopbundles .and. allocated(draft_opbundle))then
           deallocate(draft_opbundle)
           allocate(draft_opbundle(nob), stat=aerr)
           if(aerr /= 0) call memerror("density_bundles_setup 1a")
       end if
       if(.not.allocated(draft_opbundle))then
          allocate(draft_opbundle(nob), stat=aerr)
          if(aerr /= 0) call memerror("density_bundles_setup 2a")
       end if   
   else
      if(nob > nopbundles .and. allocated(opbundle))then
          deallocate(opbundle)
          allocate(opbundle(nob), stat=aerr)
          if(aerr /= 0) call memerror("density_bundles_setup 1b")
      end if
      if(.not.allocated(opbundle))then
         allocate(opbundle(nob), stat=aerr)
         if(aerr /= 0) call memerror("density_bundles_setup b2")
      end if
  end if
   if(.not.allocated(frag_draftopbundlestart)) then
      allocate(frag_draftopbundlestart(nfragments,nfragments), stat=aerr)
      if(aerr /= 0) call memerror("density_bundles_setup 3")
   end if
   nob = 0
    
   do ifrag = 1,nfragments
       do ffrag = 1,nfragments
           frag_draftopbundlestart(ifrag,ffrag) = nob+1

             call count_create_draft_opbundles_den(draft,ifrag,ffrag,.true.,'P1B','f',nob)
             call count_create_draft_opbundles_den(draft,ifrag,ffrag,.true.,'P1B','h',nob)

             call count_create_draft_opbundles_den(draft,ifrag,ffrag,.true.,'N1B','f',nob)
             call count_create_draft_opbundles_den(draft,ifrag,ffrag,.true.,'N1B','h',nob)
!             if(ifrag <= ffrag)then
	          if(ifrag < ffrag)then

              call count_create_draft_opbundles_den(draft,ffrag,ifrag,.true.,'P1B','b',nob)
              call count_create_draft_opbundles_den(draft,ffrag,ifrag,.true.,'N1B','b',nob)
           end if
       end do

   end do

   if(.not.allocated(opbundlestart)) then
      allocate( opbundlestart(0:nprocs-1), stat=aerr)
      if(aerr /= 0) call memerror("density_bundles_setup 10")
   end if
   if(.not.allocated(opbundleend)) then
      allocate( opbundleend(0:nprocs-1), stat=aerr)
      if(aerr /= 0) call memerror("density_bundles_setup 20")
   end if
   nopbundles = nob
   
   if(nprocs==1)then  ! default
	   opbundlestart(0)=1
	   opbundleend(0)  = nopbundles
   end if
   call checkbundlefragments
   return

   end subroutine density_bundles_setup

!===============================================================

! creates a draft set of bundles
!
! 
!
subroutine count_create_draft_opbundles_den(draft,ifrag,ffrag,create,optype,hchar,nob)


  use fragments
!  use subsectors
  use sectors
  use opbundles
  use operation_stats
  use jumpNbody
  use jump3body
  use basis
  use nodeinfo

  implicit none

  logical draft
  logical create
  integer nob  !   temporary for # of opbundles
  integer ifrag,ffrag
  character(3) :: optype
  character(1) :: hchar
  integer :: is,fs  ! indices for subsectors
  integer(8) :: xstartx,xendx,ystarty,yendy
  integer :: xsj, ysj, isj
  integer :: cs, ncs, ics
  integer(8) :: cstride
  integer :: ibundle
  logical okay

  type (bund), pointer :: bundle(:)
  type (bund), target :: bundle_dummy(1)
  bundle => bundle_dummy
  if(draft .and. create)bundle => draft_opbundle
  if(.not.draft .and. create)bundle => opbundle

  do is = fragmentlist(ifrag)%ssectorstart, fragmentlist(ifrag)%ssectorend
     do fs = fragmentlist(ffrag)%ssectorstart, fragmentlist(ffrag)%ssectorend

         select case(optype)

!........... ASSIGN PP .......................

         case('P1B')
!........... ASSIGN P .......................

!..........FIND START, STOP FOR P JUMPS..........
!          ASSUME SUBSECTORS SAME AS SECTORS
               xsj = -1
               do isj = 1,x1bjump(1)%nsectjumps
                   if( is == x1bjump(1)%isector(isj) .and.  & 
                        fs == x1bjump(1)%fsector(isj)  ) then
                      xsj = isj
                      exit
                   end if 

               end do

               if(xsj == -1)cycle   ! no P1B operations here
               if(hchar=='h' .and. (is/=fs .or. .not. x1bjump(1)%sjmp(xsj)%diag) )cycle
               if( ( hchar=='f' .or. hchar == 'b') .and.  & 
                   (is==fs .and. x1bjump(1)%sjmp(xsj)%diag) )      cycle
			   if(x1bjump(1)%sjmp(xsj)%njumps<1)cycle
               xstartx = x1bjump(1)%sjmp(xsj)%nstart+1
               xendx   = xstartx + x1bjump(1)%sjmp(xsj)%njumps-1

!..........FIND START, STOP FOR NEUTRON SDs.........
               cs = x1bjump(1)%csjmp(xsj)%cjump(1)  ! this is the first conjugate
                                                     ! neutron sector
               ystarty = xsd(2)%sector(cs)%xsdstart
               ncs = x1bjump(1)%csjmp(xsj)%ncjmps
               cs  = x1bjump(1)%csjmp(xsj)%cjump(ncs)  
               yendy = xsd(2)%sector(cs)%xsdend

!..........CONVERT TO STATE INDICES
               nob = nob+1

               cstride = 1

               if(create)then
                  bundle(nob)%optype     = 'P1B'

                  bundle(nob)%pxstart    = xstartx
                  bundle(nob)%pxend      = xendx
                  bundle(nob)%nxstart    = ystarty
                  bundle(nob)%nxend      = yendy
                  bundle(nob)%cstride    = cstride
                  bundle(nob  )%hchar      = hchar
				  bundle(nob)%njumps     = xendx-xstartx+1
				  bundle(nob)%nops     = (xendx-xstartx+1)*(yendy-ystarty+1)
				  
!............. 'forward'.............

                  if(hchar == 'f')then
                     bundle(nob)%ifragment  = ifrag
                     bundle(nob)%ffragment  = ffrag
                     bundle(nob)%isector = is
                     bundle(nob)%fsector = fs
                  end if

!............. 'backward'.................

                  if(hchar == 'b' .or. hchar == 'h')then
                     bundle(nob)%ifragment  = ffrag
                     bundle(nob)%ffragment  = ifrag
                     bundle(nob)%isector = fs
                     bundle(nob)%fsector = is					 
                  end if
               end if
!........... ASSIGN N .......................

         case('N1B')

!..........FIND START, STOP FOR N JUMPS..........
!          ASSUME SUBSECTORS SAME AS SECTORS
            if(is /= fs)cycle   ! cannot change proton sector

            do xsj = 1,x1bjump(2)%nsectjumps
               
                ncs =  x1bjump(2)%csjmp(xsj)%ncjmps 
                okay = .false.
                do cs = 1,ncs
                  if( is == x1bjump(2)%csjmp(xsj)%cjump(cs)) then
                      okay =.true.
                      ysj = cs
                      exit
                   end if 
                end do
                if(.not.okay)cycle
               if(hchar=='h' .and.  .not. x1bjump(2)%sjmp(xsj)%diag )cycle
               if( ( hchar=='f' .or. hchar == 'b') .and.  & 
                    x1bjump(2)%sjmp(xsj)%diag )      cycle
				if(x1bjump(2)%sjmp(xsj)%njumps<1 )cycle
                xstartx = x1bjump(2)%sjmp(xsj)%nstart+1
                xendx   = xstartx+x1bjump(2)%sjmp(xsj)%njumps -1
!..........FIND START, STOP FOR PROTON SDs.........
!          assume subsectors = sectors!
                ystarty = xsd(1)%sector(is)%xsdstart
                yendy =   xsd(1)%sector(is)%xsdend

!..........CONVERT TO STATE INDICES
                cstride = 1
                nob = nob + 1

                if(create) then  ! set up for N1b operations

                  bundle(nob)%optype     = 'N1B'
                  bundle(nob)%pxstart    = ystarty
                  bundle(nob)%pxend      = yendy
                  bundle(nob)%nxstart    = xstartx
                  bundle(nob)%nxend      = xendx
                  bundle(nob)%nsortstart    = xstartx
                  bundle(nob)%nsortend      = xendx				  
                  bundle(nob)%cstride    = cstride
                  bundle(nob  )%hchar      = hchar
				  bundle(nob)%njumps     = xendx-xstartx+1
				  bundle(nob)%nops     =  (xendx-xstartx+1)*(yendy-ystarty+1)
				  
!............. 'forward'.............

                  if(hchar=='f')then

                     bundle(nob)%ifragment  = ifrag
                     bundle(nob)%ffragment  = ffrag
                     bundle(nob)%isector = is
                     bundle(nob)%fsector = fs
                   end if
!............. 'backward'.................

                  if(hchar=='b' .or. hchar=='h')then

                     bundle(nob)%ifragment  = ffrag
                     bundle(nob)%ffragment  = ifrag
                     bundle(nob)%isector = fs
                     bundle(nob)%fsector = is
                   end if

                end if
            end do
         end select

      end do ! fs
  end do  ! is
  
  return
end subroutine count_create_draft_opbundles_den
!===============================================================
!
! routine to check ordering of fragments in bundles, checks they are contiguous
!
subroutine checkbundlefragments
    use fragments
	use opbundles
	use nodeinfo
	
	implicit none
    character(3) :: optype
	integer :: ibundle
    integer :: ifrag,ffrag
	integer :: bstart,bend
	
	
	do ifrag=1,nfragments
		do ffrag = 1,nfragments
			bstart = nopbundles+1
			bend = 0
			
			do ibundle = 1,nopbundles
				if(ifrag == opbundle(ibundle)%ifragment .and. ffrag ==  opbundle(ibundle)%ffragment)then
					bstart = min(bstart, ibundle)
					bend = max(bend,ibundle)
				end if
				
			end do
!			if(bend==0 .and. iproc==0)print*,' fragments ',ifrag,ffrag,' get no bundles '
			if(bend==0)cycle
!			if(iproc==0 )print*,' fragments ',ifrag,ffrag,' get bundles ',bstart,bend
			do ibundle = bstart,bend
				if(ifrag /= opbundle(ibundle)%ifragment .and. ffrag /=  opbundle(ibundle)%ffragment)then
					if(iproc==0)then
						print*,' Non contiguous list of bundles  for fragments ',ifrag,ffrag
						print*,' start /stop = ',bstart,bend,' but ',ibundle,' has ',opbundle(ibundle)%ifragment,opbundle(ibundle)%ffragment
					end if
					
				end if
			end do  ! ibundle
			
			
		end do ! ffrag
		
	end do ! ifrag	
!    if(iproc ==0)print*,' FINISHED check contiguiuty of opbundles and fragments '
	return
end subroutine checkbundlefragments


!===============================================================



