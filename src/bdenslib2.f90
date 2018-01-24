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
  call getfragmentstats


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
	use recruiter
	use para_main_mod
	use fragments
	use nodeinfo
	use io
	use opbundles
	use onebodypot, only: meanie
	implicit none
	
	onebodyonly = .true.
	if(.not.allocated(opbundlestart))allocate(opbundlestart(0:nprocs-1))
	if(.not.allocated(opbundleend))allocate(opbundleend(0:nprocs-1))
    call master_op_stat_den()    
	
    call density_bundles_setup()	
	call f2fjumpstats
   	call setnodal4applyonebody_DRAFT
	call distro_opbundles_over_fragments_onebody	
	
	if(auto_readin .or. meanie)then    ! reading in an established wavefunction
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
	use recruiter	
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
	if(.not.auto_readin)then
		fixedprocs = .true.
		if(.not.allocated(nodal))then
			print*,' WARNING!  nodal is not allocated '
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
				   call BMPI_ABORT(icomm,101,ierr)
				   stop
				end if
			end do			
			avgjumps = totjumps/myprocs
			
			localjumps = 0
			opbundlestart(localproc)=startopbundle
!			if(iproc==0)print*,'(ddd)',localproc,startopbundle
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
				
				if( localjumps > avgjumps .and. testproc < myprocs)then   ! reset
					opbundleend(localproc) = ibundle
					localproc = localproc+1
					
					testproc = testproc+1					
					if(ibundle < stopopbundle)then
						opbundlestart(localproc)=ibundle+1
					end if
					localjumps = 0
!.... AFTER THIS POINT, READJUST THE DISTRIBUTION....					
                    avgjumps = (totjumps-jumpssofar)/(myprocs-testproc)
				end if
			end do
			opbundleend(localproc)=stopopbundle
!............. CHECK IF NOT ALL PROCS USED...........................			
			
!            if(testproc < myprocs)then
!				if(iproc==0)print*,' For fragments ',ifrag,ffrag, myprocs -testproc,' MPI procs not used'
!			end if
			localproc = localproc+1
			if(fixedprocs)sumprocs = sumprocs + opfragstat(ifrag,ffrag)%nnodes
			if(fixedprocs .and. localproc > sumprocs)then
				if(iproc==0)print*,' somehow missed my mark too many procs ',localproc,sumprocs
			end if
				
			
		end do   ! ffrag		
	end do  ! ifrag
	
!	if(iproc==0)print*,' finished distributing, I think '

!    do testproc = 0,nproc-1
!	    if(opbundlestart(testproc)==0 .or. opbundleend(testproc)==0)then
!		   if(iproc==0)then
!		  	   print*,' ONLY ',testproc+1,' MPI PROCESSES USED '
!		   end if
!	    end if
!    end do
	return
end subroutine distro_opbundles_over_fragments_onebody

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
   subroutine density_bundles_setup

   use opbundles
   use fragments
   use verbosity
   use nodeinfo
   use operation_stats
   implicit none
   integer ifrag,ffrag
   integer nob,nob_draft
   integer(8) :: temptot
   logical draft
   integer ibundle
   integer :: aerr

   nob = 0
   draft = .false.
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
   end if
   if(nob > nopbundles .and. allocated(opbundle))then
         deallocate(opbundle)
         allocate(opbundle(nob), stat=aerr)
         if(aerr /= 0) call memerror("density_bundles_setup 1")
   end if
   if(.not.allocated(opbundle))then
      allocate(opbundle(nob), stat=aerr)
      if(aerr /= 0) call memerror("density_bundles_setup 2")
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
             if(ifrag <= ffrag)then
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
				  bundle(nob)%njumps     = yendy-ystarty+1

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



