! !===========================================================================
!  BIGSTICK
!
!  Inspired by code REDSTICK by W. E. Ormand
!  
!    History (for more, see BIGSTICK_HISTORY.txt)
!
!  VERSION 2.0: Setting up for basis -- CWJ @ SDSU -- June/July 2008
!  VERSION 3.0: BASIS DONE, hops include-- CWJ  @ SDSU -- OCT 2008
!  VERSION 4.0: jumps done                 CWJ @ SDSU  -- DEC 2008
!  VERSION 5.0   : 1-body jumps fixed in no-truncation case
!...........................................................................  
!  6.0   : Includes MPI codes from W. E. Ormand and P. G. Krastev    -- June 9, 2009
!  6.4.17: successful 3-body version                                  -- Aug 2011 (CWJ)
!  6.5.0 : Set runtime options for splitting Lanczos vector           -- February 2011 (PGK)
!  6.6.0:  merged code                                                -- August 2011 (CWJ)
!  6.6.6:  new "opbundles" for organizing application of Hamiltonian -- Sept 2011 (CWJ)
!  6.6.7: fixed strength functions; added WEO's OpenMP; Oct 2011 (CWJ)
!  6.6.9: improved usage of geneologies, chains Oct 2011 CWJ
!  6.7.6:  opbundles fixed for use with 3-body
!  6.7.7: OpenMP works with 3-body
!  6.8.0: Rebooting MPI
!  6.8.1: Single-precision Lanczos vectors; added "preconditioning"
!  6.8.2: Fixed small bug in strength functions for "full" diagonalization
!  6.8.3: Added "time-reversal" for pn matrix elements, saving some memory
!                              (but taking more time)
!  6.8.5: MPI counting of jump
!  6.8.7: MPI distribution of operations and mock serial runs
!  6.8.10: Improved instrumentation of timing;
!          Added bapplyhlib3.f90: runs an operation all at once
!  without leaving a subroutine;
!          Corrected error in estimating # of NNN operations
!  6.9.0:  Added weighting by timing
!  
!  7.0.0:  MPI Mark 1 runs correctly, including calculation of J^2, T^2 (Sept 2012, CWJ)
!  7.0.2:  Reads in J-coupled 3-body matrix elements
!  7.0.3:  restricts creation of jumps to MPI compute nodes where needed (Nov 2012, CWJ)
!  7.0.4:  storage of Lanczos across MPI compute nodes (Jan 2013, WEO)
!  7.0.6:  internal storage of Lanczos vectors option even for single-processor (Mar 2013, CWJ)
!  7.0.7:  RELEASE VERSION -- tested, although still problems in MPI
!  7.0.8:  option to write lanczos vectors to disk even when storing internally; 
!          coupled to restart option (Mar 2013, CWJ)
!	   added OpenMP to expectation values ('x'), density matrices ('d') (Apr 2013/CWJ)
!  7.0.9:  added double precision for eigenvalues (egv_prec in module precisions)
!          capability to read in uncoupled/deformed TBMEs
!  7.0.10: bugs in restart fixed
!  7.1.0 : problem with precision in MPI orthogonalization in bparallellib4.f90 fixed
!  7.1.1 : added switching order of loops in applyh if no OpenMP threading
!  7.1.2/3 : made vector "reading"/"writing"/"orthogonalization" more uniform across
!            parallelizations
!          fixed bug in MPI when sorting jumps--only sort once
!  7.1.4 : thick-restart now works when storing lanczos vectors in RAM (including MPI)
!          OpenMP for uncoupling PP/NN/PN matrix elements  (July 2013/CWJ)
!  7.1.5 : improved handling of eigenvectors for compute J^2, T^2, TRDENS output
!          now store eigenvectors internally (if lanczos vectors stored)
!  7.1.6 : switch order of loops in applyh extended to 3-body
!  7.1.7 : turned off OpenMP for computing traces
!  7.1.8 : "improved" menu for Lanczos
!  7.1.9 : improved efficiency in MPI calculation of jumps; fixed occasional bug in computing jump distribution
!  7.2.0-2: moving towards cap on storing jumps on a node
!  7.2.3 : draft routines to cap storage of jumps on nodes (need further testing)
!  7.2.4 : turned calculation of T back on even when breaking isospin
!  7.2.5 : fixed bug for density matrices in exact diagonalization; 
!          started towards p-h conjugation
!          put all interfaces into 'binterfaces.inc'
!          started towards direct storage of XX/XXX matrix elements
!          computing % of wfns in W-truncation subspaces
!  7.2.6 : p-h conjugation works now especially for just one species changed
!  7.2.8 : allowed for different W-truncation for proton, neutron spaces (Feb 2014)
!        : corrections to deformed option
!  7.2.10: calculate break-up of Lanczos vectors
!  7.2.11: RELEASE VERSION April 23 2014
!  7.2.12: Simpler (non-optimized) breakup of Lanczos vectors into fragments; distributing work bewtween fragments
!  7.2.13: Fixed an MPI bug in write_wfn_header; improved timing of jump creation and decoupling of matrix elements
!          fixed some initialization errors
!  7.2.14: More work on breaking up lanczos vectors into fragments; detailed notes on distribution 
!          (header of bparallel_lib1.f90)
!          Started work on alternative storage of PN, PPN, PNN uncoupled matrix elements NEEDS TO BE REVISED
!          Adding timing of "descendents" mostly to see where jump creation time is going.
!  7.2.15: Continued work on breaking up lanczos vectors; working draft of distribution routines
!  
!  7.3.0 : MPI distribution over fragments
!  7.3.1 : reorthgonalization routines over fragments (KSM)
!  7.3.3 : Fixing some input routines for isospin breaking (WEO)
!          Added environment variables for input files directories BIG_SPS_DIR and BIG_INT_DIR (WEO)
!  7.3.4 : Fixed some run-time information on distributing jumps
!  7.3.5 : improved reorthgonalization (KSM)
!  7.3.6 : started toward improved distribution of PN etc matrix elements (CWJ); 
!          can now write fragmented vectors to disk (KSM)
!          reinstalled particle-occupation option
!  7.3.7:  improved inline calculation of single-particle occupations
!  7.3.8:  wrote and debugged opbundles for one-body densities
!  7.3.9:  merged with updates from KSM, improved handling of MPI
!  7.4.0:  more experiments with improved PN storage, introduced "introns" to study discontinuous storage of jumps
!  7.4.1:  merged with updates from KSM, improved handling of particle occupations, MPI
!       :  routines towards improved fragments
!  7.4.2: Major upgrade of MPI wrappers (KSM) with more robust handling of data types/passing data
!         Minor bugs in MPI calls fixed (KSM)
!         More work on better distribution of jumps and deleting zeroes in jump arrays (via "introns" in bintron.f90)
!         Started work on end-of-run report (module reporter/subroutine report)
!         Fixed in part "overlap" option
!  7.4.3  Overlap, apply scalar operator, and apply 1-body operator work in serial mode (at least)
!  7.4.4  Fixing errors in running with 3-body forces; revisions to postprocessing options to provide uniformity
!  7.4.5  Eliminated "subsectors"
!  7.4.6  Returned to improved fragmentation based upon haiku blocks (in progress);
!         fixed bugs in apply one-body nonscalar operators, overlap
!  7.4.9  Improved fragmentation debugged; fixed a bug where we fail to assign processes between fragments
!  7.5.0  Small fixes; moved towards OpenMP parallelism in counting jumps
!  7.5.1  Progress on more compact jump storage (a.k.a. "intron deletion")
!  7.5.2  Revising matrix element storage and decoupling to make it easier to port to other codes
!  7.5.3  Improved storage of jumps (intron deletion) implemented in 1- and 2-body jumps
!  7.5.4  Improved storage of jumps (intron deletion) implemented in 3-body jumps
!  7.5.5  Removed "HF" and preconditioner; see older versions (7.5.4) if you want to reintroduce
!         Fixed problem (at least in serial) reading .wfn files and computing density matrices
!  7.5.6  Fixed density matrices to run in particle-hole conjugation 
!  7.5.7  Fixed minor issues in processing answers with "exact" (Householder) diagonalization
!  7.5.8 (skipped over; development of grouples for advanced storage)
!  7.5.9  Changed # of uncouple matrix elements (nmatXX, nmatpn, etc) to integer(8)
!  7.6.0  Added log file, timing of intron routines, minor fixes to p-h conjugation with truncation
!  7.6.3  Rewrote lanczos routine so pivot is introduced outside it; added options to "improve" pivot
!         by applying selected parts of the Hamiltonian
!  7.6.5  Improved OpenMP caching KSM
!  7.6.6  Fixed bugs in OpenMP caching on old compilers, compiles now on Vulcan
!         Reading in prior wavefunctions works for fragments for options 'x', 's'
!  7.6.7  Reorganized routines for 1-body densities and application; made latter
!         conform to opbundle structures 
!  7.6.8  Apply one-body should now work in MPI, with fragments 
!  7.6.9  Continued reorganization / inclusion of subroutines into modules, esp. basis routines
!  7.7.0  Inline 1-body densities work with MPI and fragments
!         more reorganization of routines into modules, esp. parallel distribution routines
!  7.7.1  1-body density matrices from old wavefunctions works in MPI now, incl. fragments
!              for more updates see file WHATSNEW.txt
!==================================================================
! IMPORTANT NOTES:  For more details see README_BIGSTICK.txt  (incomplete)
!                 also see SUBROUTINE_SUMMARY.txt (incomplete)
!
!  here protons are indexed by 1 and neutrons by 2
!  furthermore, protons with m >=0 are +1, m < 0 are -1
!               neutrons with m >0 are +2, m < 0 are -2
! (We allow for fermions without intrinsic spin, that is, j = 0, 1, 2)
!
! Truncation is via an additive quantum number W. This is the equivalent of
! "Nhw" in REDSTICK and "t" in ANTOINE. 
!
! The basis is a factorized product of proton and neutron Slater determinants.
! Each SD in turn is a factorized product of a left and right haiku.
! Left haikus are constructed from particles with m < 0, right with m >= 0
! Hence right proton haikus are indexed by +1, left proton haikus by -1 etc
! While SDs have fixed particle number, haikus do not.
!
! "Hops" are fundamental operations on haikus, either the destruction or creation
! of a particle. These change particle number
!
!  "Jumps" are operations on SDs, constructed from hops. These conserve particle number.
!  The action of the Hamiltonian is the factorized product of jumps. 
!
!===========================================================================
!
! routines in this file:
!   bigstick_main  (main program)
!   getoutputfile
!   memory_requirements 
!   setupformodeling
!
!===========================================================================
! 
! main program bigstick 
!
! SUBROUTINES CALLED:
!	clocker 	: starts internal time
!	random_seed
!	menu 		: runtime options
!       getoutputfile   : opens an output file
!       get_orbit_info  : reads in file on single-particle orbits
!       extract_sps     : inflates s.p. orbits to s.p. states
!       sort_spsqn      : sorts single-particles states by q.n.
!	define_system   : choose valence Z,N, parity, W
!       make_hspsqn     : creates subsets of s.p. states needed for constructing 'haikus'
! 	master_w_limits : truncation on W
!	basismaster     : construct basis
!	write_out_basis : (optional) write out basis to file
! 	sector_stats    : info on basis sectors (sorted by Jz, parity, W)
!	checkbasisdim	: if wfn file read in, check basis dimension consistent
! 	write_wfn_header: create header for wfn file
!	readin1bodyop	: (optional) read in 1-body operator to apply to pre-existing wfn
!       decouple1bodyop : (optional) decouple 1-body op
!       hoopmaster      : count up jumps
!	memory_requirements 
!     	timeperopmaster
! 	master_op_stat	: statistics on operations in preparation for parallelization
! 	master_para_distribute : compute parallel distribution of operations
!	threebodymaster
!	master_readin_tbmes
!	makespeme	: create arrays for applying single-particle energies
!	uncoupleXXtbme  : uncouple PP/NN matrix elements
! 	count_uncouplepn: uncouple PN matrix elements
!       masterconvert2to3bme 
!       master_read_3bmes     ! note: must come AFTER masterconvert2to3body even if no 2-body read in
!       jumpmaster      : create jumps
!     	sort_bundled1bops	: sorting for OpenMP parallelism
!	set_bundledXY_threadstart : more setup for OpenMP
!	bundle_clock	: timing
!     	lanczos_master	: main lanczos routine
!	density1b_from_oldwfn
!	expectator_p	: expecation value of Hamiltonian-like operator
!	applicator1b	: get new wfn by applying 1-body operator (transition-like)
!	applicator_h	: get new wfn by applying Hamiltonian-like operator
!	tracemaster	: compute 1st, second moments
!	clock_out	: print out timing results
!
program bigstick_main
  use program_info
  use hermit
  use io
  use verbosity
  use lanczos_info
  use menu_choices
  use basis
  use sporbit
  use coupledmatrixelements ! interaction
  use system_parameters
  use nodeinfo
  use flags3body
  use timing
  use haiku_info
  use flagger
  use wfn_mod
  use pocc_mod
  use switch_mod
  use para_main_mod
  use b3bme_mod
  use btbme_mod
  use jump_mod
  use bsector_mod
  use para_bundles_mod
  use diagh
  use apply_obs_mod
  use sectors
  use timeinfo
  use onebodypot
  use shampoo
  use sp_potter
  implicit none

  integer(4) :: ierr
  integer(8) :: nme
  integer(4) :: i
  character(1) :: ychar
!----------- FLAGS FOR RESTRICTION ON CHANGING q#s in jumps-----------------
  logical    :: change1bodyqs, change2bodyqs, change3bodyqs
  logical    :: hermflag, whermflagp, whermflagn
  logical threebody0
!--------   Thread info-----------------------------------------------------
  integer num_threads, mythread
  integer omp_get_num_threads
  integer omp_get_thread_num
  character(len=MPI_MAX_PROCESSOR_NAME) :: mpiprocname
  integer :: mpiprocnamelen
  
  integer :: is, isc, jsc
  integer(8) :: ndim
  

  ! KSM DBG: force standard in
  ! needed for pgdbg - no way to redirect stdin with MPI
  ! open(5, FILE='te123.minput', status='OLD', action='READ')
  ! open(5, FILE='xe136.minput', status='OLD', action='READ')

! UNCOMMENT IF YOU WANT BACKTRACE
  call installsignalhandlers  ! better reporting on crashes / interrupt
  call BMPI_INIT(ierr)
  icomm = MPI_COMM_WORLD
  call BMPI_COMM_RANK(icomm,iproc,ierr)
  call BMPI_COMM_SIZE(icomm,nproc,ierr)
  noisy0 = iproc == 0 .and. wantnoisy0

  call MPI_GET_PROCESSOR_NAME(mpiprocname, mpiprocnamelen, ierr)
  ! On HPC machines, the "name" will describe the location of
  ! the physical processor.   For example, on edison
  ! the routine returns a string like nid01418
  ! The node id identifies the compute node.  
  ! You can run
  !  %  xtdb2proc -f edisontopo.out
  ! to produce a table that can be used to map node ids to locations in
  ! the compute hierarchy.
  ! We can simply mpi_gather the node ids back to rank 0 to make a table.
  !
  ! print *, "iproc=", iproc, ", name=", TRIM(mpiprocname)

  call clocker('all','start')
  call random_seed

!--------- ADDED IN V 7.3.3 BY WEO --- ENVIRONMENTAL VARIABLES FOR INPUT FILES----

  if(iproc==0)then
    sps_path(1:200)=' '
    call getenv('BIG_SPS_DIR',sps_path)
    length_sps_path=index(sps_path,' ')
    sps_path(length_sps_path:length_sps_path)='/'

    int_path(1:200)=' '
    call getenv('BIG_INT_DIR',int_path)
    length_int_path=index(int_path,' ')
    int_path(length_int_path:length_int_path)='/'
  end if

!------------------INTRODUCTION---------------------------------------------
  if ( iproc == 0 ) then
     write(6,*)' '     
     write(6,*)' BIGSTICK: a CI shell-model code '
     write(6,*)' Version ',version,lastmodified
     write(6,*)' '
     write(6,*)' by C. W. Johnson, W. E. Ormand, '
     write(6,*)' K. S. McElvain, and H.Z. Shan '
     write(6,*)' For reference please cite: '
     write(6,*)' C. W. Johnson, W. E. Ormand, and P. G. Krastev '
     write(6,*)' Comp. Phys. Comm. 184, 2761-2774 (2013) '
     write(6,*)' [also found in arXiv:1303.0905] '
     write(6,*)' and report UCRL LLNL-SM-739926 '
     write(6,*)' '
	 write(6,*)' This code distributed under the MIT Open Source License '
  end if

  ! read environment and setup paths like scratch_dir
  call main_getenv()  ! needs MPI_INIT

  if ( iproc == 0 ) write(6,*) 'Number of MPI processors =',nproc
!$omp parallel shared(num_threads,iproc)  
   num_threads = omp_get_num_threads()
   mythread = omp_get_thread_num()
   if ( iproc == 0 .and. mythread == 0 ) write(6,*) 'NUM_THREADS = ', num_threads

   if(num_threads > 1 .and. .not.dosortjumps)then
       if(iproc==0.and.mythread==0)then
          write(6,*)' For OpenMP: must have flag dosortjumps = T '  
          write(6,*)' in module flagger found in bmodules_flags.f90  '
       end if
       stop
   end if
!$omp end parallel

  call menu      ! CALL MAIN MENU FOR BIGSTICK 

!---------------- SET SOME FLAGS -------------------------------------------
  hermitian =  .true.
  hermflag = .true.
  whermflagp = .true.   ! for protons
  whermflagn = .false.  ! for neutrons -- ALWAYS TURN OFF
  simulateMPI = simulateMPIdefault

  if ( op_readin ) then
     change1bodyqs = .false.
  else
     change1bodyqs = .true.
  endif
  change2bodyqs = .false.     ! NB this might change with 3-body forces

  change3bodyqs = .false.

!-----------------  Set a value defining type of Lanczos vector for MPI routines
  if ( lanc_prec == 4 ) then
     MPI_lanc = MPI_REAL
  else
     MPI_LANC = MPI_REAL8
  end if

  modifytimewts = .true.

  if(menu_char/='r' .and. menu_char /= 'c' .and. menu_char/='v' .and. menu_char/='b')call getoutputfile
  call openlogfiles
  call writelogfile('ini')
!..........................................................................
!------------------SINGLE-PARTICLE INFORMATION------------------------------
!..........................................................................

  if ( .not. auto_readin )then
       call get_orbit_info  
                     ! reads in file on single-particle orbits
!  call truncate_orbits   ! option to limit max orbit
       call extract_sps       !  inflates s.p. orbits to s.p. states
! ------------- sort single particle states by W, parity, Jz----------------
       call sort_spsqn(1)       ! protons
       call sort_spsqn(2)       ! neutrons
  endif
!------------------GET REST OF SYSTEM PARAMETERS-------------------
  call define_system
  call do_I_subsume_spes   ! ADDED 7.8.0 -- subsume s.p.e.s into 2-body
  call writelogfile('sys')
!-----------------set up hspsqn array---------------------------------
!-----------------these are s.p. states used to construct haikus------
  call make_hspsqn

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!------------------ CREATE BASIS--------------------------------------
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!----------- figure W limitations based upon truncations

  call master_w_limits

  call clocker('basis','start')

  call basismaster

  call sector_stats(1)
  call sector_stats(2)

  call clocker('basis','end')
  call writelogfile('bas')

  if(menu_char=='b')then
	  call basis_out4postprocessing
	  call BMPI_BARRIER(icomm,ierr)
	  call BMPI_FINALIZE(ierr)
	  if(iproc==0)print*,' Finished writing to file '
	  
	  stop
  end if
  if(auto_readin)call checkbasisdim

!- - - - - - - - - - - - -  IF BASIS IS VERY LARGE THEN BREAK INTO "FRAGMENTS" HERE - - - - - - -

  call fragmenter
  call sectorboundaries    ! new routine added 7.4.9 (see bbasislib6.f90) to define where sectors go
  if(nproc >1)call writelogfile('fra')
  
  if(print_basis .and. iproc==0)then
        print*,' Writing basis to file ... '
        call write_out_basis
  end if
  
  if(ham_readin .or. menu_char=='m')call threebodyquestion
  call getMPIprocs

!--------------- IMMEDIATE OPTION FOR COMPUTING OVERLAP-------

  if(menu_char == 'v')then
	  call overlap
	  call report
	  call BMPI_BARRIER(icomm,ierr)
	  call BMPI_FINALIZE(ierr)
	  if ( iproc == 0 ) write(6,*)'End of main program.'
	  stop
  end if
  if(menu_char == 'dx')then
     call BMPI_BARRIER(icomm,ierr)
     call density1b_from_oldwfn
     call report
     call BMPI_BARRIER(icomm,ierr)
     call BMPI_FINALIZE(ierr)
     if ( iproc == 0 ) write(6,*)'End of main program.'
     stop
  end if
! If we are saving eigenvectors, then write the header
! need fragments and mpi init first
  if(writeout .and. write_wfn) then
     call wfn_wopen_file(wfnfile,.false.)
     call write_wfn_header(wfnfile)
  end if

!------------------------- BASIS NOW COMPLETE --------------------------

  if(print_sectors)then
    call writeoutsectors(1)
    call writeoutsectors(2)
  endif
  if ( print4modelinfo ) then
     call writesectors4modeling
  endif

!------------------------ READ IN OPERATOR IF NEEDED 

  if ( op_readin ) then
     call readin1bodyop
     call decouple1bodyop
  end if

  if ( menu_char == 'o' .or. menu_char=='oo' ) then
     call BMPI_BARRIER(icomm,ierr)
     call applicator1b
     call report
     call BMPI_BARRIER(icomm,ierr)
     call BMPI_FINALIZE(ierr)
     if ( iproc == 0 ) write(6,*)'End of main program.'
	 stop
  end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!--------------------- SET UP TO CREATE JUMPS --------------
if (meanie)then
	call onebodysetup
else
     call clocker('ham','start')

     call hoopmaster ( hermflag, whermflagp, whermflagn, change1bodyqs, change2bodyqs, change3bodyqs)
     call clocker('ham','end')
     call clockout('ham')
     threebody0 = threebody
!----------------------- COMPUTE DISTRIBUTION OF WORK ---- 

    call timeperopmaster('set')
!	if(menu_char=='m')call memory_requirements(.false.)
    call master_para_distribute
!---------------- THE FOLLOWING ARE NOT YET NEEDED ------------------------------	
	call setsectorpossession(0)   ! routine to set flags denoting which sectors belong to which MPI process
!    call new_countupXYmes(0)   ! new routines for computing storage of PN matrix elements
!---------------------------------------------------	

   ! make sure bundles will not generate out of bounds references in
   ! actual use in bapplyhlib1.f90
    call check_opbundle_bounds
    call memory_requirements(.true.)   
    if(menu_char == 'm')then
            call clockout('bas')
            call clockout('hop')
            call clockout('jmc')
            call clockout('des')
            call clockout('mes')

            call BMPI_BARRIER(icomm,ierr)
            call BMPI_FINALIZE(ierr)
            stop
    end if
end if

!!------------------------ READ IN INTERACTION 

 if ( ham_readin ) then
     emptyHam=.true.   ! if no matrix elements, then many things aren't run

     call threebodymaster
     call master_readin_tbmes
     call pandamaster
	 if(subsume_spe)call subsume_sp_pot_into_2body

  end if

!----------------------- GENERATE JUMPS------------------
  if(menu_char /= 'p' )then   ! don't generate jumps if computing particle occupation
	                          ! don't uncouple matrix elements

   if(.not.meanie)then
      call clocker('ham','start')

      call jumpmaster(hermflag,whermflagp,whermflagn,change1bodyqs,change2bodyqs,change3bodyqs)
      call updateopbundles
    
      if(iproc==0)print*,' jumps built '
      call clocker('ham','end')
      call clockout('ham')

!- - - - - - - SET UP FOR PARALLEL RUNS - - - - - - - - - - - - - - - - - - 

! Call subroutine to sort jump arrays................................
      if ( np(1)*np(2) > 0  ) then
         call clocker('s1b','start')
         call sort_bundled_jumps
         call set_bundledXY_threadstart
         call clocker('s1b','end')
      end if
    end if   ! not meanie

!------- UNCOUPLE MATRIX ELEMENTS------------
    if ( ham_readin ) then
       call clocker('mun','start')
       if ( .not.emptyHam) then
           call makespeme(1,'H')
           call makespeme(2,'H')
           call uncoupleXXtbme(1,'H')
           call uncoupleXXtbme(2,'H') 
		   call delayedPNmatrixelements
           call count_uncouplepn('H',.true.)
        else
           if(iproc==0)print*,' Not decoupling; no 2-body matrix element file '
        end if

        if(threebody)then 
           if(.not.emptyHam)call masterconvert2to3bme
           call master3bodyuncoupler
        endif
        call clocker('mun','end')
     end if  ! ham_readin
   end if  ! menu_char not p
!......................................FINISHED MAKING JUMPS ..............

   if ( menu_char == 'n' .or. menu_char == 's' .or. menu_char == 'd' &
       .or. menu_char == 'ns' .or. menu_char=='ne' .or. menu_char == 'np') then
!.......... THIS IS A KLUGE....SHOULD PROPERLY PUT ELSEWHERE	   
       do is = 1,nsectors(1)  ! loop over proton sectors
         ndim = 0
         do isc = 1,xsd(1)%sector(is)%ncsectors
            jsc = xsd(1)%sector(is)%csector(isc)
            ndim = ndim + xsd(2)%sector(jsc)%nxsd
         end do ! isc
         xsd(1)%sector(is)%ncxsd = ndim
       end do
	   
       call BMPI_BARRIER(icomm,ierr)
       call bundle_clock(0,'set')
       call lanczos_menu(.false.)
   end if
  
   if(menu_char=='f')then
!.......... THIS IS A KLUGE....SHOULD PROPERLY PUT ELSEWHERE	   
	  do is = 1,nsectors(1)  ! loop over proton sectors
	           ndim = 0
	           do isc = 1,xsd(1)%sector(is)%ncsectors
	              jsc = xsd(1)%sector(is)%csector(isc)
	              ndim = ndim + xsd(2)%sector(jsc)%nxsd
	           end do ! isc
	           xsd(1)%sector(is)%ncxsd = ndim
	  end do	  
	  call scmf_master
   end if

   if( menu_char == 'r')then
      call BMPI_BARRIER(icomm,ierr)
      call bundle_clock(0,'set')
	  print*,' NEED TO REBUILD RESTART OPTION '
	  call lanczos_menu(.true.)   ! note: may not work correctly
   end if

   if ( menu_char == 'x' ) then
      call BMPI_BARRIER(icomm,ierr)
      call bundle_clock(0,'set')
	 
      call expectator_p
   end if

   if ( menu_char == 'p' ) then
      call BMPI_BARRIER(icomm,ierr)
      call particle_occupation_p
   end if

  if ( menu_char == 'a' ) then
     call BMPI_BARRIER(icomm,ierr)
     call applicator_h
  end if
!.............. ADDED 4/2011 by CWJ @ SDSU -- computes TRACE................
  if ( menu_char == 'c' ) then
     call BMPI_BARRIER(icomm,ierr)
     call tracemaster
  end if

  call clocker('all','end')
  call writelogfile('tim')  ! write summary information to 
  call clockout('all')
  call clockout('bas')
  call clockout('jmc')
  call clockout('des')
  call clockout('mes')
  call clockout('mun')

  call clockout('ham')

  if ( nproc > 1  )then
     call clockout('dis')
  end if
  call clockout('s1b')
  call clockout('piv')
  
  call clockout('lan')

  call clockout('hmu')
  call clockout('spe')
  if(threebody0)then  ! need to use threebody0 since threebody turned off for J,T
  call clockout('ppp')
  call clockout('ppn')
  call clockout('pnn')
  call clockout('nnn')

  else
  call clockout('pno')
  call clockout('pnb')
  call clockout('ppo')
  call clockout('ppb')

  call clockout('nno')
  endif
  call clockout('ort')
!  call clockout('rre')
!  call clockout('ror')
!  call clockout('dot')
!  call clockout('pro')
!  call clockout('swp')
  call clockout('obs')
  call clockout('aob')
  call clockout('egv')
  call clockout('eig')
  if(nproc > 1)call clockout('int')  ! ONLY DELETE INTRONS IF MPI 
  if( trdensout) call clockout('trd')

  if ( menu_char == 'n' ) call clockout('wev')
!  call bundle_clock_out  ! moved to main lanczos routine

  call report
  call BMPI_BARRIER(icomm,ierr)
  call BMPI_FINALIZE(ierr)
  if ( iproc == 0 ) write(6,*)'End of main program.'
  stop
end program bigstick_main

!!=====================================================================
! Examine environment and set up stuff:
!   paths like scratch_dir
!
!  called by main routine
!
subroutine main_getenv()
   use io
   use nodeinfo
   use bmpi_mod
   implicit none
   character(len=1024) :: buf
   integer:: estat
   integer:: ierr

   if(iproc == 0) then
      nersc_host = "none"
      scratch_dir = "."  ! use current directory
      call get_environment_variable("NERSC_HOST", buf, STATUS=estat)
      if(estat == 0) then
         nersc_host = TRIM(buf)
         if(nersc_host == "hopper" .or. nersc_host == "edison" .or. nersc_host == "cori") then
            call get_environment_variable("SCRATCH", buf, STATUS=estat)
            if(estat == 0) then
               scratch_dir = TRIM(buf)
            endif
         else
            nersc_host = "none"
         endif
      endif
      print *
   endif
   ! strings have implicit length fields that MPI won't see
   ! broadcast in common buffer, and assign after BCAST
   buf = nersc_host
   call BMPI_BCAST(buf,LEN(buf),0,icomm,ierr)
   nersc_host = TRIM(buf)
   buf = scratch_dir
   call BMPI_BCAST(buf,LEN(buf),0,icomm,ierr)
   base_scratch_dir = buf
   scratch_dir = buf   ! may add on output file as directory
   if(iproc == 0) then
      print *, "Running on NERSC_HOST: ", TRIM(nersc_host)
      print *, "scratch_dir (*.wfn,...): ", TRIM(scratch_dir)
   endif
end subroutine main_getenv


!=====================================================================
!
!  MAIN ROUTINE TO OPEN OUTPUT FILE
!
subroutine getoutputfile
  use menu_choices
  use program_info
  use io
  use nodeinfo
  use wfn_mod
  implicit none

  integer(4) :: ierr
  integer(4) :: ilast
  character(len=8) :: res_suf

1 continue
  if(iproc==0)then
     if(auto_input)then
        read(autoinputfile,'(a)')outfile
        write(6,*)outfile
     else 
        if(densityflag .or. trdensout)then
           print*,' Enter output name (required for your chosen option)!'
        else
           print*,' Enter output name (enter "none" if none)'
        endif
        read(5,'(a)')outfile
        write(autoinputfile,'(a)')outfile
	end if
	
  end if
  call BMPI_BARRIER(icomm,ierr)
  call BMPI_BCAST(outfile,20,0,icomm,ierr)
  ilast = index(outfile,' ')-1  
  if(densityflag .and. outfile(1:ilast)=='none')goto 1
  
  if(outfile(1:ilast)=='none')then
     writeout = .false.
  else
     writeout = .true.
     if ( iproc == 0 ) then
        res_suf = ".res"
        if(menu_char == 'dx') res_suf = ".dres"
        open(unit=resultfile,file=outfile(1:ilast)// TRIM(res_suf),status = 'unknown')
        write(resultfile,*)' BIGSTICK Version ',version,lastmodified
   	    select case(menu_char)
       	 case('dx')
   		 write(6,*)' Density matrices written to :', outfile(1:ilast)//".dres" 
   		 write(logfile,*)' Density matrices written to :', outfile(1:ilast)//".dres" 

       	 case('t')
   		 write(6,*)' Output written to :', outfile(1:ilast)//".res" 
   		 write(6,*)' Wfn info written to :', outfile(1:ilast)//".trwfn" 
   		 write(logfile,*)' Output written to :', outfile(1:ilast)//".res" 
   		 write(logfile,*)' Wfn info written to :', outfile(1:ilast)//".trwfn" 
		 
   		 case('d')
   		 write(6,*)' Output written to :', outfile(1:ilast)//".res" 
   		 write(logfile,*)' Output written to :', outfile(1:ilast)//".res" 
		 
   		 case default

        end select
		
     endif

     if(base_scratch_dir /= '.') then
        scratch_dir = TRIM(base_scratch_dir) // '/' // outfile(1:ilast)
        if(iproc == 0) print *, "set scratch_dir=", TRIM(scratch_dir)
        call system('mkdir -p ' // TRIM(scratch_dir))
     endif
  endif
  return
end subroutine getoutputfile

!===========================================================================
!
!  INPUT: printall = logical flag to print out all results
!       (if false, just print intermediate results)
!
! CALLED BY 
!   main routine
!
subroutine memory_requirements(printall)
   use flags3body
   use jumpNbody
   use jump3body
   use basis
   use precisions
   use nodeinfo
   use chains
   use interaction
   use interactions3body
   use jumplimits
   use jump_mod
   use fragments
   use bmpi_mod
   use butil_mod
   use para_main_mod
   use para_util_mod
   use io
   use operation_stats,only:totalops
   implicit none
   logical :: printall
   logical, parameter :: memory_all = .true. 
   integer(4) :: ierr
   integer :: nv
   integer(8) :: maxfragment
   integer :: i

  integer(kind = 8) :: n,nme

  if(iproc /= 0)return
  if(printall)then
     print*,' '
     print*,'  :   :   :   :   :  :  MEMORY ANALYSIS :   :   :   :   :   :   '
     print*,' '
  end if

   nme = 0
   print*,' (NOTE: THE COUNTS ARE APPROXIMATE AND SHOULD BE VERIFIED)'
   if(.not.threebody)then
      call write1bjmps4modeling(1,.false.,n)
      nme = nme+ n
      print*,' PN ',n
      call write2bjmps4modeling(1,.false.,n)
      nme = nme+ n
      print*,' PP ',n
      call write2bjmps4modeling(2,.false.,n)
      nme = nme+ n
      print*,' NN ',n
   else
      call writeXXXjmps4modeling(1,.false.,n)
      print*,' PPP ',n
      nme = nme+ n
      call writeXXXjmps4modeling(2,.false.,n)
      print*,' NNN ',n

      nme = nme+ n
      call writeXXYjmps4modeling(1,.false.,n)
      print*,' PPN ',n

      nme = nme+ n
      call writeXXYjmps4modeling(2,.false.,n)
      print*,' PNN ',n

      nme = nme+ n
   endif

   if(nfragments==1)then
      write(6,1)2*dimbasis*lanc_prec*1.0e-6
      write(logfile,1)2*dimbasis*lanc_prec*1.0e-6
   else
      maxfragment=0
      do i = 1,nfragments
         maxfragment = bmax(maxfragment,fragmentlist(i)%localdim)
      end do
      write(6,11)2*maxfragment*lanc_prec*1.0e-6
      write(logfile,11)2*maxfragment*lanc_prec*1.0e-6  
   end if
   ! Note - 13 is unchanged and wrong for 3 body jumps 
   write(6,2)1.0e-6*( bytesper1Bjump*(totn1bjumps+totp1bjumps)+ &
       bytesper2Bjump*(totn2bjumps+totp2bjumps) + bytesper3Bjump*(totp3bjumps + totn3bjumps) )
   write(logfile,2)1.0e-6*( bytesper1Bjump*(totn1bjumps+totp1bjumps)+ &
       bytesper2Bjump*(totn2bjumps+totp2bjumps) + bytesper3Bjump*(totp3bjumps + totn3bjumps) )
   if(nprocs > 1 .and. printall)then
	   write(6,6)maxjumpmemory
	   write(logfile,6)maxjumpmemory
   end if
   if(.not.threebody)then
      write(6,4) (float(nmatXX(1))+float(nmatXX(2))+float(nmatpn))*4e-6
      write(logfile,4) (float(nmatXX(1))+float(nmatXX(2))+float(nmatpn))*4e-6
      print*,' '
   else
      write(6,5)(nmatppp+nmatppn+nmatpnn+nmatnnn)*4e-6
      write(logfile,5)(nmatppp+nmatppn+nmatpnn+nmatnnn)*4e-6
   end if

   print*,' '

   write(6,*)' total # operations ',nme  !, totalops, totalops/nme
   write(logfile,*)' total # operations ',nme

   write(6,3)real(nme)*4.0e-9   ! this includes integers for final states
                                ! but also includes /2 for hermiticity
   write(logfile,3)real(nme)*4.0e-9
   print*,' Approximate time per iterations estimated :',totalops*1.0e-9,' sec, or ',totalops*1.0e-9/60.,' min '

   print*,' '
   print*,'  :   :   :   :   :  :  :   :   :   :   :   :   '
   print*,' '

   if(simulateMPI)call model_distribute_lanczos_pieces
   

1  format(' RAM for 2 lanczos vectors in storage         : ',f12.3,' Mb ')
11 format(' RAM for 2 lanczos vector fragments (max)         : ',f12.3,' Mb ')

2  format(' RAM for jumps in storage    (total)          : ',f15.3,' Mb ')

3  format(' Total nonzero matrix (including indices, assuming hermitian) would require ',g12.4,' Gb storage ')

4  format(' RAM for uncoupled two-body matrix elements   : ',f12.3,' Mb ')

5  format(' RAM for uncoupled three-body matrix elements : ',f12.3,' Mb ')

6  format(' Max RAM for local storage of jumps           : ',f15.3,' Mb ')

   return
end subroutine memory_requirements
!====================================================================

subroutine memerror(msg)
   use nodeinfo
   use bmpi_mod
   implicit none
   character (len=*) :: msg

   print *, "iproc=", iproc, "memerror in - ", msg
   flush(6)
   if(nproc > 1) call BMPI_ERR("memerror quit")
   stop 5
end subroutine memerror
