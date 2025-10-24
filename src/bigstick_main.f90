
!===========================================================================
!  BIGSTICK
!
!  Inspired by code REDSTICK by W. E. Ormand
!  
!  for history see file WHATSNEW.txt
!
!  A manual can be found under /docs or  arXiv:1801:08432
!
!  For references please cite
! C. W. Johnson, W. E. Ormand, and P. G. Krastev,  Comp. Phys. Comm. 184, 2761-2774 (2013)
! C. W. Johnson, W. E. Ormand, K. S. McElvain, and H.Z. Shan, arXiv:1801.08432
! (report UCRL LLNL-SM-739926)
!
! This code is distributed under MIT Open Source License. 
! Copyright (c) 2017 Lawrence Livermore National Security and San Diego State University Research Foundation
! for full licence statement see README_Licenses.txt
!
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
!  use hermit
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
  use sp_potter
  use configurations
  use threebodycentroids
  use jumpstart
  use wfn_organize
  use interaction,only:random_mscheme  ! seldom used
  implicit none

  integer(4) :: ierr
  integer(8) :: nme
  integer(4) :: i
  character(1) :: ychar
!----------- FLAGS FOR RESTRICTION ON CHANGING q#s in jumps-----------------
  logical    :: change1bodyqs, change2bodyqs, change3bodyqs
  logical    :: hermflag, whermflagp, whermflagn
  logical    :: whermflag2p, whermflag2n 
  logical threebody0
!--------   Thread info-----------------------------------------------------
  integer num_threads, mythread
  integer omp_get_num_threads
  integer omp_get_thread_num
#ifdef _MPI
  character(len=MPI_MAX_PROCESSOR_NAME) :: mpiprocname
#endif
  integer :: mpiprocnamelen
  
  integer :: is, isc, jsc
  integer(8) :: ndim
  

  ! KSM DBG: force standard in
  ! needed for pgdbg - no way to redirect stdin with MPI
  ! open(5, FILE='te123.minput', status='OLD', action='READ')
  ! open(5, FILE='xe136.minput', status='OLD', action='READ')

! UNCOMMENT IF YOU WANT BACKTRACE
  call installsignalhandlers  ! better reporting on crashes / interrupt
#ifdef _MPI  
  call BMPI_INIT(ierr)
!  icomm = MPI_COMM_WORLD
  call BMPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call BMPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
#else
  iproc=0
  nproc=1
#endif
  noisy0 = iproc == 0 .and. wantnoisy0
  
#ifdef _MPI
  call MPI_GET_PROCESSOR_NAME(mpiprocname, mpiprocnamelen, ierr)
#endif
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
     write(6,*)' BIGSTICK: a CI shell-model code, version ',version,lastmodified
     write(6,*)' '
     write(6,*)' Please cite: C. W. Johnson, W. E. Ormand, and P. G. Krastev '
     write(6,*)' Comp. Phys. Comm. 184, 2761-2774 (2013);  '
     write(6,*)' and C. W. Johnson, W. E. Ormand, K. S. McElvain, H.-Z. Shan '
     write(6,*)'  arXiv:1801.08432 and report UCRL LLNL-SM-739926 '
     write(6,*)' '
	 write(6,*)' This code distributed under the MIT Open Source License '

  end if

  ! read environment and setup paths like scratch_dir
  call main_getenv()  ! needs MPI_INIT

!  if ( iproc == 0 ) write(6,*) 'Number of MPI processors =',nproc
!$omp parallel shared(num_threads,iproc)  
   num_threads = omp_get_num_threads()
   num_threads_global=num_threads
   mythread = omp_get_thread_num()
   if ( iproc == 0 .and. mythread == 0 ) write(6,*)'Number of MPI processors =',nproc, ', NUM_THREADS = ', num_threads

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
!  hermitian =  .true.
  hermflag = .true.
  whermflagp = .true.   ! for protons
  whermflagn = .false.  ! for neutrons -- ALWAYS TURN OFF
  
  whermflag2p = .true.
  whermflag2n = .true.

  simulateMPI = simulateMPIdefault

  if ( op_readin ) then
     change1bodyqs = .false.
  else
     change1bodyqs = .true.
  endif
  change2bodyqs = .false.     ! NB this might change with 3-body forces 

  change3bodyqs = .false.
 
!--- in 7.10.4, changing this approach  
!  if(dens2bflag)then
!	  hermflag = .false.
!	  whermflagp = .false.
!	  whermflag2p =  .false.
!	  whermflag2n = .false.
!  end if

!-----------------  Set a value defining type of Lanczos vector for MPI routines
!  if ( lanc_prec == 4 ) then
!     MPI_lanc = MPI_REAL
!  else
!     MPI_LANC = MPI_REAL8
!  end if
!  modifytimewts = .true.

  if(menu_char/='r' .and. menu_char /= 'c' .and. menu_char/='v' .and. menu_char/='b' &
       .and. menu_char/='wh' .and. menu_char/='wo')call getoutputfile
  call openlogfiles
  call writelogfile('ini')

!-------------- OPEN JUMPSTART FILE (added 7.9.6) -------------------
   maxfragmentsize = 0
   if(writejumpstart .and. nproc > 1)then    ! open jupmpstart file to write to
	  call openjumpstartfile2write  
   end if
   
!   if(.not.writejumpstart .and. nproc > 1)then
!  	  call openjumpstartfile2read   !checks if there is a jumpstart file
!	  call BMPI_Barrier(icomm,ierr)
!	  if(readjumpstart)then
!		  call readbaseinfo4jumpstart
!	  end if
!   end if  
   
!..........................................................................
!------------------SINGLE-PARTICLE INFORMATION------------------------------
!..........................................................................

  if ( .not. auto_readin .and. .not. readjumpstart )then
	  
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

  if(stopafterdim) then  ! if you want to only compute the density
	  close(autoinputfile)

	  call clocker('all','end')

	  call clockout('all')
	  call clockout('bas')
	  call report
	  call writelogfile('end')
#ifdef _MPI	  
	  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
	  call BMPI_FINALIZE(ierr)
#endif	  
	  if ( iproc == 0 ) write(6,*)'End of main program.'
	  stop
  end if
  if(menu_char=='b')then
	  call basis_out4postprocessing
#ifdef _MPI	  
	  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
	  call BMPI_FINALIZE(ierr)
#endif	  
	  if(iproc==0)print*,' Finished writing to file '
	  
	  stop
  end if
  
  if(auto_readin .or. readjumpstart)call checkbasisdim

!- - - - - - - - - - - - -  IF BASIS IS VERY LARGE THEN BREAK INTO "FRAGMENTS" HERE - - - - - - -

  call fragmenter  
  call sectorboundaries    ! new routine added 7.4.9 (see bbasislib6.f90) to define where sectors go
  
  if(nproc >1)call writelogfile('fra')
  
  if(print_basis .and. iproc==0)then
        print*,' Writing basis to file ... '
        call write_out_basis
  end if
  
! - - - - - - - ADDED in 7.10.7 - - - - - OPTION TO CONVERT .trwfn TO .wfn

  if(menu_char=='tw')then
      call wfn_wopen_file(wfnfile,.false.)	 
      call write_wfn_header(wfnfile)	 	  
	  call readin_trwfn
!..... WRITE TO WFN....
	  call report
#ifdef _MPI	  
	  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
	  call BMPI_FINALIZE(ierr)
#endif
	  if ( iproc == 0 ) write(6,*)'End of main program.'	  
	  stop 
  end if  
  
  if(ham_readin .or. (menu_char=='m' .and. .not. modeldensities) .or. menu_char=='c')then
	  call threebodyquestion
	  if(.not.threebody .and. enable_centroids)call centroidquestion
  end if
  
  if(applycentroids)call bossconfigurator(.false.)

  call getMPIprocs
!
#ifdef _MPI  
  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!--------------- IMMEDIATE OPTION FOR COMPUTING OVERLAP-------

  if(menu_char == 'v')then
	  call overlap(.false.)
	  call report
#ifdef _MPI	  
	  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
	  call BMPI_FINALIZE(ierr)
#endif
	  if ( iproc == 0 ) write(6,*)'End of main program.'
	  stop 
  end if
  
!-------------- IMMEDIATE OUTPUT FOR READING IN .wfn file, WRITING AS TRDENS FILE ----
! added 7.10.7

if(menu_char=='tx')then
  call output_TRDENS(.true.)
	
  call report
#ifdef _MPI  
  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
  call BMPI_FINALIZE(ierr)
#endif
  if ( iproc == 0 ) write(6,*)'End of main program.'
  stop  
end if

!------ OPTION TO READ IN MULTIPLE FILES AND ORTHONORMALIZE-----

if(menu_char=='ro' .or. menu_char=='ru')then
    call wfn_wopen_file(wfnfile,.false.)	 
    call write_wfn_header(wfnfile)	
	if(menu_char=='ro')then
		call read_n_ortho_boss(.true.)
	else
		call read_n_ortho_boss(.false.)
	endif
    call report
#ifdef _MPI  
    call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
    call BMPI_FINALIZE(ierr)
#endif
    if ( iproc == 0 ) write(6,*)'End of main program.'
    stop  	
	
end if
  
!---------------  OPTION FOR COMPUTING DENSITIES from prior runs -------
  
  if(menu_char == 'dx')then
#ifdef _MPI 	  
     call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	 call clocker('den','sta')

     call density1b_from_oldwfn
 	 call clocker('den','end')
	 
     call report
#ifdef _MPI 
     call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	 call clockout_den1b

#ifdef _MPI 	 
     call BMPI_FINALIZE(ierr)
#endif	 
     if ( iproc == 0 ) write(6,*)'End of main program.'
     stop
  end if
  if(menu_char == 'du')then
#ifdef _MPI 
     call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	 call clocker('den','sta')

     call uncoupled_density1b_from_oldwfn
 	 call clocker('den','end')
	 
     call report
#ifdef _MPI 	 
     call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	 call clockout_den1b
#ifdef _MPI 	 
     call BMPI_FINALIZE(ierr)
#endif
     if ( iproc == 0 ) write(6,*)'End of main program.'
     stop
  end if
  if(menu_char=='m' .and. modeldensities)then
	  call onebodysetup
	  call clockout_den1b
      call report
#ifdef _MPI 
      call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
      call BMPI_FINALIZE(ierr)
#endif	  
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

!  if(writejumpstart)then
!	  call write_wfn_header(jumpstartfile)
!	  call write_info2jumpstart1
!  end if  


#ifdef _MPI
   call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
  if(print_sectors)then
    call writeoutsectors(1)
    call writeoutsectors(2)
  endif
  if ( print4modelinfo ) then
     call writesectors4modeling
  endif
  
!--------------------- SET UP CONFIGURATIONS ---------
  if(configout)then
	  if(iproc==0)then
	     print*,' '
	     print*,' Generating configurations '
	     print*,' '
      end if
	  call bossconfigurator(.false.)
  end if  

!------------------------ READ IN OPERATOR IF NEEDED 

  if ( op_readin .and. .not.uncoupled_dens) then
     call readin1bodyop
     call decouple1bodyop
  end if
  if ( op_readin .and. uncoupled_dens) then
     call readin1bodyop_uncoupled
  end if
  if ( menu_char == 'o' .or. menu_char=='oo' .or.  menu_char=='wo') then
#ifdef _MPI  
     call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	 if(printouttrans1flag)then
		 call applicator1b_print
	 else
         call applicator1b
	 end if
     call report
#ifdef _MPI  	 
     call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
     call BMPI_FINALIZE(ierr)
#endif
     if ( iproc == 0 ) write(6,*)'End of main program.'
	 stop
  end if
!--------------------- SET UP TO CREATE JUMPS --------------
if(menu_char/='cx')then
	
     call clocker('ham','start')
     call hoopmaster ( hermflag, whermflagp, whermflagn,whermflag2p,whermflag2n, change1bodyqs, change2bodyqs, change3bodyqs)
	 
     call clocker('ham','end')
     call clockout('ham')
	 
     threebody0 = threebody
!----------------------- COMPUTE DISTRIBUTION OF WORK ---- 

    call timeperopmaster('set')
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
#ifdef _MPI  
            call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
            call BMPI_FINALIZE(ierr)
#endif
            stop
    end if


!!------------------------ READ IN INTERACTION 


 if ( ham_readin .or. menu_char=='2b') then
     emptyHam=.true.   ! if no matrix elements, then many things aren't run

     call threebodymaster
     call master_readin_tbmes(.false.)
	 
	 if(applycentroids)then
		 call read_and_fill_3bodycentroidpots
		 call combinedconfigs(.true.)
	 end if
	 if(ham_readin )then
         call pandamaster
	     if(subsume_spe)call subsume_sp_pot_into_2body
	 end if
  end if
end if ! if menu_char /= cx

if(menu_char=='jp')then ! added 7.11.3
    call master_readin_tbmes(.true.)
	if(iproc==0)print*,' J^2 set up '

    call bundle_clock(0,'set')
	emptyHam = .false.
	
end if
  
!----------------------- GENERATE JUMPS------------------
  if(menu_char /= 'p' .and. menu_char/='cx')then   ! don't generate jumps if computing particle occupation
	                          ! don't uncouple matrix elements

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
         call sort_bundled_jumps  !ERROR-- JUMP ARRAYS NOT ALLOCATED
         call set_bundledXY_threadstart
         call clocker('s1b','end')
      end if

!------- UNCOUPLE MATRIX ELEMENTS------------
    if ( ham_readin .or. menu_char=='jp' ) then
       call clocker('mun','start')
       if ( .not.emptyHam .or. menu_char=='c') then
           call makespeme(1,'H')
           call makespeme(2,'H')
           call uncoupleXXtbme(1,'H')
           call uncoupleXXtbme(2,'H') 
		   call delayedPNmatrixelements
           call count_uncouplepn('H',.true.)
        else
           if(iproc==0)print*,' Not decoupling; no 2-body matrix element file '
        end if
		
		if(random_mscheme)call randomize_tbme  ! specialized call, seldom used

        if(threebody)then 
           if(.not.emptyHam)call masterconvert2to3bme
           call master3bodyuncoupler
        endif
        call clocker('mun','end')
     end if  ! ham_readin
   end if  ! menu_char not p
!......................................FINISHED MAKING JUMPS ..............

   if ( menu_char == 'n' .or. menu_char == 's' .or. menu_char == 'd' .or. menu_char=='g' &
       .or. menu_char == 'ns' .or. menu_char=='ne' .or. menu_char == 'np' .or. menu_char=='wh' & 
	   .or. menu_char == 'pv' ) then
!.......... THIS IS A KLUGE....SHOULD PROPERLY PUT ELSEWHERE	   
       do is = 1,nsectors(1)  ! loop over proton sectors
         ndim = 0
         do isc = 1,xsd(1)%sector(is)%ncsectors
            jsc = xsd(1)%sector(is)%csector(isc)
            ndim = ndim + xsd(2)%sector(jsc)%nxsd
         end do ! isc
         xsd(1)%sector(is)%ncxsd = ndim
       end do
#ifdef _MPI  	   
       call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
       call bundle_clock(0,'set')
       call lanczos_menu(.false.)
	   
   end if
   
   if(configout)call writeoutconfigocc
   if(menu_char == '2b')then    !ADDED IN 7.9.1
!.......... THIS IS A KLUGE....SHOULD PROPERLY PUT ELSEWHERE	   
       do is = 1,nsectors(1)  ! loop over proton sectors
         ndim = 0
         do isc = 1,xsd(1)%sector(is)%ncsectors
            jsc = xsd(1)%sector(is)%csector(isc)
            ndim = ndim + xsd(2)%sector(jsc)%nxsd
         end do ! isc
         xsd(1)%sector(is)%ncxsd = ndim
       end do
#ifdef _MPI  	   	   
       call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
 	   call boss_twobody_densities
	  
 	   call report
#ifdef _MPI  	   
       call BMPI_FINALIZE(ierr)
#endif
       if ( iproc == 0 ) write(6,*)'End of main program.'
       stop
	  
   end if


   if( menu_char == 'r')then
#ifdef _MPI  	   
      call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
      call bundle_clock(0,'set')
	  print*,' NEED TO REBUILD RESTART OPTION '
	  call lanczos_menu(.true.)   ! note: may not work correctly
   end if

   if ( menu_char == 'x' ) then
#ifdef _MPI  
      call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
      call bundle_clock(0,'set')
	 
      call expectator_p
   end if

   if ( menu_char == 'p' ) then
#ifdef _MPI  	   
      call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
      call particle_occupation_p
   end if

  if ( menu_char == 'a' ) then
!.......... THIS IS A KLUGE....SHOULD PROPERLY PUT ELSEWHERE	   
       do is = 1,nsectors(1)  ! loop over proton sectors
         ndim = 0
         do isc = 1,xsd(1)%sector(is)%ncsectors
            jsc = xsd(1)%sector(is)%csector(isc)
            ndim = ndim + xsd(2)%sector(jsc)%nxsd
         end do ! isc
         xsd(1)%sector(is)%ncxsd = ndim
       end do
	   
#ifdef _MPI  	   
       call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif	   
       call bundle_clock(0,'set')	  
#ifdef _MPI  	   
     call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
     call applicator_h
  end if
!.............. ADDED 4/2011 by CWJ @ SDSU -- computes TRACE................
  if ( menu_char == 'c' ) then
#ifdef _MPI  	  
     call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif	 
     call tracemaster
  end if
!....... ADDED 7/2024 by CWJ in 7.11.1 -- compute matrix elements in read-in basis

  if(menu_char=='h')then
      call bundle_clock(0,'set')
	  
	  call hmatrixelements_boss
	  
  end if
  
!.... ADDED in 7.11.3.....

if(menu_char=='jp')then
!.......... THIS IS A KLUGE....SHOULD PROPERLY PUT ELSEWHERE	   
	       do is = 1,nsectors(1)  ! loop over proton sectors
	         ndim = 0
	         do isc = 1,xsd(1)%sector(is)%ncsectors
	            jsc = xsd(1)%sector(is)%csector(isc)
	            ndim = ndim + xsd(2)%sector(jsc)%nxsd
	         end do ! isc
	         xsd(1)%sector(is)%ncxsd = ndim
	       end do
	   
#ifdef _MPI  	   
	       call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif	   
	
	call Jproject_boss	
	
	
end if  
  
  close(autoinputfile)

  call clocker('all','end')
  call detailed_hmult_timing
!  call writelogfile('tim')  ! write summary information to 
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
  if(nproc > 1)call clockout('int')  ! ONLY DELETE INTRONS IF MPI 
  
  call clockout('piv')
  
  call clockout('lan')

  call clockout('hmu')
  call clockout('spe')
  if(threebody0 )then  ! need to use threebody0 since threebody turned off for J,T
  call clockout('ppp')
  call clockout('pXf')  ! modified 7.9.2
  call clockout('pXb')
  call clockout('pYf')
  call clockout('pYb')
  call clockout('nnn')

else
  call clockout('pno')
  call clockout('pnb')
  call clockout('ppo')
  call clockout('ppb')

  call clockout('nno')
  endif
  call clockout('ort')
  call clockout('blr')
!  call clockout('rre')
!  call clockout('ror')
!  call clockout('dot')
!  call clockout('pro')
!  call clockout('swp')

  if(menu_char=='o')call clockout('aob')
  call clockout('egv')
  call clockout('eig')
  call clockout('obs')
  
  if( trdensout) call clockout('trd')

  if ( menu_char == 'n' .or. menu_char=='wh') call clockout('wev')
  if(densityflag)then
      call clockout('den')
     call clockout('obw')
  end if
  
  if(writejumpstart .and. iproc==0)then
	  call draft_bundle_clock_out
	  call writeallinfo2jumpstart
!	  call write_opwt2jumpstart
	  
	  close(jumpstartfile)  ! close the jumpstart file
  end if
  if(readjumpstart .and. iproc==0)call compare_bundle_clock_out
  call report
  call writelogfile('end')
#ifdef _MPI    
  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
  call BMPI_FINALIZE(ierr)
#endif  
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
#ifdef _MPI     
   call BMPI_BCAST(buf,LEN(buf),0,MPI_COMM_WORLD,ierr)
#endif
   nersc_host = TRIM(buf)
   buf = scratch_dir
#ifdef _MPI     
   call BMPI_BCAST(buf,LEN(buf),0,MPI_COMM_WORLD,ierr)
#endif
   base_scratch_dir = buf
   scratch_dir = buf   ! may add on output file as directory
   if(iproc == 0) then
      print *, "Running on NERSC_HOST: ", TRIM(nersc_host), ", scratch_dir (*.wfn,...): ", TRIM(scratch_dir)
   endif
end subroutine main_getenv


!=====================================================================
!
!  MAIN ROUTINE TO OPEN OUTPUT FILE
!
! called by main routine
!
subroutine getoutputfile
  use menu_choices
  use program_info
  use io
  use nodeinfo
  use wfn_mod
  use jumpstart
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
        if((densityflag.and. .not.modeldensities) .or. trdensout .or. menu_char=='2b' .or.greenflag)then
           print*,' Enter output name (required for your chosen option)!'
        else
           print*,' Enter output name (enter "none" if none)'
        endif
        read(5,'(a)')outfile		
        write(autoinputfile,'(a)')outfile
	end if
	
  end if
#ifdef _MPI  
  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(outfile,20,0,MPI_COMM_WORLD,ierr)
#endif  
  ilast = index(outfile,' ')-1  
  if(((densityflag .and. .not.modeldensities) .or. writejumpstart) .and. outfile(1:ilast)=='none')goto 1
  
  if(outfile(1:ilast)=='none')then
     writeout = .false.
  else
     writeout = .true.
     if ( iproc == 0 ) then
        res_suf = ".res"
        if(menu_char == 'dx') res_suf = ".dres"
        if(menu_char == 'du') res_suf = ".udres"
		if(menu_char == 'h')  res_suf = ".xme"

        open(unit=resultfile,file=outfile(1:ilast)// TRIM(res_suf),status = 'unknown')
        write(resultfile,*)' BIGSTICK Version ',version,lastmodified
   	    select case(menu_char)
       	 case('dx')
   		 write(6,*)' Density matrices written to :', outfile(1:ilast)//".dres" 
   		 write(logfile,*)' Density matrices written to :', outfile(1:ilast)//".dres" 
         open(unit=occresultfile,file=outfile(1:ilast)//".occres",status = 'unknown')
		 write(logfile,*)' Single-particle occupations written to : ',outfile(1:ilast)//".occres" 
		 write(resultfile,*)' '
         write(occresultfile,*)' BIGSTICK Version ',version,lastmodified
		 write(occresultfile,*)' '
		 
		 case('du')
   		 write(6,*)' Density matrices written to :', outfile(1:ilast)//".udres" 
   		 write(logfile,*)' Density matrices written to :', outfile(1:ilast)//".udres" 
		 write(resultfile,*)' '

       	 case('t')
   		 write(6,*)' Output written to :', outfile(1:ilast)//".res" 
   		 write(6,*)' Wfn info written to :', outfile(1:ilast)//".trwfn" 
   		 write(logfile,*)' Output written to :', outfile(1:ilast)//".res" 
   		 write(logfile,*)' Wfn info written to :', outfile(1:ilast)//".trwfn" 
		 
   		 case('d')
   		 write(6,*)' Output written to :', outfile(1:ilast)//".res" 
   		 write(logfile,*)' Output written to :', outfile(1:ilast)//".res" 
         open(unit=denresultfile,file=outfile(1:ilast)//".dres",status = 'unknown')
         open(unit=occresultfile,file=outfile(1:ilast)//".occres",status = 'unknown')
		 
         write(denresultfile,*)' BIGSTICK Version ',version,lastmodified
		 write(denresultfile,*)' '
         write(occresultfile,*)' BIGSTICK Version ',version,lastmodified
		 write(occresultfile,*)' '
   		 write(6,*)' Density matrices written to :', outfile(1:ilast)//".dres" 
   		 write(logfile,*)' Density matrices written to :', outfile(1:ilast)//".dres" 	 
		 write(6,*)' Single-particle occupations written to : ',outfile(1:ilast)//".occres" 
		 write(logfile,*)' Single-particle occupations written to : ',outfile(1:ilast)//".occres" 
		 
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
! SUBROUTINES CALLED:
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

  integer(kind = 8) :: n,nme,nmepp,nmenn,nmepn,nmeppp,nmeppn,nmepnn,nmennn

  if(iproc /= 0)return
  if(printall)then
!     print*,' '
     print*,'  :   :   :   :   :  :  MEMORY ANALYSIS :   :   :   :   :   :   '
     print*,' '
  end if

   nme = 0
!   print*,' (NOTE: THE COUNTS ARE APPROXIMATE AND SHOULD BE VERIFIED)'
   if(.not.threebody)then
      call write1bjmps4modeling(1,.false.,n)
      nme = nme+ n
	  nmepn = n
      call write2bjmps4modeling(1,.false.,n)
      nme = nme+ n
	  nmepp = n
      call write2bjmps4modeling(2,.false.,n)
      nme = nme+ n
	  nmenn = n
   else
      call writeXXXjmps4modeling(1,.false.,n)
      nme = nme+ n
	  nmeppp = n
	  
      call writeXXXjmps4modeling(2,.false.,n)
	  nmennn = n

      nme = nme+ n
      call writeXXYjmps4modeling(1,.false.,n)
	  nmeppn = n
	
      nme = nme+ n
      call writeXXYjmps4modeling(2,.false.,n)
	  nmepnn = n

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
      write(6,21)
      write(logfile,11)2*maxfragment*lanc_prec*1.0e-6  
	  write(logfile,21)
   end if
   ! Note - 13 is unchanged and wrong for 3 body jumps 
   write(6,2)1.0e-6*( bytesper1Bjump*(totn1bjumps+totp1bjumps)+ &
       bytesper2Bjump*(totn2bjumps+totp2bjumps) + bytesper3Bjump*(totp3bjumps + totn3bjumps) )
   write(logfile,2)1.0e-6*( bytesper1Bjump*(totn1bjumps+totp1bjumps)+ &
       bytesper2Bjump*(totn2bjumps+totp2bjumps) + bytesper3Bjump*(totp3bjumps + totn3bjumps) )
   if(nprocs > 1 .and. printall)then
	   write(6,6)maxjumpmemory*1e-6
	   write(logfile,6)maxjumpmemory*1e-6
   end if
   if(.not.threebody)then
      write(6,4) (float(nmatXX(1))+float(nmatXX(2))+float(nmatpn) )*4e-6
      write(logfile,4) (float(nmatXX(1))+float(nmatXX(2))+float(nmatpn))*4e-6
	  write(6,*)' (Breakdown of storage for uncoupled matrix elements)'
	  write(6,*)' PP: ',float(nmatXX(1))*4e-6,' Mb'
	  write(6,*)' NN: ',float(nmatXX(2))*4e-6,' Mb'
	  write(6,*)' PN: ',float(nmatpn)*4e-6,' Mb'
	  
      print*,' '
   else
      write(6,5)(nmatppp+nmatppn+nmatpnn+nmatnnn)*4e-6
      write(logfile,5)(nmatppp+nmatppn+nmatpnn+nmatnnn)*4e-6
	  write(6,*)' (Breakdown of storage for uncoupled matrix elements)'
	  write(6,*)' PPP: ',float(nmatppp)*4e-6,' Mb'
	  write(6,*)' NNN: ',float(nmatnnn)*4e-6,' Mb'
	  write(6,*)' PPN: ',float(nmatppn)*4e-6,' Mb'
	  write(6,*)' PNN: ',float(nmatpnn)*4e-6,' Mb'
	  print*,' '
   end if
   
!   print*,' '
   print*,' Counts of many-body matrix elements (operations)'
   if(.not.threebody)then

      print*,' PN ',nmepn
      print*,' PP ',nmepp
      print*,' NN ',nmenn
   else
      print*,' PPP ',nmeppp
      print*,' NNN ',nmennn
      print*,' PPN ',nmeppn
      print*,' PNN ',nmepnn
   endif
   print*,' '

   write(6,*)' total # operations ',nme ,', ~ ops/jump = ', & 
    real(nme,kind=4)/real(totn1bjumps+totp1bjumps+totn2bjumps+totp2bjumps+totp3bjumps + totn3bjumps,kind=4)
	write(6,*)' Effective storage per operation = ', real(bytesper1Bjump*(totn1bjumps+totp1bjumps)+ &
       bytesper2Bjump*(totn2bjumps+totp2bjumps) + bytesper3Bjump*(totp3bjumps + totn3bjumps),kind=8) /real(nme,kind=8),   ' bytes '
   write(logfile,*)' total # operations ',nme,' ~ ops/jump = ', & 
    real(nme,kind=4)/real(totn1bjumps+totp1bjumps+totn2bjumps+totp2bjumps+totp3bjumps + totn3bjumps,kind=4)
	write(logfile,*)' Effective storage per operation = ', real(bytesper1Bjump*(totn1bjumps+totp1bjumps)+ &
       bytesper2Bjump*(totn2bjumps+totp2bjumps) + bytesper3Bjump*(totp3bjumps + totn3bjumps),kind=8) /real(nme,kind=8),   ' bytes '

   write(6,3)real(nme)*4.0e-9   ! this includes integers for final states
                                ! but also includes /2 for hermiticity
   write(logfile,3)real(nme)*4.0e-9
   print*,' '
   
   print*,' Approximate time per iterations estimated :',totalops*1.0e-9/real(nprocs),' sec, or ',& 
      totalops*1.0e-9/60./real(nprocs),' min '
   print*,' Assuming good efficiency on ',nprocs,' MPI processes '
   print*,   '(does not include OpenMP threads, and may not be using most up-to-date timing)'

   print*,' '
   print*,'  :   :   :   :   :  :  :   :   :   :   :   :   '
   print*,' '

   if(simulateMPI)call model_distribute_lanczos_pieces
   

1  format(' RAM for 2 lanczos vectors in storage         : ',f12.3,' Mb ')
11 format(' RAM for 2 lanczos vector fragments (max)         : ',f12.3,' Mb ')
21 format(' Note: if using block lanczos, must multiply by # vectors/block ')
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
#ifdef _MPI
   if(nproc > 1) call BMPI_ERR("memerror quit")
#endif
   
   stop 5
end subroutine memerror
