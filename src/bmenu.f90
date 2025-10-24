!===================================================================
!
!  subroutine menu
!
!  master routine to set up flags to direct flow of program
!
! SUBROUTINES CALLED:
!	wfn_ropen_file  
!     	call read_wfn_header
!........................................................................
subroutine menu

  use menu_choices
  use io
  use verbosity
  use nodeinfo
  use flagger
  use wfn_mod
  use bmpi_mod
  use jumpNbody
  use onebodypot
  use densities, only : pndensities
  use coupledmatrixelements,only: dens2bflag,diag_den2b
  use configurations,only:configout
  use flags3body,only:threebodycheck
  use TRstuff
  use jumpstart
  use tracy
  implicit none

  integer(4) :: ierr
  character :: ychar

! File to log interactive inputs to so the user can create an input file
  open(unit=autoinputfile,file='autoinput.bigstick',status='unknown')  ! always do this

  auto_input = .false.
  pndensities = .false.
  dens2bflag = .false.
  diag_den2b = .false.
  useTRphase = useTRphase_def
  dotflag = .false.
  binden  = .false.
  
  if ( iproc == 0 ) then
     print*,' '
     print*,' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ' 
     print*,' *                                                                         * '
     print*,' *               OPTIONS (choose one)                                      * '
     print*,' * (i) Input automatically read from "autoinput.bigstick" file             * '
     print*,' *  (note: autoinput.bigstick file created with each nonauto run)          * '
     print*,' * (n) Compute spectrum (default); (ns) to suppress eigenvector write up   * '
!     print*,' * (j) Carry out jumpstart run to set MPI timing (MPI only)                * '
     print*,' * (d) Densities: Compute spectrum + all one-body densities (isospin fmt)  * '
     print*,' * (2) Two-body density from previous wfn (default p-n format)             * '

!     print*,' * (dx[m]) Densities: Compute one-body densities from previous run (.wfn)  * '
!     print*,' *     optional m enables mathematica output                               * '
!     print*,' * (p) Compute spectrum + single-particle occupations; (ps) to suppress wfn* '
!     print*,' * (occ) single-particle occupations (from previous wfn)                   * '
     print*,' * (x) eXpectation value of a scalar Hamiltonian (from previous wfn)       * '
     print*,' * (o) Apply a one-body (transition) operator to previous wfn and write out* '
     print*,' * (s) Strength function (using starting pivot )                           * '
!     print*,' * (r) Restart from a previous run                                         * '
!     print*,' * (a) Apply a scalar Hamiltonian to a previous wfn and write out          * '
     print*,' * (g) Apply the resolvent 1/(E-H) to a previous wfn and write out         * '  ! ADDED 7.8.8; restarted 7.9.10
	 
!     print*,' * (y) Cycle back to this menu after calculation completed (not active)    * '
!     print*,' * (v) Overlap of initial states with final states                         * '
     print*,' * (m) print information for Modeling parallel distribution                * '
!     print*,' * (b) Create file with full basis information (for postprocessing)        * '
!     print*,' * (f) Self-consistent mean-field approximation (prepare pivot)            * '    ! ADDED 7.7.4
     print*,' * (l) print license and copyright information                             * '
	 
     print*,' * (?) Print out all options                                               * '
     print*,' *                                                                         * '
     print*,' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ' 
  end if
1 continue

  if( auto_input)then
     read(autoinputfile,'(a)')menu_char
  else  if ( iproc == 0 ) then
     print*,' '
     print*,' Enter choice '
     read(5,'(a)')menu_char
     if(menu_char /= 'i' .and. menu_char /= 'I')then
         write(autoinputfile,'(a,"    ! menu choice ")')menu_char
     end if
  end if

#ifdef _MPI
  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(menu_char,3,0,MPI_COMM_WORLD,ierr)
#endif
  auto_readin = .false.
  ham_readin = .false.
  op_readin  = .false.
  print4modelinfo = .false.
  strengthflag = .false.
  greenflag = .false.
  densityflag = .false.
  trdensout   = .false.
  spoccflag   = .false.
  menu_dx_omathematica = .true.  ! flag enabling mathematica output of density matrices
  get_JT       = .true.
  onebodyonly  = .false.
  printouthamflag = .false.
  printouttrans1flag = .false.
  baseASCII=.false.
  configout = .false.
  modeldensities=.false.
  uncoupled_dens = .false.
#ifdef _BD3    
  threebodycheck=.true.   ! a beta version in 7.9.2; ! added precompile option in 7.11.1
#else
   threebodycheck=.false.   ! a beta version in 7.9.2
#endif  
  blockstrengthflag=.false.
  
  writejumpstart = .false.
  readjumpstart  = .false.
  strengthnorm = .true.
  dotflag = .false.
  stopafterdim = .false.
  complex_green = .false.
  only_odd = .false.
  
  select case (menu_char)

  case ('?')
     if ( iproc == 0 ) then
        print*,' '
        print*,' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ' 
        print*,' *                                                                         * '
        print*,' *               OPTIONS (choose 1)                                        * '
        print*,' * (i) Input automatically read from "autoinput.bigstick" file             * '
        print*,' *  (note: autoinput.bigstick file created with each nonauto run)          * '
        print*,' * (n) Compute spectrum (default); (ns) to suppress eigenvector write up   * '
        print*,' * (ne) Compute energies but NO observable (i.e. J or T)                   * '   !ADDED 7.6.8
        print*,' * (np) Compute spectrum starting from prior pivot                         * '   !ADDED 7.7.7. modified in 7.9.11
!        print*,' * (j) Carry out jumpstart run to set MPI timing (MPI only)                * '
        print*,' * (d) Densities: Compute spectrum + all one-body densities (isospin fmt)  * '
        print*,' * (dx[m]) Densities: Compute one-body densities from previous run (.wfn)  * '
        print*,' *     optional m enables mathematica output                               * '
        print*,' * (dp) Densities in proton-neutron format                                 * '
        print*,' * (dxp) Compute one-body densities from prior run (.wfn) in p-n format.   * '
		print*,' * (db) Write one-body densities to a binary file                          * ' ! added 7.10.9
		print*,' * (dxb) Compute one-body densities from prior run, write to a binary file * ' ! added 7.10.9
!		print*,' * (do) Compute proton-neutron densities only for odd-J levels             * ' ! ADDED 7.10.8; DEPRECATED 7.11/1
!		print*,' * (dxo) Compute proton-neutron densities only for odd-J from prior run    * ' ! ADDED 7.10.8
		
!        print*,' * (du) Compute uncoupled one-body densities from prior run (.wfn)         * ' ! SECRET MENU ITEM
		
        print*,' * (2) Two-body density from previous wfn (default p-n format)             * '  ! ADDED 7.9.1
        print*,' * (2d) Two-body density from previous wfn, only initial=final, Jt=0       * '  ! added 7.9.5
        print*,' * (2i) Two-body density from previous wfn ( isopin format)       * '  ! NOT YET IMIPLEMENTED

#ifdef _BD3    
        print*,' * (3) Normal spectrum but using three-body forces (beta version)          * '
#endif  		
        print*,' * (p) Compute spectrum + single-particle occupations,(ps) to supress wfn  * '
        print*,' * (occ) single-particle occupations (from previous wfn)                   * '
        print*,' * (x) eXpectation value of a scalar Hamiltonian (from previous wfn)       * '
        print*,' * (o) Apply a one-body (transition) operator to previous wfn and write out* '
!        print*,' * (oo) Apply a one-body (transition) operator with enforced orthogonality * ' SECRET MENU ITEM
!        print*,' *     (that is, forced resulting pivot to be orthogonal to starting pivot * '
!        print*,' *      to eliminate transition strength to original state)                * '
!        print*,' * (ou) Apply uncoupled one-body (transition) operator to previous wfn     * ' ! SECRET MENU ITEM
		
        print*,' * (s),(sn) Strength function (using starting pivot ) wfn out normalized   * '
        print*,' * (ss) Strength function (using starting pivot ), but no output wfn or J,T* '
		
        print*,' * (su) Strength function (using starting pivot ) wfn out unnormalized     * '
        print*,' * (sb) Strength function (using block of starting pivots )                * '
        print*,' * (sbs) Strength function (using block of starting pivots ) no wfn out    * '

        print*,' * (a) Apply a scalar Hamiltonian to a previous wfn and write out          * '
        print*,' * (h) Compute matrix elements of a scalar Hamiltonian (inputs as basis)  * '
        print*,' * (g) Apply resolvent 1/(E-H) to a previous wfn and write out             * '  ! ADDED 7.8.8; restarted 7.9.10
        print*,' * (gv) Apply resolvent 1/(E-H) to a previous wfn, then take dot prod      * '  ! ADDED 7.9.10
        print*,' * (gc) Apply resolvent 1/(E-H) to a previous wfn, E complex, and write out* '  ! ADDED 7.10.6

        print*,' * (v) Overlap of initial states with final states                         * '
!        print*,' * (vs) Relative entropy between initial states and final states                         * '
		print*,' * (pv) Read in previous vector, write out Lanczos coef, take dot prod     * '
		
        print*,' * (m) print information for Modeling parallel distribution                * '
        print*,' * (md) Modeling parallel distribution for 1-body densities                * '
        print*,' * (m3) print information for Modeling parallel distribution w/3-body force* '
        print*,' * (m0) Compute dimensions only                                            * '
		
!        print*,' * (f) Self-consistent mean-field approximation (prepare pivot)            * ' NOW OBSOLETE AND DEPRECATED
        print*,' * (t) create TRDENS-readable file for post processing                     * '
        print*,' * (tx) create TRDENS-readable file for post processing from previous wfn  * ' ! ADDED 7.10.7
        print*,' * (tw) from TRDENS-readable file create standard BIGSTICK .wfn file       * ' ! ADDED 7.10.7

		print*,' * (b) Create binary file with full basis information (for postprocessing) * '
		print*,' * (ba) Create ASCII file with full basis information (for postprocessing) * '  ! ADDED 7.9.1
		
        print*,' * (wh) write out Hamiltonian matrix to a file and stop                    * '  ! ADDED 7.8.5
        print*,' * (wo) write out one-body transition matrix to a file and stop            * '  ! ADDED 7.8.5
!        print*,' * (wm) write out uncoupled 1+2-body operator m.e.s to file and stop       * '  ! ADDED 7.9.12

        print*,' * (co) Compute configurations (partitions)                                * '   ! ADDED 7.9.1
        print*,' * (cx) Compute configurations (partitions) from prior wfn                 * '   ! ADDED 7.10.4
		
		print*,' * (jp) Project states of good J from prior wfns and normalize             * '  ! ADDED 7.11.3
		print*,' * (ro) Read in multiple files of wfns and orthonormalize                  * '  ! ADDED 7.11.3
		print*,' * (ru) Read in multiple files of wfns but DO NOT orthonormalize           * '  ! ADDED 7.11.3

        print*,' * (c) Compute traces                                                      * '
        print*,' * (l) print license and copyright information                             * '
		
!       print*,' * (l) Time Hamiltonian loops  (obsolete?)                                 * '
        print*,' * (?) Print out all options                                               * '
    
        print*,' *                                                                         * '
        print*,' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ' 
     end if
     goto 1
     
  case ('i','I')
     if(iproc == 0)print*,' Reading automatically from "autoinput.bigstick" file '
     auto_input = .true.
     goto 1
  case ('d','D')
     menu_char = 'd'
     ham_readin = .true.
     if ( iproc == 0 )print*,' Compute spectrum + all one-body densities ' 
     densityflag = .true.
     write_wfn = .true.
  case ('dp','DP')
     menu_char = 'd'
     ham_readin = .true.
     if ( iproc == 0 )print*,' Compute spectrum + all one-body densities in proton-neutron format ' 
     densityflag = .true.
     write_wfn = .true.
	 pndensities = .true.
  !
  case ('db','DB')
     menu_char = 'd'
     ham_readin = .true.
     if ( iproc == 0 )then 
		 print*,' Compute spectrum + all one-body densities in proton-neutron format ' 
		 print*,' and write to a binary file '
	 end if
     densityflag = .true.
     write_wfn = .true.
	 pndensities = .true.	 
	 binden = .true.
  case ('dx','DX')
     menu_dx_omathematica = .false.
     menu_char = 'dx'
     auto_readin = .true.
     ham_readin = .false.
     if ( iproc == 0)then
		 print*,' Compute one-body densities from previous wfn ' 
		 print*,' Output written to .dres file in isospin format '
		 print*,' NEW: You will need to select initial, final range '
		
	 end if
     densityflag = .true.
	 onebodyonly =.true.
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
	 
  case ('dxp','DXP')
     menu_dx_omathematica = .false.
     menu_char = 'dx'
     auto_readin = .true.
    ham_readin = .false.
   	 pndensities = .true.
		
        if ( iproc == 0)then
   		 print*,' Compute one-body densities from previous wfn ' 
   		 print*,' Output written to .dres file in proton-neuton format '
   		 print*,' NEW: You will need to select initial, final range '
		
   	 end if
     densityflag = .true.
   	 onebodyonly =.true.
        call wfn_ropen_file(oldwfnfile)
        call read_wfn_header(oldwfnfile,.true.)	 
  !
  case ('dxb','DXB')
     menu_dx_omathematica = .false.
     menu_char = 'dx'
     auto_readin = .true.
    ham_readin = .false.
   	 pndensities = .true.
		
        if ( iproc == 0)then
   		 print*,' Compute one-body densities from previous wfn ' 
   		 print*,' Output written to binary .dres.bin file in proton-neuton format '
   		 print*,' NEW: You will need to select initial, final range '
		
   	 end if
     densityflag = .true.
   	 onebodyonly =.true.
	 binden = .true.
        call wfn_ropen_file(oldwfnfile)
        call read_wfn_header(oldwfnfile,.true.)	 	 
  case ('dxm', 'DXM')
        menu_dx_omathematica = .true.
        menu_char = 'dx'
        auto_readin = .true.
        ham_readin = .false.
        if ( iproc == 0 )print*,' Compute one-body densities from previous wfn ' 
        if ( iproc == 0 )print*,' and write out in Mathematica-readable format ' 

        densityflag = .true.
   	 onebodyonly =.true.
        call wfn_ropen_file(oldwfnfile)
        call read_wfn_header(oldwfnfile,.true.)	 

  
  case ('du','DU')
     menu_char = 'du'
     auto_readin = .true.
     ham_readin = .false.
     if ( iproc == 0)then
		 print*,' Compute uncoupled one-body densities from previous wfn ' 
		 print*,' Output written to .du file in isospin format '
		 print*,' You will need to select initial, final range '
		
	 end if
     densityflag = .true.
	 onebodyonly =.true.
	 uncoupled_dens =.true.
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
	 
  case ('2','2i','2d')                  ! ADDED IN 7.9.1  -- TWO BODY DENSITIES FROM EXISTING WFNs
     if(menu_char=='2i')then
	    pndensities=.false.
     else
	    pndensities = .true.
     end if
	 if(menu_char=='2d' .or. menu_char=='2i')then
		 diag_den2b = .true.
	 end if
	 
     menu_char = '2b'
     auto_readin = .true.
     ham_readin = .false.
     if ( iproc == 0 )then
		 print*,' Compute two-body densities from previous wfn ' 
!		 if(nproc> 1)print*,' THIS DOES NOT YET WORK IN  MPI '
		 
		 if(diag_den2b)then
			 print*,' (Only for initial=final, Jt=0)'
			 if(pndensities)then
				 print*,' (Computed in proton-neutron formalism )'
			 else
				 print*,' (Computed in isospin formalism )'
			 end if
		 end if
	 end if
     densityflag = .true.
	 dens2bflag = .true.
	 useTRphase = .false.  ! for now
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)	 

  case ('b','B','ba','BA')         ! ADDED IN 7.7.2
     if ( iproc == 0 )print*,' Create basis (.bas) file for later postprocessing ' 

     if(menu_char =='ba' .or. menu_char=='BA') then
		 if(nproc > 1)then
			 if(iproc==0)then
				 print*,' Sorry, this option will not work in MPI '
			 end if
#ifdef _MPI
		     call BMPI_Abort(MPI_COMM_WORLD,101,ierr)
#endif
			 stop
		 end if
		 baseASCII=.true.
		 if(iproc==0)print*,' (writing out in human-readable ASCII)'
	 end if
     menu_char = 'b'
	 writeout=.true.

!  case ('r','R')
!     menu_char = 'r'
!     auto_readin = .true.
!    ham_readin = .true.
!     call open_restart_files
!     call read_wfn_header(wfnfile,.true.)
!     if ( iproc == 0 ) then
!        print*,' Restart from previous run ' 
!     end if
     
  case ('a','A')
     menu_char = 'a'
     ham_readin = .true.
     auto_readin =.true.
     write_wfn = .true.
     if ( iproc == 0 ) then
        print*,' Apply a scalar Hamiltonian to a previous wfn and write out ' 
        print*,' (You must read in a previously computed wfn) '
     end if
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
  
  case ('h','H')
     menu_char = 'h'
     ham_readin = .true.
     auto_readin =.true.
     write_wfn = .false.
!     dotflag = .true.
     writeout = .false.
	 
     if ( iproc == 0 ) then
        print*,' Compute matrix elements of a scalar Hamiltonian '
		print*,' in a basis of previously compute states ' 
        print*,' (You must read in a previously computed wfn) '
     end if
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)	 
     
  case ('o','O')
     menu_char = 'o'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.true.
     write_wfn = .true.
	 onebodyonly = .true.
	 strengthnorm = .false.

     if ( iproc == 0 )print*,' Apply a one-body operator to a previous wfn and write out ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('ou','OU')  ! added  7.10.5 -- reading in and apply uncoupled one-body operator (don't normalize result)
     menu_char = 'o'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.true.
     write_wfn = .true.
	 onebodyonly = .true.
	 uncoupled_dens = .true.
	 strengthnorm = .false.
	 
     if ( iproc == 0 )print*,' Apply uncoupled one-body operator to a previous wfn and write out (unnormalized)' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
  !
  case ('oun','OUN')  ! added  7.10.5 -- reading in and apply uncoupled one-body operator  (normalize result)
     menu_char = 'o'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.true.
     write_wfn = .true.
	 onebodyonly = .true.
	 uncoupled_dens = .true.
	 strengthnorm = .true.
	 
     if ( iproc == 0 )print*,' Apply uncoupled one-body operator to a previous wfn and normalize ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
  
  
  case ('oo','OO')
     menu_char = 'oo'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.true.
     write_wfn = .true.
	 onebodyonly = .true.
	 strengthnorm = .false.

     if ( iproc == 0 )print*,' Apply a one-body operator to a previous wfn and write out ' 
     if ( iproc == 0 )print*,' (Will force resulting wfn to be orthogonal to original) ' 

     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('wo','WO')    ! write out one-body operator matrix to file trans1.dat '
        menu_char = 'wo'
        ham_readin = .false.
        auto_readin = .true.
        op_readin  =.true.
        write_wfn = .true.
   	 onebodyonly = .true.
     printouttrans1flag = .true.
	 

        if ( iproc == 0 )then
			print*,' Writing one-body transition operator to file trans1.dat '
			print*,' Kluge: must read in existing .wfn file to set up '
		end if
		 
        call wfn_ropen_file(oldwfnfile)
        call read_wfn_header(oldwfnfile,.true.)
	

  case ('v','V','vs','VS')
     dotflag = .true.
     if(menu_char=='vs' .or. menu_char=='VS')dotflag=.false.
     menu_char = 'v'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.false.
     write_wfn = .false.
     writeout = .false.

     if ( iproc == 0 )then
	print*,' Read in a wavefunction file, select one initial state ' 
        print*,' then read in second waverfunction file and compute overlap with all final states '
        print*,' '
        print*,' INITIAL STATE '

     endif
     
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('x','X')
     menu_char = 'x'
     auto_readin = .true.
     ham_readin = .true.
     if ( iproc == 0 ) print*,' Compute expectation value of Hamiltonian operator using previous wfn ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('p','P', 'ps', 'PS')
     ! s suppresses writing wfn - takes a long time
	 if(.not. (menu_char(2:2) == 's' .or. menu_char(2:2) == 'S')) write_wfn = .true.
     menu_char = 'n'
     if ( iproc == 0 )print*,' Compute spectrum + single-particle occupations ' 
     ham_readin = .true.
     spoccflag = .true.

  case ('occ', 'OCC')
     menu_char = 'p'
     auto_readin = .true.
     ham_readin = .false.
     spoccflag = .true.
     if ( iproc == 0 ) print*,' Compute single-particle occupations using previous wfn ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('m','M')
     menu_char = 'm'
     ham_readin = .false.
     if ( iproc == 0 )print*,' Print information needed to model parallel distribution ' 
     print4modelinfo = .true.

  case ('m0','M0')
        menu_char = 'm'
        ham_readin = .false.
        if ( iproc == 0 )print*,' Compute dimensions ONLY  '
		stopafterdim = .true.
		
  case ('md','MD')
     menu_char = 'm'
	 modeldensities=.true.
	 densityflag = .true.
     ham_readin = .false.
     if ( iproc == 0 )print*,' Print information needed to model parallel distribution for 1-body densities' 
     print4modelinfo = .true.   
	 
  case ('m3','M3')
        menu_char = 'm'
        ham_readin = .false.
        if ( iproc == 0 )print*,' Print information needed to model parallel distribution ' 
        print4modelinfo = .true.
	    threebodycheck=.true.
		
  case ('s','S','sn','SN','ss','SS')
     if(menu_char=='ss' .or. menu_char=='SS')then
       write_wfn=.false.
	   get_JT     = .false.
	   
     else
	   write_wfn=.true.
     end if 
     menu_char = 's'
     ham_readin = .true.
     auto_readin = .true.
     strengthflag = .true.
	 strengthnorm = .true.
     if ( iproc == 0 )print*,' Compute strength function distribution using previous wfn ' 
	 if(.not.write_wfn .and. iproc==0)print*,' (Will not compute J, T or write wfn to file )'

     if( iproc==0 .and. write_wfn)then ! added 7.9.7
		 write(6,*)' Output wfns will  be normalized '
		 write(6,*)' To change this, use option (su) instead'
	 end if
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

     case ('su','SU')
        menu_char = 's'
		strengthnorm = .false.
        ham_readin = .true.
        auto_readin = .true.
        strengthflag = .true.
        if ( iproc == 0 )print*,' Compute strength function distribution using previous wfn ' 
	 
        if(iproc==0)then ! added 7.9.7
   		 write(6,*)' Output wfns will not be normalized, but prior normalization x strength included '
   		 write(6,*)' To change this, use option (sn) instead '
   	    end if

        call wfn_ropen_file(oldwfnfile)
        call read_wfn_header(oldwfnfile,.true.)
        write_wfn = .true.	 
  
  case ('sb','SB','sbs','SBS')
     if(menu_char=='sbs' .or. menu_char=='SBS')then
	    write_wfn=.false.
		get_JT     = .false.
	 else
		write_wfn=.true.
	 end if 
     menu_char = 's'
     ham_readin = .true.
     auto_readin = .true.
     strengthflag = .true.
	 blockstrengthflag=.true.
     if ( iproc == 0 )print*,' Compute strength function distribution using block of previous wfns ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
	 if(.not.write_wfn .and. iproc==0)print*,' (Will not compute J, T or write wfn to file )'
  	 
  case ('g','G','gv','GV','Gv','gc','GC','Gc')
     if(menu_char=='gv' .or. menu_char=='GV'.or. menu_char=='Gv')then
		dotflag=.true.
	 end if
	 if(menu_char=='GC' .or. menu_char=='gc' .or. menu_char=='Gc')then
		 complex_green=.true.
	 end if
     menu_char = 'g'
     ham_readin = .true.
     auto_readin = .true.
     greenflag = .true.
	 strengthnorm=.false.
     if ( iproc == 0 )print*,' Apply resolvent/Green function 1/(E-H) using previous wfn ' 
	 if(iproc==0 .and. dotflag)print*,' (You will have a chance to take a dot product with other wave functions at the end) '
	 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
     write_wfn = .true.
	 
  case ('t','T')
     menu_char = 'n'
     ham_readin = .true.
     if ( iproc == 0 ) print*,' Will create output usable for TRDENS post-processing '
     trdensout = .true.
     write_wfn = .true.


  case ('tx','TX')           ! ADDED 7.10.7
     menu_char = 'tx'
     ham_readin = .false.
     if ( iproc == 0 ) print*,' Will create output usable for TRDENS post-processing from previous .wfn file'
	 auto_readin = .true.
     trdensout = .true.
     op_readin  =.false.
     write_wfn = .false.
     writeout = .false.
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
	 
  case ('tw','TW')           ! ADDED 7.10.7
     menu_char = 'tw'
     ham_readin = .false.
     if ( iproc == 0 ) then
		 print*,' Will read TRDENS-compatible file, create BIGSTICK-format .wfn file'
		 print*,' NOTE: You will have to enter in the wave function description by hand '
	 end if
	 auto_readin = .false.
     trdensout = .false.
     write_wfn = .true.
	 call open_trwfn

  case ('ns','NS')
     menu_char = 'ns'
     if ( iproc == 0 ) then 
        print*,' Compute spectrum only '
        print*,' Suppress eigenvector write up '
     end if
     ham_readin = .true.
     write_wfn = .false.
	 
  case ('ne','NE')   ! ADDED 7.6.8
        menu_char = 'ne'
        if ( iproc == 0 ) then 
           print*,' Compute spectrum only '
           print*,' Suppress eigenvector write up '
		   print*,' No observables, no J, no T'
        end if
        ham_readin = .true.
        write_wfn = .false.	 
		get_JT     = .false.
		
  case('np','NP')   ! ADDED 7.7.7  -- option to read with a specified pivot 
      menu_char = 'np'
      if ( iproc == 0 )print*,' Compute spectrum starting from a prior pivot ' 
      ham_readin = .true.
      write_wfn = .true. 
      auto_readin = .true.
      call wfn_ropen_file(oldwfnfile)
      call read_wfn_header(oldwfnfile,.true.)

  !
  case('pv','PV')   ! ADDED 7.11.2  
      menu_char = 'pv'
      if ( iproc == 0 )then
		  print*,' Compute Lanczos coeff and vec starting from a prior pivot; '
		  print*,' Fixed number of Lanczos iterations; '
		  print*,' dot Lanczos vectors with existing file of vectors; no spectrum ' 
	  end if
      ham_readin = .true.
      write_wfn = .false. 
      auto_readin = .true.
      call wfn_ropen_file(oldwfnfile)
      call read_wfn_header(oldwfnfile,.true.)

case ('wh','WH')   ! ADDED 7.8.5
      menu_char = 'wh'
      if ( iproc == 0 ) then 
         print*,' writing Hamiltonian to file ham.dat'

      end if
      ham_readin = .true.
      write_wfn = .false.	 
	get_JT     = .false.
    printouthamflag = .true.

case ('jp','JP')  ! ADDED 7.11.3
    menu_char = 'jp'
	if(iproc==0)then
		print*,' Read in prior wfns, project good J and normalize '
	end if
	
	auto_readin = .true.
    write_wfn = .true.	 
    call wfn_ropen_file(oldwfnfile)
    call read_wfn_header(oldwfnfile,.true.)
	
case ('ro','RO')  ! ADDED 7.11.3
    menu_char = 'ro'
	if(iproc==0)then
		print*,' Read multiple files of wfns and orthnormalize '
		print*,' '
		print*,' Enter first file of wfns, to set up basis '
		print*,' (Later you may enter in additional files)'
		print*,' '
	end if
	
	auto_readin = .true.	
    write_wfn = .true.	 
    call wfn_ropen_file(oldwfnfile)
    call read_wfn_header(oldwfnfile,.true.)	
!
case ('ru','RU')  ! ADDED 7.11.4
    menu_char = 'ru'
	if(iproc==0)then
		print*,' Read multiple files of wfns  '
		print*,' '
		print*,' Enter first file of wfns, to set up basis '
		print*,' (Later you may enter in additional files)'
		print*,' '
	end if
	
	auto_readin = .true.	
    write_wfn = .true.	 
    call wfn_ropen_file(oldwfnfile)
    call read_wfn_header(oldwfnfile,.true.)		
	
    case ('co','CO')
       menu_char = 'n'
       if ( iproc == 0 ) then 
          print*,' Compute configurations '
       end if
       ham_readin = .true.
       write_wfn = .true.
	   configout = .true.
	   
    case ('cx','CX')
          menu_char = 'cx'
          if ( iproc == 0 ) then 
             print*,' Compute configurations from previous file '
          end if
          auto_readin = .true.
          ham_readin = .false.

          write_wfn = .false.
   	   configout = .true.   
       call wfn_ropen_file(oldwfnfile)
       call read_wfn_header(oldwfnfile,.true.)
	   
  case ('c','C')
     menu_char = 'c'
     if ( iproc == 0 ) then 
        print*,' Compute traces (moments) '
		print*,' '
		print*,' Do you want variance (2nd moment?) (y/n)'
		if(auto_input)then
			read(autoinputfile,'(a)')ychar
		else
		    read(5,'(a)')ychar
			write(autoinputfile,'(a)')ychar
		end if
		if(ychar=='y' .or.ychar=='Y')then
			getvariance=.true.
			write(6,*)' Computing both centroid and variance '
		else
			getvariance=.false.
			write(6,*)' Computing only the centroid '
		end if
     end if
#ifdef _MPI
	 call BMPI_Bcast(getvariance,1,0,MPI_COMM_WORLD,ierr)
#endif
     ham_readin = .true.
     write_wfn = .false.
	 threebodycheck=.true.

  case ('l','L')

  call printlicense
  goto 1

!  case ('f','F')
!        menu_char = 'f'
!		meanie=.true.
 !       if ( iproc == 0 )print*,' Running self-consistent mean-field ' 
  !      ham_readin = .true.
 !       write_wfn = .true.
!        get_JT    = .false.
		
  case ('3')
	menu_char = 'n'
	if ( iproc == 0 )print*,' Compute spectrum only ' 
	if ( iproc == 0 )print*,' Allowing three-body ' 
    threebodycheck=.true.
	ham_readin = .true.
	write_wfn = .true.
  case ('j','J')
    menu_char='ne'		   		
    if(iproc==0)then
		
		if(nproc==1)then
			print*,' '
			print*,' Warning! Cannot carry out jumpstart on 1 MPI process '
			print*,' You MUST run on the same number of MPI processes as for the full run '
			print*,' This is to obtain timing information '
			print*,' '
			write(logfile,*)' '
			write(logfile,*)' Warning! Cannot carry out jumpstart on 1 MPI process '
			write(logfile,*)' You MUST run on the same number of MPI processes as for the full run '
			write(logfile,*)' This is to obtain timing information '
			write(logfile,*)' '
			close(logfile)
			stop
		end if
			
		print*,'  ............................ '
		print*,'  Carrying out a jumpstart run '
		print*,'  ............................ '
		
		write(logfile,*)' ............................ '
		write(logfile,*)' Carrying out a jumpstart run '
		write(logfile,*)' ............................ '
	end if
	ham_readin = .true.
	
	writejumpstart = .true.
	readjumpstart  = .false.
	get_JT = .false.  ! suppress eigenvectors, getting J^2, T^2
  
  case default
     menu_char = 'n'
     if ( iproc == 0 )print*,' Compute spectrum only ' 
     ham_readin = .true.
     write_wfn = .true.
  end select

  return
end subroutine menu
!=====================================================================

subroutine printlicense
	
	print*,' '

    print*,'  This code is distributed under MIT Open Source License. '
	print*,' Copyright (c) 2017 Lawrence Livermore National Security '
	
	print*,' and San Diego State University Research Foundation  '
	print*,' '
			
	print*,' Permission is hereby granted, free of charge, to any person obtaining '
	print*,' a copy of this software and associated documentation files (the "Software"),'
    print*,' to deal in the Software without restriction, including without limitation'
	print*,' the rights to use, copy, modify, merge, publish, distribute, sublicense, '
	print*,' and/or sell copies of the Software, and to permit persons to whom the Software '
	print*,' is furnished to do so, subject to the following conditions: '
	print*,' '
	print*,' The above copyright notice and this permission notice shall be included '
	print*,' in all copies or substantial portions of the Software. '
	print*,' '
	print*,' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,'
	print*,' INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A '
	print*,' PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT '
	print*,' HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION '
	print*,' OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE '
	print*,' OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'

	print*,' '
	
	return
	
	
end subroutine printlicense
