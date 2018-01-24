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
  implicit none

  integer(4) :: ierr

! File to log interactive inputs to so the user can create an input file
  open(unit=autoinputfile,file='autoinput.bigstick',status='unknown')  ! always do this

  auto_input = .false.
  meanie=.false.
  pndensities = .false.
  if ( iproc == 0 ) then
     print*,' '
     print*,' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ' 
     print*,' *                                                                         * '
     print*,' *               OPTIONS (choose one)                                      * '
     print*,' * (i) Input automatically read from "autoinput.bigstick" file             * '
     print*,' *  (note: autoinput.bigstick file created with each nonauto run)          * '
     print*,' * (n) Compute spectrum (default); (ns) to suppress eigenvector write up   * '
     print*,' * (d) Densities: Compute spectrum + all one-body densities                * '
     print*,' * (dx[m]) Densities: Compute one-body densities from previous run (.wfn)  * '
     print*,' *     optional m enables mathematica output                               * '
     print*,' * (p) Compute spectrum + single-particle occupations; (ps) to suppress wfn* '
     print*,' * (occ) single-particle occupations (from previous wfn)                   * '
     print*,' * (x) eXpectation value of a scalar Hamiltonian (from previous wfn)       * '
     print*,' * (o) Apply a one-body (transition) operator to previous wfn and write out* '
     print*,' * (s) Strength function (using starting pivot )                           * '
!     print*,' * (r) Restart from a previous run                                         * '
     print*,' * (a) Apply a scalar Hamiltonian to a previous wfn and write out          * '
!     print*,' * (y) Cycle back to this menu after calculation completed (not active)    * '
     print*,' * (v) Overlap of initial states with final states                         * '
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

  call BMPI_BARRIER(icomm,ierr)
  call BMPI_BCAST(menu_char,3,0,icomm,ierr)

  auto_readin = .false.
  ham_readin = .false.
  op_readin  = .false.
  print4modelinfo = .false.
  strengthflag = .false.
  densityflag = .false.
  trdensout   = .false.
  spoccflag   = .false.
  menu_dx_omathematica = .true.  ! flag enabling mathematica output of density matrices
  get_JT       = .true.
  onebodyonly  = .false.
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
        print*,' * (np) Compute spectrum starting from prior pivot                         * '   !ADDED 7.7.7
        print*,' * (d) Densities: Compute spectrum + all one-body densities                * '
        print*,' * (dx[m]) Densities: Compute one-body densities from previous run (.wfn)  * '
        print*,' *     optional m enables mathematica output                               * '
        print*,' * (dp) Densities in proton-neutron format                                 * '
		
        print*,' * (p) Compute spectrum + single-particle occupations,(ps) to supress wfn  * '
        print*,' * (occ) single-particle occupations (from previous wfn)                   * '
        print*,' * (x) eXpectation value of a scalar Hamiltonian (from previous wfn)       * '
        print*,' * (o) Apply a one-body (transition) operator to previous wfn and write out* '
        print*,' * (oo) Apply a one-body (transition) operator with enforced orthogonality * '
        print*,' *     (that is, forced resulting pivot to be orthogonal to starting pivot * '
        print*,' *      to eliminate transition strength to original state)                * '
        print*,' * (s) Strength function (using starting pivot )                           * '
!        print*,' * (r) Restart from a previous run                                         * '  ! NOW OBSOLETE
        print*,' * (a) Apply a scalar Hamiltonian to a previous wfn and write out          * '
!       print*,' * (y) Cycle back to this menu after calculation completed (not active)    * '
        print*,' * (v) Overlap of initial states with final states                         * '
        print*,' * (m) print information for Modeling parallel distribution                * '
        print*,' * (f) Self-consistent mean-field approximation (prepare pivot)            * '
        print*,' * (t) create TRDENS-readable file for post processing                     * '
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
  case ('dx','DX','dxm', 'DXM')
     menu_dx_omathematica = .true.
     menu_char = 'dx'
     auto_readin = .true.
     ham_readin = .false.
     if ( iproc == 0 )print*,' Compute one-body densities from previous wfn ' 
     densityflag = .true.
	 onebodyonly =.true.
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('b','B')         ! ADDED IN 7.7.2
     menu_char = 'b'
     if ( iproc == 0 )print*,' Create basis (.bas) file for later postprocessing ' 
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
     
  case ('o','O')
     menu_char = 'o'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.true.
     write_wfn = .true.
	 onebodyonly = .true.

     if ( iproc == 0 )print*,' Apply a one-body operator to a previous wfn and write out ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('oo','OO')
     menu_char = 'oo'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.true.
     write_wfn = .true.
	 onebodyonly = .true.

     if ( iproc == 0 )print*,' Apply a one-body operator to a previous wfn and write out ' 
     if ( iproc == 0 )print*,' (Will force resulting wfn to be orthogonal to original) ' 

     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('v','V')
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

  case ('s','S')
     menu_char = 's'
     ham_readin = .true.
     auto_readin = .true.
     strengthflag = .true.
     if ( iproc == 0 )print*,' Compute strength function distribution using previous wfn ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
     write_wfn = .true.

  case ('t','T')
     menu_char = 'n'
     ham_readin = .true.
     if ( iproc == 0 ) print*,' Will create output usable for TRDENS post-processing '
     trdensout = .true.
     write_wfn = .true.

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
	  
  case ('c','C')
     menu_char = 'c'
     if ( iproc == 0 ) then 
        print*,' Compute traces (1st, 2nd moments) '
     end if
     ham_readin = .true.
     write_wfn = .false.

  case ('l','L')

  call printlicense
  goto 1

  case ('f','F')
        menu_char = 'f'
		meanie=.true.
        if ( iproc == 0 )print*,' Running self-consistent mean-field ' 
        ham_readin = .true.
        write_wfn = .true.
        get_JT    = .false.
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
