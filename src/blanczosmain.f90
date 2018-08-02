!====================================================================
!  LANCZOS routines for BIGSTICK
!
!  versions for 'new' parallelization scheme -- FALL 2011
!  This code uses LAPACK ROUTINES
!
!  LAPACK copyright statements and license
!
!Copyright (c) 1992-2013 The University of Tennessee and The University
!                        of Tennessee Research Foundation.  All rights
!                        reserved.
!Copyright (c) 2000-2013 The University of California Berkeley. All
!                        rights reserved.
!Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
!                        reserved.

!Additional copyrights may follow

!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are
!met:
!
!- Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!
!- Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer listed
!  in this license in the documentation and/or other materials
!  provided with the distribution.
!
!- Neither the name of the copyright holders nor the names of its
!  contributors may be used to endorse or promote products derived from
!  this software without specific prior written permission.
!
!The copyright holders provide no reassurances that the source code
!provided does not infringe any patent, copyright, or any other
!intellectual property rights of third parties.  The copyright holders
!disclaim any liability to any recipient for claims brought against
!recipient by any third party for infringement of that parties
!intellectual property rights.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===== END OF LACK COPYRIGHT STATEMENT ========================
!===========================================================================
!  SUBROUTINES IN THIS FILE
!   lanczos_menu
!	setup_localvectors       : when lanczos vector broken up
!	initialize_lanczos_vector
!       swap_vchar               : controls "direction" of mat-vec multiply
!                                   either H * vec1 = vec2  or  H * vec2 = vec1
!       lanczos_p                : main lanczos subroutine
!       initialize_final         : zeroes out next vector to be calculated
!       lanczos_output
!       find_lanczos_eigenvalues
!       density1b_output
!       exactdiag_p
!       thick_restart_sub_p
!       random_restart_p
!       open_lanczosfile
!       close_lanczosfile

!===========================================================================
!
!  subroutine lanczos_menu
!
!  more options for user-defined control of lanczos iterations and convergence
!
!  asks for # of iterations and # of vectors to keep
!
!  SUBROUTINES CALLED:
!   countlanczositerations  : if restarting
!	openlanczosfile
!	lanczos_p
!	exactdiag
!       distribute_lanczos_pieces in bparallelib4.f90
!
!   OPTIONS:
!      EX  : exact/full diagonalization
!      LD  : default lanczos with standard convergence
!      LC  : normal lanczos, but with user-defined convergence
!      LF  : normal lanczos with fixed iterations
!      TD  : thick restart with default convergence
!      TC  : thick restart with user-defined convergence
!      TF  : thick restart with fixed iterations
!      TU  : User-defined thick-restart parameters
!
!===========================================================================
subroutine lanczos_menu(restart)
  use flagger
  use io
  use system_parameters
  use basis
!  use lanczos_info
  use nodeinfo
  use localvectors
  use precisions
  use fragments
  use coupledmatrixelements
  use bmpi_mod
  use butil_mod
  use para_main_mod
  use para_util_mod
  use lanczos_util
  use apply_ham_omp
  use onebodypot, only: meanie
  implicit none
  logical    :: restart
  integer(4) :: ierr
  integer(4) :: maxiter,max_iter_restart
  integer(4) :: niter0
!  character(2) :: lanczchar
  logical :: leave  ! flag to abort
  integer :: aerr

  leave = .false.
  thick_restart = .false.

  if ( restart ) then
     call countlanczositerations
  else
     if(reorthog)call open_lanczosfile
     startiter = 0
  end if
  
!.............. IF STRENGTH FUNCTION MENU OPTION 'S' THEN ONLY FIXED LANCZOS .............

if(strengthflag)then
	lanczchar='lf'
	if(iproc==0)then
		write(6,*)' Fixed iterations ONLY: '
        write(6,*)' Enter nkeep, # iterations for lanczos  '
        write(6,*)' (nkeep = # of states printed out )'
        read(5,*)nkeep,niter
        if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
        if(niter < nkeep) niter = nkeep
        if(niter >= dimbasis) niter = int(dimbasis-1, 4)
        write(autoinputfile,*)nkeep,niter,'   ! # states to keep, # iterations '		
		
	end if
	leave = .false.
	goto 543
	
end if
!........................................................................................  

  if ( iproc == 0 ) then
11   continue

!..... THE FOLLOWING USED ONLY IN OPTION 'tx' TARGETING A SPECIFIC ENERGY; may be revised
     Etarget = 0.
     targetX = .false.
!....................................	 
     print*,' '
     print*,' / ------------------------------------------------------------------------\ ' 
     print*,' |                                                                         | '
     print*,' |    DIAGONALIZATION OPTIONS (choose one)          | '
     print*,' |                                                                         | '

     if(.not.strengthflag)then     ! cannot use strength option with exact diagonalization
     print*,' | (ex) Exact and full diagonalization (use for small dimensions only)     | '

     if(dimbasis <= fulllimit)then
     print*,'       (recommended, as dim = ',dimbasis,')'

     end if
     
     end if

     print*,' |                                                                         | '
     print*,' | (ld) Lanczos with default convergence (STANDARD)                        | '
     print*,' | (lf) Lanczos with fixed (user-chosen) iterations                        | '
     print*,' | (lc) Lanczos with user-defined convergence                              | '
     print*,' |                                                                         | '

     if(.not.strengthflag)then     ! cannot use strength option with thick-restart

     print*,' | (td) Thick-restart Lanczos with default convergence                     | '
     print*,' | (tf) Thick-restart Lanczos with fixed iterations                        | '
     print*,' | (tc) Thick-restart Lanczos with user-defined convergence                | '
!------- ADDED in 7.6.8 -------	 
     print*,' | (tx) Thick-restart Lanczos targeting states near specified energy       | '
	 

     end if
     print*,' |                                                                         | '
     print*,' | (sk) Skip Lanczos (only used for timing set up)                         | '
     print*,' |                                                                         | '
     print*,' \ ------------------------------------------------------------------------/ ' 
     print*,' '

     if(auto_input)then                
 ! ...............................AUTO INPUT..............................
        read(autoinputfile,'(a)')lanczchar
        select case (lanczchar)
           case('ex')
              read(autoinputfile,*)nkeep
              if(nkeep>dimbasis)nkeep=int(dimbasis,4)
              print*,' Full diagonalization; keeping ',nkeep,' out of ',dimbasis
           case('ld')
              read(autoinputfile,*)nkeep,maxiter
	          if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
	          if(maxiter < nkeep) maxiter = nkeep
	          if(maxiter >= dimbasis) maxiter = int(dimbasis-1, 4)
			  
              print*,' Default Lanczos convergence; keeping ',nkeep,' out of ',dimbasis
              print*,' Max ',maxiter,' iterations '
              niter = maxiter
              ncheck_ex = ncheck_ex_def
              converge_test = converge_test_def
              converge_tol  = converge_tol0_def

           case('lf')
              read(autoinputfile,*)nkeep,niter
	          if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
	          if(niter < nkeep) niter = nkeep
	          if(niter >= dimbasis) niter = int(dimbasis-1, 4)
              print*,' Fixed Lanczos iterations; keeping ',nkeep,' out of ',dimbasis
              print*,' with ',niter,' total Lanczos iterations '

           case('lc')
              print*,' User-defined Lanczos convergence '
              read(autoinputfile,*)nkeep,maxiter
	          if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
	          if(maxiter < nkeep) maxiter = nkeep
	          if(maxiter >= dimbasis) maxiter = int(dimbasis-1, 4)
              niter = maxiter
              print*,' Keeping ',nkeep,' eigenstates with max of ',maxiter,' iterations '
              read(autoinputfile,*)ncheck_ex
              print*,'( Using a total of ',nkeep+ncheck_ex,' states for convergence test ) '
              read(autoinputfile,*)converge_test
              read(autoinputfile,*)converge_tol
              select case (converge_test)
                 case (0)
              print*,'(0) Average difference in energies between one iteration and the last '
                 case (1)
              print*,'(1) Max difference in energies between one iteration and the last '
                 case (2)
              print*,'(2) Average difference in wavefunctions between one iteration and the last; '
                 case (3)
              print*,'(3) Min difference in wavefunctions between one iteration and the last; '
              end select 
              print*,' Convergence tolerance = ',converge_tol

           case('td')
              print*,' Thick-restart Lanczos with default convergence test '
              thick_restart = .true.
              fixed_iter = .false.
              read(autoinputfile,*)nkeep,niter
              max_iter_restart = bmax((nkeep+ncheck_ex_def)*thick_restart_factor_def,thick_restart_min_def)
              if(niter < max_iter_restart)then
                niter = max_iter_restart
              end if
              print*,' Keeping ',nkeep,' eigenstates with ',niter,' iterations '
              nthick_keep = max(nkeep+nthick_add_def,3*nkeep)
			  if(nthick_keep > niter - nthick_add_def)nthick_keep = niter-nthick_add_def
              read(autoinputfile,*)max_num_thick_restart
              print*,' with a maximum of ',max_num_thick_restart,' restarts '
              converge_test = 0
              converge_tol = converge_tol0_def
 
           case('tf')
              print*,' Thick-restart Lanczos with fixed iterations'
              read(autoinputfile,*)nkeep,niter
              read(autoinputfile,*)nthick_keep
              read(autoinputfile,*)max_num_thick_restart
              print*,' Printing out ',nkeep,' final eigenstates '
              print*,' Doing ',niter,' iterations before restarting '
              print*,' Keeping ',nthick_keep,' vectors after each restart '
              print*,' Doing ',max_num_thick_restart,' restarts '
              thick_restart = .true.
              fixed_iter = .true.

           case('tc')
              print*,' Thick-restart Lanczos with user-defined convergence test '

              thick_restart = .true.
              fixed_iter = .false.

              read(autoinputfile,*)nkeep,niter
              read(autoinputfile,*)nthick_keep
              read(autoinputfile,*)max_num_thick_restart
              print*,' Printing out ',nkeep,' final eigenstates '
              print*,' Doing ',niter,' iterations before restarting '
              print*,' Keeping ',nthick_keep,' vectors after each restart '
              print*,' Doing a max of ',max_num_thick_restart,' restarts '
              read(autoinputfile,*)converge_test
              read(autoinputfile,*)converge_tol
              select case (converge_test)
                 case (0)
              print*,'(0) Average difference in energies between one iteration and the last '
                 case (1)
              print*,'(1) Max difference in energies between one iteration and the last '
                 case (2)
              print*,'(2) Average difference in wavefunctions between one iteration and the last; '
                 case (3)
              print*,'(3) Min difference in wavefunctions between one iteration and the last; '
              end select 
              print*,' Convergence tolerance = ',converge_tol
			  
		   case('tx')
           print*,' Thick-restart Lanczos targeting specific energy '
           thick_restart = .true.
           fixed_iter = .false.

           read(autoinputfile,*)nkeep,niter
		   print*,' Keeping ',nkeep,' vectors,', niter, 'iterations before restarting '
           max_iter_restart = bmax((nkeep+ncheck_ex_def)*thick_restart_factor_def,thick_restart_min_def)
           if(niter < max_iter_restart)then
             niter = max_iter_restart
           end if

           read(autoinputfile,*)nthick_keep
		   print*,nthick_keep,' vectors kept after thick restart'
		  
           if(nthick_keep < nkeep)nthick_keep = nkeep

           read(autoinputfile,*)max_num_thick_restart
           print*,' Keeping ',nkeep,' eigenstates with ',niter,' iterations '
           print*,' with a maximum of ',max_num_thick_restart,' restarts '
           converge_test = 0
           converge_tol = converge_tol0_def		
			  targetX = .true.	  
           read(autoinputfile,*)Etarget
		   print*,' Targeting energy = ',Etarget

           case('sk')
              print*,' Skipping diagonalization '

        end select
     else
!................................ READ IN DIAGONALIZATION OPTIONS AND PARAMETERS......

        read(5,'(a)')lanczchar
  
        select case (lanczchar)
 
           case ('ex','EX','Ex')
              if(strengthflag)then
                 print*,' Cannot use strength function option "S" with exact diagonalization '
                 stop
              end if
              lanczchar = 'ex'
              write(autoinputfile,'(a,"    ! Lanczos menu option ")')lanczchar(1:2)
              write(6,*)' Doing full diagonalization; Enter nkeep (out of ',dimbasis,')'
              read(5,*)nkeep                 
              if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
              write(autoinputfile,*)nkeep,'    ! # of states to keep '

           case ('ld','LD','Ld')
              lanczchar = 'ld'
              write(autoinputfile,'(a,"    ! Lanczos menu option ")')lanczchar(1:2)
              write(6,*)' Enter nkeep, max # iterations for lanczos '
              write(6,*)' (nkeep = # of states printed out )'
              read(5,*)nkeep,maxiter
              if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
              if(maxiter < nkeep) maxiter = nkeep
              if(maxiter >= dimbasis) maxiter = int(dimbasis-1, 4)
              niter = maxiter
              write(autoinputfile,*)nkeep,maxiter,'     ! # states to keep, max # iterations '
              ncheck_ex = ncheck_ex_def
              converge_test = converge_test_def
              converge_tol  = converge_tol0_def

           case ('lf','LF','Lf')
              lanczchar = 'lf'
              write(autoinputfile,'(a,"    ! Lanczos menu option ")')lanczchar(1:2)

              write(6,*)' Enter nkeep, # iterations for lanczos  '
              write(6,*)' (nkeep = # of states printed out )'
              read(5,*)nkeep,niter
              if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
              if(niter < nkeep) niter = nkeep
              if(niter >= dimbasis) niter = int(dimbasis-1, 4)
              write(autoinputfile,*)nkeep,niter,'   ! # states to keep, # iterations '

           case ('lc','LC','Lc')
              lanczchar = 'lc'
              write(autoinputfile,'(a,"    ! Lanczos menu option ")')lanczchar(1:2)
              print*,' '
              write(6,*)' Enter nkeep, max # iterations for lanczos  '
              read*,nkeep,maxiter
              write(autoinputfile,*)nkeep,maxiter
              niter = maxiter
              print*,' Enter how many ADDITIONAL states for convergence test '
              print*,' ( Defaul t= ',ncheck_ex_def,'; you may choose 0 ) '
              read*,ncheck_ex
              write(autoinputfile,*)ncheck_ex

321           continue 
              print*,' '
              print*,'    Enter one of the following choices for convergence control :'
              print*,'(0) Average difference in energies between one iteration and the last; '
              print*,'(1) Max difference in energies between one iteration and the last; '
              print*,'(2) Average difference in wavefunctions between one iteration and the last; '
              print*,'(3) Min difference in wavefunctions between one iteration and the last; '
              read*,converge_test
              if(converge_test < 0 .or. converge_test > 3)goto 321

              print*,' Enter desired tolerance '
              select case (converge_test) 
                 case (0)
                 write(6,322)converge_tol0_def

                 case (1)
                 write(6,322)converge_tol1_def

                 case (2)
                 write(6,322)converge_tol2_def

                 case (3)
                 write(6,322)converge_tol3_def
              end select
322  format(" (default tol = ",E10.3," ) ")

              read*,converge_tol
              write(autoinputfile,*)converge_test
              write(autoinputfile,*)converge_tol
              print*,' '

!.......... NOTES ON THICK RESTART......................
!   WE TARGET nkeep EIGENVECTORS
!   niter IS MAX # OF LANCZOS ITERATIONS BEFORE RESTARTING
!   nthick_keep IS HOW MANY "LANCZOS" VECTORS WE RESTART WITH
!   max_num_thick_restart = max # of times we restart

           case ('td','TD','Td')
              lanczchar = 'td'
              write(autoinputfile,'(a)')lanczchar(1:2)

!................ SET DEFAULT THICK-RESTART PARAMETERS................

              if(strengthflag)then
                 print*,' Cannot use strength function option "S" with thick restart Lanczos '
                 leave=.true.
                 goto 543
                 
              end if
              thick_restart = .true.
              fixed_iter    = .false.
              print*,' Enter # of states to keep, number of iterations before restarting '
              read*,nkeep,niter
              write(autoinputfile,*)nkeep,niter
              max_iter_restart = bmax((nkeep+ncheck_ex_def)*thick_restart_factor_def,thick_restart_min_def)
              if(niter < max_iter_restart)then
                niter = max_iter_restart
              end if
              nthick_keep = max(nkeep+nthick_add_def,3*nkeep)
			  if(nthick_keep > niter - nthick_add_def)nthick_keep = niter-nthick_add_def
              write(6,*)' Enter max # of restarts '
              read(5,*)max_num_thick_restart
              write(autoinputfile,*)max_num_thick_restart
              converge_test = 0
              converge_tol = converge_tol0_def
			  
           case('tx','TX','Tx')
		      lanczchar='tx'
              write(autoinputfile,'(a)')lanczchar(1:2)
              print*,' Thick-restart Lanczos targeting specific energy '
              thick_restart = .true.
              fixed_iter = .false.
              print*,' Enter # of states to keep, number of iterations before restarting '
              read*,nkeep,niter
              write(autoinputfile,*)nkeep,niter
              max_iter_restart = bmax((nkeep+ncheck_ex_def)*thick_restart_factor_def,thick_restart_min_def)
              if(niter < max_iter_restart)then
                niter = max_iter_restart
              end if
              print*,' Enter # of vectors to keep after thick-restart '
              print*,' (Typical value would be ',  & 
                 bmax((nkeep+ncheck_ex_def)*thick_restart_factor_def,thick_restart_min_def), ')'
              print*,' ( Must be between ',nkeep,' and ',niter,')'
              read(5,*)nthick_keep
              write(autoinputfile,*)nthick_keep
			  
              if(nthick_keep < nkeep)nthick_keep = nkeep
              write(6,*)' Enter max # of restarts '
              read(5,*)max_num_thick_restart
              write(autoinputfile,*)max_num_thick_restart
              print*,' Keeping ',nkeep,' eigenstates with ',niter,' iterations '
              print*,' with a maximum of ',max_num_thick_restart,' restarts '
              converge_test = 0
              converge_tol = converge_tol0_def		
   			  targetX = .true.
   			  print*,' Enter your target energy '
   			  read*,Etarget	  
              write(autoinputfile,*)Etarget
			  

           case ('tf','TF','Tf')
              lanczchar = 'tf'
              write(autoinputfile,'(a)')lanczchar(1:2)

              if(strengthflag)then
                 print*,' Cannot use strength function option "S" with thick restart Lanczos '
                 leave = .true.
                 goto 543

              end if
              print*,' Enter # of states to keep, # of iterations before thick-restart '
              read(5,*)nkeep,niter
              print*,' Enter # of vectors to keep after thick-restart '
              print*,' (Typical value would be ',  & 
                 bmax((nkeep+ncheck_ex_def)*thick_restart_factor_def,thick_restart_min_def), ')'
              print*,' ( Must be between ',nkeep,' and ',niter,')'
              read(5,*)nthick_keep
              if(nthick_keep < nkeep)nthick_keep = nkeep
              thick_restart = .true.
              fixed_iter    = .true.
              write(6,*)' Enter # of restarts '
              read(5,*)max_num_thick_restart
              write(autoinputfile,*)nkeep,niter
              write(autoinputfile,*)nthick_keep
              write(autoinputfile,*)max_num_thick_restart

           case ('tc','TC','Tc')
              lanczchar = 'tc'
              write(autoinputfile,'(a)')lanczchar(1:2)

!................ USER SETS THICK-RESTART PARAMETERS................

              if(strengthflag)then
                 print*,' Cannot use strength function option "S" with thick restart Lanczos '
                 leave = .true.
                 goto 543

              end if
              thick_restart = .true.
              fixed_iter    = .false.

              print*,' Enter # of states to keep, number of iterations before restarting '
              read*,nkeep,niter
              write(autoinputfile,*)nkeep,niter
              print*,' Enter # of vectors to keep after thick-restart '
              print*,' (Typical value would be ',  & 
                 bmax((nkeep+ncheck_ex_def)*thick_restart_factor_def,thick_restart_min_def), ')'
              print*,' ( Must be > ',nkeep,')'
              read(5,*)nthick_keep
              write(autoinputfile,*)nthick_keep
              write(6,*)' Enter max # of restarts '
              read(5,*)max_num_thick_restart
              write(autoinputfile,*)max_num_thick_restart

421           continue 
              print*,' '
              print*,'    Enter one of the following choices for convergence control :'
              print*,'(0) Average difference in energies between one iteration and the last; '
              print*,'(1) Max difference in energies between one iteration and the last; '
              print*,'(2) Average difference in wavefunctions between one iteration and the last; '
              print*,'(3) Min difference in wavefunctions between one iteration and the last; '
              read*,converge_test
              if(converge_test < 0 .or. converge_test > 3)goto 421

              print*,' Enter desired tolerance '
              select case (converge_test) 
                 case (0)
                 write(6,322)converge_tol0_def

                 case (1)
                 write(6,322)converge_tol1_def

                 case (2)
                 write(6,322)converge_tol2_def

                 case (3)
                 write(6,322)converge_tol3_def
              end select

              read*,converge_tol
              write(autoinputfile,*)converge_test
              write(autoinputfile,*)converge_tol
              print*,' '

           case ('sk','SK','q','Q')
              print*,'skipping diagonalization '
              lanczchar = 'sk'
              write(autoinputfile,'(a)')lanczchar(1:2)

           case default
              print*,' Not one of the options '
              goto 11

        end select

     end if
  end if

  call BMPI_BCAST(leave,1,0,icomm,ierr)
543 continue
  if(leave)then
     call BMPI_ABORT(icomm,120,ierr)
     stop
  end if
  call BMPI_BCAST(lanczchar,2,0,icomm,ierr)
  call BMPI_BCAST(niter,1,0,icomm,ierr)
  call BMPI_BCAST(nkeep,1,0,icomm,ierr)
  call BMPI_BCAST(maxiter,1,0,icomm,ierr)
  call BMPI_BCAST(thick_restart,1,0,icomm,ierr)
  call BMPI_BCAST(nthick_keep,1,0,icomm,ierr)
  call BMPI_BCAST(max_num_thick_restart,1,0,icomm,ierr)


  call BMPI_BCAST(converge_test,1,0,icomm,ierr)
  call BMPI_BCAST(ncheck_ex,1,0,icomm,ierr)
  call BMPI_BCAST(converge_tol,1,0,icomm,ierr)
  call BMPI_BCAST(Etarget,1,0,icomm,ierr)
  call BMPI_BCAST(targetX,1,0,icomm,ierr)
  nthick_add = nthick_keep - nkeep     ! how may to add on

  if((nkeep==0 .and. niter==0) .or. lanczchar=='sk')then  ! DO NOT TRY TO DIAGONALIZE
       close(autoinputfile)
       return       
  end if

  if(lanczchar=='lf' .or. lanczchar=='tf')then
     fixed_iter = .true.
  else
     fixed_iter = .false.
  end if

  if(lanczchar /= 'sk' .and. emptyHam)then
     if(iproc==0)then
         print*,' '
         print*,' ** WARNING WARNING NO INTERACTION READ IN WARNING WARNING ** '
		 print*,' *** YOU ARE LIKELY TO GET AN ERROR  ********  '
         print*,' '
         write(logfile,*)' ** WARNING WARNING NO INTERACTION READ IN WARNING WARNING ** '
		 write(logfile,*)' *** YOU ARE LIKELY TO GET AN ERROR  ********  '
     end if
  end if

!............. OPTIONAL STORAGE OF VECTORS IN RAM ..............

  if(enablelanczosincore1 .and. nproc == 1 .and. & 
     (useNewReorthog .or. dimbasis*(niter)*lanc_prec/2*1.0e-9 < maxlanczosstorage1)) then
          storelanczosincore1 = .true.
          if(.not. useNewReorthog)then
             allocate( lvec(dimbasis, niter), stat=aerr )
             if(aerr /= 0) call memerror("lanczos_menu 1");
          end if
          print*,' Internal storage of all lanczos vectors = ', & 
               dimbasis*(niter)*lanc_prec/2*1.0e-9,' Gb '
          print*,dimbasis,niter,lanc_prec
  else
          storelanczosincore1 = .false.
          if(nproc==1)then
             print*,' Internal storage of all lanczos vectors would require ', & 
               dimbasis*(niter)*lanc_prec/2*1.0e-9,' Gb '
             print*,' You can change variable maxlanczosstorage1 = ',maxlanczosstorage1
	     print*,' in module flagger in file bmodules_flags.f90 '

          end if
  end if

!................................................
  if(enablelanczosincoreMPI .and. nproc > 1)then
          storelanczosincoreMPI = .true.
  else
          storelanczosincoreMPI = .false.
  end if
  call distribute_lanczos_pieces


  if( lanczchar=='ex')then
 	  if(nproc > 1)then
 		  call exactdiag_MPI
 	  else
        call exactdiag_p
     end if 
  else
	  call setup_for_lanczos   ! ADDED IN 7.6.3 
      if (useHZSomp .and. .not. meanie) then
         if (iproc == 0 .and. verbose_ompham) print *, "MPI Time (before getLoadDist) : ", BMPI_Wtime() - mpiStartTime
         call getLoadDist()
      end if
      if (iproc == 0 .and. verbose_ompham) print *, "MPI Time (before lanczos_p) : ", BMPI_Wtime() - mpiStartTime
      call lanczos_p
      if (iproc == 0 .and. verbose_ompham) print *, "MPI Time (after lanczos_p) : ", BMPI_Wtime() - mpiStartTime
	  
  end if
  return
end subroutine lanczos_menu


!=================================================================
! SUBROUTINE LANCZOS_P
!
! Main Lanczos subroutine for BIGSTICK
!
!  -- modified for "new" parallel regime (Sept 2011 by CWJ)
!
! Can read in a pivot and generate a strength function;
! Can compute one-body densities matrices from eigenvectors
!
! Thick-restart added May 2009 by CWJ
!
!
! SUBROUTINES CALLED:
!   initialize_lanczos_vector
!   setup_localvectors
!   br_grab_vec1
!   br_normalize
!   br_restore_vec1
!   br_add2hist
!   print_vec_head
!   write_lanczos_vector_a
!   initialize_final
!   applyHbundled_g
!   swap_vchar
!   br_grab_vec2
!   br_remove_prev_overlap
!   br_orthogonalize
!   reorthogonalize_a
!   dnormvec_p
!   random_restart_p
!   find_lanczos_eigenvalues
!   thick_restart_sub_p
!   lanczos_output
!
!=================================================================
subroutine lanczos_p

  use system_parameters
  use sporbit
  use io
  use menu_choices
  use basis
  use obs
  use densities
  use haiku_info
  use verbosity
  use nodeinfo
  use localvectors
  use flagger
  use fragments
  use mod_reorthog
  use bmpi_mod
  use butil_mod
  use bvectorlib_mod
  use lanczos_util
  use apply_ham
  use apply_ham_omp

  implicit none
  integer(4) :: ierr

  integer(4) :: iter
  real(kind=8) :: da,db
  integer(4) :: i,j,k
  real(kind=8) :: dnorm, dnorm0
  real(4) :: xj,xt

!--------------------- CONVERGENCE -----------------------------------------

  real    :: ediff0,ediff1

!-------------- ALTERNATE CONVERGENCE ---------------------------

  real(8)     :: vdiff2,vdiff3,vdot

  logical smallflag,zeroflag
  
  integer :: i_thick_restart
  integer(4) :: i_thick

  integer :: nread,jvec

  integer inode
  character(1) :: vchar
  integer :: aerr

  if(noisy0) print *, "Entering lanczos_p"
  
  !.......... NOTES ON THICK RESTART......................
!   WE TARGET nkeep EIGENVECTORS
!   niter IS MAX # OF LANCZOS ITERATIONS BEFORE RESTARTING
!   nthick_keep = nkeep+nthick_add  IS HOW MANY "LANCZOS" VECTORS WE RESTART WITH
!   max_num_thick_restart = max # of times we restart
!
!....... CONVERGENCE......
  if (converge_test == 2 .or. converge_test == 3)then
     altconverge = .true.
  else
    altconverge = .false.
  end if

  if ( .not. fixed_iter ) then
	 if(targetX)then
	   ncheck = int(bmin(dimbasis, nkeep), 4) 
	 else
       ncheck = int(bmin(dimbasis, nkeep+ncheck_ex), 4) ! check extra states to prevent diving
     end if
     if(.not.allocated(eold)) allocate( eold(ncheck), stat=aerr )
     if(aerr /= 0) call memerror("lanczos_p 1")
  else
     ncheck = 0
  end if

!-----------------------------------------------------------------
  if(.not.allocated(alpha)) allocate( alpha(niter),beta(niter), stat=aerr )
  alpha = 0.0
  beta  = 0.0
  if(aerr /= 0) then
     call memerror("lanczos_p 2")
     stop 5
  end if
  if(.not.allocated(e)) allocate( eiglvec(niter,niter),e(niter), stat=aerr )
  if(aerr /= 0) then
     call memerror("lanczos_p 3")
     stop 5
  end if
  if(.not.allocated(eiglvecold)) allocate( eiglvecold(niter,niter), stat=aerr )
  if(aerr /= 0) then
     call memerror("lanczos_p 4")
     stop 5
  end if
   
  if(.not.initializing) call clocker('lan','sta')

  finished = .false.
  iter = 0
  i_thick_restart = 0
  i_thick = 0

!----------- RESTART OPTION -------------
  if ( startiter > 0 ) then
     iter = startiter-1
     if(iproc==0)then
          rewind(coef_file)
          do i = 1,startiter-1
             read(coef_file,*)j,alpha(i),beta(i)
          end do
     endif
     do i = 1,startiter-1
         call BMPI_BCAST(alpha(i),1,0,icomm,ierr)
         call BMPI_BCAST(beta(i),1,0,icomm,ierr)
     end do
  endif

  vchar='n'
!........ INITIAL NORMALIZATION.......

  if(noisy0) print *, "Adding first vector to history"
  call br_grab_vec1()          ! push vec1 -> br_reg
  call br_normalize(dnorm0)
  call br_restore_vec1()    !    pull vec1 <- br_reg
  call br_add2hist(1)       !   adds br_reg to history stack

  if(iproc==0 .and. strengthflag)then 
      print*,dnorm0*dnorm0,'   = total input strength '
     if(writeout)write(resultfile,*)dnorm0*dnorm0,' = total strength '
  endif  

  call   procOP_clock(0,'set','all')

! START LANCZOS ITERATIONS......................................
  actual_iterations = 0
  call BMPI_BARRIER(icomm, ierr)
  if(iproc == 0) then
     print *, "Before loop in lanczos_p"
     print *, "nfragments=", nfragments  ! This nfragments reference is ok
     print *, "vchar=", vchar
     print *, "storelanczosincoreMPI=", storelanczosincoreMPI
     print *, "alpha_before_orthog=", alpha_before_orthog
  end if
  if(wantnoisy0) call BMPI_BARRIER(icomm, ierr)

  do while ( .not. finished )
     iter = iter + 1
     if(noisy0) then
        print *, "Starting Lanczos iteration ", iter
        flush(6)
     end if
     if(wantnoisy0) call BMPI_BARRIER(icomm, ierr)
     if(vchar == 'n') then
        call print_vec_head(vec1, "vec1")
     else 
        call print_vec_head(vec2, "vec2")
     end if
     if (.not. useNewReorthog) then   ! write vector to disk
        if ((startiter == 0 .or. iter > startiter ) .and. reorthog ) then
           ! exception restart
           call write_lanczos_vector_a(vchar,'i',iter,lvec_file)
        end if
     end if
     if(noisy0) then
        print *, "initialize_final"
        flush(6)
     end if
     if(wantnoisy0) call BMPI_BARRIER(icomm, ierr)
     ! zeros output vector and thread specific output vectors if enabled
     call initialize_final(vchar)
     if(noisy0) then
        print *, "clocker"
        flush(6)
     end if
     if(wantnoisy0) call BMPI_BARRIER(icomm, ierr)	 
      if(.not.initializing) call clocker('hmu','sta')

     if(noisy0) then
        print *, "Calling applyHbundled_g"
        flush(6)
     end if
     if(wantnoisy0) call BMPI_BARRIER(icomm, ierr)


     if (useHZSomp) then
           call applyHbundled_omp(vchar)    ! THE CENTRAL ENGINE: APPLY MAT-VEC MULTIPLY
     else
            call applyHbundled_g(vchar)    ! THE CENTRAL ENGINE: APPLY MAT-VEC MULTIPLY
     end if
     if(noisy0) then
        print *, "After Calling applyHbundled_g"
        flush(6)
     end if
     if(wantnoisy0) call BMPI_BARRIER(icomm, ierr)

     if(.not.initializing) call clocker('hmu','end')
	  
     if(.not. useNewReorthog) call swap_vchar(vchar)
     actual_iterations = actual_iterations + 1

!.......... CAN COMPUTE ALPHA(ITER) EITHER BEFORE OR AFTER REORTHOGONALIZATION
!           while for exact arithmetic does not make a difference, it can 
!           be different for finite arithmetic.  Tests so far suggest 
!           converged energies agree to less than 1 keV and converged 
!           wavefunctions agree in overlap to less than 10^4 but be wary
!
!           USE DEFAULT-- ALPHA BEFORE REORTHOGONALIZTION
!.................... STANDARD REORTHOGONALIZATION from KSM...................
     if(useNewReorthog) then
        if(noisy0) print *, "Orthogonalizing and normalizing new Lanczos vector"
		if(iproc==0)call clocker('ort','sta')
        call br_grab_vec2()   ! push vec2 -> br_reg
        call br_remove_prev_overlap(da)
        call br_orthogonalize(1) ! 1 says ignore last history entry (prev)

        call br_normalize(db)
        ! Restore back to vec1 for next iteration
        ! This is instead of swap_vchar in older code
        call br_restore_vec1()
        alpha(iter) = real(da,kind(0.0e0))
        beta(iter) = real(db, kind(0.0e0))
        ! Finally, add br_reg to history
        if(iter < niter) call br_add2hist(iter+1)
		if(iproc==0)call clocker('ort','end')
		

        ! KSM TODO:  Have to implement restart - see below
     else
		 
!............... OLD REORTHOGONALIZATION... KEPT FOR BACK UP..............		 
!............... OLD REORTHOGONALIZATION... KEPT FOR BACK UP..............		 
!............... OLD REORTHOGONALIZATION... KEPT FOR BACK UP..............		 

        if(noisy0) print *, "Orthogonalizing and normalizing new Lanczos vector - old way"
        call dvecproj_p(vchar,da)   
        alpha(iter) = real(da,kind(0.0e0))
!----------- REORTHOGONALIZE ------------------------------------
        ! da is the overlap of v_i with v_{i-1}
        call reorthogonalize_a(iter,vchar,da)  
        if ( iter <= niter ) then
           if(iproc == 0 .or. storelanczosincoreMPI)then
              call dnormvec_p(vchar,'i',db,smallflag)     
              ! Calculate BETA.................................
              beta(iter) = real(db,kind(0.0e0))      
!------------- check for restart ------------------------------
              if ( smallflag .and. iter < niter) then
                 if(strengthflag)then
                     if(iproc ==0)write(6,*)' Exhausted space, exiting '
                     go to 3
                 end if
                 if ( iproc == 0 ) then
                    write(6,*) ' Restarting with random vector at iteration ',iter,beta(iter)
                 end if
                 call random_restart_p(iter,vchar)
              end if
            end if
        end if
        if(.not. storelanczosincoreMPI .and. nproc > 1)then
!......... BROADCAST RESULT............................
!          MUST BREAK UP FOR LARGE DIMENSION
           if(vchar=='n')then
! old              call Bcast_vec(vec1,dimbasis,MPI_lanc,0,icomm,ierr)    ! Changed to account for large vectors
              call BMPI_Bcast(vec1,dimbasis,0,icomm,ierr)   !  Now supports chunking
           else
! old              call Bcast_vec(vec2,dimbasis,MPI_lanc,0,icomm,ierr) 
              call BMPI_Bcast(vec2,dimbasis,0,icomm,ierr)  ! Now supports chunking
           end if
        end if
!.................. END OF OLD REORTHOGONALIZATION....................
!.................. END OF OLD REORTHOGONALIZATION....................
!.................. END OF OLD REORTHOGONALIZATION....................
		
     end if

!............ Write out Lanczos coefficients..............................

     if ( iproc ==0 ) write(coef_file,'(i10,2(2x,f20.15))')iter,alpha(iter), beta(iter)
     if ( iter == dimbasis ) cycle

!---------- compute and print intermediate eigenvalues.....................
! V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V 
     flyeig = (iter ==(iter/itskipeig)*itskipeig)
     if ( iproc==0 .and. .not.initializing .and. (( iter > nkeep .and. flyeig .and. checkeig) .or.      & 
          (iter >= ncheck .and. .not.fixed_iter .and. .not. skipeig .and. checkeig)  .or. & 
		     ( thick_restart .and. .not. finished .and. iter+1 == niter )))then

        if(noisy0) print *, "Finding lanczos eigenvalues"
        call find_lanczos_eigenvalues(iter,altconverge,  i_thick)

        if ( flyeig .or. thick_restart .and. iproc == 0 ) then       ! printout
           write(6,*)' '
           do j = 1, nkeep
              write(6,'(i5,f20.5)')j,e(j)
           end do
           write(6,*)' '
           write(6,*)iter,' iterations ' ! ' eigenvalues '
           if ( i_thick > 0 .and. iproc == 0 ) then
              write(6,*)' ( total = ',iter + niter*i_thick_restart,' ) '
           end if
        end if
        flush(6)
     end if
	 
! ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^
!----------------- end of compute and print intermediate eigenvalues..................

!---------------- CONVERGENCE CONTROL.........................

     if(noisy0) print *, "Checking Convergence"
	 
!--------------- CHECK FOR INITIALIZTION ----------------
!                IF FINISHED, THEN USE THICK_RESTART ROUTINE TO
!                CREATE INITIAL VECTOR(S)

     if(initializing .and. iter == initial_maxiter)then
		  call thick_restart_sub_p(nsave,initial_maxiter+1,0,niter,vchar)
		  call procOP_clock_out
		  return
	  end if
	 
!----------------- NORMAL CHECK -------------------------	 	 

     if ( iter == niter ) then
        finished = .true.
        
        if ( .not. fixed_iter .and. iproc == 0 ) then
           write(6,*)' Did not converge fully '
        end if
     end if
     if ( .not. fixed_iter ) then
		 
        if ( iter > ncheck  .and. (.not.skipeig .or. (skipeig .and. flyeig))) then
!................ CHECK DIFFERENCES IN ENERGY..........................

           ediff0 = 0.0
           ediff1 = 0.0
           do j = 1,ncheck
              ediff0 = ediff0 + abs(eold(j) - e(j) )
              ediff1 = bmax(ediff1, abs(  real( eold(j) -e(j),kind=4) ) )
           end do
           ediff0 = ediff0/sqrt(float(nkeep+1))*(1+i_thick_restart)

!................ CHECK DIFFERENCES IN VECTORS......................
           if(altconverge)then
              vdiff2 = 0.0d0
              vdiff3 = 0.0d0

              do j = 1,ncheck
                  vdot = 0.d0

                  do k = 1,iter-1
                      vdot = vdot + real(eiglvec(k,j),kind=8)*real(eiglvecold(k,j),kind=8)
                  end do
                  vdiff2 = vdiff2 + (abs(vdot)-1.d0)**2
                  vdiff3 = bmax( vdiff3, dabs(abs(vdot)-1.d0) )

              end do
              vdiff2 = dsqrt(vdiff2)/float(ncheck)
           endif 

!................. CONVERGE CHECK.............
           if ( flyeig .and. iproc == 0 ) then
                select case (converge_test)
                case (0)
		write(6,777)ediff0,converge_tol 
                case (1)
		write(6,777)ediff1,converge_tol 
                case (2)
		write(6,778)vdiff2,converge_tol 
                case (3)
		write(6,778)vdiff3,converge_tol 
                end select
           end if
777        format('      (energy convergence  ',f8.5,' > criterion ',f8.5,')')
778        format('      (wavefn convergence  ',f8.5,' > criterion ',f8.5,')')
           select case (converge_test)
           case ( 0)
              if ( ediff0 < converge_tol ) finished = .true.
           case ( 1)
              if ( ediff1 < converge_tol ) finished = .true.
           case ( 2)
              if ( vdiff2 < converge_tol ) finished = .true.
           case ( 3)
              if ( vdiff3 < converge_tol ) finished = .true.
           end select
        end if
		
        do j = 1,ncheck
           eold(j) = e(j)
           if(altconverge)then
                do k = 1,iter
                  eiglvecold(k,j) = eiglvec(k,j)
                end do
           end if
        end do
		
     end if
     call BMPI_BCAST(finished,1,0,icomm,ierr) 
 
!--------------- CHECK FOR THICK-RESTART --------------------------

     if ( thick_restart .and. .not. finished .and. iter+1 == niter ) then           
        if ( i_thick_restart > max_num_thick_restart ) then
           if ( iproc == 0 ) write(6,*)' WARNING could not finish '
           finished = .true.
           cycle
        end if
        
        if ( iproc == 0 ) write(6,*)' Doing a thick restart ',i_thick_restart

        call thick_restart_sub_p(nkeep+nthick_add,iter,i_thick_restart,niter,vchar)
        i_thick_restart = i_thick_restart + 1
        i_thick = nkeep + nthick_add   !ncheck
        iter = nkeep+nthick_add !  ncheck
     end if
  end do  ! iter

  if(noisy0) print *, "End of Lanczos interations"

!----------------------- LANCZOS ITERATIONS FINISHED ---------
   if(.not.initializing) call clocker('lan','end')

  if ( iproc == 0 ) then
     print*,iter,' iterations '
	 write(logfile,*)
     print*,' '
  end if
  if ( i_thick > 0 .and. iproc == 0) then
     print*,' ( total = ',iter + niter*i_thick_restart,' ) '
  end if
  if ( writeout .and. iproc == 0 ) then
     write(resultfile,*)iter,' iterations '
	 write(logfile,*)iter,' iterations '
     write(resultfile,*)' '
     if ( i_thick > 0 .and. iproc == 0 ) then
        write(resultfile,*)' ( total = ',iter + niter*i_thick_restart,' ) '
        write(logfile,*)' ( total = ',iter + niter*i_thick_restart,' ) '

     end if
  end if
3 continue

!  call proc_clock_out
!  call procOP_clock_out
!  call bundle_clock_out
  call all_clocks_out
  call timeperopmaster('out')  ! MUST BE LAST

  if(reorthog)then
     call  lanczos_output(iter,i_thick,dnorm0)
  end if

  return

end subroutine lanczos_p

!---- NOTE: density1b_output moved to bdenslib1.f90
!=================================================================
!
!  if in a small space, then diagonalize fully.
!  this helps to eliminate problems with Lanczos
!  the criterion is fulllimit found in module lanczos_info
!
! SUBROUTINES CALLED:
!   setup_localvectors
!
!=================================================================
subroutine exactdiag_p

  use system_parameters
  use flagger
  use sporbit
  use io
  use menu_choices
  use basis
!  use lanczos_info
!  use precisions
  use obs 
  use densities
  use haiku_info
  use nodeinfo
  use localvectors
  use coupledmatrixelements
  use fragments
  use mod_reorthog
  use wfn_mod
  use pocc_mod
  use bmpi_mod
  use butil_mod
  use btbme_mod
  use bvectorlib_mod
  use lanczos_util
  use jump_mod
  use apply_ham
  use apply_obs_mod
  implicit none
  
  integer(4) :: ierr
  real(kind=egv_prec), allocatable :: h(:,:)
  real(kind=lanc_prec), allocatable :: vi(:), vf(:)
  real(kind=egv_prec), allocatable :: u(:,:), work(:)
!  real(kind=egv_prec), allocatable :: e(:), u(:,:), work(:),eiglvec(:,:)

  integer(kind=basis_prec) :: i, j, k
  integer :: vidx
  real :: xj,xt

!---------------- VARIABLES FOR DENSITY MATRICES -----------------------
  integer ji,ti,jf,tf
  integer jmin,jmax,jt,tmin,tmax,tt
  real, allocatable :: denmat(:,:)
  real xjj,xtt,ei,ef
  logical smallflag,zeroflag
  integer n,m
  logical numberflag  ! flag to check on total number of particles; used for debugging
  logical printouthamflag  ! flag for printing out nonzeroes
  real :: nparticles ! temporary

  real(kind=8) :: dnorm0,dsclrprod,da
  character(1) :: vchar
  integer jvec
  integer info
  integer iorb
  character  :: outchoice  !used for steering output options

  printouthamflag = .false.
!---------- SET UP ARRAY FOR SINGLE-PARTICLE OCCUPATIONS --
  outchoice ='d'  ! default option
  if (spoccflag .and. .not.densityflag)then
       if(npeff(1)>0 .and. npeff(2) > 0)outchoice='b'
       if(npeff(1)>0 .and. npeff(2) == 0)outchoice='p'
       if(npeff(1)==0 .and. npeff(2) > 0)outchoice='n'
  end if

  if(nfragments > 1) then ! This nfragments reference is ok
     if(iproc == 0) print *, "exactdiag_p:  Not suported with nfragments > 1" ! this nfragments reference is ok
     stop 1
  end if
  niter = dimbasis
  call setup_localvectors
  vchar = 'n'
  allocate(h(dimbasis,dimbasis), stat=aerr)
  if(aerr /= 0) then
     call memerror("exactdiag_p 1")
     stop 5
  end if
  call   procOP_clock(0,'set','all')

  call clocker('lan','sta')
!----------------- CREATE FULL HAMILTONIAN MATRIX -----------
  if(iproc==0)open(unit=60,file='ham.dat',status='unknown')
  if(iproc==0 .and.printouthamflag)then
	  write(60,*)dimbasis
	  print*,' '
	  print*,' Writing hamiltonian to ham.dat '
  end if
  do i = 1, dimbasis
     vec1(:) = 0.d0
     vec1(i) = 1.d0
     vec2(:) = 0.d0
     if(useVec2Thread) vec2threadflat(:) = 0.0
     call applyHbundled_g(vchar)
     if(nproc > 1) then
         ! leaves data only on isfragroot nodes
         ! In normal use, would move to vector space (breorthog.f90)
         ! need to broadcast
         if(vchar == 'n') then
            ! output in vec2 with input 'n'
            call BMPI_BCAST(vec2, size(vec2), 0, fcomm2, ierr)
         else
            call BMPI_BCAST(vec1, size(vec1), 0, fcomm1, ierr)
         end if
     end if
     do j = 1,dimbasis
        h(i,j) = real(vec2(j),kind(egv_prec))   ! store as real(4); this can be changed

     end do  !J 
     if(dimbasis < 100 .and. iproc==0 .and. .not. printouthamflag)write(60,1)(h(i,j),j=1,dimbasis)
	 if(printouthamflag)then
		 do j = 1,i
			 if(h(i,j)/=0.0)write(60,*)i,j,h(i,j)
		 end do
	 end if
1  format(9f8.4)
  end do  !i
  print*,' Finished with creating the Hamiltonian '
  close(60)
  if(printouthamflag)stop
!------------------ DIAGONALIZE VIA HOUSEHOLDER -------------------

  allocate( e(dimbasis), work(dimbasis*3), stat=aerr)
  if(aerr /= 0) call memerror("exactdiag_p 2")
  if(egv_prec==4)then
  call SSYEV( 'V','U', dimbasis, h, dimbasis, e, WORK, 3*dimbasis, INFO )
  else
  call DSYEV( 'V','U', dimbasis, h, dimbasis, e, WORK, 3*dimbasis, INFO )

  end if
!  allocate( e(dimbasis), u(dimbasis,dimbasis), work(dimbasis), stat=aerr)
!  if(aerr /= 0) call memerror("exactdiag_p 3")
!  call eig(h,ndim,ndim,e,u,work)
!  call eigsrt(e,u,ndim,ndim)
!   call eig(h,dimbasis,dimbasis,e,u,work)
!   call eigsrt(e,u,dimbasis,dimbasis)

  if(iproc==0)then
     do i = 1,nkeep
        write(6,'(i4,f12.4)')i,e(i)
     end do
  end if
  call clocker('lan','end')

  call clocker('obs','start')

  call wfn_write_nkeep(nkeep) ! write number of vectors to wfn file

!------------  COMPUTE J^2, T^2
!------------ SET UP FOR ANG MOM
  call setup4obsmaster('J')
  call setup4obsmaster('T')
  if ( iproc == 0 ) print*,' '

  call reset_jumps4obs

!...................................................................
!----------- COMPUTE SINGLE PARTICLE OCCUPATIONS---------- added in 7.3.7
if(spoccflag .and. .not.densityflag)then
  call pocc_init_spoccx() ! allocate and initialize array
  call pocc_write_orbits()
end if
  
  if ( iproc == 0 ) print*,' '

88 format(2f10.5)
     if ( iproc == 0 ) then
        write(6,*)' State      E        Ex         J       T '
        if(writeout)write(resultfile,*)' State      E        Ex         J      T '
     end if
  
! Construct nkeep eigenvectors of H from  Lanczos eigenvectors....

  if(densityflag .or. trdensout)then  ! must store in lvec
     if(allocated(lvec))deallocate(lvec)
     allocate( lvec(dimbasis, nkeep), stat=aerr )
     if(aerr /= 0) call memerror("exactdiag_p 4")
     allocate(energy(nkeep), xjlist(nkeep),xtlist(nkeep), stat=aerr )
     if(aerr /= 0) call memerror("exactdiag_p 5")
  end if

  do vidx = 1, nkeep

        do k = 1, dimbasis
           vec1(k) = real(h(k,vidx), kind=lanc_prec)  ! u(k,vidx)
        end do !k
        call br_grab_vec1()
        call br_load_vec2_from_vec1()
		call br_add2hist(vidx)

     if(densityflag .or. trdensout)call write_lanczos_vector_a('n','i',vidx,lvec_file)
!------------ COMPUTE J2, T2................................................
     call clocker('aob','sta')
     twoobsflag = .true.  ! compute both J^2 and T^2
     xj2 = 0.0
     xt2 = 0.0
     call applyobsbundled(1)
     call clocker('aob','end')
     xj = real( -0.5 + sqrt(xj2 + 0.25), kind=4)
     xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
	 
     if(densityflag .or. trdensout)then
        energy(vidx) = real(e(vidx), kind(energy(1)))
        xjlist(vidx) = xj
        xtlist(vidx) = xt
     end if
!----------------WRITE OUT WFN..............................................
     call clocker('wev','sta')
     if ( writeout .and. write_wfn) then
        call wfn_writeeigenvec(wfnfile, frag1, vec1, vidx, real(e(vidx),kind=4), xj, real(xt2,kind=4))
     end if
     call clocker('wev','end')
	 
!----------- COMPUTE SINGLE PARTICLE OCCUPATIONS---------- added in 7.3.7 
	  call pocc_compute_spocc(vidx, .true.)  ! true for restore J/T setup

!----------------- WRITE OUT RESULTS.............................
	         if ( iproc == 0 ) then
	            select case (outchoice)

	                case('d')
	                  call pocc_write_ejt(vidx, e, xj, xt)

	                case('b')
	                  ! write(6,12)vidx,e(vidx), e(vidx) - e(1),xj,xt,(spoccx(1,vidx,iorb),iorb=1,numorb(1))
	                  call pocc_write_ejt(vidx, e, xj, xt)
	                  call pocc_write_occvec(6, spoccx, vidx, 1, "    p occ:")
	                  call pocc_write_occvec(6, spoccx, vidx, 2, "    n occ:")
	                  if ( writeout ) then
	                     call pocc_write_occvec(resultfile, spoccx, vidx, 1, "    p occ:")
	                     call pocc_write_occvec(resultfile, spoccx, vidx, 2, "    n occ:")
	                  end if
	                case('p')
	                  call pocc_write_ejt(vidx, e, xj, xt)
	                  call pocc_write_occvec(6, spoccx, vidx, 1, "    p occ:")

	                  if ( writeout ) then
	                     call pocc_write_occvec(resultfile, spoccx, vidx, 1, "    p occ:")
	                  end if

	                case('n')
	                  call pocc_write_ejt(vidx, e, xj, xt)
	                  call pocc_write_occvec(6, spoccx, vidx, 2, "    n occ:")

	                  if ( writeout ) then
	                     call pocc_write_ejt(vidx, e, xj, xt)
	                     call pocc_write_occvec(resultfile, spoccx, vidx, 2, "    n occ:")
	                  end if
	             end select
	         end if	 
	 
  end do
11 format(i5,3x,2f10.5,2x,2f8.3)

  if(densityflag)call density1b_output !(e,eiglvec)

  if ( trdensout  ) then
     call output_TRDENS
  end if
  call close_lanczosfile

  call clocker('obs','end')
  return
end subroutine exactdiag_p

!=================================================================
!
!  if in a small space, then diagonalize fully.
!  this helps to eliminate problems with Lanczos
!  the criterion is fulllimit found in module lanczos_info
!
! SUBROUTINES CALLED:
!   setup_localvectors
!
!=================================================================
subroutine exactdiag_MPI

  use system_parameters
  use sporbit
  use io
  use menu_choices
  use basis
!  use lanczos_info
!  use precisions
  use obs
! use results   
  use densities
  use haiku_info
  use nodeinfo
  use localvectors
!  use interaction
  use coupledmatrixelements
  use fragments
  use mod_reorthog
  use wfn_mod
  use pocc_mod
  use bmpi_mod
  use butil_mod
  use btbme_mod
  use bvectorlib_mod
  use lanczos_util
  use jump_mod
  use apply_ham
  use apply_obs_mod
  implicit none
!  include 'binterfaces.inc'
  
  integer(4) :: ierr
  real(kind=egv_prec), allocatable :: h(:,:),hvec(:)
  real(kind=egv_prec) :: htmp
  real(kind=lanc_prec), allocatable :: vi(:), vf(:)
  real(kind=egv_prec), allocatable ::  u(:,:), work(:)
!  real(kind=egv_prec), allocatable :: e(:), u(:,:), work(:),eiglvec(:,:)

  integer(kind=basis_prec) :: i,j,k
  real :: xj,xt

!---------------- VARIABLES FOR DENSITY MATRICES -----------------------
  integer ji,ti,jf,tf
  integer jmin,jmax,jt,tmin,tmax,tt
  real, allocatable :: denmat(:,:)
  real xjj,xtt,ei,ef
  logical smallflag,zeroflag
  integer n,m
  logical numberflag  ! flag to check on total number of particles; used for debugging
  real :: nparticles ! temporary

  real(kind=8) :: dnorm0,dsclrprod,da
  character(1) :: vchar
  integer jvec
  integer info
  integer :: ki

  call setup_localvectors
  vchar = 'n'
  k = 1;
  if(iproc==0) k = dimbasis;
  allocate(h(k,k), stat=aerr) ! always allocate, suppresses compiler warning
  if(aerr /= 0) then
     call memerror("exactdiag_p 1")
     stop 5
  end if
  if(iproc==0)h=0.d0
  allocate(hvec(dimbasis*dimbasis), stat=aerr)
  if(aerr /= 0) then
     call memerror("exactdiag_p 1.5")
     stop 5
  end if
  hvec= 0.d0
  call   procOP_clock(0,'set','all')

  call clocker('lan','sta')
!----------------- CREATE FULL HAMILTONIAN MATRIX -----------
  if(iproc==0)open(unit=60,file='ham.dat',status='unknown')
!  print*,' sub dimensions ',v1s,v1e,v2s,v2e
  do i = 1, dimbasis
     vec1(:) = 0.d0
     if(i >= v1s .and. i <= v1e)vec1(i) = 1.d0
     vec2(:) = 0.d0
     if(useVec2Thread) vec2threadflat(:) = 0.0
     call applyHbundled_g(vchar)
	 
!     if(nproc > 1) then
         ! leaves data only on isfragroot nodes
         ! In normal use, would move to vector space (breorthog.f90)
         ! need to broadcast
!         if(vchar == 'n') then
            ! output in vec2 with input 'n'
!            call BMPI_BCAST(vec2, size(vec2), 0, fcomm2, ierr)
!         else
!            call BMPI_BCAST(vec1, size(vec1), 0, fcomm1, ierr)
!         end if
!     end if
     !      if(i ==1)write(46,*)vf
	 if( i>=v1s .and. i <=v1e)then
!		 print*,i,iproc,vec2
     do j = v2s,v2e
		 
        hvec((i-1)*dimbasis+j) = real(vec2(j),kind(egv_prec))   ! store as real(4); this can be changed

     end do  !J 
 end if


  end do  !i
  print*,' all done ',iproc,icomm
!  print*,iproc,' hvec ',hvec
  call BMPI_ALLREDUCE(hvec,size(hvec),MPI_SUM,icomm,ierr)
  do i = 1,dimbasis
 	 do j = 1,dimbasis
!		 call MPI_REDUCE(hvec((i-1)*dimbasis+j),htmp,1,MPI_REAL8,MPI_SUM,0,icomm,ierr)
 		 if(iproc==0)h(i,j)=hvec((i-1)*dimbasis+j)*nfragments**2/nproc
 	 end do
	 
  end do
  if(iproc /=0)return
  
 1  format(9f8.4)
  
  if(iproc==0)then 
	  do i = 1,dimbasis
	      if(dimbasis < 100 .and. iproc==0)write(60,1)(h(i,j),j=1,dimbasis)
	  end do
	  close(60)
	  print*,' WROTE OUT HAMILTONIAN '
  end if
!------------------ DIAGONALIZE VIA HOUSEHOLDER -------------------

  allocate( e(dimbasis), work(dimbasis*3), stat=aerr)
  if(aerr /= 0) call memerror("exactdiag_p 2")
  if(egv_prec==4)then
  call SSYEV( 'V','U', dimbasis, h, dimbasis, e, WORK, 3*dimbasis, INFO )
  else
  call DSYEV( 'V','U', dimbasis, h, dimbasis, e, WORK, 3*dimbasis, INFO )

  end if
!  allocate( e(dimbasis), u(dimbasis,dimbasis), work(dimbasis), stat=aerr)
!  if(aerr /= 0) call memerror("exactdiag_p 3")
!  call eig(h,ndim,ndim,e,u,work)
!  call eigsrt(e,u,ndim,ndim)
!   call eig(h,dimbasis,dimbasis,e,u,work)
!   call eigsrt(e,u,dimbasis,dimbasis)

  if(iproc==0)then
     do ki = 1,nkeep
        write(6,'(i4,f12.4)')ki,e(ki)
     end do
  end if
  call clocker('lan','end')

  call clocker('obs','start')

  call wfn_write_nkeep(nkeep) ! write number of vectors to wfn file

!------------  COMPUTE J^2, T^2
!------------ SET UP FOR ANG MOM
  call setup4obsmaster('J')
  call setup4obsmaster('T')
  if ( iproc == 0 ) print*,' '

  call reset_jumps4obs

!...................................................................

  if ( iproc == 0 ) print*,' '

88 format(2f10.5)
     if ( iproc == 0 ) then
        write(6,*)' State      E        Ex         J       T '
        if(writeout)write(resultfile,*)' State      E        Ex         J      T '
     end if
  
! Construct nkeep eigenvectors of H from  Lanczos eigenvectors....

  if(densityflag)then  ! must store in lvec
     if(useNewReorthog) then
         print *, "lvec not supported when new reorthog is used"
         stop 1
     end if
     if(allocated(lvec))deallocate(lvec)
     allocate( lvec(dimbasis, nkeep), stat=aerr )
     if(aerr /= 0) call memerror("exactdiag_p 4")
     allocate(energy(nkeep), xjlist(nkeep),xtlist(nkeep), stat=aerr )
     if(aerr /= 0) call memerror("exactdiag_p 5")
  end if

  do ki = 1, nkeep

     do k = 1, dimbasis
           vec1(k) = real( h(k,ki), kind=lanc_prec)  ! u(k,ki)
     end do !k
     call br_load_vec2_from_vec1()

     if(densityflag)call write_lanczos_vector_a('n','i',ki,lvec_file)
!------------ COMPUTE J2, T2................................................
     call clocker('aob','sta')
     twoobsflag = .true.  ! compute both J^2 and T^2
     xj2 = 0.0
     xt2 = 0.0
     call applyobsbundled(1)
     call clocker('aob','end')
     xj = real( -0.5 + sqrt(xj2 + 0.25), kind=4)
     xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)

!----------------- WRITE OUT RESULTS........................................
     call pocc_write_ejt(ki, e, xj, xt)

     if(densityflag)then
        energy(ki) = real(e(ki), kind(energy(1)))
        xjlist(ki) = xj
        xtlist(ki) = xt
     end if
!----------------WRITE OUT WFN..............................................
     call clocker('wev','sta')
     if ( writeout .and. write_wfn) then
        call wfn_writeeigenvec(wfnfile, frag1, vec1, ki, real(e(ki),kind=4), xj, real(xt2,kind=4))
     end if
     call clocker('wev','end')
  end do
11 format(i5,3x,2f10.5,2x,2f8.3)

  if(densityflag)call density1b_output !(e,eiglvec)

  if ( trdensout  ) then
     call output_TRDENS
  end if
  call close_lanczosfile

  call clocker('obs','end')
end subroutine exactdiag_MPI

