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
  use nodeinfo
  use localvectors
  use localblocks
  use precisions
  use fragments
  use coupledmatrixelements
  use bmpi_mod
  use butil_mod
  use para_main_mod
  use para_util_mod
  use lanczos_util
  use apply_ham_omp
  use menu_choices
  use convergence
  use mod_balg
  use jumpstart
  
  implicit none
  logical    :: restart
  integer(4) :: ierr
  integer(4) :: maxiter,max_iter_restart
  integer(4) :: niter0
!  character(2) :: lanczchar
  logical :: leave  ! flag to abort
  integer :: aerr
  integer(8) :: maxfragment
  integer  :: i
  double precision :: get_Wtime

  kpow = 1  ! default exponent
  coef_only=.false.
  leave = .false.
  thick_restart = .false.
  ! Block Lanczos main loop switch
  block_flag = 0

  if ( restart ) then
 !    call countlanczositerations
     print*,' RESTARTING NOT CURRENTLY IMPLEMENTED '
  else
     if(reorthog)call open_lanczosfile
     startiter = 0
  end if
  
!............. MENU CHOICE TO WRITE OUT HAMILTONIAN TO FILE......

if(menu_char=='wh')then
	
	if(iproc==0)print*,' Writing out Hamilonian nonzeroes to ham.dat '
    if(enablelanczosincoreMPI .and. nproc > 1)then
            storelanczosincoreMPI = .true.
    else
            storelanczosincoreMPI = .false.
    end if
    call distribute_lanczos_pieces
	call exactdiag_p
	return
end if  

!................ IF JUMPSTART THEN FIXED NUMBER OF ITERATIONS.........
! ADDED 7.9.6
! DELETE 7.11.4

  
!.............. IF STRENGTH FUNCTION MENU OPTION 'S' THEN ONLY FIXED LANCZOS .............

if(strengthflag)then
	
	if(blockstrengthflag)then
		lanczchar='bf'
		if(iproc==0)then
			if(auto_input)then
				read(autoinputfile,*)nkeep,dimblock,maxblockiter
				if(nkeep < dimblock*(maxblockiter+1))then
					nkeep = dimblock*(maxblockiter+1)
				end if
				write(6,*)' parameters of block Lanczos: '
				write(6,*)' # states to keep, dimension of block, # block iterations '	
					
		        write(6,*)nkeep,dimblock, niter
			else
			    write(6,*)' Fixed block iterations ONLY: '
	            write(6,*)' Enter nkeep, dim of block,  # block Lanczos iterations   '
	            write(6,*)' (nkeep = # of states printed out, should be ~ dim block x (# block iters +1) )'
				if(nproc > 1)then
					write(6,*)' Max block dimeension is ',2e9/maxfragsize_all,'; if you need more '
					write(6,*)' decrease the fragment size and possibly add more MPI processes '
				end if
				
	            read(5,*)nkeep,dimblock, maxblockiter
				if(nkeep < dimblock*(maxblockiter+1))then
					nkeep = dimblock*(maxblockiter+1)
					print*,' Change Nkeep to ',nkeep
				end if
	    	end if
			
			if(nproc > 1 .and. dimblock*maxfragsize_all > 2e9)then
				if(iproc==0)then
					write(6,*)' Warning! Your block dimension is too large to work with MPI '
					write(6,*)' Max block dimeension is ',2e9/maxfragsize_all,'; if you need more '
					write(6,*)' decrease the fragment size and possibly add more MPI processes '
					write(logfile,*)' Warning! Your block dimension is too large to work with MPI '
					write(logfile,*)' Max block dimeension is ',2e9/maxfragsize_all,'; if you need more '
					write(logfile,*)' decrease the fragment size and possibly add more MPI processes '
					write(resultfile,*)' Warning! Your block dimension is too large to work with MPI '
					write(resultfile,*)' Max block dimeension is ',2e9/maxfragsize_all,'; if you need more '
					write(resultfile,*)' decrease the fragment size and possibly add more MPI processes '
				end if
			end if
			
			
	        niter = maxblockiter*dimblock + dimblock ! We can leave br_histmax alone
			
	        if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
	        if(niter < nkeep) then
				maxblockiter = nkeep/niter
				niter = maxblockiter*dimblock + dimblock
				
				niter = nkeep
			end if
	        if(niter >= dimbasis) niter = int(dimbasis-1, 4)
	        if(.not. auto_input)then
				write(autoinputfile,*)nkeep,dimblock,maxblockiter, & 
				 '   ! # states to keep, dim of block, # block iters '	
			end if
			block_flag = 1	
		
		end if
		
	else
	
	lanczchar='lf'
	if(iproc==0)then
		if(auto_input)then
			read(autoinputfile,*)nkeep,niter
	        write(6,*)nkeep,niter,'   ! # states to keep, # iterations '		
		else
		    write(6,*)' Fixed iterations ONLY: '
            write(6,*)' Enter nkeep, # iterations for lanczos  '
            write(6,*)' (nkeep = # of states printed out )'
            read(5,*)nkeep,niter
    	end if
		
        if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
        if(niter < nkeep) then
			niter = nkeep
		end if
        if(niter >= dimbasis) niter = int(dimbasis-1, 4)
        if(.not. auto_input)write(autoinputfile,*)nkeep,niter,'   ! # states to keep, # iterations '		
		
	end if
    end if
	leave = .false.
	goto 543
	
end if
!..................... RESOLVENT/ GREEN'S FUNCTION.......... added 7.9.10............
!                      MODIFIED FOR COMPLEX ENERGY 7.10.6
Egf = 0.0
Egfi= 0.0

if(greenflag)then

		if(iproc==0)then
		    print*,' '
			print*,' + + + + + + + + + + + + + + + + + + + + + + + + + + + '
			print*,' + Calculating resolvent/Green function on an initial vector '
			print*,' '
			if(auto_input)then
				if(complex_green)then
  				   read(autoinputfile,*)Egf,Egfi
				   print*,' Complex Green function Energy = ',Egf,' + i ',Egfi
					
				else
				
 				   read(autoinputfile,*)Egf
				   print*,' Green function Energy = ',Egf
			   end if
			   read(autoinputfile,*)green_tol
				print*,' Tolerance for convergence = ',green_tol
			else
			  print*,' '
			  print*,' Enter energy E in resolvent/Green function 1/(E-H)'
			  if(complex_green)then
				  print*,' Enter real, imaginary parts of E '
				  read*,Egf,Egfi
			  	  print*,' Energy = ',Egf,' + i',Egfi
			      write(autoinputfile,*)Egf,Egfi,'  ! complex energy for resolvent/Green function'
			  else
				  read*,Egf
			  	  print*,' Energy = ',Egf
			      write(autoinputfile,*)Egf,'  ! energy for resolvent/Green function'
			  end if
			  print*,' Enter tolerance for convergence of Green function (typical = 0.00001) '
			  read*,green_tol
			  write(autoinputfile,*)green_tol,' ! tolerance for Green function convergence '
		    end if
		    if(auto_input)then
			  read(autoinputfile,*)niter
			  print*,' Will do a maximum ',niter,' iterations '
		    else
			  print*,' Enter number of iterations '
			  print*,' (If number of iterations >= dimension, will do full matrix )'
			  read*,niter
			  write(autoinputfile,*)niter,' ! # of iterations'
		    end if
		  
		end if
#ifdef _MPI		
	    call BMPI_BCAST(Egf,1,0,MPI_COMM_WORLD,ierr)
		call BMPI_BCAST(niter,1,0,MPI_COMM_WORLD,ierr)
#endif
		if(complex_green)then
			nkeep=2
		else
		   nkeep = 1
	    end if
	    niter = min(niter,dimbasis)
		maxiter = niter
		fixed_iter =.true.
	
		allocate(gvec(niter+2),gvec_old(niter+2))
		if(complex_green)allocate(gvec_i(niter+2),gvec_i_old(niter+2))
	leave = .false.
	
	goto 543
end if
!........................................................................................  
! ... related to Green's function/resolvent: option 'pv'
!     This option does Lanczos on a previously computed vector and writes out Lanczos coefficients
!     but does not compute a spectrum; instead it computes the overlaps of Lanczos vectors (not eigenvectors)
!     on a read-in file of previously computed vectors. This is for computing Green's function without fixed energy.


if(menu_char=='pv')then
	if(iproc==0)then
	    print*,' '
		print*,' + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +'
		print*,' + Calculating elements of resolvent/Green function on an initial vector '
		print*,' + (output Lanczos coefficients to .lcoef file, plus '
		print*,' +  overlaps of Lanczos vectors with precomputed wave vectors)'
		print*,' '
	    if(auto_input)then
		  read(autoinputfile,*)niter
		  print*,' Will do ',niter,' iterations (fixed) '
	    else
		  print*,' Enter number of fixed iterations '
		  print*,' (If number of iterations >= dimension, will do full matrix )'
		  read*,niter
		  write(autoinputfile,*)niter,' ! # of iterations'
	    end if		
		
		
	end if
#ifdef _MPI		
			call BMPI_BCAST(niter,1,0,MPI_COMM_WORLD,ierr)
#endif
	niter = min(niter,dimbasis-1)
	maxiter = niter
	fixed_iter =.true.
	lanczchar='lf'
	coef_only=.true.
	goto 543
	
end if

!........................................................................................  

  if ( iproc == 0 ) then
11   continue

!..... THE FOLLOWING USED ONLY IN OPTION 'tx' TARGETING A SPECIFIC ENERGY; may be revised
     Etarget = 0.d0
     targetX = .false.
!....................................	 
     print*,' '
     print*,' / ------------------------------------------------------------------------\ ' 
     print*,' |                                                                         | '
     print*,' |    DIAGONALIZATION OPTIONS (choose one)                                 | '
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
!---- ADDED in 7.9.9
!     print*,' | (lp) Lanczos with (H-E0)^k acceleration/interior eigenpairs.            | '
     print*,' |                                                                         | '
     print*,' | (bd) Block Lanczos with default convergence (STANDARD)                  | '
     print*,' | (bf) Block Lanczos with fixed (user-chosen) iterations                  | '
     print*,' | (bc) Block Lanczos with user-defined convergence                        | '

     print*,' |                                                                         | '
!     print*,' |                                                                         | '

 if(.not.strengthflag)then     ! cannot use strength option with thick-restart

     print*,' | (td) Thick-restart Lanczos with default convergence                     | '
     print*,' | (tf) Thick-restart Lanczos with fixed iterations                        | '
     print*,' | (tc) Thick-restart Lanczos with user-defined convergence                | '
!------- ADDED in 7.6.8 -------	 
     print*,' | (tx) Thick-restart Lanczos targeting states near specified energy       | '
!---- ADDED in 7.9.9
!     print*,' | (tp) Thick-restart power-filter Lanczos targeting specified energy      | '
     !---- ADDED in 7.9.12     
     print*,' | (tb) Thick-restart block Lanczos with default convergence               | '
end if
     print*,' |                                                                         | '
     print*,' | (sk) Skip Lanczos (only used for timing set up)                         | '
     print*,' | (li) Lanczos iterations only, no further eigensolutions                 | '
	 
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
			  
		   case ('li')  ! do LANCZOS ITERATIONS ONLY-- HIDDEN OPTION----
   		      coef_only=.true.		   

		      write(6,*)' Lanczos iterations ONLY to get Lanczos coefficients '
		      read(autoinputfile,*)niter
              if(niter >= dimbasis) niter = int(dimbasis-1, 4)
              print*,' Fixed Lanczos iterations with # iter = ',niter		
		   
           case('lp')
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
			  read(autoinputfile,*)kpow
			  if(mod(kpow,2)==0)read(autoinputfile,*)Etarget
			  
			  write(6,*)' Filtering with exponent k = ',kpow,', target energy of ',Etarget 
			  
			  
		   
		   case('bd')
		      print*,' Block Lanczos with default convergence'
		      read(autoinputfile,*)nkeep,dimblock,maxblockiter
              print*,' Keeping ',nkeep,' eigenstates with maximum ',maxblockiter,' blocks iterations ',&
		          ' of blocks of size ',dimblock
		      if(menu_char=='np')then
			     print*,' You have read in pivot wave functions (recommended) '
			     piv_flag='y'
		      else
			     print*,' You did not read in block of pivot wave functions  '
			     print*,' Starting from random pivots (not recommended) '	
			     print*,' (to read in pivots, chose "np" in starting menu )'  
			     piv_flag='n'
		      end if

!              if((piv_flag .eq. 'y') .and. (nkeep .lt. dimblock))nkeep = dimblock
		   
              niter = maxblockiter*dimblock + dimblock ! We can leave br_histmax alone
		      block_flag = 1
              ncheck_ex = 2
              converge_test = converge_test_def
              converge_tol  = converge_tol0_def
              fixed_iter = .false.
			    
		   case('bf')
		      print*,' Block Lanczos with fixed iterations '
		      read(autoinputfile,*)nkeep,dimblock,maxblockiter
              print*,' Keeping ',nkeep,' eigenstates with ',maxblockiter,' blocks iterations ',&
		          ' of blocks of size ',dimblock
		      if(menu_char=='np')then
			     print*,' You have read in pivot wave functions (recommended) '
			     piv_flag='y'
		      else
			     print*,' You did not read in block of pivot wave functions  '
			     print*,' Starting from random pivots (not recommended) '	
			     print*,' (to read in pivots, chose "np" in starting menu )'  
			     piv_flag='n'
		      end if

 !             if((piv_flag .eq. 'y') .and. (nkeep .lt. dimblock))nkeep = dimblock
		   
              niter = maxblockiter*dimblock + dimblock ! We can leave br_histmax alone
		      block_flag = 1
              fixed_iter = .true.
		   
           case('bc')
              print*,' User-defined block Lanczos convergence '
		      read(autoinputfile,*)nkeep,dimblock,maxblockiter
              print*,' Keeping ',nkeep,' eigenstates with maximum ',maxblockiter,' blocks iterations ',&
		          ' of blocks of size ',dimblock
		      if(menu_char=='np')then
			     print*,' You have read in pivot wave functions (recommended) '
			     piv_flag='y'
		      else
			     print*,' You did not read in block of pivot wave functions  '
			     print*,' Starting from random pivots (not recommended) '	
			     print*,' (to read in pivots, chose "np" in starting menu )'  
			     piv_flag='n'
		      end if

!              if((piv_flag .eq. 'y') .and. (nkeep .lt. dimblock))nkeep = dimblock
		   
              niter = maxblockiter*dimblock + dimblock ! We can leave br_histmax alone
		      block_flag = 1
			  
	          if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
	          if(maxblockiter*dimblock < nkeep) maxiter = nkeep/dimblock
	          if(maxiter*dimblock >= dimbasis) maxiter = int(dimbasis/dimblock-1, 4)
              niter = maxblockiter*dimblock + dimblock ! We can leave br_histmax alone
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
		   
		   case('tp')
           print*,' Thick-restart power-filter Lanczos targeting specific energy '
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
		   read(autoinputfile,*)kpow
		   write(6,*)' Power k for (H-E0)^k is ',kpow
           read(autoinputfile,*)Etarget
		   print*,' Targeting energy = ',Etarget

           case('sk')
              print*,' Skipping diagonalization '
			  
           case('tb')
              print*,' Thick-restart block Lanczos with default convergence '
              maxblockiter = 1
		      read(autoinputfile,*)nkeep,dimblock,maxblockiter
              print*,' Keeping ',nkeep,' eigenstates with maximum ',maxblockiter,' blocks iterations ',&
		     &' of blocks of size ',dimblock
           
           niter = maxblockiter*dimblock + dimblock
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
	
		      if(menu_char=='np')then
			     print*,' You have read in pivot wave functions (recommended) '
			     piv_flag='y'
		      else
			     print*,' You did not read in block of pivot wave functions  '
			     print*,' Starting from random pivots (not recommended) '	
			     print*,' (to read in pivots, chose "np" in starting menu )'  
			     piv_flag='n'
		      end if

		        block_flag = 1
              ncheck_ex = 2
              converge_test = converge_test_def
              converge_tol  = converge_tol0_def
              thick_restart = .true.
              fixed_iter = .false.

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
			  
			  if(nkeep > dimbasis .or. maxiter >= dimbasis .or. maxiter < nkeep-1)then
				  print*,' '
				  print*,' Too many final states/iterations requested '
				  print*,' (You may find the "ex" option more appropriate)'
				  print*,' '
			  end if
              if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
              if(maxiter < nkeep) maxiter = nkeep
              if(maxiter >= dimbasis) maxiter = int(dimbasis-1, 4)
              niter = maxiter
              write(autoinputfile,*)nkeep,maxiter,'     ! # states to keep, max # iterations '
              ncheck_ex = ncheck_ex_def
              converge_test = converge_test_def
              converge_tol  = converge_tol0_def
			  
		   case ('li','LI','Li')
		      lanczchar = 'li'
              write(autoinputfile,'(a,"    ! Lanczos menu option ")')lanczchar(1:2)
			  write(6,*)' Lanczos iterations ONLY to get Lanczos coefficients '
              write(6,*)' Enter # iterations for lanczos '
			  read(5,*)niter
              if(niter >= dimbasis) niter = int(dimbasis-1, 4)
              write(autoinputfile,*)niter,'     ! # iterations '
			  coef_only=.true.
		   !
           case('lp','LP','Lp')
  		      lanczchar='lp'
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
			  write(6,*)' Enter power k for (H-E0)^k (typically 2, 3, 4,...)'
			  read(5,*)kpow  
			  write(autoinputfile,'(i3," ! power for Lanczos polynomial filtering ")')kpow
			  if(mod(kpow,2)==0)then
				  targetX = .true.
				  write(6,*)' Enter target energy E0 (absolute, not relative to g.s.) '
				  read(5,*)Etarget
				  write(autoinputfile,'(f10.5," ! target energy ")')Etarget
				  
			  end if
			  
           case ('bd','BD','Bd')
              lanczchar = 'bd'
              write(autoinputfile,'(a,"    ! Block Lanczos menu option ")')lanczchar(1:2)
              write(6,*)' Enter nkeep, dimension of blocks, max # of block iterations '
              write(6,*)' (nkeep = # of states printed out; typically ~ dim block )'
              read(5,*) nkeep,dimblock,maxblockiter
              if(nkeep > dimbasis)nkeep = int(dimbasis, 4)

              write(autoinputfile,*)nkeep,dimblock,maxblockiter,'     ! # states to keep, dim of block,# block iters '		
			  if(menu_char=='np')then
				  print*,' You have read in pivot wave functions (recommended) '
				  piv_flag='y'
			  else
				  print*,' You did not read in block of pivot wave functions  '
				  print*,' Starting from random pivots (not recommended) '	
				  print*,' (to read in pivots, chose "np" in starting menu )'  
				  piv_flag='n'
			  end if
!              if((piv_flag .eq. 'y') .and. (nkeep .lt. dimblock))nkeep = dimblock
              niter = maxblockiter*dimblock + dimblock ! We can leave br_histmax alone
              block_flag = 1
			  
              ncheck_ex = 0
              converge_test = converge_test_def
              converge_tol  = converge_tol0_def
              fixed_iter = .false.
			  
           case ('bf','BF','Bf')
              lanczchar = 'bf'
              write(autoinputfile,'(a,"    ! Block Lanczos menu option ")')lanczchar(1:2)
              write(6,*)' Enter nkeep, dimension of blocks, # of block iterations '
              write(6,*)' (nkeep = # of states printed out; typically ~ dim block )'
              read(5,*) nkeep,dimblock,maxblockiter
              if(nkeep > dimbasis)nkeep = int(dimbasis, 4)

              write(autoinputfile,*)nkeep,dimblock,maxblockiter,'     ! # states to keep, dim of block,# block iters '		
			  if(menu_char=='np')then
				  print*,' You have read in pivot wave functions (recommended) '
				  piv_flag='y'
			  else
				  print*,' You did not read in block of pivot wave functions  '
				  print*,' Starting from random pivots (not recommended) '	
				  print*,' (to read in pivots, chose "np" in starting menu )'  
				  piv_flag='n'
			  end if
!              if((piv_flag .eq. 'y') .and. (nkeep .lt. dimblock))nkeep = dimblock
              niter = maxblockiter*dimblock + dimblock ! We can leave br_histmax alone
              block_flag = 1
			  converge_test=0
              fixed_iter = .true.
			  
           case ('bc','BC','Bc')
              lanczchar = 'bc'
              write(autoinputfile,'(a,"    ! Lanczos menu option ")')lanczchar(1:2)
              print*,' '
              write(6,*)' Enter nkeep, dimension of blocks, max # of block iterations '
              write(6,*)' (nkeep = # of states printed out; typically ~ dim block )'
              read(5,*) nkeep,dimblock,maxblockiter
              if(nkeep > dimbasis)nkeep = int(dimbasis, 4)

              write(autoinputfile,*)nkeep,dimblock,maxblockiter,'     ! # states to keep, dim of block,# block iters '		
			  if(menu_char=='np')then
				  print*,' You have read in pivot wave functions (recommended) '
				  piv_flag='y'
			  else
				  print*,' You did not read in block of pivot wave functions  '
				  print*,' Starting from random pivots (not recommended) '	
				  print*,' (to read in pivots, chose "np" in starting menu )'  
				  piv_flag='n'
			  end if
!              if((piv_flag .eq. 'y') .and. (nkeep .lt. dimblock))nkeep = dimblock
              niter = maxblockiter*dimblock + dimblock ! We can leave br_histmax alone
              block_flag = 1

              print*,' Enter how many ADDITIONAL states for convergence test '
              print*,' ( Default= ',ncheck_ex_def,'; you may choose 0 ) '
              read*,ncheck_ex
              write(autoinputfile,*)ncheck_ex

3321           continue 
              print*,' '
              print*,'    Enter one of the following choices for convergence control :'
              print*,'(0) Average difference in energies between one iteration and the last; '
              print*,'(1) Max difference in energies between one iteration and the last; '
              print*,'(2) Average difference in wavefunctions between one iteration and the last; '
              print*,'(3) Min difference in wavefunctions between one iteration and the last; '
              read*,converge_test
              if(converge_test < 0 .or. converge_test > 3)goto 3321

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
		   !
           case('tp','TP','Tp')
		      lanczchar='tp'
              write(autoinputfile,'(a)')lanczchar(1:2)
              print*,' Thick-restart Lanczos targeting specific energy using power filter '
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
			  print*,' Enter power k for (H-E0)^k'
			  read*,kpow
			  write(autoinputfile,*)kpow
   			  print*,' Enter your target energy E0 '
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
           case('tb','TB','Tb')
		      lanczchar='tb'
              write(autoinputfile,'(a)')lanczchar(1:2)
              write(6,*)' Enter nkeep, dimension of blocks, max # of block iterations before restarting '
              write(6,*)' (nkeep = # of states printed out; typically ~ dim block )'
              maxblockiter = 1
              read(5,*) nkeep,dimblock,maxblockiter
              if(nkeep > dimbasis)nkeep = int(dimbasis, 4)
              write(autoinputfile,*)nkeep,dimblock,maxblockiter,'     ! # states to keep, dim of block,# block iters '		
			  if(menu_char=='np')then
				  print*,' You have read in pivot wave functions (recommended) '
				  piv_flag='y'
			  else
				  print*,' You did not read in block of pivot wave functions  '
				  print*,' Starting from random pivots (not recommended) '	
				  print*,' (to read in pivots, chose "np" in starting menu )'  
				  piv_flag='n'
			  end if
              niter = maxblockiter*dimblock + dimblock ! We can leave br_histmax alone
              print*,' Keeping ',nkeep,' vectors,', niter, 'iterations before restarting '

              max_iter_restart = bmax((nkeep+ncheck_ex_def)*thick_restart_factor_def,thick_restart_min_def)
              if(niter < max_iter_restart)then
                  niter = max_iter_restart
              end if
              print*,"enter number of vectors kept after thick restart (nthick)"
              read(5,*)nthick_keep

              ncheck_ex=2 ! number of additional vectors to check for convergence testing
              if(nthick_keep < (dimblock + ncheck_ex)) nthick_keep = dimblock + ncheck_ex
              print*,"enter max number of restarts"
              read(5,*)max_num_thick_restart
           write(autoinputfile,*)max_num_thick_restart
           print*,' Keeping ',nkeep,' eigenstates with ',niter,' iterations '
           print*,' with a maximum of ',max_num_thick_restart,' restarts '
           
              block_flag = 1
              converge_test = 0
              converge_tol  = converge_tol0_def
              thick_restart = .true.
              fixed_iter = .false.
           case default
              print*,' Not one of the options '
              goto 11

        end select

     end if
  end if

#ifdef _MPI
  call BMPI_BCAST(leave,1,0,MPI_COMM_WORLD,ierr)
#endif
543 continue
  if(leave)then
#ifdef _MPI
     call BMPI_ABORT(MPI_COMM_WORLD,120,ierr)
#endif
     stop
  end if
#ifdef _MPI
  call BMPI_BCAST(lanczchar,2,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(niter,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(nkeep,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(maxiter,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(maxblockiter,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(dimblock,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(block_flag,1,0,MPI_COMM_WORLD,ierr)

  call BMPI_BCAST(thick_restart,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(nthick_keep,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(max_num_thick_restart,1,0,MPI_COMM_WORLD,ierr)

  call BMPI_BCAST(piv_flag,1,0,MPI_COMM_WORLD,ierr)

  call BMPI_BCAST(converge_test,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(ncheck_ex,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(converge_tol,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(Etarget,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(targetX,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(kpow,1,0,MPI_COMM_WORLD,ierr)
#endif
  nthick_add = nthick_keep - nkeep     ! how may to add on

  if((nkeep==0 .and. niter==0) .or. lanczchar=='sk')then  ! DO NOT TRY TO DIAGONALIZE
       close(autoinputfile)
       return       
  end if

  if(lanczchar=='lf' .or. lanczchar=='tf' .or. lanczchar=='li' .or. writejumpstart .or. greenflag)then
     fixed_iter = .true.
  else
     fixed_iter = .false.
  end if

  if(lanczchar /= 'sk' .and. emptyHam .and. .not. writejumpstart)then
     if(iproc==0)then
         print*,' '
         print*,' ** WARNING WARNING NO INTERACTION READ IN WARNING WARNING ** '
		 print*,' *** YOU ARE LIKELY TO GET AN ERROR  ********  '
         print*,' '
         write(logfile,*)' ** WARNING WARNING NO INTERACTION READ IN WARNING WARNING ** '
		 write(logfile,*)' *** YOU ARE LIKELY TO GET AN ERROR  ********  '
     end if
  end if
  
if(blockflag .and. nproc > 1 .and. dimblock*maxfragsize_all > 2e9)then
	if(iproc==0)then
		write(6,*)' Warning! Your block dimension is too large to work with MPI '
		write(6,*)' Max block dimeension is ',2e9/maxfragsize_all,'; if you need more '
		write(6,*)' decrease the fragment size and possibly add more MPI processes '
		write(logfile,*)' Warning! Your block dimension is too large to work with MPI '
		write(logfile,*)' Max block dimeension is ',2e9/maxfragsize_all,'; if you need more '
		write(logfile,*)' decrease the fragment size and possibly add more MPI processes '
		write(resultfile,*)' Warning! Your block dimension is too large to work with MPI '
		write(resultfile,*)' Max block dimeension is ',2e9/maxfragsize_all,'; if you need more '
		write(resultfile,*)' decrease the fragment size and possibly add more MPI processes '
	end if
end if
  

!............. OPTIONAL STORAGE OF VECTORS IN RAM ..............

  if(enablelanczosincore1 .and. nproc == 1 .and. & 
     (useNewReorthog .or. dimbasis*(niter)*lanc_prec/2*1.0e-9 < maxlanczosstorage1)) then
          storelanczosincore1 = .true.
!          if(.not. useNewReorthog)then
!             allocate( lvec(dimbasis, niter), stat=aerr )
!             if(aerr /= 0) call memerror("lanczos_menu 1");
!          end if
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
!........ WRITE OUT STORAGE FOR BLOCK LANCZOS  
  if(dimblock > 1)then
	  
      if(nfragments==1)then
         maxfragment=dimbasis
      else
         maxfragment=0
         do i = 1,nfragments
            maxfragment = bmax(maxfragment,fragmentlist(i)%localdim)
         end do
 
      end if
	 if(iproc==0)then 
	  write(6,*)' '
	  write(6,*)' Storage of two Lanczos blocks requires ',2*maxfragment*lanc_prec*1.0e-9 *dimblock,' Gb' 
	  write(6,*)' '
	  write(logfile,*)' Storage of two Lanczos blocks requires ',2*maxfragment*lanc_prec*1.0e-8 *dimblock,' Gb ' 
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
      if (useHZSomp ) then
         if (iproc == 0 .and. verbose_ompham) print *, "MPI Time (before getLoadDist) : ", get_Wtime() - mpiStartTime
         call getLoadDist()
      end if
      if (iproc == 0 .and. verbose_ompham) print *, "MPI Time (before lanczos_p) : ", get_Wtime() - mpiStartTime
	  if(block_flag==0)then
         call lanczos_p
	  else
		  call setup_mpi_block_reg
		 call lanczos_block
	  end if
      if (iproc == 0 .and. verbose_ompham) print *, "MPI Time (after lanczos_p) : ", get_Wtime() - mpiStartTime
  end if
  
  call writelogfile('lan')
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
! SUBROUTINES CALLED:
!   initialize_lanczos_vector
!   br_grab_vec1
!   br_normalize
!   br_restore_vec1
!   br_add2hist
!   print_vec_head
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
  use localblocks
  use flagger
  use fragments
  use mod_reorthog
  use mod_balg! block alg
  use bmpi_mod
  use butil_mod
  use bvectorlib_mod
  use lanczos_util
  use apply_ham
  use apply_ham_omp
  use apply_ham_block
  use convergence
  use jumpstart,only:writejumpstart
  implicit none
  integer(4) :: ierr
  integer(4) :: mrk
  integer(4) :: iter
  real(kind=8) :: da,db,dbb
  integer(4) :: i,j,k
  real(kind=8) :: dnorm, dnorm0
  real(4) :: xj,xt

  logical smallflag,zeroflag
  
  integer :: i_thick_restart
  integer(4) :: i_thick

  integer :: nread,jvec

  integer inode
  character(1) :: vchar
  integer :: aerr
  
  real :: smallbetatol = 0.0001 ! tolerance for a small beta
  logical :: loss
  integer :: ipow
  
!..... modified skipping of diagonalization, added 7.10.7
  
  integer :: itskipeig
  integer :: nprint_e
  
  itskipeig = itskipeig_def
!  if(large_skip_enabled)itskipeig = max(itskipeig,int(frac_large_skip*niter))

  if(large_skip_enabled)itskipeig = max(itskipeig,int(sqrt(float(nkeep))))

  
  if(itskipeig > itskipeig_def .and. iproc==0)then
	  print*,' Diagonalizing only every ',itskipeig,' iterations '
	  write(logfile,*)' Diagonalizing only every ',itskipeig,' iterations '
  end if

  
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
  use_restricted_esolver=.false. ! default
  
  if ( .not. fixed_iter ) then
	 if(targetX)then
	   ncheck = int(bmin(dimbasis, nkeep), 4) 
	 else
       ncheck = int(bmin(dimbasis, nkeep+ncheck_ex), 4) ! check extra states to prevent diving
	   
! added 7.11.2 to speed up convergence for a large # of target levels
	   
	   if(use_test_converge .and. .not.thick_restart)then
		   if(nkeep > nswitch_esolver)then
			   use_restricted_esolver = .true.
			   ntest_converge = int(sqrt(float(ncheck))+1)
		   else
			   use_restricted_esolver=.false.
			   ntest_converge = ncheck
		   end if
	   else
		   use_restricted_esolver=.false.
		   ntest_converge = ncheck
	   end if
	   
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
  if(lanczchar=='tx')then       ! ADDED 7.9.0
	  allocate(lanczos_mat(niter,niter),stat=aerr)
	  allocate(etx(niter),eiglvectx(niter,niter))
      if(aerr /= 0) then
         call memerror("lanczos_p 3.5")
         stop 55
      end if
	  lanczos_mat=0.d0
  end if
  if(.not.allocated(eiglvecold)) allocate( eiglvecold(niter,niter), stat=aerr )
  if(aerr /= 0) then
     call memerror("lanczos_p 4")
     stop 5
  end if
   
  call clocker('lan','sta')
  finished = .false.
  iter = 0
  i_thick_restart = 0
  i_thick = 0

  vchar='n'
!........ INITIAL NORMALIZATION.......

  call br_grab_vec1()          ! push vec1 -> br_reg
  call br_normalize(dnorm0)
  call br_restore_vec1()    !    pull vec1 <- br_reg
  call br_add2hist(1)       !   adds br_reg to history stack

  if(iproc==0 .and. (strengthflag .or. greenflag))then 
      print*,dnorm0*dnorm0,'   = total input strength '
     if(writeout)write(resultfile,*)dnorm0*dnorm0,' = total strength '
  endif  

  call   procOP_clock(0,'set','all')

! START LANCZOS ITERATIONS......................................
  actual_iterations = 0
#ifdef _MPI
  call BMPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
!  if(iproc == 0) then
!     print *, "Starting Lanczos with nfragments=", nfragments  ! This nfragments reference is ok
!  end if

  do while ( .not. finished )
     iter = iter + 1
     ! zeros output vector and thread specific output vectors if enabled
     call initialize_final(vchar)

     if (useHZSomp .and. num_threads_global>1) then
           call applyHbundled_omp(vchar)    ! THE CENTRAL ENGINE: APPLY MAT-VEC MULTIPLY
     else
            call applyHbundled_g(vchar)    ! THE CENTRAL ENGINE: APPLY MAT-VEC MULTIPLY
     end if
	 
	 if(kpow > 1)then
		 do ipow = 2,kpow
			 call br_grab_vec2()
			 call br_restore_vec1() ! this puts vec2 back into vec1
			 
		     call initialize_final(vchar)
			 
		     if (useHZSomp .and. num_threads_global>1) then
		           call applyHbundled_omp(vchar)    ! THE CENTRAL ENGINE: APPLY MAT-VEC MULTIPLY
		     else
		            call applyHbundled_g(vchar)    ! THE CENTRAL ENGINE: APPLY MAT-VEC MULTIPLY
		     end if			 			 
		 end do
	 end if

     actual_iterations = actual_iterations + 1

!.......... CAN COMPUTE ALPHA(ITER) EITHER BEFORE OR AFTER REORTHOGONALIZATION
!           while for exact arithmetic does not make a difference, it can 
!           be different for finite arithmetic.  Tests so far suggest 
!           converged energies agree to less than 1 keV and converged 
!           wavefunctions agree in overlap to less than 10^4 but be wary
!
!           USE DEFAULT-- ALPHA BEFORE REORTHOGONALIZTION
!.................... STANDARD REORTHOGONALIZATION from KSM...................


		if(iproc==0)call clocker('ort','sta')
        call br_grab_vec2()   ! push vec2 -> br_reg
        call br_remove_prev_overlap(da)
        call br_orthogonalize(1) ! 1 says ignore last history entry (prev)

        call br_normalize(db)
		
!....... CHECK FOR LOST OF ORTHOGONALITY.......SIGNALLED BY SMALL db....
!   added in 7.9.4
!   Note: very tiny db generally requires a random restart
!        however if db is just "smaller" than previous beta, but not tiny
!        check for loss of orthogonality
!

        if(iter > 1 .and. lanczchar /='tx')then
			if(db < beta(iter-1)*smallbetatol)then
				if(iproc==0)then
					print*,' Random restart on iteration ',iter
					write(logfile,*)' Random restart on iteration ',iter
					
				end if
				call random_restart_with_history
				go to 4334

			end if
			if(confirm_orthogonality .or. db < 0.1*beta(iter-1))then
				if(iproc==0)then
					print*,' Testing for loss of orthogonality on iteration ',iter
					write(logfile,*)' Testing for loss of orthogonality on iteration ',iter
				end if
				call br_check_orthogonal(loss)
			    if(loss)then
				   if(iproc ==0)then
					   print*, ' loss of orthogonality  ',da,db,beta(iter-1)
					   write(logfile,*)' loss of orthogonality  ',da,db,beta(iter-1)
				   end if
		            call br_orthogonalize(0)
				end if
			end if

		end if
	
4334 continue
		
        ! Restore back to vec1 for next iteration
        ! This is instead of swap_vchar in older code
        call br_restore_vec1()  ! vec1 <-- br_reg
        alpha(iter) = da !real(da,kind(0.0e0))
        beta(iter) = db !real(db, kind(0.0e0))
		if(lanczchar=='tx')then
			lanczos_mat(iter,iter)=alpha(iter)
			lanczos_mat(iter,iter+1)=beta(iter)
			lanczos_mat(iter+1,iter)=beta(iter)		
			do i = 1,iter-1
				lanczos_mat(iter,i)= br_orthogdotf(i)
				lanczos_mat(i,iter)= br_orthogdotf(i)
				
			end do	
		end if
        ! Finally, add br_reg to history
        if(iter < niter) call br_add2hist(iter+1)
		if(iproc==0)call clocker('ort','end')
		
!............ Write out Lanczos coefficients..............................

     if ( iproc ==0 ) write(coef_file,'(i10,2(2x,f20.15))')iter,alpha(iter), beta(iter)
     if ( iter == dimbasis ) exit

!---------- compute and print intermediate eigenvalues.....................
! V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V 
     flyeig = (iter ==(iter/itskipeig)*itskipeig)
     if ( iproc==0 .and. .not. coef_only .and.  & 
	      (( iter > nkeep .and. flyeig .and. checkeig) .or.      & 
          (iter >= ncheck .and. .not.fixed_iter .and. .not. skipeig .and. checkeig)  .or. & 
		     ( thick_restart .and. .not. finished .and. iter+1 == niter )).and..not.greenflag)then

        if(noisy0) print *, "Finding lanczos eigenvalues"
        call find_lanczos_eigenvalues(iter,altconverge,  i_thick)

        if ( flyeig .or. thick_restart .and. iproc == 0 ) then       ! printout
           write(6,*)' '
		   if(ntest_converge==ncheck)then
			   nprint_e = nkeep
		   else
			   nprint_e = ntest_converge
		   end if
           do j = 1,nprint_e
!           do j = 1, nkeep
			   if(lanczchar=='tx'.and.txchoice==2)then
                   write(6,'(i5,f20.5)')j,etx(j)
			   else
                  write(6,'(i5,f20.5)')j,e(j)
			  end if
           end do
           write(6,*)' '
           write(6,*)iter,' iterations ' ! ' eigenvalues '
           if ( i_thick > 0 .and. iproc == 0 ) then
              write(6,*)' ( total = ',iter + niter*i_thick_restart,' ) '
           end if
        end if
        flush(6)
     end if
	 
	 if(flyeig .and. coef_only .and. iproc==0)write(6,*)iter,' iterations ' ! ' eigenvalues '
	 
! ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^
!----------------- end of compute and print intermediate eigenvalues..................

!---------------- CONVERGENCE CONTROL.........................
!--------------- CHECK FOR INITIALIZTION ----------------
!                IF FINISHED, THEN USE THICK_RESTART ROUTINE TO
!                CREATE INITIAL VECTOR(S)
	 
	  call convergence_test(iter,i_thick_restart,finished)
 
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

!----------------------- LANCZOS ITERATIONS FINISHED ---------
  call clocker('lan','end')

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

  call all_clocks_out  
  call timeperopmaster('out')  ! MUST BE LAST
    if(reorthog .and. .not.coef_only .and. .not. writejumpstart)then
     call  lanczos_output(iter,i_thick,dnorm0)
  end if
  
  if(coef_only .and. iproc==0)then
	  print*,' '
	  print*,' Lanczos coefficients only '
	  print*,' '
  end if
  if(menu_char=='pv')then
	  call overlap_lanczos
  end if
  return
end subroutine lanczos_p

!=================================================================
!
!  Block Lanczos routines primarily by RZ 2018, added in 7.8.6-7.8.7
!  separate routine created in 7.8.8 by CWJ
!
!  Modified in 3/2022 RMZ - added HZOMP algorithm capabilities for Hvec
!  Modified in 5/2022 RMZ - added thick-restart 
!  SUBROUTINES CALLED:
!
subroutine lanczos_block
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
  use localblocks
  use flagger
  use fragments
  use mod_reorthog
  use mod_balg! block alg
  use bmpi_mod
  use butil_mod
  use bvectorlib_mod
  use lanczos_util
  use apply_ham
  use apply_ham_omp
  use applyhamblockomp
  use apply_ham_block
  use convergence
  use jumpstart,only:writejumpstart
  
  implicit none
  integer(4) :: ierr
  integer(4) :: mrk
  integer(4) :: iter
  real(kind=8) :: da,db
  integer(4) :: i,j,k
  real(kind=8) :: dnorm, dnorm0
  real(4) :: xj,xt

  logical smallflag,zeroflag
  
  integer :: i_thick_restart
  integer(4) :: i_thick

  integer :: nread,jvec

  integer inode
  character(1) :: vchar
  integer :: aerr

		
dimlanc = 0
dimlanc = niter  ! dimension of lanczos matrix, and used to control iterations.
!dd = dimbasis ! debugging purposes

if (converge_test == 2 .or. converge_test == 3)then
   altconverge = .true.
else
  altconverge = .false.
end if

if ( .not. fixed_iter ) then
     ncheck = int(bmin(dimbasis, nkeep+ncheck_ex), 4) ! check extra states to prevent diving
   if(.not.allocated(eold)) allocate( eold(ncheck), stat=aerr )
   if(aerr /= 0) call memerror("lanczos_p 1")
   ntest_converge=ncheck
else
   ncheck = 0
end if

! Allocate data structures
!-----------------------------------------------------------------
  if(.not.allocated(alpha_block)) allocate( alpha_block(dimblock,dimblock),beta_block(dimblock,dimblock), stat=aerr )
  alpha_block = 0.d0
  beta_block  = 0.d0 
  
  if(aerr /= 0) then
     call memerror("lanczos_p 2")
     stop 5
  end if
  if(.not.allocated(ev)) allocate( lanczos_mat(dimlanc,dimlanc),ev(dimblock), stat=aerr )
  
  !... THESE ARE ALLOCATED IN TWO PLACES... NEED TO FIX
  if(.not.allocated(e)) allocate( eiglvec(niter,niter),e(niter), stat=aerr )
  if(aerr /= 0) then
     call memerror("lanczos_p 5")
     stop 5
  end if
  if(.not.allocated(eiglvecold)) allocate( eiglvecold(niter,niter), stat=aerr )
  if(aerr /= 0) then
     call memerror("lanczos_p 6")
     stop 5
  end if  
  
  lanczos_mat = 0.d0
  ev = 0.d0
  if(aerr /= 0) then
     call memerror("lanczos_p 3")
     stop 5
  end if

  !.... IMPORTANT NEED TO ALLOCATE CORRECTLY)

!    allocate( block_1(dimblock, dimbasis ), stat=aerr )
    allocate( blockvec_1(1+dimblock*(v1s-1):dimblock*v1e), stat=aerr )
	
!	print*,iproc,' blockvec 1 ',frag1,int(1+dimblock*(v1s-1),kind=4),int(dimblock*v1e,kind=4)

    if(aerr /= 0) call memerror("setup_block1")
!    allocate( block_2(dimblock, dimbasis ), stat=aerr )
    allocate( blockvec_2(1+dimblock*(v2s-1):dimblock*v2e), stat=aerr )
	
    if(aerr /= 0) call memerror("setup_block2")
!    block_1 = 0
!    block_2 = 0
blockvec_1 = 0.0
blockvec_2 = 0.0
  !if(.not.allocated(blk_s)) allocate( eiglvecold(niter,niter), stat=aerr )
 !if(aerr /= 0) then
  !   call memerror("lanczos_p 4")
  !   stop 5
  !end if

 if(.not.allocated(blk_s)) allocate( blk_s(dimlanc),blk_e(dimlanc), stat=aerr )
 if(aerr /= 0) then
    call memerror("lanczos_p 4")
    stop 5
 end if
 blk_s = 0
 blk_e = 0
   call clocker('lan','sta')
  finished = .false.
  iter = 0
  i_thick_restart = 0
  i_thick = 0
  vchar = 'n'	
!SETUP HISTORY BLOCK BOOKKEEPING

call setup_hist_index(dimlanc) !set up labeling vectors for start/end of blocks in hist.

if(iproc==0)print*,' dimlanc = ',dimlanc

!........ BLOCK LANCZOS INITIAL NORMALIZATION.......
if(menu_char .ne.'np' .and. .not.blockstrengthflag) then
    call initialize_pivot_block ! random pivot
else
	allocate(orthonorm_pivot(dimblock,dimblock))
	if(iproc==0)print*,' Orthonormalizing pivot block '
	call spectral_decomp(0)
	orthonorm_pivot = beta_block
	beta_block=0.d0
endif

 call   procOP_clock(0,'set','all')

! START LANCZOS ITERATIONS......................................
  if(iproc==0)print*, "starting lanczos main loop"
  actual_iterations = 0
  
  iter = 0

  if(iproc == 0) then
  end if
!  if(wantnoisy0) call BMPI_BARRIER(MPI_COMM_WORLD, ierr)
   !change main loop control structure for iter to be updated by thick-restart
  !do iter = 1,  niter/dimblock-1 ! dimlanc/dimblock !- 1 
    do while ( .not. finished )
     iter = iter + 1
	  alpha_block = 0.
	  beta_block  = 0.  	  
     if(noisy0) then
        print *, "Starting Block Lanczos iteration ", iter
        flush(6)
     end if
!     if(wantnoisy0) call BMPI_BARRIER(MPI_COMM_WORLD, ierr)

     if(noisy0) then
        print *, "initialize_final"
        flush(6)
     end if
!     if(wantnoisy0) call BMPI_BARRIER(MPI_COMM_WORLD, ierr)
     if(noisy0) then
        print *, "clocker"
        flush(6)
     end if
	 
!     if(wantnoisy0) call BMPI_BARRIER(MPI_COMM_WORLD, ierr)
     if(noisy0) then
        print *, "Calling applyHbundled_block"
        flush(6)
     end if
!     if(wantnoisy0) call BMPI_BARRIER(MPI_COMM_WORLD, ierr)

     call br_pull_blreg_from_hist(blk_s(iter))

     call br_pull_block1_from_reg()
     if (useHZSomp .and. num_threads_global>1) then
      ! Shan's omp distribution algorithm applied to block Lanczos
         call applyHbundled_block_omp   ! THE CENTRAL ENGINE: APPLY MAT-VEC MULTIPLY
     else
      call applyHbundled_block
     endif

	 actual_iterations=actual_iterations+dimblock
     if(noisy0) then
        print *, "After Calling applyHbundled_block"
        flush(6)
     end if
!     if(wantnoisy0) call BMPI_BARRIER(MPI_COMM_WORLD, ierr)

	if(iproc==0)call clocker('ort','sta')

	 call br_push_block2_to_reg()   ! COULD BE PROBLEMATIC
	 
	 call br_push_blreg_to_hist(blk_s(iter+1))


!compute alpha
	 
     call grab_alpha_blockv2(iter)

! add alpha to Tmat
     !print*," adding alpha to lanczoz matrix"
     call add_alpha_Tmat(iter)
     !print*," adding alpha to lanczoz matrix done"

!compute W'
     if(iter .gt. 1) then
        !print*,"computing W'"
        call compute_W(iter)
         !print*,"computing W' done"
	 end if
     !print*,"computing spectral decompositon"
     call spectral_decomp(iter)
     !print*,"computing spectral decompositon done"
	 

     !print*,"adding beta to lanczos matirx"
     call add_beta_Tmat(iter)
     !print*,"adding beta to lanczos matirx done"
! endif	 
if(iproc==0)call clocker('ort','end')
 

! write alpha and beta to coeff file
     call report_alpha_beta(iter)
	 
!print itermediate eigenvalues
!     call find_lblock_eigenvalues(nkeep,.false.,iter)
if(iter*dimblock >= nkeep)then
   call find_lanczos_eigenvalues(iter*dimblock,.false.,i_thick)

if(iproc==0)then
    write(6,*)' '
    write(6,*) iter," block iterations = ",iter*dimblock,' mat vec iterations '
	
    do j = 1, nkeep
       write(6,'(i5,f20.5)')j,e(j)
    end do
	end if 

	!call convergence_test(iter*dimblock,0,finished)
   call convergence_test(iter*dimblock,i_thick_restart,finished) ! for thick restart
	 !     print*,"----------------------------------------------------"
    !  print*,lanczos_mat(1:nthick_keep+6,1:nthick_keep+6)
    !  print*,"----------------------------------------------------"

   !--------------- CHECK FOR THICK-RESTART --------------------------
   ! added V7.9.12 -RMZ SDSU 2022
   !(iter+1)*dimblock >= niter in case the block dimension is not a multiple of nkeep
     if ( thick_restart .and. .not. finished .and. (iter+1)*dimblock .ge. niter ) then           
        if ( i_thick_restart > max_num_thick_restart ) then
           if ( iproc == 0 ) write(6,*)' WARNING could not finish '
           finished = .true.
           cycle
        end if
        
        if ( iproc == 0 ) write(6,*)' Doing a thick restart ',i_thick_restart

        call block_thick_restart_sub_p(nkeep+nthick_add,iter*dimblock,i_thick_restart,niter,vchar)
        i_thick_restart = i_thick_restart + 1
        i_thick = nkeep + nthick_add   !ncheck
        !iter = nkeep+nthick_add !  ncheck
        iter = (nkeep+nthick_add)/dimblock !  ncheck
   end if






	if(finished .and. iproc==0)then
		print*,' '
		print*,' Convergence satisfied! '
		print*,' '
	end if
! if stopping via fixed iterations deals with do while loop for thick-restart.
if(lanczchar .ne. 'tb' .and. (iter+1 .gt. niter/dimblock-1))then
   if(.not. finished) then !converges at last iter.
      if(iproc==0)then
         print*,' '
         print*,' Out of block iterations, Convergence not satisfied! '
         print*,' setting finished to true'
         print*,' '
      end if
      finished = .true.
   endif
end if 
	if(finished)exit


end if

blockvec_1 = 0.
blockvec_2 = 0.
alpha_block = 0.
beta_block = 0.

  end do  ! iter


call find_lanczos_eigenvalues(iter*dimblock,.true.,i_thick)
	 
! check for final iteration
call all_clocks_out
call timeperopmaster('out')  ! MUST BE LAST

! let's record the number of iterations for thick-restart
 if ( iproc == 0 ) then
     print*,iter*dimblock,' iterations '
	 write(logfile,*)
     print*,' '
  end if
  if ( i_thick > 0 .and. iproc == 0) then
     print*,' ( total = ',iter*dimblock + niter*i_thick_restart,' ) '
  end if
  if ( writeout .and. iproc == 0 ) then
     write(resultfile,*)iter*dimblock,' iterations '
	 write(logfile,*)iter*dimblock,' iterations '
     write(resultfile,*)' '
     if ( i_thick > 0 .and. iproc == 0 ) then
        write(resultfile,*)' ( total = ',iter*dimblock + niter*i_thick_restart,' ) '
        write(logfile,*)' ( total = ',iter*dimblock + niter*i_thick_restart,' ) '

     end if
  end if

   !if(.not.writejumpstart)call  lanczos_output(iter*dimblock,0,dnorm0) i_thick
   if(.not.writejumpstart)call  lanczos_output(iter*dimblock,i_thick,dnorm0)
   return

end subroutine lanczos_block

!---- NOTE: density1b_output moved to bdenslib1.f90
!=================================================================
!
!  if in a small space, then diagonalize fully.
!  this helps to eliminate problems with Lanczos
!  the criterion is fulllimit found in module lanczos_info
!
! SUBROUTINES CALLED:
!   setup_localvectors
!  applyhbundled_g
!  ssyev/dsyev  LAPACK routines
!   wfn_write_nkeep
!   setup4obsmaster
!   reset_jumps4obs
!   write_lanczos_vector_a   SHOULD BE REPLACED
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
  real :: nparticles ! temporary

  real(kind=8) :: dnorm0,dsclrprod,da
  character(1) :: vchar
  integer jvec
  integer info
  integer iorb
  integer :: xpar
  character  :: outchoice  !used for steering output options

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
#ifdef _MPI
            call BMPI_BCAST(vec2, size(vec2), 0, fcomm2, ierr)
#endif
        else
#ifdef _MPI
            call BMPI_BCAST(vec1, size(vec1), 0, fcomm1, ierr)
#endif
         end if
     end if
     do j = 1,dimbasis
        h(i,j) = real(vec2(j),kind(egv_prec))   ! store as real(4); this can be changed

     end do  !J 
     if(dimbasis < 400 .and. iproc==0 .and. .not. printouthamflag)write(60,1)(h(i,j),j=1,dimbasis)
	 if(printouthamflag)then
		 do j = i,dimbasis
			 if(h(i,j)/=0.0)write(60,*)i,j,h(i,j)
		 end do
	 end if
1  format(9f13.4)
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
  
  if(.not.skip_T2)then
     call setup4obsmaster('T')
  end if
  if ( iproc == 0 ) print*,' '

  call reset_jumps4obs

!...................................................................
!----------- COMPUTE SINGLE PARTICLE OCCUPATIONS---------- added in 7.3.7
if(spoccflag .and. .not.densityflag)then
  call pocc_init_spoccx() ! allocate and initialize array
  call pocc_write_orbits()
end if


!  move later to be consistent
!if(densityflag) then
!	if(iproc==0)write(denresultfile,*)' '
!	 call pocc_write_orbits_alt(denresultfile)
! end if
  
  if ( iproc == 0 ) print*,' '

88 format(2f10.5)
     if ( iproc == 0 ) then
		 if(print_parity)then
        write(6,*)' State      E        Ex         J       T    par'
        if(writeout)write(resultfile,*)' State      E        Ex         J      T     par'
		if(densityflag)then
			write(denresultfile,*)' State      E        Ex         J      T     par'
			write(occresultfile,*)' State      E        Ex         J      T     par'
			
		end if
		
	    else
        write(6,*)' State      E        Ex         J       T      par'
        if(writeout)write(resultfile,*)' State      E        Ex         J      T      par'
		if(densityflag)then
			write(denresultfile,*)' State      E        Ex         J      T      par'
			write(occresultfile,*)' State      E        Ex         J      T      par'
		end if
			
	   endif
     end if
  
! Construct nkeep eigenvectors of H from  Lanczos eigenvectors....

  if(densityflag .or. trdensout)then  ! must store in lvec
     if(allocated(lvec))deallocate(lvec)
     allocate( lvec(dimbasis, nkeep), stat=aerr )
     if(aerr /= 0) call memerror("exactdiag_p 4")
     allocate(energy(nkeep), xjlist(nkeep),xtlist(nkeep), stat=aerr )
     if(aerr /= 0) call memerror("exactdiag_p 5")
  end if
  allocate(parity(nkeep))
  do vidx = 1, nkeep

        do k = 1, dimbasis
           vec1(k) = real(h(k,vidx), kind=lanc_prec)  ! u(k,vidx)
        end do !k
        call br_grab_vec1()
        call br_load_vec2_from_vec1()
		call br_add2hist(vidx)
!... THIS SHOULD BE REWRITTEN TO REMOVE OBSOLETE ROUTINE...
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
	 if(print_parity .or.densityflag)call find_parity(xpar,1)
	 
     if(densityflag .or. trdensout)then
        energy(vidx) = real(e(vidx), kind(energy(1)))
        xjlist(vidx) = xj
        xtlist(vidx) = xt
		parity(vidx) = xpar
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
	                  call pocc_write_ejt(vidx, e, xj, xt,xpar)

	                case('b')
	                  ! write(6,12)vidx,e(vidx), e(vidx) - e(1),xj,xt,(spoccx(1,vidx,iorb),iorb=1,numorb(1))
	                  call pocc_write_ejt(vidx, e, xj, xt,xpar)
	                  call pocc_write_occvec(6, spoccx, vidx, 1, "    p occ:")
	                  call pocc_write_occvec(6, spoccx, vidx, 2, "    n occ:")
	                  if ( writeout ) then
	                     call pocc_write_occvec(resultfile, spoccx, vidx, 1, "    p occ:")
	                     call pocc_write_occvec(resultfile, spoccx, vidx, 2, "    n occ:")
	                  end if
	                case('p')
	                  call pocc_write_ejt(vidx, e, xj, xt,xpar)
	                  call pocc_write_occvec(6, spoccx, vidx, 1, "    p occ:")

	                  if ( writeout ) then
	                     call pocc_write_occvec(resultfile, spoccx, vidx, 1, "    p occ:")
	                  end if

	                case('n')
	                  call pocc_write_ejt(vidx, e, xj, xt,xpar)
	                  call pocc_write_occvec(6, spoccx, vidx, 2, "    n occ:")

	                  if ( writeout ) then
	                     call pocc_write_ejt(vidx, e, xj, xt,xpar)
	                     call pocc_write_occvec(resultfile, spoccx, vidx, 2, "    n occ:")
	                  end if
	             end select
	         end if	 
	 
  end do
  
  if(densityflag) then
  	if(iproc==0)then
		write(denresultfile,*)' '
		write(occresultfile,*)' '
	end if
  	 call pocc_write_orbits_alt(denresultfile)
  	 call pocc_write_orbits_alt(occresultfile)
	 
   end if

11 format(i5,3x,2f10.5,2x,2f8.3)

  if(densityflag)call density1b_output !(e,eiglvec)

  if ( trdensout  ) then
     call output_TRDENS(.false.)
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
!   write_lanczos_vector_a
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
  integer :: xpar

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
!  print*,iproc,' hvec ',hvec
#ifdef _MPI
  call BMPI_ALLREDUCE(hvec,size(hvec),MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
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
	  if(dimbasis < 100 .and. iproc==0)print*,' Wrote out Hamiltonian '
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
  if(.not.skip_T2)then
     call setup4obsmaster('T')
  end if
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
     call pocc_write_ejt(ki, e, xj, xt,xpar)

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
     call output_TRDENS(.false.)
  end if
  call close_lanczosfile

  call clocker('obs','end')
end subroutine exactdiag_MPI

