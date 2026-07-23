! file BDENSITYLIB5.f90
!
! routines for two-body densities, added 7.9.1   Feb 2019 by CWJ
!
! Two-body densities always generated from existing wfn files
!
! Steps for two-body densities:
!  -- Read in wfn header and create basis (done in main)
!  -- set up two-body jump arrays; worry about hermiticity
!        --> need to adjust indexing
!  -- set up uncoupled density arrays <--> similar to uncoupled two-body matrix element arrays
!        however must allow for non-Hermitian values
!        (by which I mean  the value < a^+ b^+ d c >  not same as  < c^+ d^+ b a >  
!         and must be stored independently )
!  -- set up coupled density arrays; also must allow for non-Hermitian structure
!  
!    Loop over wave functions:
!  -- compute uncoupled densities (routines in bdenslib4.f90)
!  -- couple densities and write to file




!........... BOSS ROUTINE FOR 2-BODY DENSITIES......
!
!=========================================================================
!
! subroutine BOSS_TWOBODY_DENSITIES
!
! based upon subroutine DENSITY1B_FROM_OLDWFN
!
! computes density from pre-existing wavefunction file
!
!  CALLED BY
!       main routine (bigstick_main.f90)
!  SUBROUTINES CALLED:
!       setup_localvectors
!       set1bsectorjumps
!       masterfindconjumps1b
!       master1bodyjumps
!       master_op_stat()
!       defaultopbundles()
!       setnodaltribution
!       density_bundles_setup()
!       setMPIlimits4densities()
!       wfn_read_nkeep
!       pocc_write_orbits()
!       wfn_readeigenvec
!       master_density1b_bundled
!       coupled_densities
!       dm_add_tab              ! ADDED BY KSM 7.4.3
!       density1b_write_math    ! ADDED BY KSM 7.4.3
!       call dm_reset()         ! ADDED BY KSM 7.4.3
!
!  FUNCTIONS CALLED
!       closest2J
!

subroutine boss_twobody_densities
	use applytwobodydens
	use twobodydens_util
	
  use menu_choices
  use system_parameters
  use sporbit
  use spstate
  use localvectors
  use nodeinfo
  use io
  use basis
  use densities
  use haiku_info
  use precisions
  use fragments
  use mod_reorthog
  use wfn_mod 
  use pocc_mod
  use menu_choices
  use dm_mod
  use para_main_mod
  use hoppy
  use jump_mod
  use bvectorlib_mod
  use lanczos_util
  use densities
  use interaction
  implicit none

  real(4) :: xj,xt,ei,ef,xtt,xjj
  real(4) :: xt2
!  integer :: nkeep
  integer(4) :: i,j,m,n
  integer(4) :: ji,ti,jf,tf,jt,tt,jmin,jmax
  integer(kind=basis_prec) :: k
  logical :: evenAJ,evenAT    ! needed to force "correct" J, T
   
  logical :: zeroflag
  
  real, allocatable :: denmat(:,:,:)
  real :: nparticles
  logical :: numberflag  ! flag to check on total number of particles; used for debugging

  real(kind=4), allocatable :: stateE(:), stateJ(:), stateT(:)
  real(kind=4) :: me ! matrix element
  
  integer :: istatestart,istatestop,fstatestart,fstatestop  ! options to start and stop
  
  integer :: ierr

!  if(iproc==0)print*,' WARNING SUBROUTINE density1b_from_oldwfn NEEDS TO BE VALIDATED '

  numberflag = .true.

!........ DETERMINE IF EVEN OR ODD SYSTEM............

  if( mod(np(1)+np(2),2) == 1 )then

! if "spinless" then J is always integer (because it is really L)
! however T can be half-integer (if "spinless" then this is really S)
     if(spinless)then
       evenAJ = .true.
     else
       evenAJ = .false.
     end if
     evenAT = .false.
  else
     evenAJ = .true.
     evenAT = .true.
  end if

!...................

  if ( iproc == 0 ) then
     print*,' '
     print*,' Computing two-body density matrices '
	 if(diag_den2b)then
		 print*,' Only computing for initial=final, Jt = 0 '
	 end if
     print*,' '
  end if
  
!  print*,' BEFORE SETTING UP LOCAL VECTORS, NEED TO CREATE JUMPS '
call setup_localvectors  
!...................LOOP OVER WFNS..............................
  call wfn_read_nkeep(oldwfnfile, nkeep)

  allocate(stateE(1:nkeep), stateJ(1:nkeep), stateT(1:nkeep), stat=aerr)
  if(aerr /= 0) then
     call memerror(" boss_twobody_densities 10")
     stop 5
  end if
  if(iproc==0)then
     print*,' '
    print*,nkeep,' states in file '
	print*,' '
	
	if(auto_input)then
		read(autoinputfile,*)istatestart,istatestop
		if(diag_den2b)then
 		   fstatestart=istatestart
 		   fstatestop =istatestop
		else
			 read(autoinputfile,*)fstatestart,fstatestop
		 endif
	else
	   if(diag_den2b)then
		   print*,' Enter start, stop for states '
	   else
		   print*,' Enter start, stop for initial states '
	   end if
		   	
	   print*,' (This is because two-body densities are large)'
	   print*,' (Enter 0,0  to read all )'
	   read*,istatestart,istatestop
	   
	   if(diag_den2b)then
		   fstatestart=istatestart
		   fstatestop =istatestop
	   else

   	       print*,' Enter start, stop for final states '
	       print*,' (Enter 0,0  to read all )'	
	       read*,fstatestart,fstatestop
	   end if
	
    end if
	
	if(istatestart <=0 .or. istatestop<=0)then
		istatestart =1
		istatestop = nkeep
	end if
	istatestop = min(istatestop,nkeep)
	if(istatestop < istatestart)istatestart = istatestop
	if(fstatestart <=0 .or. fstatestop<=0)then
		fstatestart =1
		fstatestop = nkeep
	end if
	fstatestop = min(fstatestop,nkeep)
	if(fstatestop < fstatestart)fstatestart = fstatestop
	
	if(auto_input)then
		if(diag_den2b)then
			print*,' Densities for:  states from, to ',istatestart,istatestop
		else
		     print*,' Densities for: initial states from, to ',istatestart,istatestop
		     print*,' Densities for: final states from, to ',fstatestart,fstatestop
		 end if
	else
		if(diag_den2b)then
 		   write(autoinputfile,*)istatestart,istatestop," ! start/stop for states for 2b dens "
			
		else
			
		   write(autoinputfile,*)istatestart,istatestop," ! start/stop for initial states for 2b dens "
	   	   write(autoinputfile,*)fstatestart,fstatestop," ! start/stop for final states for 2b dens "
	   endif
	end if
	
  end if
  
#ifdef _MPI
  call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! call OLDMPI_BCAST(ychar,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(istatestart,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(istatestop,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(fstatestart,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(fstatestop,1,0,MPI_COMM_WORLD,ierr)
#endif
  ! allocate storage to store density matrices
    
  allocate(dmatpp(nmatXX(1)),dmatpphc(nmatXX(1)), stat=aerr)
  if(aerr /= 0) call memerror("allocating dmatpp in boss_twobody_densities ")

  allocate(dmatnn(nmatXX(2)),dmatnnhc(nmatXX(2)), stat=aerr)
  if(aerr /= 0) call memerror("allocating dmatnn in boss_twobody_densities ")

  call delayedPNmatrixelements  
  call count_uncouplepn('2',.true.)  ! BUT WILL NOT ACTUALLY UNCOUPLE MATRIX ELEMENTS 
    
  call bundle_clock(0,'set')
  call   procOP_clock(0,'set','all')
  call setup_PNarrays
  call setup_2bdens_output
  
  do i = istatestart,istatestop ! 1,nkeep
     ! new interface, we say which vec to read.  It checks
     call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,ei,xj,xtt)

     ji = closest2J(evenAJ,xj)
     ! if !isoflag, then what?  orig code didn't calculate ti
     ti = closest2J(evenAT,-0.5 + sqrt(xtt+0.25))
     stateE(i) = ei
     stateJ(i) = ji
     stateT(i) = ti

     if(diag_den2b)then
	    fstatestart=i
	    fstatestop=i
	 end if
     do j = fstatestart,fstatestop ! 1,nkeep		 
           ! new interface, we say which vec to read.  It seeks and checks
           call wfn_readeigenvec(oldwfnfile,frag2, fcomm2_index, vec2, j, ef,xjj,xtt)

           jf = closest2J(evenAJ,xjj) ! nint(2*xjj)          
           tf = closest2J(evenAT,-0.5 + sqrt(xtt+0.25)) ! nint(-1 + sqrt(4*xtt+1))
           jmax = (jf + ji)/2
           jmin = abs( jf -ji)/2
		   
		   if(diag_den2b)then
			   jmin=0
			   jmax=0
		   else
	           jmax = (jf + ji)/2
	           jmin = abs( jf -ji)/2
		   end if

333        format(' Initial state # ',i4,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
334        format(' Final state   # ',i4,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
433        format(' Initial state # ',i4,' E = ',f10.5,' 2xJ   = ',i4) 
434        format(' Final state   # ',i4,' E = ',f10.5,' 2xJ   = ',i4) 

             if(np(1)>1)dmatpp= 0.0
             if(np(1)>1)dmatpphc= 0.0
			 
			 if(np(2)>1)dmatnn=0.0
			 if(np(2)>1)dmatnnhc=0.0
			 
			 if(np(1)*np(2)>0)dmatpn =0.0
			 if(np(1)*np(2)>0)dmatpnhc =0.0
			 
             call boss_den2b_bundled
			 
!...................... NOW COUPLE DENSITY MATRICES.............
    
             if(np(1)> 1)then
				 call  couple_2bdensXX(1,Ji,Jf,.true.)  ! zero out			 
				 call  couple_2bdensXX(1,Ji,Jf,.false.)
			 end if			 
             if(np(2)> 1)then
				 call  couple_2bdensXX(2,Ji,Jf,.true.) ! zero out			 
				 call  couple_2bdensXX(2,Ji,Jf,.false.)
			 endif
			 
             if(np(2)*np(1) > 0)then
				 call  couple_2bdensPN(Ji,Jf,.true.) ! zero out			 
				 call  couple_2bdensPN(Ji,Jf,.false.)
			 endif
			 if(iproc==0)then
                 call print_out_2bdens(i,j,Ei,Ef,0.5*Ji,0.5*Jf,0.5*Ti,0.5*Tf)
			 end if

     enddo   !j
  enddo  ! i
  if(iproc == 0)then
     write(dens2bfile,*)'!# +++++++++++++++++++++++++++++++++++++++++++++'
     write(dens2bfile,*)' '
	 write(dens2bfile,*)'!# (Note: value of -999 means: Clebsch-Gordan vanishes)'
     write(dens2bfile,*)' '

     write(dens2bfile,*)' '
     write(dens2bfile,*)'!# +++++++++++++++++++++++++++++++++++++++++++++'
     write(dens2bfile,*)' '
     write(resultfile,*)'!!# +++++++++++++++++++++++++++++++++++++++++++++'
     write(resultfile,*)' '
     write(resultfile,*)'!#  The above are expectation values of scalar two-body operators '
	 write(resultfile,*)'!#  that correspond to terms in the Hamiltonian. '
     write(resultfile,*)'!#  See BIGSTICK Manual for precise definitions. '
     write(resultfile,*)' '
     write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'

  end if

  return
		
end subroutine boss_twobody_densities


	
