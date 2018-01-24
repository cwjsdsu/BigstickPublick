!=================================================================
!
! file BDENSLIB1.f90
!
! master routines for computing density matrices 
!                     and applying one-body operators
!  (reorganized in 7.6.7)
!=================================================================
!
! these routines read in a one-body operator and apply it to a vector
! operators are read in as reduced matrix elements, using the conventions
! from Edmonds "Angular momentum in quantum mechanics," Eq. (5.4.1)
!
!(Jf || O_K || Ji) = [Jf] ( Jf Mf, KM | Ji Mi)^-1 ( Jf Mf | O_KM | Ji Mi)
!
!Special case: number operator:  One can easily see that 
!rho(a^+a) = sqrt(2 j_a + 1)/ sqrt( 2 Jf +1 ) n(a)
!If one includes isospin then
!rho(a^+a) = sqrt(2*(2 j_a +1 ))/ [ sqrt(2 Jf +1) sqrt(2 Tf +1)] * n(a)
!
!
!=====================================================
!
!  master subroutine for densities
!
! CALLED BY:
!   lanczos_output
!   exactdiag_p
!
! SUBROUTINES CALLED:
!    set1bsectorjumps
!    masterfindconjumps1b
!    master1bodyjumps
!    density_bundles_setup
!    master_op_stat_den
!    setMPIlimits4densities
!    master_density1b_bundled
!    coupled_densities

subroutine density1b_output 
  use system_parameters
  use sporbit
  use spstate
  use localvectors
  use nodeinfo
  use io
  use basis
  use obs
  use lanczos_info
  use localvectors
  use densities
  use haiku_info
  use precisions
  use fragments
  use mod_reorthog
  use wfn_mod
  use butil_mod
  use bvectorlib_mod
  use para_main_mod
  use jump_mod
  use lanczos_util
  use densities
  implicit none
!  include 'binterfaces.inc'

  real(4) :: xj,xt,ei,ef,xtt,xjj
  integer(4) :: i,j,m,n
  integer(4) :: jdummy
  integer(4) :: ji,ti,jf,tf,jt,tt,jmin,jmax
  integer(kind=basis_prec) :: k
  logical :: evenAJ,evenAT    ! needed to force "correct" J, T

  integer :: ierr
  integer :: aerr
  logical, parameter :: usesymmetry = .false.   ! signals to only do i->f for i >= f
!--------------------- CONVERGENCE -----------------------------

  logical smallflag,zeroflag
  
  real, allocatable :: denmat(:,:,:)
  real :: nparticles
  logical numberflag  ! flag to check on total number of particles; used for debugging
  integer inode
  character(1) :: vchar

  numberflag = .true.
  vchar = 'n'

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

  numorbmax=bmax(numorb(1),numorb(2) )

     allocate( denmat( numorbmax, numorbmax, 0:1 ), stat=aerr )
     if(aerr /= 0) then
        call memerror("density1b_output 1")
        stop 5
     end if
     allocate( p1bopme(-nhspsmax:nhspsmax, -nhspsmax:nhspsmax), stat=aerr )
     if(aerr /= 0) then
        call memerror("density1b_output 2")
        stop 5
     end if
     allocate( n1bopme(-nhspsmax:nhspsmax, -nhspsmax:nhspsmax), stat=aerr )
     if(aerr /= 0) then
        call memerror("density1b_output 3")
        stop 5
     end if

  if ( iproc == 0 ) then
        print*,' '
        print*,' Computing density matrices '
        print*,' '
  end if

  if(np(1) > 0)then
        call set1bsectorjumps(1, .false. , .false. , .false., .false. )
        call set1bsectorjumps(1, .true.  , .false. , .false. , .false.)
  endif
  if(np(2) > 0)then
        call set1bsectorjumps(2, .false. , .false. , .false., .false. )
        call set1bsectorjumps(2, .true.  , .false. , .false., .false. )
  endif

  if(np(1) > 0)then
!--------------- find conjugate jumps.......................................
     call masterfindconjumps1b(1,.false.)
     call master1bodyjumps(1,.false.)

  endif

  if(np(2) > 0)then

!--------------- find conjugate jumps..............................
        call masterfindconjumps1b(2,.false.)
        call master1bodyjumps(2,.false.)

  endif
  
  call master_para_distribute_onebody     
  if(np(1) > 0)then
     call master1bodyjumps(1,.true.)
  endif

  if(np(2) > 0)then
        call master1bodyjumps(2,.true.)
  endif
!...................LOOP OVER WFNS..............................

  if(.not.storelanczosincore1 .and. .not.storelanczosincoreMPI)then
     call wfn_rewind(wfnfile)   
     call read_wfn_header(wfnfile,.false.)
     call wfn_read_nkeep(wfnfile, j)  ! dummy reading in nkeep
  endif
  vchar='n'
  
  do i = 1,nkeep

     if(storelanczosincore1 .or. storelanczosincoreMPI)then
        if(useNewReorthog) then
           call br_retrieve_hist(i)
           ! It turns out that we need the vector loaded into both vec1 and vec2
           ! of course, they are different slices on each node
           call br_restore_vec1()
        else
           vec1 = 0.0 ! in all ranks
           call read_lanczos_vector_a(vchar,'i',i,lvec_file) ! only in rank=0
           if(storelanczosincoreMPI) then
            ! each mpi process reads one slice.  allreduce is overkill but works
            call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, icomm, ierr) ! in place
           end if
        end if
        ei = energy(i)
        xj = xjlist(i)   ! fix
        xtt= xtlist(i)   ! fix
        ji = closest2J(evenAJ,xjlist(i))
        if(isoflag)ti = closest2J(evenAT,xtlist(i))

     else
        do j = 1,i  ! temp
           ! new interface, we say which vector to read and it checks
           call wfn_readeigenvec(wfnfile, frag1, fcomm1, vec1, j, ei, xj, xtt)
        end do  ! temp
        call wfn_rewind(wfnfile)
        call read_wfn_header(wfnfile,.false.)
        call wfn_read_nkeep(wfnfile, j)    ! dummy reading in nkeep
        ji = closest2J(evenAJ,xj)
        if(isoflag)ti = closest2J(evenAT,-0.5 + sqrt(xtt+0.25))
     end if

     do j = 1,nkeep
		  if(usesymmetry .and. j > i)cycle                   ! only compute for f <= i
        if(storelanczosincore1 .or. storelanczosincoreMPI)then
           if(useNewReorthog) then
              call br_retrieve_hist(j)
              call br_restore_vec2()
           else
              vec2 = 0.0
              call read_lanczos_vector_a(vchar,'f',j,lvec_file)
              if(storelanczosincoreMPI) then
                 ! each mpi process reads one slice.  allreduce is overkill but works
                 call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, icomm, ierr) ! in place
              end if
           end if

           ef = energy(j)

           jf = closest2J(evenAJ,xjlist(j)) !         
           if(isoflag)then
              tf = closest2J(evenAT,xtlist(j)) ! 
           endif
        else
           ! new interface, we say which vec (j) to read and it checks
           call wfn_readeigenvec(wfnfile, frag2, fcomm2, vec2, j, ef, xjj, xtt)
           jf = closest2J(evenAJ,xjj) ! nint(2*xjj)          
           if(isoflag)then
              tf = closest2J(evenAT,-0.5 + sqrt(xtt+0.25)) ! nint(-1 + sqrt(4*xtt+1))
           endif
        end if 
        jmax = (jf + ji)/2
        jmin = abs( jf -ji)/2
        if ( iproc == 0 .and. isoflag) then
              write(resultfile,*)' '
              write(resultfile,333)i,ei,Ji,Ti 
              write(resultfile,334)j,ef,Jf,Tf
        end if
        if ( iproc == 0 .and. .not.isoflag) then
              write(resultfile,*)' '
              write(resultfile,433)i,ei,Ji
              write(resultfile,434)j,ef,Jf
        end if
333     format(' Initial state #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
334     format(' Final state   #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
433     format(' Initial state #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
434     format(' Final state   #',i5,' E = ',f10.5,' 2xJ   = ',i4) 

       denmat(:,:,:) = 0.0  
       call master_density1b_bundled
        do jt = jmin,jmax
              call coupled_densities(Jt,Ji,Ti,Jf,Tf,Jz,npeff(1)-npeff(2),  & 
                   zeroflag, numorbmax ,denmat)
!--------------------- CHECK PARTICLE OCCUPATIONS------------------
              if(iproc == 0)then

              if(numberflag .and. i == j .and. Jt == 0 )then
                 if(isoflag .and. .not.pndensities)then
                    nparticles = 0.0
                    do n = 1,numorbmax
                       xjj = 0.5* orbqn(1,n)%j
                       nparticles = nparticles + (denmat(n,n,0))*sqrt( 2*(2*xjj+1))
                       
                    end do  ! n
                    nparticles = nparticles /sqrt( float(jf+1)*float(tf+1))
!......... MODIFY FOR p-h CONJUGATION ....
                     if( abs( nparticles -npeff(1) -npeff(2)) > 0.001 .and. iproc==0)then
                         print*,nparticles,' particles total for state ',i
						       write(logfile,*)nparticles,' particles for state ',i,' when expecting ',npeff(1)+npeff(2)
                     end if
                 else
                    nparticles = 0.0
                    do n = 1,numorb(1)
                       xjj = 0.5* orbqn(1,n)%j
                       nparticles = nparticles + (denmat(n,n,0))*sqrt( (2*xjj+1))
                       
                    end do  ! n
                    nparticles = nparticles /sqrt( float(jf+1))
                    if( abs( nparticles -npeff(1) ) > 0.001 .and. iproc==0)then
                       print*,nparticles,' protons total '
 				    end if
                    nparticles = 0.0
                    do n = 1,numorb(2)
                       xjj = 0.5* orbqn(2,n)%j
                       nparticles = nparticles + (denmat(n,n,1))*sqrt( (2*xjj+1))
                       
                    end do  ! n
                    nparticles = nparticles /sqrt( float(jf+1))
                    if( abs( nparticles -npeff(2) ) > 0.001 .and. iproc==0)then
                       print*,nparticles,' neutrons total '
				    end if

                 endif
              endif

!---------------------------------------------------------------

              if(zeroflag)cycle
!---------------------- WRITE OUT DENSITY MATRIX ELEMENTS
              if(isoflag .and. .not.pndensities)then
                 write(resultfile,31)Jt
31               format(' Jt = ',i3,', Tt = 0        1 ' )
              else
                 write(resultfile,32)Jt
32               format(' Jt = ',i3,', proton      neutron ')

              endif
              do m = 1,numorbmax
                 do n = 1,numorbmax
                    if ( (denmat(m,n,0) /=  0.0 .or. denmat(m,n,1) /=  0.0) & 
                         .and. iproc == 0 )then
                       write(resultfile,36)m,n,(denmat(m,n,tt),tt=0,1)
                    end if
36                  format(2i5, 2f10.5)
                 end do ! b
              end do   ! a
		  end if ! iproc == 0
           end do ! jt

        end do   !j
!..... TEMPORARY
     if(.not.storelanczosincore1 .and. .not.storelanczosincoreMPI)then
        call wfn_rewind(wfnfile)     
        call read_wfn_header(wfnfile,.false.)
        call wfn_read_nkeep(wfnfile, j)  ! dummy reading in nkeep
     end if
!......END TEMPORARY
  end do  ! i
     
  if(iproc == 0)then
     write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
     write(resultfile,*)' '
     write(resultfile,*)' Definition of density matrices : '
     write(resultfile,*)' rho_K(a^+b) =   (Jf || (a^+ b)_K || Ji) / sqrt(2K+1) '
     write(resultfile,*)'  where reduced matrix element is convention of Edmonds '
     write(resultfile,*)' (note: if isospin is good symmetry, then '
     write(resultfile,*)'  doubly reduced/ divided by sqrt(2T+1) as well'
     write(resultfile,*)' '
     write(resultfile,*)' Note time-reversal symmetry relation: '
     write(resultfile,*)' rho_K(a^+b, Ji->Jf) = (-1)^(ja-jb + Ji-Jf) rho_K(b^+a, Jf->Ji) '
     write(resultfile,*)' For isospin, add in factor (-1)^(Ti-Tf)'
     write(resultfile,*)' '
     write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
     write(resultfile,*)' '

  end if
  call clocker('obs','end')

end subroutine density1b_output

!=========================================================================
!
! subroutine DENSITY1B_FROM_OLDWFN
!
! computes density from pre-existing wavefunction file
! modified/fixed by KSM in V7.4.3
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
subroutine density1b_from_oldwfn
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

  if(iproc==0)print*,' WARNING SUBROUTINE density1b_from_oldwfn NEEDS TO BE VALIDATED '

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

  numorbmax=MAX(numorb(1),numorb(2) )
 

     allocate( denmat( numorbmax, numorbmax, 0:1 ), stat=aerr )
     if(aerr /= 0) then
        call memerror("density1b_from_oldwfn 2")
        stop 5
     end if
     allocate( p1bopme(-nhspsmax:nhspsmax, -nhspsmax:nhspsmax), stat=aerr )
     if(aerr /= 0) then
        call memerror("density1b_from_oldwfn 3")
        stop 5
     end if
     allocate( n1bopme(-nhspsmax:nhspsmax, -nhspsmax:nhspsmax), stat=aerr )
     if(aerr /= 0) then
        call memerror("density1b_from_oldwfn 4")
        stop 5
     end if

  if ( iproc == 0 ) then
     print*,' '
     print*,' Computing density matrices '
     print*,' '
  end if

  call onebodysetup

!...................LOOP OVER WFNS..............................
  call wfn_read_nkeep(oldwfnfile, nkeep)

  allocate(stateE(1:nkeep), stateJ(1:nkeep), stateT(1:nkeep), stat=aerr)
  if(aerr /= 0) then
     call memerror("density1b_from_oldwfn 10")
     stop 5
  end if
  if(iproc==0)then
     print*,' '
    print*,nkeep,' states '
  end if
  ! allocate storage to store density matrices

  call pocc_write_orbits()

  do i = 1,nkeep
     ! new interface, we say which vec to read.  It checks
     call wfn_readeigenvec(oldwfnfile, frag1, fcomm1, vec1,i,ei,xj,xtt)

     ji = closest2J(evenAJ,xj)
     ! if !isoflag, then what?  orig code didn't calculate ti
     ti = closest2J(evenAT,-0.5 + sqrt(xtt+0.25))
     stateE(i) = ei
     stateJ(i) = ji
     stateT(i) = ti
	 
     do j = 1,nkeep		 
           ! new interface, we say which vec to read.  It seeks and checks
           call wfn_readeigenvec(oldwfnfile,frag2, fcomm2, vec2, j, ef,xjj,xtt)
		   
           jf = closest2J(evenAJ,xjj) ! nint(2*xjj)          
           ! if ! isoflag, then what
           tf = closest2J(evenAT,-0.5 + sqrt(xtt+0.25)) ! nint(-1 + sqrt(4*xtt+1))
           jmax = (jf + ji)/2
           jmin = abs( jf -ji)/2
           if ( iproc == 0 .and. isoflag) then
              write(resultfile,*)' '
              write(resultfile,333)i,ei,Ji,Ti 
              write(resultfile,334)j,ef,Jf,Tf
           end if
           if ( iproc == 0 .and. .not.isoflag) then
              write(resultfile,*)' '
              write(resultfile,433)i,ei,Ji
              write(resultfile,434)j,ef,Jf
           end if

333        format(' Initial state # ',i4,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
334        format(' Final state   # ',i4,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
433        format(' Initial state # ',i4,' E = ',f10.5,' 2xJ   = ',i4) 
434        format(' Final state   # ',i4,' E = ',f10.5,' 2xJ   = ',i4) 
              if(np(1) > 0) p1bopme(:,:) = 0.0
              if(np(2) > 0) n1bopme(:,:) = 0.0
              denmat(:,:,:) = 0.0
              ! compute density matrix			  
              call master_density1b_bundled
			 if(iproc==0)then
              do jt = jmin,jmax
                 call coupled_densities(Jt,Ji,Ti,Jf,Tf,Jz,npeff(1)-npeff(2),  & 
                   zeroflag, numorbmax ,denmat)
!--------------------- CHECK PARTICLE OCCUPATIONS------------------
                 if(numberflag .and. i == j .and. Jt == 0 )then
                    if(isoflag .and. .not.pndensities)then
                       nparticles = 0.0
                       do n = 1,numorbmax
                          xjj = 0.5* orbqn(1,n)%j
                          nparticles = nparticles + (denmat(n,n,0))*sqrt( 2*(2*xjj+1))
                          print *, "n=", n, ",  xjj=", xjj, ", denmat(n,n,0)=", denmat(n,n,0), ", npart=", nparticles
                          
                       end do  ! n
                       nparticles = nparticles /sqrt( float(jf+1)*float(tf+1))
                       print *, "KSM: jf=", jf,", jt=", jt,", tf=", tf, ", nparticles=", nparticles
                       if( abs( nparticles -npeff(1) -npeff(2)) > 0.001 .and. iproc==0)then
                             print*,nparticles,' particles total for state ',i
                       end if
                    else
                       nparticles = 0.0
                       do n = 1,numorb(1)
                          xjj = 0.5* orbqn(1,n)%j
                          nparticles = nparticles + (denmat(n,n,0))*sqrt( (2*xjj+1))
                          
                       end do  ! n
                       nparticles = nparticles /sqrt( float(jf+1))
                       if ( iproc == 0 )print*,nparticles,' protons total '
                       nparticles = 0.0
                       do n = 1,numorb(2)
                          xjj = 0.5* orbqn(2,n)%j
                          nparticles = nparticles + (denmat(n,n,1))*sqrt( (2*xjj+1))
                          
                       end do  ! n
                       nparticles = nparticles /sqrt( float(jf+1))
                       if ( iproc == 0 )print*,nparticles,' neutrons total '
                    endif
                 endif

!---------------------------------------------------------------

                 if(zeroflag)cycle
!---------------------- WRITE OUT DENSITY MATRIX ELEMENTS
                 if(isoflag)then
                       call dm_add_tab(i, j, Jt, numorbmax)
                       write(resultfile,31)Jt
31                     format(' Jt = ',i3,', Tt = 0        1 ' )
                 else
                       call dm_add_tab(i, j, Jt, numorbmax)
                       write(resultfile,32)Jt
32                     format(' Jt = ',i3,', proton      neutron ')

                 endif
                 do m = 1,numorbmax
                    do n = 1,numorbmax
                       if ( (denmat(m,n,0) /=  0.0 .or. denmat(m,n,1) /=  0.0) & 
                                .and. iproc == 0 )then
                          ! save for printing to Mathematica file
                          dmhead%tab(m,n,0) = denmat(m,n,0)
                          dmhead%tab(m,n,1) = denmat(m,n,1)
                          ! write to standard result file
                          write(resultfile,36)m,n,(denmat(m,n,tt),tt=0,1)
                       endif
36                     format(2i5, 2f10.5)
                    end do ! b
                 end do   ! a
              enddo ! jt
		 end if ! iproc == 0

     enddo   !j
  enddo  ! i
  if(iproc == 0)then
     write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
     write(resultfile,*)' '
     write(resultfile,*)' Definition of density matrices : '
     write(resultfile,*)' rho_K(a^+b) =   (Jf || (a^+ b)_K || Ji) / sqrt(2K+1) '
     write(resultfile,*)'  where reduced matrix element is convention of Edmonds '
     write(resultfile,*)' (note: if isospin is good symmetry, then '
     write(resultfile,*)'  doubly reduced/ divided by sqrt(2T+1) as well'
     write(resultfile,*)' '
     write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
     write(resultfile,*)' '

  end if
  call density1b_write_math(nkeep, stateE, stateJ, stateT)
  call dm_reset()
  return
end subroutine density1b_from_oldwfn

!=============================================================
!
! master routine
!
! subroutines called:
!   applyP1Bdenbundled_g
!   applyN1Bdenbundled_g
!   phconjdensity
!
  subroutine master_density1b_bundled

  use system_parameters
  use precisions
  use densities
  use nodeinfo
  use opbundles
  use bmpi_mod
  use haiku_info
  implicit none
  integer :: asps,ia,ath,it,asgn,bsps,ib,bth,bsgn,ierr
  real(kind=obs_prec) :: xme
  

  p1bopme(:,:) = 0.0
  n1bopme(:,:) = 0.0

  call applyP1Bdenbundled_g(opbundlestart(iproc),opbundleend(iproc))
  call applyN1Bdenbundled_g(opbundlestart(iproc),opbundleend(iproc))

!............ CHECK FOR PARTICLE-HOLE CONJUGATION AND SWAP IF NECESSARY.....
  if(phconj(1))call phconjdensity(1)
  if(phconj(2))call phconjdensity(2)
  
!............... NOW REDUCE .........................  
if(nproc > 1)then	
	do it = 1,2
    do asps = 1, nhsps(it)+nhsps(-it)
       if(asps > nhsps(-it))then
           ia = asps - nhsps(-it)
           ath = it
           asgn = 1
       else
           ia = asps
           ath = -it
           asgn = -1
       endif

       do bsps = 1, nhsps(it)+nhsps(-it)
          if(bsps > nhsps(-it))then
             ib = bsps - nhsps(-it)
              bth = it
             bsgn = 1
          else
             ib = bsps
             bth = -it
             bsgn = -1
          endif
		  
		  if(np(it)> 0)then
			  if(it==1)call BMPI_ALLREDUCE(p1bopme(ia*asgn,ib*bsgn),xme,1,MPI_SUM,icomm,ierr)
			  if(it==2)call BMPI_ALLREDUCE(n1bopme(ia*asgn,ib*bsgn),xme,1,MPI_SUM,icomm,ierr)
              if(it==1)p1bopme(ia*asgn,ib*bsgn)=xme
              if(it==2)n1bopme(ia*asgn,ib*bsgn)=xme

		  end if
		  
	  end do  ! bsps
	  
  end do  ! asps
end do   ! it
end if
  
  return

  end subroutine master_density1b_bundled

!=====================================================

! these routines read in a one-body operator and apply it to a vector
! operators are read in as reduced matrix elements, using the conventions
! from Edmonds "Angular momentum in quantum mechanics," Eq. (5.4.1)
!
!(Jf || O_K || Ji) = [Jf] ( Jf Mf, KM | Ji Mi)^-1 ( Jf Mf | O_KM | Ji Mi)
!
!Special case: number operator:  One can easily see that 
!rho(a^+a) = sqrt(2 j_a + 1)/ sqrt( 2 Jf +1 ) n(a)
!If one includes isospin then
!rho(a^+a) = sqrt(2*(2 j_a +1 ))/ [ sqrt(2 Jf +1) sqrt(2 Tf +1)] * n(a)
!
!
!=====================================================!
!   subroutine applicator1b
!       master subroutine to control application of a 1-body operator
!   revised in 7.4.4 to better use current constructs
!
!  CALLED BY:  main routine
!
!  SUBROUTINES CALLED:
!   set1bsectorjumps
!   masterfindconjumps1b
!    master1bodyjumps
!    master_op_stat_den
!  setnodaltribution  
!  setup_localvectors
!  density_bundles_setup
!   setMPIlimits4densities
!  wfn_read_nkeep
!  wfn_write_nkeep
!   wfn_readeigenvec
!  master_apply1bop
!  wfn_writeeigenvec
!
   subroutine applicator1b
   use system_parameters
   use precisions
   use basis
   use io
   use nodeinfo
   use menu_choices
   use wfn_mod
   use localvectors
   use jumpNbody
   use para_main_mod
   use hoppy
   use jump_mod
   use bvectorlib_mod
   implicit none

   integer nkeep
   integer i,j,n
   real :: e,xj,xt2
   integer :: aerr,ierr

   call onebodysetup

   call wfn_read_nkeep(oldwfnfile, nkeep)

   if(iproc==0)then
      print *, ' '
      print *, nkeep, ' states '
   end if
   call wfn_write_nkeep(nkeep) ! write number of vectors to wfn file
   do i = 1,nkeep
      call wfn_readeigenvec(oldwfnfile, frag1, fcomm1, vec1,i,e,xj,xt2) ! KSM: updated
      call master_apply1bop        ! new version added in 7.6.7
!......... OPTION TO ENFORCE ORTHOGONALITY.........CURRENTLY TURNED OFF

      if(writeout) call wfn_writeeigenvec(wfnfile,frag2, vec2,i,e,xj,xt2)
   end do

   return
   end subroutine applicator1b
   
   
!============================================================
!
!  
!   

subroutine onebodysetup
    use system_parameters
    use hoppy
    use jump_mod
    use bvectorlib_mod
	use btbme_mod
	use onebodypot, only: meanie
	use ntuple_info
	use nodeinfo
	use fragments, only: nodal
	implicit none
	integer :: aerr
    call hopmaker
	if(meanie)then
!----------------- RESTRICTIONS
! In order to save space, one needs to restrict the uncoupled matrix elements 
! created.  In order to do this, one needs the hops beforehand
		call master_q_ntuple(2,1)
	    call master_q_ntuple(2,2)
	    call master_cross_ntuple(1,1)
	    call prepare_to_uncoupleXXtbme(1)	
	    call prepare_to_uncoupleXXtbme(2)
	    call prepare_to_uncouplePNtbme

	end if
	
  !-------------- SET UP JUMPS
     if(np(1) > 0)then
          call set1bsectorjumps(1, .false. , .false. , .false., .false. )
          call set1bsectorjumps(1, .true.  , .false. , .false. , .false.)
     endif
     if(np(2) > 0)then
          call set1bsectorjumps(2, .false. , .false. , .false., .false. )
          call set1bsectorjumps(2, .true.  , .false. , .false., .false. )
     endif

     if(np(1) > 0)then
  !--------------- find conjugate jumps
          call masterfindconjumps1b(1,.false.)
          call master1bodyjumps(1,.false.)
     endif
     if(np(2) > 0)then

  !--------------- find conjugate jumps
          call masterfindconjumps1b(2,.false.)
          call master1bodyjumps(2,.false.)
     endif

   
     call master_para_distribute_onebody   
     if(.not.meanie)call setup_localvectors  

   
     if(np(1) > 0)then
        call master1bodyjumps(1,.true.)
     endif

     if(np(2) > 0) then
           call master1bodyjumps(2,.true.)
     endif
	 
	 return
	
end subroutine onebodysetup

!============================================================
!
!  CALLED BY main routine 
!

  subroutine readin1bodyop
  use opmatrixelements
  use sporbit
  use nodeinfo
  use bmpi_mod
  use system_parameters
  implicit none

  logical success
  character*15 :: appfile
  character*40 :: title
  integer ilastap
  integer i,i1
  integer a,b
  real xj
  real :: xxx,yyy
 
!-------------- DUMMIES FOR CONSISTENCY CHECK --------------
  integer n,l,j
  integer :: aerr,ierr
  

!---------------- ALLOCATED ARRAY FOR COUPLED OPERATOR MATRIX ELEMENTS

!  if(.not.isoflag)then
!    print*,' Sorry, cannot yet handle pn for densities '
!    stop
!  endif

!------------- READ IN OPERATOR FILE AND DECOUPLE

if(iproc == 0)then
  success = .false.
  do while(.not.success)  
     write(6,*)' Enter name of .opme  file '
     read(5,'(a)')appfile
     ilastap = index(appfile,' ') -1
     open(unit = 23,file=appfile(1:ilastap)//'.opme',status = 'old',err=301)
     success = .true.
     cycle
301  continue
     print*,'File ',appfile(1:ilastap),'.opme does not exist '
   enddo
!------- READ IN HEADER???
   read(23,'(a)')title
   print*,title
   
   select case (title(1:3))
      case ('iso')
	  pnoperators = .false.
  
	  case ('pns')
	  pnoperators = .true.
	  
	  case ('xpn')
	  print*,' Currently BIGSTICK does not accept .opme files in xpn format '
	  print*,' Instead for proton-neutron breaking it must be in pns format, with two columns '
	  print*,' See section 4.6.2 in the BIGSTICK Manual/Overview '
	  stop
	
      case default
	  print*,'  Something wrong with first line of .opme file '
	  print*,title(1:3)
	  print*,title
	  stop
   end select
end if

!.......... SET UP OPERATOR ARRAYS.................

call BMPI_BCAST(pnoperators,1,0,icomm,ierr)
if(pnoperators)then
  if(np(1)> 0)then
	  
     allocate ( pop1bod( numorb(1), numorb(1) ), stat=aerr)
     if(aerr /= 0) call memerror("readin1bodyop -- p")
     pop1bod = 0.0
  end if
  if(np(2)> 0)then
	  
     allocate ( nop1bod( numorb(2), numorb(2) ), stat=aerr)
     if(aerr /= 0) call memerror("readin1bodyop -- n")
     nop1bod = 0.0
  end if   	
	
else
  allocate ( op1bod( numorb(1), numorb(1) ), stat=aerr)
  if(aerr /= 0) call memerror("readin1bodyop -- iso")
  op1bod = 0.0	
end if

if(iproc==0)then
   
   
!------------- CHECK FOR CONSISTENCY----------------
   read(23,*)n
   if(n /= numorb(1))then
     print*,' Mismatch in single particle orbits ',n,numorb(1) 
     stop
   endif

   do i = 1, numorb(1)
     read(23,*)i1,n,l,xj
     if(i1 /= i)then
         print*,' ooopsie ',i1,i
         stop
     endif
     j = nint(2*xj)
     if( n/= orbqn(1,i)%nr .or. l /= orbqn(1,i)%l .or. j /= orbqn(1,i)%j)then
        print*,i,' Mismatch in orbits ',n,l,j
        print*, orbqn(1,i)%nr, orbqn(1,i)%l, orbqn(1,i)%j
        stop
     endif

   enddo
!---------- READ IN J,T for operator----------------------
  if(.not.pnoperators)then
    read(23,*)jop,top      ! J and T of operator
    write(6,*)' Operator has J, T = ',jop,top
  else
      read(23,*)jop     ! J and T of operator
      write(6,*)' Operator has J ',jop
	  top = 2
  end if

!--- MKGK option - only for r^2 Y00 ----
!--- for strength only to excited states and not back to the gs
  if (jop == 0 .and. subtract_enabled) then 
     subtract = .true.
  else
     subtract = .false.
  end if

!------------ READ IN REDUCED MATRIX ELEMENTS ----------------
  do i = 1, 1000
	  if(pnoperators)then
		  read(23,*,end=348)a,b,xxx,yyy
	  else
        read(23,*,end =348)a,b,xxx
	  end if
!------------- ERROR TRAP---------------------
    if( orbqn(1,a)%j+orbqn(1,b)%j < 2*jop .or. abs( orbqn(1,a)%j-orbqn(1,b)%j) > 2*jop) then
       if(iproc == 0)then
         write(6,*)' Mismatch in angular momentum of operator '
         write(6,*)' one-body states : ',a,b
         write(6,*)orbqn(1,a)%j/2.,orbqn(1,b)%j/2. , jop

       endif
       stop
    endif
	if(pnoperators)then
		if(np(1)>0)pop1bod(a,b)=xxx
		if(np(2)>0)nop1bod(a,b)=yyy
	else
        op1bod(a,b) = xxx
    end if
  enddo

348      continue

  close(unit=23)
end if
!--------------- NOW BROADCAST -------------------------------
!    a bit of a kludge
call BMPI_BCAST(jop,1,0,icomm,ierr)
call BMPI_BCAST(top,1,0,icomm,ierr)
do a = 1,numorb(1)
	do b = 1,numorb(1)
		if(pnoperators)then
            if(np(1)> 0)call BMPI_BCAST(pop1bod(a,b),1,0,icomm,ierr)
            if(np(2)> 0)call BMPI_BCAST(nop1bod(a,b),1,0,icomm,ierr)
			
			
		else
           call BMPI_BCAST(op1bod(a,b),1,0,icomm,ierr)
	   
       end if
    end do
end do
!print*,iproc,' op1bod mes ',op1bod
  return
  end subroutine readin1bodyop
!=====================================================
!
!  Routine to read in one-body transition operators 
!  allowing for 3 formats  (added in 7.5.2); not yet fully implemented
!
!  CALLED BY not yet implemented
!
  subroutine readin1bodyop_gen
  use densities
  use opmatrixelements
  use sporbit
  use nodeinfo

  implicit none

  logical success
  character*15 :: appfile
  character*40 :: title
  integer ilastap
  integer i,i1
  integer a,b
  real xj
  real :: xxx
!-------------- DUMMIES FOR CONSISTENCY CHECK --------------
  integer n,l,j
  integer :: aerr

!---------------- ALLOCATED ARRAY FOR COUPLED OPERATOR MATRIX ELEMENTS

  if(.not.isoflag)then
    print*,' Sorry, cannot yet handle pn for densities '
    stop
  endif
  allocate ( op1bod( numorb(1), numorb(1) ), stat=aerr)
  if(aerr /= 0) call memerror("readin1bodyop")
  op1bod = 0.0
!------------- READ IN OPERATOR FILE AND DECOUPLE

  success = .false.
  do while(.not.success)  
     write(6,*)' Enter name of .opme  file '
     read(5,'(a)')appfile
     ilastap = index(appfile,' ') -1
     open(unit = 23,file=appfile(1:ilastap)//'.opme',status = 'old',err=301)
     success = .true.
     cycle
301  continue
     print*,'File ',appfile(1:ilastap),'.opme does not exist '
   enddo
!------- READ IN HEADER???
   read(23,'(a)')title
   print*,title
!------------- CHECK FOR CONSISTENCY----------------
   read(23,*)n
   if(n /= numorb(1))then
     print*,' Mismatch in single particle orbits ',n,numorb(1) 
     stop
   endif

   do i = 1, numorb(1)
     read(23,*)i1,n,l,xj
     if(i1 /= i)then
         print*,' ooopsie ',i1,i
         stop
     endif
     j = nint(2*xj)
     if( n/= orbqn(1,i)%nr .or. l /= orbqn(1,i)%l .or. j /= orbqn(1,i)%j)then
        print*,i,' Mismatch in orbits ',n,l,j
        print*, orbqn(1,i)%nr, orbqn(1,i)%l, orbqn(1,i)%j
        stop
     endif

   enddo
!---------- READ IN J,T for operator----------------------

  read(23,*)jop,top      ! J and T of operator
  write(6,*)' Operator has J, T = ',jop,top

!--- MKGK option - only for r^2 Y00 ----
!--- for strength only to excited states and not back to the gs
  if (jop == 0 .and. subtract_enabled) then 
     subtract = .true.
  else
     subtract = .false.
  end if

!------------ READ IN REDUCED MATRIX ELEMENTS ----------------
  do i = 1, 1000
    read(23,*,end =348)a,b,xxx
!------------- ERROR TRAP---------------------
    if( orbqn(1,a)%j+orbqn(1,b)%j < 2*jop .or. abs( orbqn(1,a)%j-orbqn(1,b)%j) > 2*jop) then
       if(iproc == 0)then
         write(6,*)' Mismatch in angular momentum of operator '
         write(6,*)' one-body states : ',a,b
         write(6,*)orbqn(1,a)%j/2.,orbqn(1,b)%j/2. , jop

       endif
       stop
    endif

    op1bod(a,b) = xxx
  enddo

348      continue

  close(unit=23)
  return
  end subroutine readin1bodyop_gen
!=====================================================
!
!
! CALLED BY main routine
!
!

  subroutine decouple1bodyop

  use opmatrixelements
  use onebodypot
  use sporbit
  use spstate
  use haiku_info
  use nodeinfo
  use system_parameters

  implicit none

  integer a,b
  integer ja,jb
  integer ma,mb
  integer ia,ib
  integer asps, bsps
  integer asgn, bsgn, ath, bth
  integer iphase
  real   cleb
  integer :: aerr

!  if(.not.isoflag)then
!    print*,' Sorry, cannot yet handle pn for densities '
!    stop
!  endif
!------------- ALLOCATE UNCOUPLED ARRAYS ------------------

  allocate(pPOT_h(-nhsps(-1):nhsps(1), -nhsps(-1):nhsps(1)), stat=aerr)
  if(aerr /= 0) call memerror("decouple1bodyop 1")
  allocate(nPOT_h(-nhsps(-2):nhsps(2), -nhsps(-2):nhsps(2)), stat=aerr)
  if(aerr /= 0) call memerror("decouple1bodyop 2")

  pPOT_h = 0.0
  nPOT_h = 0.0

!--------------- PROTONS ------------------------
  do asps = 1, nhsps(1)+nhsps(-1)
     if(asps > nhsps(-1))then
         ia = asps - nhsps(-1)
         ath = 1
         asgn = 1
     else
         ia = asps
         ath = -1
         asgn = -1
     endif

     a = hspsqn(ath,ia)%orb
     ja= hspsqn(ath,ia)%j
     ma= hspsqn(ath,ia)%m


     do bsps = 1, nhsps(1)+nhsps(-1)
        if(bsps > nhsps(-1))then
           ib = bsps - nhsps(-1)
            bth = 1
            bsgn = 1
        else
           ib = bsps
           bth = -1
           bsgn = -1
        endif
        b = hspsqn(bth,ib)%orb
        jb= hspsqn(bth,ib)%j
        mb= hspsqn(bth,ib)%m

        if(pnoperators )then
			if(np(1)==0)cycle
            if(pop1bod(a,b) == 0.0)cycle
			
            if( ma /= mb) cycle
            if(jop > (ja + jb)/2 )cycle
            if(jop < abs(ja - jb)/2 ) cycle
            iphase = (-1)**( (jb -mb)/2)

            pPOT_h(ia*asgn,ib*bsgn) = iphase*cleb(ja,ma,jb,-mb,2*jop,ma-mb)* & 
                   pop1bod(a,b)/sqrt(2.*jop+1.)
			
		else
           if(op1bod(a,b) == 0.0)cycle
           if( ma /= mb) cycle
           if(jop > (ja + jb)/2 )cycle
           if(jop < abs(ja - jb)/2 ) cycle
           iphase = (-1)**( (jb -mb)/2)

           pPOT_h(ia*asgn,ib*bsgn) = iphase*cleb(ja,ma,jb,-mb,2*jop,ma-mb)* & 
                  op1bod(a,b)/sqrt(2.*top+1)/sqrt(2.*jop+1.)*cleb(1,1,1,-1,2*top,0)
				  
	    end if
     enddo  ! bsps
  enddo   !asps
!-------------------- NEUTRONS------------------
  do asps = 1, nhsps(2)+nhsps(-2)
     if(asps > nhsps(-2))then
         ia = asps - nhsps(-2)
         ath = 2
         asgn = 1
     else
         ia = asps
         ath = -2
         asgn = -1
     endif

     a = hspsqn(ath,ia)%orb
     ja= hspsqn(ath,ia)%j
     ma= hspsqn(ath,ia)%m


     do bsps = 1, nhsps(2)+nhsps(-2)
        if(bsps > nhsps(-2))then
           ib = bsps - nhsps(-2)
            bth = 2
            bsgn = 1
        else
           ib = bsps
           bth = -2
           bsgn = -1
        endif
        b = hspsqn(bth,ib)%orb
        jb= hspsqn(bth,ib)%j
        mb= hspsqn(bth,ib)%m
        if(pnoperators)then
			if(np(2)==0)cycle
            if(nop1bod(a,b) == 0.0)cycle
			
            if( ma /= mb) cycle
            if(jop > (ja + jb)/2 )cycle
            if(jop < abs(ja - jb)/2 ) cycle
            iphase = (-1)**( (jb -mb)/2)

            nPOT_h(ia*asgn,ib*bsgn) = iphase*cleb(ja,ma,jb,-mb,2*jop,ma-mb)* & 
                   nop1bod(a,b)/sqrt(2.*jop+1.)
			
		else
           if(op1bod(a,b) == 0.0)cycle
           if( ma /= mb) cycle
           if(jop > (ja + jb)/2 )cycle
           if(jop < abs(ja - jb)/2 ) cycle
!           iphase = (-1)**( (jb -mb)/2)
           iphase = (-1)**( (jb -mb)/2)*(-1) ! for neutrons

           nPOT_h(ia*asgn,ib*bsgn) = iphase*cleb(ja,ma,jb,-mb,2*jop,ma-mb)* & 
                 op1bod(a,b)/sqrt(2.*top+1)/sqrt(2.*jop+1.)*cleb(1,-1,1,1,2*top,0)				  
	    end if


     enddo  ! bsps
  enddo   !asps
  
!  print*,iproc,' nPOT_h mes : ',nPOT_h(6,6)

  return
  end subroutine decouple1bodyop

!=====================================================
!
!  orchestrates application of a one-body operator on the wavefunction
!  REVISION 7.6.7 (July 2016) to better conform to opbundles
!  in preparation to fully functionality with MPI w/fragments
!
!  CALLED BY
!    applicator1b
!  CALLS
!    applyP1bopbundled_g
!    applyN1bopbundled_g
!
!
     subroutine master_apply1bop

    use precisions
    use basis
    use sectors
    use jumpNbody
    use system_parameters
    use opbundles
    use opmatrixelements, only : subtract
    use nodeinfo, only : iproc,icomm,nproc
    use localvectors
    use fragments
    use bmpi_mod
	use mod_reorthog
	use lanczos_util
	use bvectorlib_mod
    implicit none

    integer :: ps,ns
    integer :: ips,ins
    integer :: fps,fns
    integer cs
    integer (kind=basis_prec) :: i
    real(kind=obs_prec) :: dot
    real(kind=8) :: dot8
    integer:: ierr
	
! TO BE DELETED
    real(kind=lanc_prec) :: vin(v1s:v1e), vout(v2s:v2e)
    ! locals
    real(kind=lanc_prec) :: vin2(v2s:v2e)

    call initialize_final('n')
	
  !---------------------------------- PROTONS --------------------
    if(np(1) > 0  )then
		call applyP1bopbundled_g(	opbundlestart(iproc),opbundleend(iproc))
    endif
  !-------------------------------- NEUTRONS ---------------------------
    if(np(2) > 0 )then
		call applyN1bopbundled_g(	opbundlestart(iproc),opbundleend(iproc))
    endif
	
    if(useNewReorthog) then
           ! call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
           ! Do reduce only onto root node.  We will be sending this data to 
           ! the slices from the isfragroot nodes (rank=0 in fcomm1, fcomm2, hcomm),
           ! so other nodes don't need it.
!		   print*,iproc,' before ',vec2(28)
           call BMPI_REDUCE(vec2, size(vec2), MPI_SUM, 0,fcomm2, ierr) ! in place
		   call br_grab_vec2()
!           if(nproc > 1)call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
!		   print*,iproc,' after ',vec2(28)
		   
    else
        if(nproc > 1 )then
           call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
        endif
    end if
!    if(iproc==0)print*,'(D5)'
!    call BMPI_BARRIER(icomm,ierr)  !just for debugging
    if(useNewReorthog) then
    	! have to reduce, fragment by fragment
!  	    call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place reduce
  	! copy vin to vin2.  Non trival
!      	if(isfragroot) vin2 = vin
!  	    call BMPI_BCAST(vin2, size(vin2), 0, fcomm2, ierr)

!--- Special MKGK option 
  	if(subtract) then
          if (iproc == 0) write(*,*) 'subtract in effect'
		  if(iproc == 0)write(*,*)' Sorry, MKGK, subtract currently unavailable'
!          dot=0.d0
!          do i = v2s,v2e
!             dot = dot + vout(i)*vin2(i)
!          enddo
!  		! note, only one node per fragment can contribute to reduction
!  		dot8 = 0
!  		if(isfragroot) dot8 = dot
!  		call BMPI_ALLREDUCE(dot8, 1,  MPI_SUM, icomm, ierr) ! in place reduce
!  		dot = dot8
!          vout(:) = vout(:) - real(dot*vin2(:), kind=lanc_prec)
  	end if
    else
!--- Special MKGK option 
    if (subtract) then
       if (iproc == 0) write(*,*) 'subtract in effect'
	  if(iproc == 0)write(*,*)' Sorry, MKGK, subtract currently unavailable'

!       dot=0.d0
!       do i = 1, dimbasis
!          dot = dot + vout(i)*vin(i)
!       enddo
!       vout(:) = vout(:) - real(dot*vin(:), kind=lanc_prec)
    end if
    end if

    return

    end subroutine master_apply1bop

!======================================================

!
! routine to compute overlaps
! reads in one (1) initial wavefunction from one file, then 
! computes overlaps with all wavefunctions in a second file
!
! added 11/2010 by CWJ @ SDSU
! fixed in 7.4.3 by CWJ
! Fix some issues with wfn interface  5 Mar 2015 by KSM
!
!  CALLED BY:
!    main routine
!
!  SUBROUTINES CALLED:
!    defaultopstat
!    setnodaltribution   
!    setup_localvectors
!    wfn_read_nkeep
!    wfn_readeigenvec
!    wfn_close_file
!    wfn_ropen_file
!    read_wfn_header
!    doverlapvec
!
subroutine overlap
   use system_parameters
   use sporbit
   use haiku_info
   use menu_choices
   use precisions
   use basis
   use io
   use densities
   use nodeinfo
   use wfn_mod
   use localvectors
   use bmpi_mod
   use para_main_mod
   use bvectorlib_mod
   use lanczos_util
   implicit none

   real (kind = lanc_prec), allocatable :: v(:),w(:)
   real (kind = 8) :: dsclrprod
!   integer nkeep
   integer ikeep
   integer i,j,n, tidx
   real :: xe,xj,xt2
   logical zeroflag
   integer(kind=basis_prec) :: oldbasis
   integer ierr
   integer :: aerr
   integer, parameter :: overlapfn = 81

   if(iproc==0)then
      print*,' in overlap '
      open(unit = overlapfn, file='overlap.dat',status='unknown')
   end if
   dimbasis = dimbasischeck
   oldbasis = dimbasischeck
   ! allocate(v(dimbasis),w(dimbasis))
   ! We are doing vector ops, so both v and w should match frag1
 
   call overlaptribution
   call setup_localvectors

   call wfn_read_nkeep(oldwfnfile, nkeep) ! does BCAST
   ikeep = nkeep

   if(iproc==0) print*,' There are ', nkeep,' wavefunctions '
   do i = 1,nkeep
      ! new interface - we say which vec to read, it checks
      ! KSM:  This will be very slow, only need to read head part of each vector
      call wfn_readeigenvec(oldwfnfile,frag1, fcomm1, vec1,i,xe,xj,xt2)  
      if(iproc==0)print*,i,xe,xj,xt2
   enddo

   if(iproc == 0) then
      tidx = -1
      do while(tidx < 1 .or. tidx > nkeep)
         print*,' Which do you want as initial state? '
         if(auto_input)then
            read(autoinputfile,*) tidx
            if(tidx < 1 .or. tidx > nkeep) then
               print *, "vector selection out of range: ", tidx
               stop 1
            end if
         else
            read(5,*) tidx
            if(tidx < 1 .or. tidx > nkeep) then
               print *, "vector selection out of range: ", tidx, ", please try again"
            else
               write(autoinputfile,*) tidx
            end if
         end if
      end do
      write(overlapfn,*)' Initial state = ',tidx
   end if
   
   call BMPI_BCAST(tidx,1,0,icomm,ierr)
   ! new interface - we say which vec to read, it seeks and reads
   ! no need to rewind and read forward
   call wfn_readeigenvec(oldwfnfile, frag1, fcomm1, vec1,tidx,xe,xj,xt2)
   call wfn_close_file(oldwfnfile)
   oldbasis = dimbasis
   if(iproc==0)then
      print*,' '
      print*,' FINAL STATES '
   end if
   call wfn_ropen_file(wfnfile)
   call read_wfn_header(wfnfile,.false.)
   if(dimbasischeck /= oldbasis)then
      if(iproc == 0)then
         print*,' mismatch in basis dimensions '
         print*,' initial ',oldbasis
         print*,' final ',dimbasischeck
      endif
      stop
   endif
   call wfn_read_nkeep(wfnfile, nkeep)

   if(iproc ==0)write(overlapfn,*)' state      E      J    T       <i|f>     |<i|f>|^2 '
   do i = 1,nkeep
       ! new interface, we say which vec to read.  It checks
       call wfn_readeigenvec(wfnfile, frag2, fcomm2, vec2,i,xe,xj,xt2)
       dsclrprod = 0.d0
	   call doverlapvec(dsclrprod)  ! modified in 7.7.0 so as to work with fragments
       if(iproc == 0)then
          write(6,101)i,xe,xj,xt2,dsclrprod,dsclrprod**2
          write(overlapfn,101)i,xe,xj,xt2,dsclrprod,dsclrprod**2
       end if
101 format(i4,3x,f10.5,2f5.1,2f12.5)
   enddo
   if(iproc==0)close(unit = overlapfn)
   if(iproc ==0)write(6,*)' results written to overlap.dat '
   stop
  end subroutine overlap
!==============================================================
!  called when carrying out overlaps
!
!  called by:
!    overlap
! 
subroutine defaultopstat 
   use fragments
   use nodeinfo
   implicit none
   integer :: ifrag,jfrag
   integer :: aerr
	
   allocate(opfragstat(nfragments,nfragments),stat=aerr)
   if(aerr /= 0) call memerror("defaultopfrag")
   do ifrag=1,nfragments
      do jfrag=1,nfragments
         opfragstat(ifrag,jfrag)%nnodes=1
      end do   ! jfrag
   end do   ! ifrag	
end subroutine defaultopstat 

!======================================================
!
!  subroutine to set distribution of fragments and processes
!  for overlap calculation
!
subroutine overlaptribution
    use fragments
    use nodeinfo
    use io
	integer :: i
    integer :: aerr
    
	if(allocated(nodal))then
		if(iproc==0)print*,' ERROR nodal already allocated before we got to overlaptribution '
		stop
	end if  
    allocate(nodal(0:nprocs-1), stat=aerr )
	nodal(:)%ifragment = -1
	nodal(:)%ffragment = -1
	nodal(:)%ifirst    = .false.
    if(aerr /= 0) call memerror("overlaptribution")

	if(nproc < nfragments)then
		if(iproc==0)print*,' PROBLEM with fragments and MPI processes ',nfragments, nproc
	    stop
	end if
	do i = 1,nfragments
		nodal(i-1)%ifragment = i
		nodal(i-1)%ffragment = i
		nodal(i-1)%ifirst    =.true.
	end do
!......... DUMMY INDICES......
    do i = nfragments+1,nproc
		nodal(i-1)%ifragment=nfragments
		nodal(i-1)%ffragment=nfragments
		nodal(i-1)%ifirst = .false.
		
	end do	
	return
	
end subroutine overlaptribution

!======================================================

!
!  swaps indices (+phase) if there is particle-hole conjugation
!  (Originally added 7.2.5 by CWJ @ SDSU )
!
!  NOTES added in version 7.5.6 by CWJ 9/2015
!  a hole creation operator of angular momentum j,m is equivalent to 
!  a particle annihilation operator j,-m with phase (-1)^(j+m)
!
!  The arrays x1bopme(asps,bsps) is empty unless m_asps = m_bsps
!  To properly map this, we must find the conjugates for asps and bsps
!  (and the associated phase)
!  Note: the time-reverse state is already stored in hspsqn(it,i)%tr
!
!  NOTE: MIGHT NOT WORK FOR "SPINLESS" SYSTEMS
!
subroutine phconjdensity(it)
   use system_parameters
   use spstate
   use densities
   use precisions
   use haiku_info
   implicit none
!........ INPUT...............
   integer :: it   ! species 
!....... INTERNAL.......
   real(kind=obs_prec), pointer :: x1bopme(:,:)
   integer a,b,ja,ma,jb,mb
   integer ia,ib
   integer asps,ath,bsps,bth,asgn,bsgn
   integer :: iatr, ibtr       ! time-reversed states
   real(kind=obs_prec) :: diag
   real(kind=obs_prec) :: tmpden

   if(.not. phconj(it) .or. np(it) < 1 ) return
   if(it==1)then
      x1bopme=> p1bopme
   else
      x1bopme=> n1bopme
   end if
   do asps = 1, nhsps(it)   ! only loop over states with m >= 0
      ia = asps 
      ath = it
      asgn = 1
      a = hspsqn(ath,ia)%orb
      ja= hspsqn(ath,ia)%j
      ma= hspsqn(ath,ia)%m
  	  iatr = hspsqn(ath,ia)%tr

      do bsps = 1, nhsps(it)
         ib = bsps 
         bth = it
         bsgn = 1

         b = hspsqn(bth,ib)%orb
         jb= hspsqn(bth,ib)%j
         mb= hspsqn(bth,ib)%m
     	 ibtr = hspsqn(bth,ib)%tr		 

         if( ma /= mb) cycle       ! note: one-body density should be empty unless ma == mb
!         if(asps /=bsps)then    !simply swap
         if(asps==bsps)then
   	         diag = 1.0
         else
	          diag = 0.0
          end if
             tmpden = x1bopme(ia*asgn,ib*bsgn)
             x1bopme(ia*asgn,ib*bsgn)= diag-x1bopme(-ibtr*bsgn,-iatr*asgn) *(-1)**( (ja-jb)/2 )
             x1bopme(-ibtr*bsgn,-iatr*asgn)=diag- tmpden *(-1)**( (ja-jb)/2 )

      end do ! b
   end do  ! a
   return
end subroutine phconjdensity
!======================================================
