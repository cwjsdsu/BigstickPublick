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
!  MODIFIED IN 7.8.2
!  
!
! CALLED BY:
!   lanczos_output
!   exactdiag_p
!
subroutine density1b_output
	use densities
	use nodeinfo
	use io
	use writebinden1b
	use lanczos_info,only:nkeep
	implicit none
	
	call clocker('all','end')
	call clockout('tmp')
	call clocker('den','sta')
	
		call density1b_compute
		
		if(iproc==0)write(6,*)' One-body densities computed '
		call clocker('all','end')
		call clockout('tmp')
		call clockout_den1b

		if(iproc==0)then
		   write(6,*)' '
		   write(6,*)' Writing densities to file '
	    end if
		if(binden)then
			call binary_density_boss(1,nkeep)
		else
		   call density1b_writeout(1,nkeep,1,nkeep)
	    end if
	
	return
end subroutine density1b_output
	
!=========================================================================
!
!  master subroutine for densities
!  NEW VERSION 7.8.2:
!     uses time-reversal symmetry and stores values in a derived type, 
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

subroutine density1b_compute
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
  implicit none

  real(4) :: xj,xt,ei,ef,xtt,xjj
  integer(4) :: i,j,m,n
  integer(4) :: jdummy
  integer(4) :: ji,ti,jf,tf,jt,tt,jmin,jmax
  integer(4) :: phase
  integer(4) :: ja,jb
  integer(kind=basis_prec) :: k
  logical :: evenAJ,evenAT    ! needed to force "correct" J, T

  integer :: ierr
  integer :: aerr
!--------------------- CONVERGENCE -----------------------------

  logical smallflag,zeroflag
  
  real(8), allocatable :: denmat(:,:,:)
  real :: nparticles
  logical numberflag  ! flag to check on total number of particles; used for debugging
  integer inode
  character(1) :: vchar
  logical :: sameif

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
        call memerror("density1b_compute 1")
        stop 5
     end if
     allocate( p1bopme(-nhspsmax:nhspsmax, -nhspsmax:nhspsmax), stat=aerr )
     if(aerr /= 0) then
        call memerror("density1b_compute 2")
        stop 5
     end if
     allocate( n1bopme(-nhspsmax:nhspsmax, -nhspsmax:nhspsmax), stat=aerr )
     if(aerr /= 0) then
        call memerror("density1b_compute 3")
        stop 5
     end if
	 allocate( densitybag(nkeep,nkeep), stat=aerr)
     if(aerr /= 0) then
        call memerror("density1b_compute 4")
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
     call masterfindconjumps1b(1,.false.,.true.)
     call master1bodyjumps(1,.false.)

  endif

  if(np(2) > 0)then

!--------------- find conjugate jumps..............................
        call masterfindconjumps1b(2,.false.,.true.)
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
#ifdef _MPI
            call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
#endif
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
           call wfn_readeigenvec(wfnfile, frag1, fcomm1_index, vec1, j, ei, xj, xtt)
        end do  ! temp
        call wfn_rewind(wfnfile)
        call read_wfn_header(wfnfile,.false.)
        call wfn_read_nkeep(wfnfile, j)    ! dummy reading in nkeep
        ji = closest2J(evenAJ,xj)
        if(isoflag)ti = closest2J(evenAT,-0.5 + sqrt(xtt+0.25))
     end if
	   if(only_odd .and. (ji/4)*4==ji)then
		   cycle
	   end if

     do j = 1,nkeep
!		  if(usesymmetry .and. j > i)cycle                   ! only compute for f <= i
		  if( j > i)cycle                   ! only compute for f <= i
        if(storelanczosincore1 .or. storelanczosincoreMPI)then
           if(useNewReorthog) then
              call br_retrieve_hist(j)
              call br_restore_vec2()
           else
              vec2 = 0.0
              call read_lanczos_vector_a(vchar,'f',j,lvec_file)
              if(storelanczosincoreMPI) then
                 ! each mpi process reads one slice.  allreduce is overkill but works
#ifdef _MPI
                 call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
#endif
              end if
           end if

           ef = energy(j)

           jf = closest2J(evenAJ,xjlist(j)) !         
           if(isoflag)then
              tf = closest2J(evenAT,xtlist(j)) ! 
           endif
        else
           ! new interface, we say which vec (j) to read and it checks
           call wfn_readeigenvec(wfnfile, frag2, fcomm2_index, vec2, j, ef, xjj, xtt)
           jf = closest2J(evenAJ,xjj) ! nint(2*xjj)          
           if(isoflag)then
              tf = closest2J(evenAT,-0.5 + sqrt(xtt+0.25)) ! nint(-1 + sqrt(4*xtt+1))
           endif
        end if 
  	   if(only_odd .and. (jf/4)*4==jf)cycle
		
        jmax = (jf + ji)/2
        jmin = abs( jf -ji)/2
!.......... NOW SET UP DERIVED TYPE FOR STORING DENSITY MATRICES..		
		densitybag(i,j)%jmin=jmin
		densitybag(i,j)%jmax=jmax
		densitybag(j,i)%jmin=jmin
		densitybag(j,i)%jmax=jmax
		
		if(.not.allocated( densitybag(i,j)%denmat ))then
			allocate(densitybag(i,j)%denmat(jmin:jmax,numorbmax,numorbmax,0:1), stat=aerr )
            if(aerr /= 0) then
               call memerror("density1b_compute 7")
               stop 5
            end if
		end if
		densitybag(i,j)%denmat=0.d0
		allocate(densitybag(i,j)%zeroflag(jmin:jmax), stat=aerr )
        if(aerr /= 0) then
            call memerror("density1b_compute 7b")
            stop 5
        end if
		
		if(i/=j)then
			allocate(densitybag(j,i)%denmat(jmin:jmax,numorbmax,numorbmax,0:1), stat=aerr )
	        if(aerr /= 0) then
	            call memerror("density1b_compute 8")
	            stop 5
	        end if
			densitybag(j,i)%denmat=0.d0
			allocate(densitybag(j,i)%zeroflag(jmin:jmax), stat=aerr )
	        if(aerr /= 0) then
	            call memerror("density1b_compute 9")
	            stop 5
	        end if
			
		
		end if
!.................... FINISHED ALLOCATING MEMORY......................................		

333     format(' Initial state #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
334     format(' Final state   #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
433     format(' Initial state #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
434     format(' Final state   #',i5,' E = ',f10.5,' 2xJ   = ',i4) 

       denmat(:,:,:) = 0.d0  
	   if(i==j)then
		   sameif=.true.
	   else
		   sameif=.false.
	   end if
       call master_density1b_bundled(sameif)
        do jt = jmin,jmax
              call coupled_densities(Jt,Ji,Ti,Jf,Tf,Jz,npeff(1)-npeff(2),  & 
                   zeroflag, numorbmax ,denmat)
				   densitybag(i,j)%zeroflag(jt)=zeroflag
				   densitybag(j,i)%zeroflag(jt)=zeroflag
				   
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
!                 write(resultfile,31)Jt
31               format(' Jt = ',i3,', Tt = 0        1 ' )
              else
!                 write(resultfile,32)Jt
32               format(' Jt = ',i3,', proton      neutron ')

              endif
              do m = 1,numorbmax
                 do n = 1,numorbmax
                    if ( (denmat(m,n,0) /=  0.0 .or. denmat(m,n,1) /=  0.0) & 
                         .and. iproc == 0 )then
						 densitybag(i,j)%denmat(jt,m,n,0)=denmat(m,n,0)
						 densitybag(i,j)%denmat(jt,m,n,1)=denmat(m,n,1)
!............. USE SYMMETRY.............................................
                         if(i/=j)then
							 ja = (orbqn(1,m)%j+1)/2
							 jb = (orbqn(1,n)%j+1)/2
							 
							 phase=(-1)**((Ji-Jf)/2 +ja-jb)
							 if(isoflag.and..not.pndensities) phase=phase*(-1)**((Ti-Tf)/2)
							 densitybag(j,i)%denmat(jt,n,m,0)=phase* denmat(m,n,0)
							 densitybag(j,i)%denmat(jt,n,m,1)=phase* denmat(m,n,1)
							 
						 end if						 
						 
!                       write(resultfile,36)m,n,(denmat(m,n,tt),tt=0,1)
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
     
  call clocker('den','end')

end subroutine density1b_compute

!=========================================================================
!
!  write out one-body densities at the end
!  NEW VERSION 7.8.2:
!     
!
! CALLED BY:
!   lanczos_output
!   exactdiag_p
!

subroutine density1b_writeout(istart,istop,fstart,fstop)
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
  use pocc_mod
  
  implicit none
  integer(4),intent(in) :: istart,istop,fstart,fstop ! start, stop for initial, final states

  real(4) :: xj,xt,ei,ef,xtt,xjj
  integer(4) :: i,j,m,n
  integer(4) :: ji,ti,jf,tf,jt,tt,jmin,jmax
  logical :: evenAJ,evenAT    ! needed to force "correct" J, T

  integer :: ierr
!  integer :: aerr
!--------------------- CONVERGENCE -----------------------------

  logical smallflag,zeroflag
  
  real(8), allocatable :: denmat(:,:,:)
  real :: dentol = 1.0e-6
  real :: nparticles
  logical numberflag  ! flag to check on total number of particles; used for debugging
  integer inode
 ! character(1) :: vchar
  real(4) :: nxprot,nxneut
  real(4) :: xtz

  if(iproc/=0)return
    
  call clocker('obw','sta')

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
  do i = istart,istop
        ei = energy(i)
        xj = xjlist(i)  
        xtt= xtlist(i)   
        ji = closest2J(evenAJ,xjlist(i))
        if(isoflag)ti = closest2J(evenAT,xtlist(i))
  	   if(only_odd .and. (ji/4)*4==ji)cycle		

     do j = fstart,fstop

           ef = energy(j)

           jf = closest2J(evenAJ,xjlist(j)) !         
           if(isoflag)then
              tf = closest2J(evenAT,xtlist(j)) ! 
           endif
	  	   if(only_odd .and. (jf/4)*4==jf)cycle
		   
        jmax = (jf + ji)/2
        jmin = abs( jf -ji)/2
			
        if (  isoflag) then
              write(denresultfile,*)' '
              write(denresultfile,333)i,ei,Ji,Ti 
              write(denresultfile,334)j,ef,Jf,Tf
			  if(i==j)then
				  write(occresultfile,*)' '
				  write(occresultfile,335)i,ei,Ji,Ti 
			  end if
        end if
        if ( .not.isoflag) then
              write(denresultfile,*)' '
              write(denresultfile,433)i,ei,Ji
              write(denresultfile,434)j,ef,Jf
			  if(i==j)then
				  write(occresultfile,*)' '
				  write(occresultfile,336)i,ei,Ji
			  end if
			  
        end if
333     format(' Initial state #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
334     format(' Final state   #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
335     format(' State #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
336     format(' State   #',i5,' E = ',f10.5,' 2xJ = ',i4) 

433     format(' Initial state #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
434     format(' Final state   #',i5,' E = ',f10.5,' 2xJ   = ',i4) 


!---------------------------------------------------------------

!---------------------- WRITE OUT DENSITY MATRIX ELEMENTS
             do jt = jmin,jmax
                 if(densitybag(i,j)%zeroflag(jt))cycle

              if(isoflag .and. .not.pndensities)then
                 write(denresultfile,31)Jt
31               format(' Jt = ',i3,', Tt = 0        1 ' )
              else
                 write(denresultfile,32)Jt
32               format(' Jt = ',i3,', proton      neutron ')

              endif
              do m = 1,numorbmax
                 do n = 1,numorbmax
					 
					if(orbqn(1,m)%j + orbqn(1,n)%j < 2*Jt .or. abs(orbqn(1,m)%j - orbqn(1,n)%j) > 2*Jt)cycle
                    if ( abs(densitybag(i,j)%denmat(jt,m,n,0)) > dentol .or.  & 
					     abs(densitybag(i,j)%denmat(jt,m,n,1)) > dentol)then					 
                       write(denresultfile,36)m,n,(densitybag(i,j)%denmat(jt,m,n,tt),tt=0,1)
                    end if
36                  format(2i5, 2f12.7)
                 end do ! b
				 if(i==j .and. jt ==0)then
                     xjj = 0.5* orbqn(1,m)%j
					 
					 if(isoflag .and. .not. pndensities)then
						 xtz = 0.5*(npeff(1)-npeff(2))
						 nxprot = densitybag(i,j)%denmat(jt,m,m,0) 
						 if(ti >0 )nxprot = nxprot + sqrt(3.)*xtz/sqrt(0.5*Ti*(0.5*Ti+1.0))* densitybag(i,j)%denmat(jt,m,m,1)
						 nxneut = densitybag(i,j)%denmat(jt,m,m,0)
						 if(ti> 0)nxneut=nxneut - sqrt(3.)*xtz/sqrt(0.5*Ti*(0.5*Ti+1.0))* densitybag(i,j)%denmat(jt,m,m,1)
						 
						 nxprot = nxprot* sqrt(2.)*sqrt(2*xjj+1)*0.5/sqrt((Ji+1.0)*(Ti+1.0))
						 nxneut = nxneut* sqrt(2.)*sqrt(2*xjj+1)*0.5/sqrt((Ji+1.0)*(Ti+1.0))
						 write(occresultfile,701)m,nxprot,nxneut
					 else
					 
					 write(occresultfile,701)m,(densitybag(i,j)%denmat(jt,m,m,tt)*sqrt(2*xjj+1.0)/sqrt(jf+1.0),tt=0,1)
				     end if
				 end if
701 format(i4,2f10.4)
				 
				 
              end do   ! a
           end do ! jt

        end do   !j

  end do  ! i
     
  if(iproc == 0)then
     write(denresultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
     write(denresultfile,*)' '
     write(denresultfile,*)' Definition of density matrices : '
     write(denresultfile,*)' rho_K(a^+b) =   (Jf || (a^+ b)_K || Ji) / sqrt(2K+1) '
     write(denresultfile,*)'  where reduced matrix element is convention of Edmonds '
     write(denresultfile,*)' (note: if isospin is good symmetry, then '
     write(denresultfile,*)'  doubly reduced/ divided by sqrt(2T+1) as well'
     write(denresultfile,*)' '
     write(denresultfile,*)' Note time-reversal symmetry relation: '
     write(denresultfile,*)' rho_K(a^+b, Ji->Jf) = (-1)^(ja-jb + Ji-Jf) rho_K(b^+a, Jf->Ji) '
     write(denresultfile,*)' For isospin, add in factor (-1)^(Ti-Tf)'
     write(denresultfile,*)' '
     write(denresultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
     write(denresultfile,*)' '
	 write(occresultfile,*)' '

  end if
  call clocker('obw','end')

end subroutine density1b_writeout


!
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
  use obs,only:parity,energy,xjlist,xtlist
  use writebinden1b
  implicit none

  real(4) :: xj,xt,ei,ef,xtt,xjj
  real(4) :: xt2
!  integer :: nkeep
  integer(4) :: i,j,m,n
  integer(4) :: ji,ti,jf,tf,jt,tt,jmin,jmax
  integer(kind=basis_prec) :: k
  logical :: evenAJ,evenAT    ! needed to force "correct" J, T
   
  logical :: zeroflag
  
  real(8), allocatable :: denmat(:,:,:)
  real :: nparticles,protonsum,neutronsum
  logical :: numberflag  ! flag to check on total number of particles; used for debugging

  real(kind=4), allocatable :: stateE(:), stateJ(:), stateT(:)
  real(kind=4) :: me ! matrix element
  logical :: sameif
  integer :: istatestart,istatestop,fstatestart,fstatestop  ! options to start and stop
  integer :: ierr
  logical, allocatable :: write_select(:)
  integer :: ipar
  integer :: ja,jb,phase
  
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
 allocate(energy(nkeep))
 allocate(xjlist(nkeep),xtlist(nkeep))
!--- ADDED in 7.9.: choose a range

if(iproc==0)then
   if(auto_input)then
	  read(autoinputfile,*)istatestart,istatestop
      if(.not.binden)then
		  read(autoinputfile,*)fstatestart,fstatestop
	  else
		  fstatestart = istatestart
		  fstatestop  = istatestop
	  end if
   else
	   
	  if(binden)then
	      print*,' Enter start, stop for states '	   	
	      print*,' (Enter 0,0  to read all )'
	      read*,istatestart,istatestop		  
		  fstatestart = istatestart
		  fstatestop  = istatestop		  
		  
	  else 
         print*,' Enter start, stop for initial states '	   	
         print*,' (Enter 0,0  to read all )'
         read*,istatestart,istatestop
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

	     print*,' Densities for: initial states from, to ',istatestart,istatestop
	     print*,' Densities for: final states from, to ',fstatestart,fstatestop
   else		
	     write(autoinputfile,*)istatestart,istatestop," ! start/stop for initial states for 1b dens "
   	     write(autoinputfile,*)fstatestart,fstatestop," ! start/stop for final states for 1b dens "
   end if

end if 

#ifdef _MPI
call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
! call OLDMPI_BCAST(ychar,1,MPI_CHARACTER,0,icomm,ierr)
call BMPI_BCAST(istatestart,1,0,MPI_COMM_WORLD,ierr)
call BMPI_BCAST(istatestop,1,0,MPI_COMM_WORLD,ierr)
call BMPI_BCAST(fstatestart,1,0,MPI_COMM_WORLD,ierr)
call BMPI_BCAST(fstatestop,1,0,MPI_COMM_WORLD,ierr)
#endif

allocate( densitybag(nkeep,nkeep), stat=aerr)
if(aerr /= 0) then
   call memerror("density1b_compute 4")
   stop 5
end if
densitybag(:,:)%filled = .false.

!------------- ADDED IN 7.9.8 need to write out energies

allocate(write_select(nkeep),e(nkeep))
write_select = .false.

do i = istatestart,istatestop
	write_select(i)=.true.
end do
do i = fstatestart,fstatestop
	write_select(i)=.true.
end do

if(iproc==0)write(resultfile,*)' '

call pocc_write_table_header() !  E   Ex ...


if(.not.write_select(1))then
	call wfn_readeigenvec(oldwfnfile,frag2, fcomm2_index, vec2, 1, ei,xjj,xtt)
	e(1)=ei
end if

!if(print_parity)then
	allocate(parity(nkeep))
!end if

do i = 1, nkeep
	
	if(write_select(i))then
        call wfn_readeigenvec(oldwfnfile,frag2, fcomm2_index, vec2, i, ei,xjj,xtt)
		e(i)=ei
		energy(i)=ei
		xjlist(i)=xjj
		xtlist(i)= 0.5*closest2J(evenAT,-0.5 + sqrt(xtt+0.25))
        xt = real( -0.5 + sqrt(xtt + 0.25), kind=4)
		call find_parity(ipar,2)
		parity(i)=ipar
		if(print_parity)then
!			write(6,111)i,e(i),e(i)-e(1),xjj,xt,ipar
 		   write(resultfile,111)i,e(i),e(i)-e(1),xjj,xt,ipar
 		   write(occresultfile,111)i,e(i),e(i)-e(1),xjj,xt,ipar
			
		else
!     	  write(6,11)i,e(i),e(i)-e(1),xjj,xt
   	      write(resultfile,11)i,e(i),e(i)-e(1),xjj,xt
 	      write(occresultfile,11)i,e(i),e(i)-e(1),xjj,xt
		   
	    end if
		
	end if
	
end do

if(iproc==0)then
	write(resultfile,*)' '
	write(occresultfile,*)' '	
	call pocc_write_orbits_alt(resultfile)
	call pocc_write_orbits_alt(occresultfile)	
end if
	  
  do i = istatestart,istatestop

     ! new interface, we say which vec to read.  It checks
     call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,ei,xj,xtt)


     ji = closest2J(evenAJ,xj)
	   if(only_odd .and. (ji/4)*4==ji)then
		   cycle
	   end if
     ! if !isoflag, then what?  orig code didn't calculate ti
     ti = closest2J(evenAT,-0.5 + sqrt(xtt+0.25))
     stateE(i) = ei
     stateJ(i) = ji
     stateT(i) = ti
     do j = fstatestart,fstatestop		
		   if(densitybag(i,j)%filled)cycle
           ! new interface, we say which vec to read.  It seeks and checks		   
           call wfn_readeigenvec(oldwfnfile,frag2, fcomm2_index, vec2, j, ef,xjj,xtt)		
           jf = closest2J(evenAJ,xjj) ! nint(2*xjj) 
	  	   if(only_odd .and. (jf/4)*4==jf)cycle
		            
           ! if ! isoflag, then what
           tf = closest2J(evenAT,-0.5 + sqrt(xtt+0.25)) ! nint(-1 + sqrt(4*xtt+1))
           jmax = (jf + ji)/2
           jmin = abs( jf -ji)/2
	
		   !.......... NOW SET UP DERIVED TYPE FOR STORING DENSITY MATRICES..		
		   		densitybag(i,j)%jmin=jmin
		   		densitybag(i,j)%jmax=jmax
		   		densitybag(j,i)%jmin=jmin
		   		densitybag(j,i)%jmax=jmax
  			
		   		if(.not.allocated( densitybag(i,j)%denmat )) then
					allocate(densitybag(i,j)%denmat(jmin:jmax,numorbmax,numorbmax,0:1), stat=aerr )
 		           if(aerr /= 0) then
 					   print*,' mem error in den1b from oldwfn, allocation denbag ',i,j,jmin,jmax,numorbmax
 		               call memerror("density1b_compute 701")
 		               stop 5
 		           end if
			    end if
				if(i/=j .and. .not.allocated( densitybag(j,i)%denmat))then
					 allocate(densitybag(j,i)%denmat(jmin:jmax,numorbmax,numorbmax,0:1), stat=aerr )
		            if(aerr /= 0) then
					   print*,' mem error in den1b from oldwfn, allocation denbag ',i,j,jmin,jmax,numorbmax
		               call memerror("density1b_compute 70")
		               stop 5
		            end if
				end if

		   		densitybag(i,j)%denmat=0.d0
		   		densitybag(j,i)%denmat=0.d0

		   		if(.not.allocated( densitybag(i,j)%zeroflag ))then
					allocate(densitybag(i,j)%zeroflag(jmin:jmax), stat=aerr )
 		           if(aerr /= 0) then
 		               call memerror("density1b_compute 70bb")
 		               stop 5
 		           end if
 			   end if
			   
		   		if(i/=j .and. .not.allocated( densitybag(j,i)%zeroflag ))then
					allocate(densitybag(j,i)%zeroflag(jmin:jmax), stat=aerr )
				
		           if(aerr /= 0) then
		               call memerror("density1b_compute 70bb")
		               stop 5
		           end if
			   end if
           if(np(1) > 0) p1bopme(:,:) = 0.0
           if(np(2) > 0) n1bopme(:,:) = 0.0
           denmat(:,:,:) = 0.d0
              ! compute density matrix			  
	   	   if(i==j)then
	   		   sameif=.true.
	   	   else
	   		   sameif=.false.
	   	   end if

           call master_density1b_bundled(sameif)
		   
!		   if(iproc==0)then
              do jt = jmin,jmax
                 call coupled_densities(Jt,Ji,Ti,Jf,Tf,Jz,npeff(1)-npeff(2),  & 
                   zeroflag, numorbmax ,denmat)
				   densitybag(i,j)%zeroflag(jt)=zeroflag
				   densitybag(j,i)%zeroflag(jt)=zeroflag
				   
!--------------------- CHECK PARTICLE OCCUPATIONS------------------
                 if(iproc==0)then
                 if(numberflag .and. i == j .and. Jt == 0 )then
                    if(isoflag .and. .not.pndensities)then
                       nparticles = 0.0
					   protonsum  = 0.0
					   neutronsum = 0.0
                       do n = 1,numorbmax
                          xjj = 0.5* orbqn(1,n)%j
                          nparticles = nparticles + (denmat(n,n,0))*sqrt( 2*(2*xjj+1))
 !                         print *, "n=", n, ",  xjj=", xjj, ", denmat(n,n,0)=", denmat(n,n,0), ", npart=", nparticles
                          
                       end do  ! n
                       nparticles = nparticles /sqrt( float(jf+1)*float(tf+1))
!                       print *, "KSM: jf=", jf,", jt=", jt,", tf=", tf, ", nparticles=", nparticles
                       if( abs( nparticles -npeff(1) -npeff(2)) > 0.001 .and. iproc==0)then
                             print*,nparticles,' particles total for state ',i
                       end if
                    else
						
                       nparticles = 0.0
					   protonsum  = 0.0
					   neutronsum = 0.0
                       do n = 1,numorbmax
                          xjj = 0.5* orbqn(1,n)%j
                          protonsum = protonsum + (denmat(n,n,0))*sqrt( (2*xjj+1))
                          neutronsum = neutronsum + (denmat(n,n,1))*sqrt( (2*xjj+1))
!                          write(occresultfile,701)n,denmat(n,n,0)*sqrt(2*xjj+1.0)/sqrt(jf+1.0) & 
!						   , denmat(n,n,1)*sqrt(2*xjj+1.0)/sqrt( real(jf+1,8))
                       end do  ! n
                       protonsum = protonsum /sqrt( float(jf+1))
                       neutronsum = neutronsum /sqrt( float(jf+1))
					   
                       if ( iproc == 0 )print*,protonsum,' protons total, ',neutronsum,' neutrons total'

                    endif
					
                 endif
			 endif ! iproc==0
				 
				 701 format(i4,2f10.4)
!---------------------------------------------------------------

                 if(zeroflag)cycle
			     densitybag(i,j)%filled = .true.
			     densitybag(j,i)%filled = .true.				 
                 if(iproc == 0)then			 
!---------------------- WRITE OUT DENSITY MATRIX ELEMENTS
                 if(isoflag .and. .not. pndensities)then
                       call dm_add_tab(i, j, Jt, numorbmax)
!                       write(resultfile,31)Jt
31                     format(' Jt = ',i3,', Tt = 0        1 ' )
                 else
                       call dm_add_tab(i, j, Jt, numorbmax)
!                       write(resultfile,32)Jt
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
!                          write(resultfile,36)m,n,(denmat(m,n,tt),tt=0,1)
						  densitybag(i,j)%denmat(Jt,m,n,0)=denmat(m,n,0)
						  densitybag(i,j)%denmat(Jt,m,n,1)=denmat(m,n,1)
						  
!
!............. USE SYMMETRY.............................................
                         if(i/=j )then
							 
							 
							 ja = (orbqn(1,m)%j+1)/2
							 jb = (orbqn(1,n)%j+1)/2
							 
							 phase=(-1)**((Ji-Jf)/2 +ja-jb)
							 if(isoflag.and..not.pndensities) phase=phase*(-1)**((Ti-Tf)/2)
							 densitybag(j,i)%denmat(jt,n,m,0)=phase* denmat(m,n,0)
							 densitybag(j,i)%denmat(jt,n,m,1)=phase* denmat(m,n,1)
							 
						 end if								  						  
                       endif
36                     format(2i5, 2f12.7)
                    end do ! b
                 end do   ! a
			     end if ! iproc = 0
              enddo ! jt
!		 end if ! iproc == 0
     enddo   !j
  enddo  ! i
  if(binden)then
	  call binary_density_boss(istatestart,istatestop)
	  
  else
	  denresultfile =resultfile
	  call density1b_writeout(istatestart,istatestop,fstatestart,fstatestop)
	  
  end if
  if( menu_dx_omathematica )then
     call density1b_write_math(nkeep, stateE, stateJ, stateT)
     call dm_reset()
   end if
  return
11 format(i5,3x,2f10.5,2x,2f8.3)
111 format(i5,3x,2f10.5,2x,2f8.3,3x,i2)
  
end subroutine density1b_from_oldwfn

!====================================================
! old version--to be deleted ----

!subroutine density1b_from_oldwfn_old DELETED

!=========================================================================
!
! subroutine UNCOUPLED DENSITY1B_FROM_OLDWFN
!  added 7.10.5
!  modification of routine density1b_from_oldwfn
!
! computes uncoupled density from pre-existing wavefunction file
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
subroutine uncoupled_density1b_from_oldwfn
  use menu_choices
  use system_parameters
  use sporbit
  use spstate
  use localvectors
  use nodeinfo
  use io
  use basis
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
  
  real(8), allocatable :: denmat(:,:,:)
  real :: nparticles
  logical :: numberflag  ! flag to check on total number of particles; used for debugging

  real(kind=4), allocatable :: stateE(:), stateJ(:), stateT(:)
  real(kind=4) :: me ! matrix element
  logical :: sameif
  integer :: istatestart,istatestop,fstatestart,fstatestop  ! options to start and stop
  integer :: ierr
  logical, allocatable :: write_select(:)
  integer :: a,b
  integer :: asps,bsps,asgn,bsgn,ath,bth,ia,ib
  integer :: it
  integer :: ja,ma,mb
  
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

!--- ADDED in 7.9.: choose a range

if(iproc==0)then
   if(auto_input)then
	  read(autoinputfile,*)istatestart,istatestop
      read(autoinputfile,*)fstatestart,fstatestop
   else

      print*,' Enter start, stop for initial states '	   	
      print*,' (Enter 0,0  to read all )'
      read*,istatestart,istatestop

      print*,' Enter start, stop for final states '
      print*,' (Enter 0,0  to read all )'	
      read*,fstatestart,fstatestop

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

	     print*,' Densities for: initial states from, to ',istatestart,istatestop
	     print*,' Densities for: final states from, to ',fstatestart,fstatestop
else		
	     write(autoinputfile,*)istatestart,istatestop," ! start/stop for initial states for 1b dens "
   	     write(autoinputfile,*)fstatestart,fstatestop," ! start/stop for final states for 1b dens "
   end if

end if 

#ifdef _MPI
call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
call BMPI_BCAST(istatestart,1,0,MPI_COMM_WORLD,ierr)
call BMPI_BCAST(istatestop,1,0,MPI_COMM_WORLD,ierr)
call BMPI_BCAST(fstatestart,1,0,MPI_COMM_WORLD,ierr)
call BMPI_BCAST(fstatestop,1,0,MPI_COMM_WORLD,ierr)
#endif

!------------- ADDED IN 7.9.8 need to write out energies

allocate(write_select(nkeep),e(nkeep))
write_select = .false.

do i = istatestart,istatestop
	write_select(i)=.true.
end do
do i = fstatestart,fstatestop
	write_select(i)=.true.
end do

if(iproc==0)write(resultfile,*)' '

call pocc_write_table_header() !  E   Ex ...

if(.not.write_select(1))then
	call wfn_readeigenvec(oldwfnfile,frag2, fcomm2_index, vec2, 1, ei,xjj,xtt)
	e(1)=ei
end if


do i = 1, nkeep
	if(write_select(i))then
        call wfn_readeigenvec(oldwfnfile,frag2, fcomm2_index, vec2, i, ei,xjj,xtt)
		e(i)=ei
        xt = real( -0.5 + sqrt(xtt + 0.25), kind=4)
		write(resultfile,11)i,e(i),e(i)-e(1),xjj,xt
	end if
end do

!........ write out uncoupled s.p. states
if(iproc==0)then
	write(resultfile,'("#Single particle states: Protons ")')
	it = 1
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

       a = hspsqn(ath,ia)%orb
       ja= hspsqn(ath,ia)%j
       ma= hspsqn(ath,ia)%m
	   write(resultfile,*)asps,a,ja,ma
   end do
write(resultfile,'("#Single particle states: Neutrons ")')
it = 2
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

      a = hspsqn(ath,ia)%orb
      ja= hspsqn(ath,ia)%j
      ma= hspsqn(ath,ia)%m
   write(resultfile,*)asps,a,ja,ma
  end do	
	
	
end if

if(iproc==0)then
	write(resultfile,'("# Range of initial")')
	write(resultfile,*)istatestart,istatestop
	write(resultfile,'("# Range of final")')
	write(resultfile,*)fstatestart,fstatestop
	write(resultfile,*)' '
	write(resultfile,'("# Initial   Final ")')

end if
	  
  do i = istatestart,istatestop
     ! new interface, we say which vec to read.  It checks
	 	 
     call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,ei,xj,xtt)

     ji = closest2J(evenAJ,xj)
     ! if !isoflag, then what?  orig code didn't calculate ti
     ti = closest2J(evenAT,-0.5 + sqrt(xtt+0.25))
     stateE(i) = ei
     stateJ(i) = ji
     stateT(i) = ti
	 
     do j = fstatestart,fstatestop		
           ! new interface, we say which vec to read.  It seeks and checks
           call wfn_readeigenvec(oldwfnfile,frag2, fcomm2_index, vec2, j, ef,xjj,xtt)
		   
           jf = closest2J(evenAJ,xjj) ! nint(2*xjj)   
		          
           ! if ! isoflag, then what
           tf = closest2J(evenAT,-0.5 + sqrt(xtt+0.25)) ! nint(-1 + sqrt(4*xtt+1))
           jmax = (jf + ji)/2
           jmin = abs( jf -ji)/2
           if ( iproc == 0) then
              write(resultfile,*)' '
			  write(resultfile,*)i,j
           end if

           if(np(1) > 0) p1bopme(:,:) = 0.d0
           if(np(2) > 0) n1bopme(:,:) = 0.d0
              ! compute density matrix			  
	   	   if(i==j)then
	   		   sameif=.true.
	   	   else
	   		   sameif=.false.
	   	   end if
		   
           call master_density1b_bundled(sameif)
		   if(iproc==0)then
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

			      a = hspsqn(ath,ia)%orb
			      ja= hspsqn(ath,ia)%j
			      ma= hspsqn(ath,ia)%m
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

			         mb= hspsqn(bth,ib)%m
			         if( ma /= mb) cycle
					 if(it==1)then
	  				   write(resultfile,'(2i4,f12.7)')asps,bsps,p1bopme(ia,ib)
						 
					 else
  	  				   write(resultfile,'(2i4,f12.7)')asps,bsps,n1bopme(ia,ib)
					 end if
					 
					 
				 end do
			 end do
			 write(resultfile,'(2i4,f12.6)')-1,-1,0.0
			 
		     end do


		   end if ! iproc == 0

     enddo   !j
  enddo  ! i
  if(iproc == 0)then
	  print*,' ++++++++++++++++++++++++++++'
	  print*,' Output written to .udres file '
	  print*,' +++++++++++++++++++++++++++'
	  
	  
 

  end if

  return
11 format(i5,3x,2f10.5,2x,2f8.3)
  
end subroutine uncoupled_density1b_from_oldwfn

!=============================================================
!
! master routine
!
! subroutines called:
!   applyP1Bdenbundled_g
!   applyN1Bdenbundled_g
!   phconjdensity
!
  subroutine master_density1b_bundled(sameif)

  use system_parameters
  use precisions
  use densities
  use nodeinfo
  use opbundles
  use bmpi_mod
  use haiku_info
  use spstate,only:nhspsmax
  implicit none
  logical :: sameif
  integer :: asps,ia,ath,it,asgn,bsps,ib,bth,bsgn,ierr
  real(kind=obs_prec) :: xme
  
!........ ADDED IN 7.8.4 to reduce matrix element
!         KLUGEY because I couldn't get reduce to work correctly
!
  real(kind=obs_prec), allocatable :: x1bopmevec(:),y1bopmevec(:)
  integer :: sizex1bopmevec,ivec
   
  sizex1bopmevec = (2*nhspsmax+1)**2
  if(.not.allocated(x1bopmevec))allocate(x1bopmevec(sizex1bopmevec))
  if(.not.allocated(y1bopmevec))allocate(y1bopmevec(sizex1bopmevec))

  p1bopme(:,:) = 0.d0
  n1bopme(:,:) = 0.d0

  call applyP1Bdenbundled_g(opbundlestart(iproc),opbundleend(iproc))
  call applyN1Bdenbundled_g(opbundlestart(iproc),opbundleend(iproc))

!............ CHECK FOR PARTICLE-HOLE CONJUGATION AND SWAP IF NECESSARY.....
!   FOR MPI, MOVED AFTER WARDS
!  if(phconj(1).and. sameif)call phconjdensity(1)
!  if(phconj(2).and.sameif)call phconjdensity(2)
!if(phconj(1))call phconjdensity(1,sameif)
!if(phconj(2))call phconjdensity(2,sameif) 
!............... NOW REDUCE .........................  
if(nproc > 1)then	
	do it = 1,2
		if(np(it)< 1)cycle
	    x1bopmevec(:)=0.0
		y1bopmevec(:)=0.0
		
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
			  
			  ivec = (ia*asgn+nhspsmax)*(2*nhspsmax+1)+ib*bsgn+nhspsmax+1
			  if(it==1)x1bopmevec(ivec)=p1bopme(ia*asgn,ib*bsgn)
			  if(it==2)x1bopmevec(ivec)=n1bopme(ia*asgn,ib*bsgn)

		  
	  end do  ! bsps
	  
  end do  ! asps

#ifdef _MPI
    call BMPI_REDUCE_REAL8_VECIP(x1bopmevec, sizex1bopmevec, MPI_SUM,0, MPI_COMM_WORLD, ierr) ! in place
#endif
!
! NOTE if obs_prec = 4 then must use the followin	
	
!    call BMPI_REDUCE_REAL4_VECIP(x1bopmevec, sizex1bopmevec, MPI_SUM,0, icomm, ierr) ! in place
		!end if
	if(iproc==0)then
	do ia = -nhspsmax,nhspsmax
		
		do ib = -nhspsmax,nhspsmax
		  ivec = (ia+nhspsmax)*(2*nhspsmax+1)+ib+nhspsmax+1
		  if(it==1)  p1bopme(ia,ib)=x1bopmevec(ivec)
		  if(it==2)  n1bopme(ia,ib)=x1bopmevec(ivec)

			
		end do 
	end do
   end if
end do   ! it
end if
if(phconj(1))call phconjdensity(1,sameif)
if(phconj(2))call phconjdensity(2,sameif) 
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
   real(8) :: dnorm
   logical :: smallflag

   call onebodysetup

   call wfn_read_nkeep(oldwfnfile, nkeep)

   if(iproc==0)then
      print *, ' '
      print *, nkeep, ' states '
	  if(strengthnorm)print*,' Normalizing output vectors '
   end if
   call wfn_write_nkeep(nkeep) ! write number of vectors to wfn file
   
   
   do i = 1,nkeep
      call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,e,xj,xt2) ! KSM: updated
      call master_apply1bop        ! new version added in 7.6.7
	  
	  if(strengthnorm)then  ! normalize
	      call dnormvec_p('n','f',dnorm,smallflag)
		  
	  end if

      if(writeout) call wfn_writeeigenvec(wfnfile,frag2, vec2,i,e,xj,xt2)
   end do

   return
   end subroutine applicator1b
!============================================================
   subroutine applicator1b_print
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
   
   real, allocatable :: onebdenmat(:,:)
   
   print*,' '

   if(.not.printouttrans1flag)then
	   print*,' Bad flag printouttrans1bflag  TRUE '
	   return
   end if
   call onebodysetup
   
   allocate(onebdenmat(dimbasis,dimbasis))
   if(iproc==0)open(unit=61,file='trans1b.dat',status='unknown')
   
   write(61,*)dimbasis

   do i = 1,dimbasis
!      call wfn_readeigenvec(oldwfnfile, frag1, fcomm1, vec1,i,e,xj,xt2) ! KSM: updated
      vec1(:) = 0.0
      vec1(i)=1.0
	  vec2(:)=0.0
      call master_apply1bop        ! new version added in 7.6.7
      do j = 1,dimbasis
         onebdenmat(i,j) = real(vec2(j),kind(egv_prec))   ! store as real(4); this can be changed

      end do  !J 
	 do j = 1,i
		 if(onebdenmat(i,j)/=0.0)write(61,*)i,j,onebdenmat(i,j)
	 end do

   end do
   close(61)
   
   print*,' '
   print*,' Written to trans1b.dat '
   print*,' '

   return
   end subroutine applicator1b_print

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
	use ntuple_info
	use nodeinfo
	use fragments, only: nodal
	use io,only:modeldensities
	use menu_choices,only:menu_char
	implicit none
	integer :: aerr
    call hopmaker
	
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
          call masterfindconjumps1b(1,.false.,.true.)
          call master1bodyjumps(1,.false.)
     endif
     if(np(2) > 0)then

  !--------------- find conjugate jumps
          call masterfindconjumps1b(2,.false.,.true.)
          call master1bodyjumps(2,.false.)
     endif

   
     call master_para_distribute_onebody   
	 
	 if(modeldensities)return
     call setup_localvectors  

   
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
!  MODIFIED in 7.8.4 to read in 'xpn' format for .opme files
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
  integer i,i1,i2
  integer a,b
  real xj
  real :: xxx,yyy
  logical :: xpnformat,nvals,pvals    ! ADDED in 7.8.4 to allow reading in xpn format opme
  integer :: neutron_offset ! ADDED in 7.8.4 to allow reading in xpn format opme
 
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
   
   xpnformat= .false.
   select case (title(1:3))
      case ('iso')
	  pnoperators = .false.
  
	  case ('pns')
	  pnoperators = .true.
	  xpnformat = .false.
	  
	  case ('xpn')
	  pnoperators = .true.
	  xpnformat = .true.
!	  print*,' Currently BIGSTICK does not accept .opme files in xpn format '
!	  print*,' Instead for proton-neutron breaking it must be in pns format, with two columns '
!	  print*,' See section 4.6.2 in the BIGSTICK Manual/Overview '
!	  stop
	
      case default
	  print*,'  Something wrong with first line of .opme file '
	  print*,title(1:3)
	  print*,title
	  stop
   end select
end if

!.......... SET UP OPERATOR ARRAYS.................

#ifdef _MPI
call BMPI_BCAST(pnoperators,1,0,MPI_COMM_WORLD,ierr)
#endif
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

   if(xpnformat)then
	   neutron_offset=numorb(1)
	   do i = 1, numorb(1)
	     read(23,*)i1,i2,n,l,xj
	     if(i1 /= i)then
	         print*,' ooopsie ',i1,i
	         stop
	     endif
	     if(i2 /= i+neutron_offset)then
	         print*,' ooopsie ',i2,i+neutron_offset
	         stop
	     endif		 
	     j = nint(2*xj)
	     if( n/= orbqn(1,i)%nr .or. l /= orbqn(1,i)%l .or. j /= orbqn(1,i)%j)then
	        print*,i,' Mismatch in orbits ',n,l,j
	        print*, orbqn(1,i)%nr, orbqn(1,i)%l, orbqn(1,i)%j
	        stop
	     endif

	   enddo	   
	   
   else
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
   end if
   
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
  do i = 1, 10000
	  if(pnoperators .and. .not. xpnformat)then
		  read(23,*,end=348)a,b,xxx,yyy
	  else
        read(23,*,end =348)a,b,xxx
	  end if
	  if(xpnformat)then
		  if(a <= neutron_offset .and. b<= neutron_offset)then
			  nvals = .false.
			  pvals = .true.
			  yyy = 0.
 
		  end if
		  
		  if((a > neutron_offset .and. b<= neutron_offset) .or. & 
		      (a<=neutron_offset .and. b > neutron_offset))then
			  print*,' A charge-changing matrix elements, skipped ',a,b
			  cycle
	  
		  end if
		  if(a > neutron_offset .and. b> neutron_offset)then
			  nvals = .true.
			  pvals = .false.
			  a = a -neutron_offset
			  b = b -neutron_offset	  
			  yyy = xxx
			  xxx = 0.0
		  end if

		  
	  else
		  nvals = .true.
		  pvals = .true.
		  
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
	
	!end if


	if(pnoperators)then
		if(np(1)>0 .and. pvals )pop1bod(a,b)=xxx
		if(np(2)>0 .and. nvals)nop1bod(a,b)=yyy
	else
        op1bod(a,b) = xxx
    end if
  enddo

348      continue

  close(unit=23)
end if   ! iproc == 0
!--------------- NOW BROADCAST -------------------------------
!    a bit of a kludge
#ifdef _MPI
call BMPI_BCAST(jop,1,0,MPI_COMM_WORLD,ierr)
call BMPI_BCAST(top,1,0,MPI_COMM_WORLD,ierr)
do a = 1,numorb(1)
	do b = 1,numorb(1)
		if(pnoperators)then
            if(np(1)> 0)call BMPI_BCAST(pop1bod(a,b),1,0,MPI_COMM_WORLD,ierr)
            if(np(2)> 0)call BMPI_BCAST(nop1bod(a,b),1,0,MPI_COMM_WORLD,ierr)
			
			
		else
           call BMPI_BCAST(op1bod(a,b),1,0,MPI_COMM_WORLD,ierr)
	   
       end if
    end do
end do
#endif
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
  
!================================================
!============================================================
!
!  CALLED BY main routine 
!
!  added in 7.10.5
!

  subroutine readin1bodyop_uncoupled
  use opmatrixelements
  use sporbit
  use nodeinfo
  use bmpi_mod
  use system_parameters
  use haiku_info
  use spstate
  use densities
  use onebodypot
  implicit none

  logical success
  character*15 :: appfile
  character*40 :: title
  integer ilastap
  integer i,i1,i2,j1,istate,istatemax
  integer :: istatestart,istatestop,fstatestart,fstatestop  ! options to start and stop
  
  real xj
  real :: xxx,yyy
  logical :: xpnformat,nvals,pvals    ! ADDED in 7.8.4 to allow reading in xpn format opme
  integer :: neutron_offset ! ADDED in 7.8.4 to allow reading in xpn format opme
  integer :: it
 
!-------------- DUMMIES FOR CONSISTENCY CHECK --------------
  integer n,l,j
  integer :: aerr,ierr  
  integer :: a,asps,asgn,ath,ja,ma,ia
  integer :: b,bsps,bsgn,bth,jb,mb,ib

!------------- READ IN OPERATOR FILE AND DECOUPLE

if(iproc == 0)then
  success = .false.
  do while(.not.success)  
     write(6,*)' Enter name of .udres  file '
     read(5,'(a)')appfile
     ilastap = index(appfile,' ') -1
     open(unit = 23,file=appfile(1:ilastap)//'.udres',status = 'old',err=301)
     success = .true.
     cycle
301  continue
     print*,'File ',appfile(1:ilastap),'.udres does not exist '
   enddo
!------- READ IN HEADER???
   read(23,'(a)')title
   print*,title
!......................... CHECK Z, N
   read(23,*)i1,i2
   if(i1/=np(1) .or. i2 /=np(2))then
	   print*,' Mismatch of Z, N '
	   print*,np(1),np(2),i1,i2
	   stop
   end if
   read(23,'(a)')title
   print*,title
   
   success=.false.
   do while (.not.success)
	   read(23,'(a)')title
	   if(title(3:5)=='Sta')success=.true.
   end do
   
   print*,' States '
   do i = 1, 10000
	   read(23,*,err=111)i1,xxx,yyy,xj
	   write(6,'(i5,2f12.5,f10.3)')i1,xxx,yyy,xj
	   istatemax = i1
   end do
   
!.................................. CHECK SPS   
111 continue   
   backspace(23)
   do it = 1,2
       success=.false.
	   do while (.not.success)
		   read(23,'(a)')title
		   if(title(2:5)=='Sing')success=.true.
	   end do 
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

          a = hspsqn(ath,ia)%orb
          ja= hspsqn(ath,ia)%j
          ma= hspsqn(ath,ia)%m
   	   read(23,*)i,i1,j1,i2
	   if(i/=asps .or. i1/=a .or. j1/=ja .or. i2/=ma)then
		   print*,' mismatch in s.p. states'
		   stop
	   end if
      end do	   
   end do
   
!....... READ IN LIMITS FOR INITIAL, FINAL....
   
   success=.false.
   do while (.not.success)
	   read(23,'(a)')title
	   if(title(3:5)=='Ran')success=.true.
   end do
   read(23,*)istatestart,istatestop
   read(23,'(a)')title
   read(23,*)fstatestart,fstatestop
   
!.........................      
   success=.false.
   do while (.not.success)
	   read(23,'(a)')title
	   if(title(3:5)=='Ini')success=.true.
   end do 
!.......... NOW TO CHOOSE INITIAL, FINAL STATES......

   success = .false.
   
   print*,' Enter the labels for initial and final states '
   print*,' Initial in the range ',istatestart,istatestop
   print*,' Final in the range   ',fstatestart,fstatestop
   
   
   do while(.not.success)
      read*,i1,i2
	  
	  if((i1 < istatestart .or. i1 > istatestop) .or. & 
	  (i2 < fstatestart .or. i2 > fstatestop))then
		  print*,' initial must be in range ',istatestart,istatestop
		  print*,' initial must be in range ',istatestart,istatestop
	  else
		  success=.true.
	  end if
   end do
   
   success = .false.
!............... FIND MATCHING SET OF INITIAL, FINAL STATES..............   
   do while (.not. success)
	   read(23,*)a,b
	   if(a==i1 .and. b==i2)then
		   success=.true.
	   else  ! read all the matrix elements
		   do it =1,2
			   do i = 1,10000
				   read(23,*,err=222)a,b,xxx
				   if(a==-1 .and. b==-1)exit
			   end do
		   end do
	   end if
   end do   
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
   allocate(pPOT_h(-nhsps(-1):nhsps(1), -nhsps(-1):nhsps(1)), stat=aerr)
   if(aerr /= 0) call memerror("decouple1bodyop 1")
   allocate(nPOT_h(-nhsps(-2):nhsps(2), -nhsps(-2):nhsps(2)), stat=aerr)
   if(aerr /= 0) call memerror("decouple1bodyop 2")

   pPOT_h = 0.0
   nPOT_h = 0.0
   
   
!................. NOW TO READ IN THE UNCOUPLED MATRIX ELEMENTS............  
p1bopme = 0.d0
n1bopme = 0.d0
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

     a = hspsqn(ath,ia)%orb
     ja= hspsqn(ath,ia)%j
     ma= hspsqn(ath,ia)%m
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

        mb= hspsqn(bth,ib)%m
        if( ma /= mb) cycle
		read(23,*)a,b,xxx
!		print*,a,b,xxx,asps,bsps,ma,mb
		if(a/=asps .or. b/=bsps)then
			print*,' Mismatch in reading matrix elements ',a,b,xxx,asps,bsps
			stop
		end if
		
	 if(it==1)then
	!   p1bopme(ia,ib)=xxx
	   pPOT_h(ia,ib)=xxx
		 
	 else
!  	   n1bopme(ia,ib)=xxx
  	   nPOT_h(ia,ib)=xxx

	 end if
	 
	 
     end do  ! bsps

 end do ! asps
 read(23,*)a,b
 if(a/=-1 .or.b/=-1)then
	 print*,' did not get to the end ',a,b
	 stop
 end if
 
 
end do
   
   
end if
close(23)

print*,' Note: MPI not yet implemented for readin1bodyop_uncoupled '
  return
  
222 continue

   print*,' some problem in reading '
   stop  
  end subroutine readin1bodyop_uncoupled  
  
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
    use nodeinfo, only : iproc,nproc
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
#ifdef _MPI
           call BMPI_REDUCE(vec2, size(vec2), MPI_SUM, 0,fcomm2, ierr) ! in place
#endif
		   call br_grab_vec2()
!           if(nproc > 1)call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
		   
    else
        if(nproc > 1 )then
#ifdef _MPI
           call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
#endif
        endif
    end if

    if(useNewReorthog) then

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
subroutine overlap(final_only)
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

   logical :: final_only ! added 7.9.10, to allow to take dot products with result of a green function
   
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
	  if(dotflag)then
         open(unit = overlapfn, file='overlap.dat',status='unknown')
      else
          open(unit = overlapfn, file='relentropy.dat',status='unknown')
	  end if
   end if
   dimbasis = dimbasischeck
   oldbasis = dimbasischeck
   
   if(final_only)go to 345  ! added 7.9.10
   
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
      call wfn_readeigenvec(oldwfnfile,frag1, fcomm1_index, vec1,i,xe,xj,xt2)  
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
   
#ifdef _MPI
   call BMPI_BCAST(tidx,1,0,MPI_COMM_WORLD,ierr)
#endif
   ! new interface - we say which vec to read, it seeks and reads
   ! no need to rewind and read forward
   call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,tidx,xe,xj,xt2)
   call wfn_close_file(oldwfnfile)
   
345 continue      ! skip here if final_only=T, used for green function
   
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

   if(iproc ==0 .and. dotflag)write(overlapfn,*)' state      E      J    T       <i|f>      |<i|f>|^2 '
   if(iproc ==0 .and. .not. dotflag)write(overlapfn,*)' state      E      J    T       S_rel(i|f)     '
   
   do i = 1,nkeep
       ! new interface, we say which vec to read.  It checks
       call wfn_readeigenvec(wfnfile, frag2, fcomm2_index, vec2,i,xe,xj,xt2)
       dsclrprod = 0.d0
	   call doverlapvec(dsclrprod,dotflag)  ! modified in 7.7.0 so as to work with fragments
       if(iproc == 0 .and. dotflag)then
          write(6,101)i,xe,xj,xt2,dsclrprod,dsclrprod**2
          write(overlapfn,101)i,xe,xj,xt2,dsclrprod,dsclrprod**2
       end if
       if(iproc == 0 .and. .not. dotflag)then
          write(6,101)i,xe,xj,xt2,dsclrprod
          write(overlapfn,101)i,xe,xj,xt2,dsclrprod
       end if
101 format(i4,3x,f10.5,2f5.1,2x,2f12.5)
   enddo
   if(iproc==0)close(unit = overlapfn)
   if(iproc ==0 .and. dotflag)write(6,*)' results written to overlap.dat '
   if(iproc ==0 .and. .not. dotflag)write(6,*)' results written to relentropy.dat '

   stop
  end subroutine overlap
!
!======================================================

!
! routine to compute overlaps of lanczos vectors with a file of vecteors

!
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
subroutine overlap_lanczos
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
   use mod_reorthog
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
         open(unit = overlapfn, file='overlap_lanczos.dat',status='unknown')
   end if
   dimbasis = dimbasischeck
   oldbasis = dimbasischeck
 
!   call overlaptribution
!   call setup_localvectors
   
   oldbasis = dimbasis
   if(iproc==0)then
      print*,' '
      print*,' Enter file of FINAL STATES '
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

   if(iproc ==0)write(overlapfn,*)' Lanczvec     Finalstate       overlap    '
!   if(iproc ==0 .and. .not. dotflag)write(overlapfn,*)' state      E      J    T       S_rel(i|f)     '
   
   do i = 1,nkeep
       ! new interface, we say which vec to read.  It checks
       call wfn_readeigenvec(wfnfile, frag2, fcomm2_index, vec2,i,xe,xj,xt2)
	   
!..........now loop over lanczos vectors
       do j = 1,niter	   
		   vec1=0.d0
   		  call br_retrieve_hist(j)
   		  call br_restore_vec1()
          dsclrprod = 0.d0
	      call doverlapvec(dsclrprod,.true.)  ! modified in 7.7.0 so as to work with fragments
          if(iproc == 0 )then
             write(6,101)j,i,dsclrprod
             write(overlapfn,101)j,i,dsclrprod
		  end if
       end do ! j
101 format(i9,5x,i5,5x,f12.5)
   enddo
   if(iproc==0)close(unit = overlapfn)
   if(iproc ==0 )write(6,*)' results written to overlap_lanczos.dat '

   stop
  end subroutine overlap_lanczos  
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
subroutine phconjdensity(it,sameif)
   use system_parameters
   use spstate
   use sporbit
   use densities
   use precisions
   use haiku_info
   implicit none
!........ INPUT...............
   integer :: it   ! species 
   logical :: sameif
!....... INTERNAL.......
   real(kind=obs_prec), pointer :: x1bopme(:,:)
   integer a,b,ja,ma,jb,mb
   integer ia,ib
   integer asps,ath,bsps,bth,asgn,bsgn
   integer :: iatr, ibtr       ! time-reversed states
   real(kind=obs_prec) :: diag
   real(kind=obs_prec) :: tmpden

   if(.not. phconj(it) .or. npeff(it) < 1 ) return
   
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
	  
	  if(orbqn(it,a)%w==99)cycle 

      do bsps = 1, nhsps(it)
         ib = bsps 
         bth = it
         bsgn = 1

         b = hspsqn(bth,ib)%orb
   	     if(orbqn(it,b)%w==99)cycle 
		 
         jb= hspsqn(bth,ib)%j
         mb= hspsqn(bth,ib)%m
     	 ibtr = hspsqn(bth,ib)%tr		 

         if( ma /= mb) cycle       ! note: one-body density should be empty unless ma == mb
!         if(asps /=bsps)then    !simply swap
         if(asps==bsps .and. sameif)then   ! only pick up a 1 if initial and final states are the same, and a=b
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
