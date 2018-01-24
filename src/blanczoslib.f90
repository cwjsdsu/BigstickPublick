!====================================================================
!  LANCZOS routines for BIGSTICK
!
!  versions for 'new' parallelization scheme -- FALL 2011
!
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

module lanczos_util
	
	use lanczos_info
	use precisions
	
	implicit none
	
    integer startiter   ! # of starting iterations already

    real(kind=4), allocatable :: alpha(:),beta(:)  ! lanczos coefficients
    real(kind=egv_prec), allocatable ::  eiglvec(:,:),e(:)
    logical :: flyeig

  !--------------------- CONVERGENCE -----------------------------------------
  
    logical :: finished
    integer :: ncheck
    real(kind=egv_prec), allocatable :: eold(:)  ! old energies
!    real    :: ediff0,ediff1

  !-------------- ALTERNATE CONVERGENCE ---------------------------
  !               based on wavefunctions
    logical :: altconverge

    real(kind=egv_prec), allocatable :: eiglvecold(:,:)
	
contains
	
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
!================================================
!
!  routine to find eigenvalues of lanczos matrix
!  if vecflag = .true. then find eigenvectors too
!  if nthick > 0  then arrange according to thick restart
!  nwindow = # of eigenvalues desired; only needed for targeted thick-restart
!
!  MODIFIED 7.6.8 so one can choose a target energy; REMOVED 7.7.4
!  

subroutine find_lanczos_eigenvalues(n,vecflag,nthick)
!	subroutine find_lanczos_eigenvalues(n,alpha,beta,e,vec,vecflag,nthick)

  use precisions
  use lanczos_info,only:Etarget,lanczchar
  implicit none

  integer(4):: n  ! actual dimension
!  integer :: lwork

!  real(kind = 4) :: alpha(np), beta(np)
!  real(kind = egv_prec) :: e(np), vec(np,np)
  logical :: vecflag  ! flag to get eigenvectors or not
  integer :: nthick
  real(kind = egv_prec) ::  work3(3*niter)  ! work(np)
  integer i,j,k
!  integer(4) :: nwindow, thickshift
  real(kind = egv_prec) :: etmp
  real(kind=egv_prec) ::  ediff(n),ediff0,xswap
  
  integer info

  eiglvec(:,:) = 0.0
  do i = 1,n
     eiglvec(i,i) = real(alpha(i),kind=egv_prec)
  enddo ! i
!-------------- put in column of thick-restart matrix elements

  if(nthick > 0)then
      do i = 1,nthick
         eiglvec(i,nthick+1) = real(beta(i),kind=egv_prec)
         eiglvec(nthick+1,i) = real(beta(i),kind=egv_prec)
      enddo
  endif
!-------- what about beta(nthick+1?)
  do i = nthick+1,n-1
     eiglvec(i,i+1) = real(beta(i),kind=egv_prec)
     eiglvec(i+1,i) = real(beta(i),kind=egv_prec)
  enddo
  call clocker('egv','sta')
  if(vecflag)then
      if(egv_prec==4)then
        call SSYEV( 'V','U', N, eiglvec, Niter, e, WORK3, 3*niter, INFO )
      else
        call DSYEV( 'V','U', N, eiglvec, Niter, e, WORK3, 3*niter, INFO )
      end if
  else
      if(egv_prec==4)then

         call SSYEV( 'N','U', N, eiglvec, Niter, e, WORK3, 3*niter, INFO )
      else
         call DSYEV( 'N','U', N, eiglvec, Niter, e, WORK3, 3*niter, INFO )

      end if
  endif
  call clocker('egv','end')
  
!--------- OPTION TO CHOOSE EXCITED STATES --------------------
!  added in 7.7.8 based in part upon work by R. Zbikowski
!  SORT EIGENVALUES, EIGENVECTORS BASED UPON PROXIMITY TO Etarget
!  
  if(lanczchar=='tx')then
!................  SET UP DIFFERENCES OF ENERGIES .............	
  	do i = 1,n
  		ediff(i)= abs(e(i)-etarget)
  	end do
!........... NOW SIMPLE BUBBLE SORT .......................	
      do i = 1,N-1
  		k = i
  		ediff0 = ediff(i)
  		do j = i+1,N
  			if(ediff(j) < ediff0)then
  				k=j
  				ediff0 = ediff(j)
  			end if
  		end do
  		if(k /= i)then
  			ediff(k)=ediff(i)
  			ediff(i)=ediff0
  			xswap   =e(k)
  			e(k)    =e(i)
  			e(i)    =xswap
  			do j = 1,N
  				xswap = eiglvec(j,k)
  				eiglvec(j,k)=eiglvec(j,i)
  				eiglvec(j,i)=xswap
  			end do
  		end if
      end do
	  
!.......... DO A SECOND SORT JUST ON THE Nkeep "LOWEST"........
!          THIS PUTS EVERYTHING INTO A "PRETTY" ORDER
do i = 1,Nkeep-1
   k = i
   ediff0 = e(i)
   do j = i+1,Nkeep
	   	if(e(j) < ediff0)then
		   k=j
		   ediff0 = e(j)
	    end if
   end do
   if(k /= i)then
 	  e(k)=e(i)
	  e(i)=ediff0
 	  do j = 1,N
		xswap = eiglvec(j,k)
		eiglvec(j,k)=eiglvec(j,i)
		eiglvec(j,i)=xswap
	  end do
   end if
end do
	  

  end if    
  
  return
end subroutine find_lanczos_eigenvalues

!====================================================================
!
! routine to control construction and "improvement" of initial vector
! added in 7.6.3 by CWJ @ SDSU
!
!  CALLED BY:  lanczos_menu
!  SUBROUTINES CALLED:
!  br_open_u
!  setup_localvectors  
!  intialize_lanczos_vector
!
subroutine setup_for_lanczos
	use flagger
    use mod_reorthog
	use io
	use bvectorlib_mod
	use onebodypot, only: meanie
    implicit none
!	integer :: startiter   ! number of iterations from previous run; only > 0 if restarting

    if(noisy0) print *, "setup_localvectors"
    call setup_localvectors
    if(noisy0) print *, "initialize_lanczos_vector"
    call initialize_lanczos_vector
	if(startiter == 0 .and. .not.strengthflag .and. .not.meanie) then
	   call dynamic_vector_initializer
    end if
    return
		
end subroutine setup_for_lanczos

!====================================================================
!
! routine to initialize Lanczos vector
!
!  NOTE: currently very kludgy!
!
!  INPUT:
!     startiter:  where iterations start
!              needed for restart option
!
!  OUTPUT:
!     dnorm0  = if read in pivot, what is magnitue
!     vchar   = keeps track of which vectors is initial, which final
!             'n' (normal) initial vec1, final vec2
!             'r' (reverse) initial vec2, final vec1
!
!  CALLED BY:
!    setup_for_lanczos
!  SUBROUTINES CALLED:
!    readpivot (old version)  reads in pivot from file, normalizes
!    dnormvec_p  : normalizes a vector (or a fragment)
!
subroutine initialize_lanczos_vector
   use flagger
   use localvectors
!   use precisions
   use basis
   use fragments
   use io
!   use lanczos_info
   use onebodypot
   use nodeinfo
   use wfn_mod
   use bvectorlib_mod
   use menu_choices
   implicit none 
   character(1) :: vchar
   real(kind=8) :: dnorm, dnorm0,da
   real :: rv
   integer startiter
   integer (kind = basis_prec) :: i
   integer (kind=4) :: iter
   logical smallflag
   integer :: pivotoption   ! different choices

   pivotoption = 3

   vchar = 'n'
   if ( strengthflag  .or. menu_char=='np') then

      call readpivot   ! FIXED so it works with fragments

	  return

   else if ( startiter == 0 ) then

      dnorm = dsqrt(real(dimbasis,kind=8))

      select case (pivotoption)

      case (0)
	  do i = v1s, v1e

!....... TO PREVENT ROUNDOFF ERRORS AND ALSO TO MINIMIZE
!        ACCIDENTAL STOPS, HAVE DECREASING WEIGHT
          vec1(i) = real(1.d0 /dnorm,kind=lanc_prec)
      end do
      case (1)
	  do i = v1s, v1e

!....... TO PREVENT ROUNDOFF ERRORS AND ALSO TO MINIMIZE
!        ACCIDENTAL STOPS, HAVE DECREASING WEIGHT
          call random_number(rv)
          vec1(i) = real( (1.d0+ 0.1*rv) /dnorm, kind=lanc_prec)
      end do

      case (2)
	  do i = v1s, v1e

!....... TO PREVENT ROUNDOFF ERRORS AND ALSO TO MINIMIZE
!        ACCIDENTAL STOPS, HAVE DECREASING WEIGHT

          vec1(i) = real(1.d0 /dsqrt(real(i,kind=8)+dnorm ), kind=lanc_prec)

       end do

       case (3)
	  do i = v1s, v1e

!....... TO PREVENT ROUNDOFF ERRORS AND ALSO TO MINIMIZE
!        ACCIDENTAL STOPS, HAVE DECREASING WEIGHT

          vec1(i) = real( 1.d0 * (-1)**i /dnorm, kind=lanc_prec)

       end do

       end select

   end if

!----------- RESTART OPTION -------------
!          NEEDS WORK!

  if ( startiter > 0 ) then

     if(iproc==0)rewind(lvec_file)
     if(.not. useNewReorthog) call swap_vchar(vchar)
     do iter = 1,startiter
		  if(.not. useNewReorthog) call swap_vchar(vchar)
        call read_lanczos_vector_restart_a(vchar,'i',iter,lvec_file)
     end do
  endif

   return
end subroutine initialize_lanczos_vector
!====================================================================

!  added in 7.6.3 March 2016
!  options to generate initial pivot vector; 
!  also useful in pointing towards a preconditioner in LOBPCG
!
! OPTIONS: 
!  1) PP, NN, + trace of PN
!  2) diagonal parts of PP, NN + trace PN only
!  3) diagonal parts of PP,NN, PN only
!  3) Hartree-Fock like initial vector
!

! SUBROUTINES CALLED:
!   masterpnsectortrace (optional)

subroutine dynamic_vector_initializer
!	use lanczos_info
	use nodeinfo
	use bmpi_mod
	use io
	use mod_reorthog
	use wfn_mod
	implicit none
	integer :: choice_initializing
	integer :: ierr
	integer :: tmpkeep
    integer i,j,n, tidx
    real :: e,xj,xt2
	
	nsave =  1  ! could be different, especially if creating vectors
	            ! to construct a preconditioner
	if(.not.initializing_enabled)then
		initializing=.false.
		applypntrace =.false.
		applyXXonly  =.false.
		diagonalsectorsonly = .false.
		write(autoinputfile,'(" ! Not optimizing initial pivot vector ")')
		return
	end if
	if(iproc==0)then
		print*,' '
		print*,' You get to dynamically prepare the pivot (initial vector) '
		print*,' Choose one of the following options: '
		print*,' (0) No further initializing '
		print*,' (1) Read in pivot from prior calculation '   ! ADDED 7.7.4
		
		print*,' (2) Use PP, NN, and trace of PN '
		print*,' (3) Sector-diagonal parts of PP, NN, and PN only '
		
		if(auto_input)then
			read(autoinputfile,*)choice_initializing
			write(6,*)' Choice = ',choice_initializing
		else
		    read(5,*)choice_initializing
			write(autoinputfile,'(i3,"      ! pivot preparation option")')choice_initializing
		end if
		if(choice_initializing > 1)then
		   print*,' How many iterations to generate ? (must be <= ',niter,')'
		   if(auto_input)then
			   read(autoinputfile,*)initial_maxiter
			   initial_maxiter = min(initial_maxiter,niter)
			   write(6,*)initial_maxiter,' iterations for preparing pivot '
		   else
		       read(5,*)initial_maxiter
		        initial_maxiter = min(initial_maxiter,niter)
				write(autoinputfile,'(i3,"    ! # of initial iterations on pivot ")')initial_maxiter
			end if
	    end if
	end if
    call BMPI_BCAST(choice_initializing,1,0,icomm,ierr)
    call BMPI_BCAST(initial_maxiter,1,0,icomm,ierr)
	
	select case (choice_initializing)
	
   	   case (0)
  	      initializing=.false.
 	      applypntrace =.false.
	      applyXXonly  =.false.
	      diagonalsectorsonly = .false. 
		  return
		  
	   case (2)
          initializing=.true.
          applypntrace =.true.
          applyXXonly  =.true.
          diagonalsectorsonly = .false. 	
		  call masterpnsectortrace   
	
	   case (3)
          initializing=.true.
          applypntrace =.false.
          applyXXonly  =.false.
          diagonalsectorsonly = .true. 	   
!		  call mastersectortrace
      case (1)
          initializing=.false.
          applypntrace =.false.
          applyXXonly  =.false.
          diagonalsectorsonly = .false. 	
		  
!............. READ IN PIVOT .............................  ADDED 7.7.4..........		  
          call wfn_ropen_file(oldwfnfile)
          call read_wfn_header(oldwfnfile,.false.)
          if(dimbasischeck /= dimbasis)then
             if(iproc == 0)then
                print*,' mismatch in basis dimensions '
                print*,' expecting ',dimbasis
                print*,' final ',dimbasischeck
             endif
             stop
           endif
           call wfn_read_nkeep(oldwfnfile, tmpkeep)
		   if(iproc==0) print*,' There are ', tmpkeep,' wavefunctions '
		   do i = 1,tmpkeep
		      ! new interface - we say which vec to read, it checks
		      ! KSM:  This will be very slow, only need to read head part of each vector
		      call wfn_readeigenvec(oldwfnfile,frag1, fcomm1, vec1,i,e,xj,xt2)  
		      if(iproc==0)print*,i,e,xj,xt2
		   enddo

		   if(iproc == 0) then
		      tidx = -1
		      do while(tidx < 1 .or. tidx > tmpkeep)
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
		      write(logfile,*)' Initial state = ',tidx
		   end if
   
		   call BMPI_BCAST(tidx,1,0,icomm,ierr)
		   ! new interface - we say which vec to read, it seeks and reads
		   ! no need to rewind and read forward
		   call wfn_readeigenvec(oldwfnfile, frag1, fcomm1, vec1,tidx,e,xj,xt2)
		   call wfn_close_file(oldwfnfile)
		   return
    end select
	
	if(initializing)then
		call clocker('piv','sta')
		write(logfile,*)' Initializing pivot vector with ',initial_maxiter,' iterations '
		write(logfile,*)' Choosing option ',choice_initializing
		write(logfile,*)' Flags: applypntrace = ',applypntrace,', applyXXonly = ',applyXXonly, & 
		',  diagonalsectorsonly = ',diagonalsectorsonly
		call lanczos_p(0)
		initializing=.false.    ! return values to defaults
		applypntrace =.false.
		applyXXonly  =.false.
		diagonalsectorsonly = .false.
		br_histpos = 0   !IMPORTANT -- need to return to zero
		if(iproc==0)print*,' Finished initializing pivot vector '
		call clocker('piv','end')
		call clockout('piv')
	end if
	return
end subroutine dynamic_vector_initializer
!====================================================================
subroutine swap_vchar(vchar)
  use fragments
  implicit none
  character(1) :: vchar

  if(useNewReorthog) then
     print *, "swap_vchar: should not call with new reorthog"
	 stop 1
  end if

  if(vchar == 'n')then
      vchar = 'r'
   else
      vchar = 'n'
   endif
   return
end subroutine swap_vchar

subroutine print_vec_head(vec, msg)
   use nodeinfo
   use precisions
   implicit none
   integer :: i
   real(kind=lanc_prec) :: vec(*)
   character (len=*) :: msg
   if(iproc /= 0) return

   ! print *, msg, ": ", vec(1), vec(2), vec(3), vec(4), vec(5)
   return
end subroutine print_vec_head

!================================================
!
!  WRITE OUT RESULTS
!
!  CALLS:
!   find_lanczos_eigenvalues
!   setup4obsmaster
!   reset_jumps4obs
!   read_lanczos_vector_a
!   block_reduce
!   applyobsbundled
!   clocker
!   wfn_writeeigenvec (formerly writeeigenvec_p)
!   close_lanczosfile
!   Wcounter
!   density1b_output
!   output_TRDENS
!   makespeme
!
subroutine lanczos_output(iter,i_thick,dnorm0)

  use sporbit
  use localvectors
  use nodeinfo
  use io
  use basis
  use obs
!  use lanczos_info
!  use precisions
  use fragments
  use mod_reorthog
  use flagger
  use coupledmatrixelements
  use system_parameters
  use wfn_mod
  use pocc_mod 
  use butil_mod
  use btbme_mod
  use bvectorlib_mod
  use jump_mod
  use diagh
  use apply_obs_mod
  implicit none

  integer(4) :: iter,i_thick
!  real    :: alpha(niter), beta(niter)
!  real(kind=egv_prec)    :: e(niter), eiglvec(niter,niter)
  real (kind =8 )::dnorm0,da    ! magnitude for strengths

  real(4) :: xj,xt
  logical :: flyeig
  integer(4) :: i,j
  integer(kind=basis_prec) :: k
  integer(kind=basis_prec) :: istate
!--------------------- CONVERGENCE -----------------------------
  
!  logical :: finished
  integer :: ncheck
  real, allocatable :: eold(:)  ! old energies
  real    :: ediff

  logical smallflag,zeroflag
  
  integer :: nread,jvec

  integer inode
  character(1) :: vchar
  integer :: ierr
  real(kind = 8) :: dtmp
  real(kind = 8) :: tmpvec(nkeep)
  ! real(kind=4),allocatable   :: spoccx(:,:,:)  ! single-particle occupations
                                   ! spoccx(it,istate,iorb)  species(p/n), eigenstate, orbit
  integer(4) :: iorb
  character  :: outchoice  !used for steering output options

  vchar = 'n'
  call find_lanczos_eigenvalues(iter,get_JT,i_thick)
  call clocker('tt8','end')

  call clocker('lan','end')
  call clocker('obs','start')

  call pocc_init_spoccx() ! allocate and initialize array
  
  call wfn_write_nkeep(nkeep) ! write number of vectors to wfn file
!---------- SET UP ARRAY FOR SINGLE-PARTICLE OCCUPATIONS --
  outchoice ='d'  ! default option
  if (spoccflag .and. .not.densityflag)then
     if(npeff(1)>0 .and. npeff(2) > 0)outchoice='b'
     if(npeff(1)>0 .and. npeff(2) == 0)outchoice='p'
     if(npeff(1)==0 .and. npeff(2) > 0)outchoice='n'
  end if

!------------ RECONSTRUCT EIGENVECTORS AND COMPUTE J^2, T^2.........
  
  if(get_JT)then
	  
!------------ SET UP FOR ANG MOM................................
  
  call setup4obsmaster('J')
  call setup4obsmaster('T')
!.................. ADDED DEC 2010 by CWJ .........................
! IF USING 3-BODY (OR POSSIBLY SKIPPING ZEROES)
! NEED TO RESET ALL THE JUMPS

  call reset_jumps4obs

  end if
!...................................................................

  if ( iproc == 0 ) print*,' '
  if ( strengthflag ) then
     if ( iproc == 0 ) then
        write(resultfile,*)' Energy   Strength '
        write(resultfile,*)' ______   ________ '
        write(6,*)' Energy   Strength '
        write(6,*)' ______   ________ '

        do j = 1,iter
           write(6,88)e(j),eiglvec(1,j)**2*dnorm0*dnorm0
           if(writeout)write(resultfile,88)e(j),eiglvec(1,j)**2*dnorm0*dnorm0
        enddo
        write(resultfile,*)' ______   ________ '
        write(resultfile,*)' '
        write(6,*)' ______   ________ '
        write(6,*)' '       
     end if
!     return
  endif
88 format(2f10.5)
!---------------- WRITE OUT SINGLE PARTICLE STATES with occupation report --------
  if(spoccflag .and. .not.densityflag)  call pocc_write_orbits()

  call pocc_write_table_header() !  E   Ex ...
  
  allocate(energy(nkeep), xjlist(nkeep),xtlist(nkeep), stat=aerr )
  if(aerr /= 0) call memerror("lanczos_output 1");

!----- Construct nkeep eigenvectors of H from  Lanczos eigenvectors....
!---- If lanczos vectors stored in core, then construct in core

if(get_JT)then                                           !added 7.6.8 
  if(storelanczosincoreMPI .or. storelanczosincore1)then
     call clocker('eig','sta')
     ! KSM lvec(bitpos, iter) is the old way of storing old vectors
     ! It will go away with useNewReorthog
     if(useNewReorthog) then
        ! This applys the same transform as the else clause, but to the 
        ! basis stored in br_histbuf

        call br_transform_basis(eiglvec, nkeep, iter)

     else

        ! KSM 17Aug2014 - my interpretation of the following code
        ! We have a sub basis stored in lvec(:, i)
        ! We extracted the \alpha and \beta coefficients against this basis
        ! This represents a projection of H onto this much smaller basis.
        ! The eigenvectors of the \alpha,\beta matrix can be used to construct
        ! approximate eigenvectors in the untruncated space by taking the linear
        ! combination of lvec(:, i) with weights being the component of the eigenvector.
        !
!$omp parallel do private(tmpvec,i,j,dtmp), shared(eiglvec, Lvec)
        do k = 1,Ldim     ! each row (labeled by kl) can be transformed independently
           dtmp = 0.d0
           do i = 1,nkeep
               dtmp = 0.d0
               do j = 1,iter
                   dtmp = dtmp +   real(lvec(k,j),kind=8)* real(eiglvec(j,i),kind=8)
               end do ! j
               tmpvec(i) = dtmp
           end do  ! i
           do i = 1,nkeep
               lvec(k,i) =real(tmpvec(i),kind=lanc_prec)
           end do
        end do ! k
!$omp end parallel do
     end if
     call clocker('eig','end')
     vchar = 'n'
!---- COMPUTE J,T AND WRITE OUT VECTORS-------------------------------
     do i = 1,nkeep
        ! KSM - the idea here is to restore a history vector into vec1
        if(useNewReorthog) then
           call br_retrieve_hist(i)
           ! It turns out that we need the vector loaded into both vec1 and vec2
           ! of course, they are different slices on each node
           call br_restore_vec1()
        else
           vec1 = 0.0
           call read_lanczos_vector_a(vchar,'i',i,lvec_file)
           if(storelanczosincoreMPI) then
              ! call block_reduce(dimbasis,vec1) 
              ! each mpi process reads one slice.  allreduce is overkill but works
              call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, icomm, ierr) ! in place
           end if
        end if
        call br_load_vec2_from_vec1()

        energy(i) = real(e(i), kind(energy(1)))
!------------ COMPUTE J2, T2.................................
        call clocker('aob','sta')
        xj2 = 0.0e0_obs_prec
        xt2 = 0.0e0_obs_prec
        twoobsflag = .true.  ! compute both J^2 and T^2
        call applyobsbundled(1)
        call clocker('aob','end')
        xj = real( -0.5 + sqrt(xj2 + 0.25), kind=4)
        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)

        xjlist(i) = xj
        xtlist(i) = xt
!!----------------WRITE OUT WFN..............................................
!!  KSM: do this before computing occupations because data like xt2 will be overwritten
!!       not clear why we are saving xt2 instead of xt, but I'm not changing it.
        call clocker('wev','sta')

        if ( writeout .and. write_wfn) then
           call wfn_writeeigenvec(wfnfile, frag1, vec1, i, real(e(i),kind=4), xj, real(xt2,kind=4))
        end if
        call clocker('wev','end')

!----------- COMPUTE SINGLE PARTICLE OCCUPATIONS---------- added in 7.3.7 
        call pocc_compute_spocc(i, .true.)  ! true for restore J/T setup

!----------------- WRITE OUT RESULTS.............................
           if ( iproc == 0 ) then
           select case (outchoice)

               case('d')
                 call pocc_write_ejt(i, e, xj, xt)

               case('b')
                 ! write(6,12)i,e(i), e(i) - e(1),xj,xt,(spoccx(1,i,iorb),iorb=1,numorb(1))
                 call pocc_write_ejt(i, e, xj, xt)
                 call pocc_write_occvec(6, spoccx, i, 1, "    p occ:")
                 call pocc_write_occvec(6, spoccx, i, 2, "    n occ:")
                 if ( writeout ) then
                    call pocc_write_occvec(resultfile, spoccx, i, 1, "    p occ:")
                    call pocc_write_occvec(resultfile, spoccx, i, 2, "    n_occ:")
                 end if
               case('p')
                 ! write(6,12)i,e(i), e(i) - e(1),xj,xt,(spoccx(1,i,iorb),iorb=1,numorb(1))
                 call pocc_write_ejt(i, e, xj, xt)
                 call pocc_write_occvec(6, spoccx, i, 1, "    p occ:")

                 if ( writeout ) then
                    !  write(resultfile,12)i,e(i),e(i) - e(1),xj,xt,(spoccx(1,i,iorb),iorb=1,numorb(1))
                    call pocc_write_occvec(resultfile, spoccx, i, 1, "    p occ:")
                 end if

               case('n')
                 ! write(6,13)i,e(i), e(i) - e(1),xj,xt,(spoccx(2,i,iorb),iorb=1,numorb(2))
                 call pocc_write_ejt(i, e, xj, xt)
                 call pocc_write_occvec(6, spoccx, i, 2, "    n occ:")

                 if ( writeout ) then
                    ! write(resultfile,13)i,e(i),e(i) - e(1),xj,xt,(spoccx(1,i,iorb),iorb=1,numorb(2))
                    call pocc_write_ejt(i, e, xj, xt)
                    call pocc_write_occvec(resultfile, spoccx, i, 2, "    n occ:")
                 end if
            end select
        end if
     end do  ! i

!.....  else construct by writing to disk  (slow).....
  else
   if(useNewReorthog) then
      if(iproc == 0) print *, "useNewReorthog  oops"
      stop 1
   end if
   do i = 1, nkeep

     do k = 1, dimbasis
        vec1(k) = 0.0e0_lanc_prec
     end do !k

     call clocker('eig','sta')
     do jvec = 1,iter
        da = real(eiglvec(jvec,i),kind=8)
        if(.not.storelanczosincoreMPI)then
           call read_lanczos_vector_a(vchar,'f',jvec,lvec_file)
           do k = 1, dimbasis
              vec1(k) = vec1(k) + real(vec2(k),kind=8)*da
           end do  !k
           ! need data on both sides:  note that nodes don't have the same slices on both sides
           call br_load_vec2_from_vec1()
        else
!$omp parallel do private(istate, k)                               &
!$omp              shared(Lstart, Ldim, vec1, Lvec, jvec, da)
           do istate = 1, Ldim
              k = Lstart + istate - 1
              vec1(k) = vec1(k) + real(Lvec(istate,jvec),kind=8)*da
           end do  !istate
!$omp end parallel do
        end if
     end do  !j

     if(storelanczosincoreMPI) then
        ! call block_reduce(dimbasis,vec1)
        ! each mpi process reads one slice.  allreduce is overkill but works
        call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, icomm, ierr) ! in place
     end if
     call clocker('eig','end')

!------------ COMPUTE J2, T2.................................
     call clocker('aob','sta')
     xj2 = 0.0e0_obs_prec
     xt2 = 0.0e0_obs_prec
     twoobsflag = .true.  ! compute both J^2 and T^2
     call applyobsbundled(1)
     call clocker('aob','end')
     xj = real( -0.5 + sqrt(xj2 + 0.25), kind=4)
     xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)

!----------- COMPUTE SINGLE PARTICLE OCCUPATIONS---------- added in 7.3.7
     if(spoccflag .and. .not.densityflag)then
        do iorb = 1,bmax(numorb(1),numorb(2) )
           if(numorb(1) > 0 .and. iorb < numorb(1))then
              pspe = 0.0
              pspe(iorb) = 1.0
              call makespeme(1,'H')
        end if
           if(numorb(2) > 0 .and. iorb < numorb(2))then
              nspe = 0.0
              nspe(iorb) = 1.0
              call makespeme(2,'H')
           end if
           xj2 = 0.0e0_obs_prec
           xt2 = 0.0e0_obs_prec
           call applyspoccbundled(1)
           if(np(1) > 0)spoccx(1,i,iorb)= xj2
           if(np(2) > 0)spoccx(2,i,iorb)= xt2
        end do
     end if
!----------------- WRITE OUT RESULTS.............................
!     if ( isoflag ) then
     if ( iproc == 0 ) then
        select case (outchoice)

            case('d')
               write(6,11)i,e(i), e(i) - e(1),xj,xt
               if ( writeout ) write(resultfile,11)i,e(i),e(i) - e(1),xj,xt

            case('b')
               write(6,121)i,e(i), e(i) - e(1),xj,xt
               call pocc_write_occvec(6, spoccx, i, 1, "    p occ:")
               call pocc_write_occvec(6, spoccx, i, 2, "    n occ:")
               if ( writeout ) then
                  write(resultfile,121)i,e(i), e(i) - e(1),xj,xt
                  call pocc_write_occvec(resultfile, spoccx, i, 1, "    p occ:")
                  call pocc_write_occvec(resultfile, spoccx, i, 2, "    n occ:")
               end if
            case('p')

            case('n')
         end select
      end if
      call pocc_write_orbits()

!----------------WRITE OUT WFN..............................................

     call clocker('wev','sta')
     if ( writeout .and. write_wfn) then
        call wfn_writeeigenvec(wfnfile, frag1, vec1, i, real(e(i),kind=4), xj, real(xt2,kind=4))
     end if
     call clocker('wev','end')
  end do
  end if
  
  else    ! just print out energies for get_JT=.false.
      do i = 1,nkeep
          write(6,11)i,e(i), e(i) - e(1)
          if ( writeout ) write(resultfile,11)i,e(i),e(i) - e(1)
	  end do
	  
  end if
  ! add blank line for formatting
  if(iproc == 0) then
     write(6,*) " "
     if(writeout) write(resultfile, *) " "
  end if

11 format(i5,3x,2f10.5,2x,2f8.3)
12 format(i5,3x,2f10.5,2x,2f8.3,' p occ: ',20f7.3)
121 format(i5,3x,2f10.5,2x,2f8.3)
13 format(i5,3x,2f10.5,2x,2f8.3,' n occ: ',20f7.3)
14 format(46x,' n occ: ',12f7.3)

  if(.not.allsamew .and. get_JT)call Wcounter    ! added in 7.2.5; computes and prints out truncation components

  call close_lanczosfile

  if(densityflag)call density1b_output  !(e,eiglvec)

  if ( trdensout  ) then
     call clocker('trd','sta')
     call output_TRDENS
     call clocker('trd','end')
  end if

  call pocc_cleanup()

  return
end subroutine lanczos_output
!================================================
!
!  function to force conversion of unconverged xJ to integer J
!  that is, odd 2 x J for odd A, and even for even A
!
  function closest2J(evenA,xj)

  implicit none
  integer closest2J
  real xj
  logical evenA

  if(evenA)then
     closest2J = 2*nint(xj)
     if(closest2J < 0)closest2J = 0
  else
     closest2J = 2*nint(xj-0.5)+1
     if(closest2J < 1)closest2J = 1
  end if

  return
  end function closest2J

!---- NOTE: density1b_output moved to bdenslib1.f90
!=============================================================================
!
!  thick-restart lanczos:
!  after some iterations, take nthick lowest eigenvectors
!  and reintroduce them as lanczos vectors
!  then start up iterations again
! 
!  NB: in normal runs, up to this point have done iter = nthick+1 iterations
!
! nthick = # of vectors to keep for thick-restart
! iter = # of lanczos iterations in total
! iter_thick = which iteration of thick_restart (starting at 0)
! maxiter = dimension of arrays
! e(:), eiglev(:,:) : eigenvectors
! alpha(),beta()
! vchar (no longer in significant use)
!
! ALSO USED but not passed in arguments: niter, dimension of truncated Hamiltonian to diagonalize
! NOT USED here: nkeep, the # of final eigenvectors
!
! CALLS:
!  find_lanczos_eigenvalues
!  swap_vchar
!  read_lanczos_vector_a
!  write_lanczos_vector_a
!  block_reduce
!
subroutine thick_restart_sub_p(nthick,iter,iter_thick,maxiter,vchar)
!	subroutine thick_restart_sub_p(nthick,iter,iter_thick,maxiter,e,eiglvec,alpha,beta,vchar)

  use localvectors
  use precisions
  use basis
  use lanczos_info
  use nodeinfo
  use fragments
  use bmpi_mod
  use butil_mod
  use bvectorlib_mod
  use mod_reorthog
  use flagger
  implicit none

  integer  :: nthick  ! # of vectors to keep for thick-restart
  integer  :: iter    ! # of lanczos iterations so far
  integer  :: iter_thick   ! which thick_restart iteration this is, starting with 0
  integer :: maxiter     ! dimension 
!  real(kind=egv_prec) :: e(maxiter), eiglvec(maxiter,maxiter)
  real(kind=8)   :: tmpvec(maxiter)
!  real(kind=4) :: alpha(maxiter), beta(maxiter)
  character(1) :: vchar

  real (kind = lanc_prec),pointer :: v1(:), v2(:)

  integer nthicktmp
  integer(4):: i,j,k
  integer(4):: thickshift  ! for targeted thick-restart
  integer(4):: iclosest
  real(kind=8) :: eclosest
  logical :: verbosetarget = .true.  ! used for debugging excited state thick-restart
  integer, parameter :: tmpfile = 58
  real (kind = lanc_prec) :: da
  integer(kind=basis_prec) :: kl
  integer(4) :: file_mode,file_info,ierror
  integer(4) :: inode

  real(kind = 8) :: dtmp

  integer :: ierr
  real(kind=8)  :: thickconverge
  
  
  if(iproc==0)print*,' In thick restart niter, nthick = ',niter,nthick

!-------------------SOLVE FOR TRUNCATED EIGENVALUES -----------

  if(iter_thick == 0)then
    nthicktmp = 0
  else
    nthicktmp = nthick
  endif
  call find_lanczos_eigenvalues(iter-1,.true.,nthicktmp)

!--------------------CREATE NEW COEFICIENTS -------------------
  do i = 1,nthick
     alpha(i) = e(i)
     beta(i) = beta(iter-1)*eiglvec(iter-1,i) ! correction 6/2011 by WEO
  enddo
  
!------------ CHECK SIZE OF COUPLING -------------------------

  thickconverge = 0.d0
  do i = 1,nkeep
	  thickconverge = thickconverge + beta(i)**2/real(nkeep,kind=8)
  end do
  if(iproc==0)print*,' Alternate thick-restart convergence : ',thickconverge

!.................. REWRITE TO FILE......
  if(iproc ==0)then
     rewind(coef_file)
     print*,' rewriting to file ',coef_file,nthick
     do i = 1,nthick
       write(coef_file,'(i10,2(2x,f20.15))')-i,alpha(i), beta(i)
     end do
  end if

  if(useNewReorthog) then
     ! use eigenvectors of tridiagonal alpha-beta matrix to
     ! mix saved lanczos vectors into approximate eigenstates of H.
     ! we will use these as the first nthick saved lanczos vectors (the history)
     ! going forward.
     ! I understand alpha(i) = e(i) above
     ! I think the beta issue is related to the last vector produced in
     ! the lanczos process.   We have dropped the last beta and associated vector
     ! by truncating the alpha-beta matrix at nthick.
     if(iproc .eq. 0) then
        print *, "New Reorthog Thick Restart"
        print *, "iter=", iter, ", nthick=", nthick
        print *, "br_histpos=", br_histpos
        print *, "nthick_add=", nthick_add
     end if

     call br_transform_basis(eiglvec, nthick, iter-1)

     ! init for next iteration. 
     ! load vec1 from the top of the history.
     call br_retrieve_hist(iter)
     call br_restore_vec1()

     ! ends with nthick as the most recent vector.
     call br_set_histpos(nthick)  ! set last written position
     call br_add2hist(iter+1) ! br_retrieve_hist left in br_reg, put at new head
     return
  end if

!................IF STORING LANCZOS IN CORE .....
  if(storelanczosincore1 .or. storelanczosincoreMPI)then

!$omp parallel do private(tmpvec,i,j,dtmp), shared(eiglvec, Lvec,nthick,iter,Ldim)
     do kl = 1,Ldim     ! each row (labeled by kl) can be transformed independently
        dtmp = 0.d0
        do i = 1,nthick
            dtmp = 0.d0
            do j = 1,iter-1
                dtmp = dtmp +   real(lvec(kl,j),kind=8)* real(eiglvec(j,i),kind=8)
            end do ! j
            tmpvec(i) = dtmp
        end do  ! i
        do i = 1,nthick
            lvec(kl,i) =real(tmpvec(i),kind=lanc_prec)
        end do
     end do ! kl
!$omp end parallel do

!........ SET VCHAR...
     vchar = 'n'
     do i = 1,nthick
        call swap_vchar(vchar)
     end do
     if(vchar=='n')then
        do kl = 1,dimbasis
           vec1(kl) = 0.0
        end do
     end if
     if(vchar=='r')then
        do kl = 1,dimbasis
           vec2(kl) = 0.0
        end do
     end if
     if(useNewReorthog .and. iproc == 0) print *, "read_lanczos_vector in thick_restart_sub_p"
     call read_lanczos_vector_a(vchar,'i',iter,tmpfile)
     if(storelanczosincoreMPI) then
        if(vchar == 'n') call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, icomm, ierr) ! in place
        if(vchar == 'r') call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, icomm, ierr) ! in place
     end if
     return
  end if
!................ IF STORING LANCZOS VECTORS ON DISK.................  
  vchar = 'n'
  open(unit=tmpfile,status ='scratch',form='unformatted')

  do i = 1,nthick
     if(vchar == 'n')then
        v1 => vec1
        v2 => vec2
     else
        v1 => vec2
        v2 => vec1
     end if
     do kl = 1,dimbasis
        v1(kl) = 0.0
     enddo ! j
     do j = 1,iter-1
         if(i == 1)then   
              call read_lanczos_vector_a(vchar,'f',j,lvec_file)
              call write_lanczos_vector_a(vchar,'f',j,tmpfile)
         else
              call read_lanczos_vector_a(vchar,'f',j,tmpfile)
         endif
         da = real(eiglvec(j,i), kind(da))
         do kl = 1,dimbasis
            v1(kl) = v1(kl) + v2(kl)*da
         enddo  !kl
     enddo  !j
     if(i == 1)then
         call read_lanczos_vector_a(vchar,'f',iter,lvec_file)
         call write_lanczos_vector_a(vchar,'f',iter,tmpfile)
     endif
     call write_lanczos_vector_a(vchar,'i',i,lvec_file)
     call swap_vchar(vchar)

  enddo  ! i

  call read_lanczos_vector_a(vchar,'i',iter,tmpfile)
333 format(i4,4f10.6)
  if(nproc ==1)then
     close(tmpfile)
  else
     call BMPI_FILE_CLOSE(tmpfile,ierror)
  end if

  return
end subroutine thick_restart_sub_p

!================================================
!
!  caution: not thoroughly tested
!
!  CALLS:
!  dnormvec_p
!  reorthogonalize_a
!
subroutine random_restart_p(iter,vchar)

   use localvectors
   use fragments
!   use precisions
   use basis
!   use lanczos_info
   use nodeinfo
   use bmpi_mod
   use butil_mod
   use bvectorlib_mod
   implicit none

   character(1) :: vchar
   integer(4) :: iter
   real(kind=lanc_prec), pointer :: w(:)
   real(kind=8):: da
   integer i,j
   real(kind=8) :: dnorm
   logical :: smallflag
   real :: rv
   integer(kind=basis_prec) :: jl,vstart,vstop

!   integer iseed =1959   ! to be used for restarting lanczos if it stops

   if(vchar == 'n')then
      w => vec1
      vstart = basestart(frag1)
      vstop  = basestop (frag1)
   else
      w => vec2
      vstart = basestart(frag2)
      vstop  = basestop (frag2)
   endif

   dnorm = 1.d0 /dsqrt(real(dimbasis,kind=8))

   do jl = vstart,vstop
      call random_number(rv)
      w(jl) = dnorm*(rv-0.5)    !ran(iseed)       ! RANDOM FUNCTION -- may need another
   end do
   call dnormvec_p(vchar,'i',dnorm,smallflag)

!----------- orthogonalize against previous vectors

   call reorthogonalize_a(iter,vchar,dnorm)
   call dnormvec_p(vchar,'i',dnorm,smallflag)

   if(smallflag)then
       if(iproc == 0)print*,' Ooops, zero vec in reorthogonalization of restart '
       stop
   endif

   return
end subroutine random_restart_p

!===========================================================================
! Open Lanczos file
! revision 6/2012 by CWJ
! requires further modification when vectors stored in memory
!
! CALLED BY:
!	lanczos_master in BLANCZOS.f90
!
!  NOTE: in 7.8.1, opening .lvec file turned OFF
!
!===========================================================================
subroutine open_lanczosfile
  use nodeinfo
  use lanczos_info
  use precisions
  use io
  use bmpi_mod
  use butil_mod
  implicit none
  character (len=25) :: filename
  character (len=4)  :: proc_name
  integer(4)         :: ilast
  integer(4)         :: file_mode
  integer(4)         :: file_info
  integer(4)         :: ierror
  writetodisk = .true.

  if ( writeout .and. iproc == 0) then
     ilast = index(outfile,' ') - 1
!     open(unit=lvec_file,file=outfile(1:ilast)//'.lvec',status = 'unknown',form='unformatted')
     open(unit=coef_file,file=outfile(1:ilast)//'.lcoef',status = 'unknown',form='formatted')

  elseif(iproc==0)then
     filename = 'lanczosvec'
     ilast = index(filename,' ')-1
!     open(unit = lvec_file,file=filename(1:ilast)//'.lvec',status='unknown',form ='unformatted')
     open(unit = coef_file,file=filename(1:ilast)//'.lcoef',status='unknown',form ='formatted')

  end if
  return
end subroutine open_lanczosfile

!===========================================================================
!  subroutine close_lanczosfile
! revision 6/2012 by CWJ
!===========================================================================
subroutine close_lanczosfile
  use lanczos_info
  use io
  use nodeinfo
  use bmpi_mod
  use butil_mod
  implicit none

  integer(4) :: ierror
  if ( writetodisk ) then
     if ( .not. restart_enabled ) then
        if ( iproc == 0 ) then ! use normal I/O
           rewind(lvec_file)
           write(lvec_file)0
        end if
     end if
     if ( iproc == 0 ) then
        close(lvec_file)
     end if
     if (iproc == 0) close(coef_file)
  end if
  return
end subroutine close_lanczosfile

!======================================================================
!
!  routine to compute distribution of W values in final wavefunctions
!  added in 7.2.5 by CWJ SDSU  11/2013
!  
subroutine Wcounter
   use W_info
   use sporbit
   use basis
   use sectors
   use io
   use nodeinfo
   use lanczos_info
   use localvectors
   use flagger
   use fragments
!   use  tribution
   use mod_reorthog
   use wfn_mod
   use bmpi_mod
   use butil_mod
   use bvectorlib_mod
   implicit none

   integer(4) :: ierr
   integer is,isc,jsc
   integer(8) :: ip,in
   integer(kind=basis_prec) :: ibasis
   integer Wp, Wn,W
   logical, allocatable :: Wallowed(:)
   integer :: nWvals,n
   integer, allocatable :: Wlist(:)
   real(8), allocatable :: Wfrac(:)
   real(8) ::  ftmp,dv
   integer istate
   integer :: idummy
   real(4) :: xj,xt,ei,ef,xtt,xjj
   integer :: aerr

   if(allsameW)return

   allocate(Wallowed(minWtot:maxWtot), stat=aerr )
   if(aerr /= 0) then
      call memerror("Wcounter 1")
      stop 5
   end if
   Wallowed(:) = .false.
!................ FIND POSSIBLE W IN BASIS.................
   do is = 1,nsectors(1)
      wp = xsd(1)%sector(is)%Wx

      do isc = 1,xsd(1)%sector(is)%ncsectors
         jsc = xsd(1)%sector(is)%csector(isc)
         wn  = xsd(2)%sector(jsc)%Wx
!................ ERROR TRAP..............
         if( Wp + Wn < minWtot .or. Wp+Wn > maxWtot)then
             print*,' wrong ws ',wp,wn, minwtot,maxwtot
             stop
         end if
         Wallowed(wp+wn)=.true.
      end do  ! isc
   end do   ! is
   nWvals = 0
!.............. COUNT HOW MANY................................
   do w = minWtot,maxWtot
      if(Wallowed(w))nWvals = nWvals+1
   end do
   if(nWvals==1)then
       if(iproc==0)then
           print*,' '
           print*,' Only one value of W in truncation encountered '
           print*,' '
           if(writeout)then
                write(resultfile,*)' '
                write(resultfile,*)' Only one value of W in truncation encountered '
                write(resultfile,*)' '
           end if
       end if
       return
   end if
   allocate(wlist(nWvals), stat=aerr)
   if(aerr /= 0) call memerror("Wcounter 2")
   n = 0
!........... CREATE LIST OF Ws............
   do w = minWtot,maxWtot
      if(Wallowed(w))then
            n=n+1
            Wlist(n) = W
      end if
   end do  
   allocate(Wfrac(minWtot:maxWtot), stat=aerr )
   if(aerr /= 0) call memerror("Wcounter 3")
!................ GO THROUGH WAVEFUNCTIONS AND COUNT UP W.....
   if(iproc==0)then
       write(6,*)' '
       write(6,*)'       % occupation of W-truncation subspaces '
       write(6,'(a8,10i7)')'State/W=',(wlist(n)-minwtot,n=1,nWvals)
       if(writeout)then 
           write(resultfile,*)' '
           write(resultfile,*)'       % occupation of W-truncation subspaces '
           write(resultfile,'(a8,10i7)')'State/W=',(wlist(n)-minWtot,n=1,nWvals)
       end if
   end if

!...................LOOP OVER WFNS..............................

  if(.not.storelanczosincore1 .and. .not.storelanczosincoreMPI)then
     call wfn_rewind(wfnfile)   
     call read_wfn_header(wfnfile,.false.)
     call wfn_read_nkeep(wfnfile, n)  ! dummy reading in nkeep
     ! read(wfnfile)n   ! dummy reading in nkeep
  endif
  do istate = 1,nkeep
     if(storelanczosincore1 .or. storelanczosincoreMPI)then
        if(useNewReorthog) then
           call br_retrieve_hist(istate)
           ! It turns out that we need the vector loaded into both vec1 and vec2
           ! of course, they are different slices on each node
           call br_restore_vec1()
        else
           vec1 = 0.0 ! all ranks
           call read_lanczos_vector_a('n','i',istate,lvec_file) ! read in rank=0
           if(storelanczosincoreMPI) then
              ! call block_reduce(dimbasis,vec1)
              ! each mpi process reads one slice.  allreduce is overkill but works
              call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, icomm, ierr) ! in place
           end if
        end if
     else
        ! new interface, we say which vector to read and it checks
        call wfn_readeigenvec(wfnfile,frag1, fcomm1, vec1,istate,ei,xj,xtt)
     end if

     wfrac(:) = 0.d0

     ! figure out if this node is the first node in ifragment.
     if(nodal(iproc)%ifirst) then
        do is = 1,nsectors(1)
           wp = xsd(1)%sector(is)%Wx
   
           do isc = 1,xsd(1)%sector(is)%ncsectors
              jsc = xsd(1)%sector(is)%csector(isc)
              wn  = xsd(2)%sector(jsc)%Wx
              ftmp = 0.d0
              do ip = xsd(1)%sector(is)%xsdstart,xsd(1)%sector(is)%xsdend
                 do in = xsd(2)%sector(jsc)%xsdstart,xsd(2)%sector(jsc)%xsdend
                    ibasis = pstart(ip)+nstart(in)
                    if(ibasis .ge. v1s .and. ibasis .le. v1e) then
                       dv = vec1(ibasis)
                       ftmp = ftmp + dv*dv
                    end if
                 end do
              end do
              wfrac(wp+wn) = wfrac(wp+wn)+ftmp
           end do ! isc
        end do   ! is
     end if
     if(useNewReorthog) then
        ! have to reduce.  Note that condition %ifirst above suppresses nodes that 
        ! are not "first" in their fragment.   This is not very efficient, but it works.
        call BMPI_ALLREDUCE(wfrac(:), SIZE(wfrac), MPI_SUM, icomm, ierr) ! in place reduce
     end if
     if(iproc==0)then
        write(6,'(i4,4x,10f7.2)')istate,(Wfrac(wlist(n))*100.,n=1,nWvals)
        if(writeout)write(resultfile,'(i4,4x,10f7.2)')istate,(Wfrac(wlist(n))*100.,n=1,nWvals)
     end if
  end do  ! istate

  return
end subroutine Wcounter

!======================================================================

! counts up # of lanczos iterations so far
! used in restart option
!
! NOTE: not yet fully parallelized

subroutine countlanczositerations
  use lanczos_info
  use basis
  use precisions
  use localvectors
  use nodeinfo
  use flagger
  use bmpi_mod
  use butil_mod
  implicit none
  integer niter0
 
  integer(4)               :: ierr
  integer(4)               :: i,j,k
  integer(4)               :: iunit
  integer(4)               :: iread
  real(kind=lanc_prec)     :: v
  integer(kind=basis_prec) :: jl,vstart,vstop
  real(kind=4) :: a,b

  if(iproc==0)then
     startiter = 0
     do i = 1,10000
        read( lvec_file,end=101) ( v, jl = 1, dimbasis )
        startiter = startiter+1
     end do
101  continue
     rewind(lvec_file)
!........... CHECKS THAT THIS AGREES WITH lanczos coefficients file
     k = 0
     thick_restart = .false.
     do i = 1,startiter
        read(coef_file,*,end=102)j,a,b
        if(j < 0)then
           thick_restart = .true.
           nkeep = -j-nthick_add
        end if
        k = k+1
     end do
     if(thick_restart)then
         print*,' Sorry, restart not fully working with thick-restart lanczos '
         stop
!     print*,' Keeping ',nkeep,' states '
     end if

     return
102 continue
    if(k < startiter -1)then
        if(iproc == 0)print*,' Problem with restarting ',startiter,k
        call BMPI_ABORT(icomm,101,ierr)
        stop
     endif
     startiter = bmin(startiter,k)

  end if
!.... BROADCAST VARIABLES
  call BMPI_BCAST(startiter,1,0,icomm,ierr)
  call BMPI_BCAST(nkeep,1,0,icomm,ierr)
  call BMPI_BCAST(thick_restart,1,0,icomm,ierr)
  return
end subroutine countlanczositerations

end module lanczos_util