!===================================================================2
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

    real(kind=8), allocatable :: alpha(:),beta(:)  ! lanczos coefficients
    real(kind=egv_prec), allocatable ::  eiglvec(:,:)
    real(kind=egv_prec), allocatable ::  eiglvectx(:,:)
    real(kind=egv_prec), allocatable,target ::  e(:),etx(:),j2tx(:),t2tx(:)
	integer, allocatable :: partx(:)
	
    integer :: txchoice = 1   ! ADDED in 7.9.0
	
    logical :: flyeig
    integer :: nthick_add   ! # of extra states to add for thick-restart; must be >= ncheck_ex
	

  !--------------------- CONVERGENCE -----------------------------------------
  
    logical :: finished
    integer :: ncheck
    real(kind=egv_prec), allocatable :: eold(:)  ! old energies
    real    :: ediff0,ediff1

  !-------------- ALTERNATE CONVERGENCE ---------------------------
  !               based on wavefunctions
    logical :: altconverge

    real(kind=egv_prec), allocatable :: eiglvecold(:,:)
!------------- restricted convergence test, added 7.11.2-----
!            to speed up convergence tests, compute only some eigenvalues

   integer :: ntest_converge ! = sqrt( nkeep + ncheck)
   logical :: use_restricted_esolver,use_switch_esolver
   	
!--------------------------- GREEN'S FUNCTION-------------	
	
	real :: Egf,Egfi   !Energy in Green function 1/(Egf-H); Egfi = imaginary part
	real :: green_tol ! tolerance for convergence of green's function
	real(kind=8), allocatable :: gvec(:),gvec_old(:)
	real(kind=8),allocatable :: gvec_i(:),gvec_i_old(:)
!.... FOR IMPROVED LAPACK ROUTINES, added 7.11.2....
    real(8), allocatable :: zzvec(:,:),workx(:)
	integer, allocatable :: iwork(:),isuppz(:)
	integer :: lwork,liwork	
    real(8), allocatable :: diag(:),subdiag(:)
	
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
!  MODIFIED 7.8.7 to allow for block lanczos
!
! CALLED BY:
!    lanczos_p
!    lanczos_block
!    lanczos_output
!    thick_restart_sub_p
!    block_thick_restart_sub_p
!

subroutine find_lanczos_eigenvalues(n,vecflag,nthick)

  use precisions
  use lanczos_info,only:Etarget,lanczchar
  use localblocks
  use mod_balg
  use flagger
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
  integer :: itarget
  real(8) :: abstol
  integer :: mmm ! needed for DSYEVR
  
  integer info
  real(8) :: DLAMCH

  if(use_switch_esolver .and. .not.thick_restart)abstol = 2*DLAMCH('S')
  
  if(use_switch_esolver .and. .not.thick_restart .and. block_flag == 0)then
!  if(use_switch_esolver  .and. block_flag == 0)then

	  if(.not.allocated(isuppz))allocate(isuppz(2*Niter))
	  lwork = 26*niter
	  if(.not.allocated(workx))allocate(workx(lwork))
	  liwork= 10*niter
	  if(.not.allocated(iwork))allocate(iwork(liwork))
	  if(.not.allocated(zzvec))allocate(zzvec(niter,niter))
	  
  end if
  
 if(use_restricted_esolver .and. n > ntest_converge .and. nthick==0)then
     if(.not.allocated(diag))allocate(diag(niter))
     if(.not.allocated(subdiag))allocate(subdiag(niter))
	 diag = 0.d0
	 subdiag = 0.d0
	 do i = 1,n
		 diag(i)=real(alpha(i),kind=8)
		 subdiag(i)=real(beta(i),kind=8)
	 end do
 end if 
if(block_flag == 0 .and. lanczchar/='tx')then
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
else   ! READ IN BLOCK STRUCTURE
	do i = 1,n
		do j = 1,n
			eiglvec(i,j)=lanczos_mat(i,j)
!			if(lanczchar=='tx')eiglvectx(i,j)=lanczos_mat(i,j)
		end do
	end do
	
	
end if  
  
  
  call clocker('egv','sta')
  if(vecflag)then
      if(egv_prec==4)then
        call SSYEV( 'V','U', N, eiglvec, Niter, e, WORK3, 3*niter, INFO )
      else
		  
		if(use_switch_esolver .and. .not.thick_restart .and. block_flag ==0)then
			! CANNOT USE THIS FOR THICK RESTART, NEED MORE THAN nkeep VECTORS
			! also cannot use for block lanczos
			call DSYEVR('V','I','U',N,eiglvec,Niter,0.d0,0.d0,1,nkeep,abstol,mmm,e,zzvec,Niter, & 
			    isuppz,WORKx,lwork,IWORK,liwork,info)	  
				eiglvec = zzvec
		else
           call DSYEV( 'V','U', N, eiglvec, Niter, e, WORK3, 3*niter, INFO )
		   
       end if
    end if
  else
      if(egv_prec==4)then

         call SSYEV( 'N','U', N, eiglvec, Niter, e, WORK3, 3*niter, INFO )
      else
  		if(use_restricted_esolver .and. n > ntest_converge .and. nthick==0  .and. .not. thick_restart .and. block_flag == 0)then
!  			call DSYEVR('N','I','U',N,eiglvec,Niter,0.d0,0.d0,1+ncheck-ntest_converge,ncheck,abstol,mmm,e,zzvec,Niter, & 
 ! 			    isuppz,WORKx,lwork,IWORK,liwork,info)	  
			call dstegr('N','I',n,diag,subdiag,0.d0,0.d0,1+ncheck-ntest_converge,ncheck,abstol,Mmm,e,zzvec,niter, & 
			   isuppz,workx,lwork,iwork,liwork,info)
				
		else
			
               call DSYEV( 'N','U', N, eiglvec, Niter, e, WORK3, 3*niter, INFO )
		   end if
      end if
  endif
  call clocker('egv','end')
  
  
  if(kpow > 1 .and. .not. vecflag)then   ! only do for convergence test
	  do i =1,n
		  if(mod(kpow,2)==0)then
			  e(i)= e(i)**(1./float(kpow))
		  else
			  xswap=e(i) ! HERE xswap IS JUST A CONVENIENT DUMMY
			  e(i)=(abs(e(i)))**(1./float(kpow))
			  if(xswap < 0)e(i)=-e(i)
		  end if
		  
	  end do
	  
  end if
!--------- OPTION TO CHOOSE EXCITED STATES --------------------
!  added in 7.7.8 based in part upon work by R. Zbikowski
!  SORT EIGENVALUES, EIGENVECTORS BASED UPON PROXIMITY TO Etarget
!  
  if(lanczchar=='tx')then
	  
if(txchoice==1)then	  
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
!.............. instead, solve for (H-E0)^2
!   ADDED 7.9.9, 
else  ! txchoice = 2
	do i = 1,n
		lanczos_mat(i,i)=lanczos_mat(i,i)-etarget
	end do
	eiglvectx= 0.0
   do i = 1,n
	   do j = 1,n
		   do k = 1,n
			   eiglvectx(i,j)=eiglvectx(i,j)+ lanczos_mat(i,k)*lanczos_mat(k,j)
			   
			   
		   end do
	   end do
   end do	
   do i = 1,n
	   lanczos_mat(i,i)=lanczos_mat(i,i)+etarget
   end do
   if(vecflag)then
       if(egv_prec==4)then
         call SSYEV( 'V','U', N-1, eiglvectx, Niter, etx, WORK3, 3*niter, INFO )
       else
         call DSYEV( 'V','U', N-1, eiglvectx, Niter, etx, WORK3, 3*niter, INFO )
       end if
   else
       if(egv_prec==4)then

          call SSYEV( 'N','U', N-1, eiglvectx, Niter, etx, WORK3, 3*niter, INFO )
       else
          call DSYEV( 'N','U', N-1, eiglvectx, Niter, etx, WORK3, 3*niter, INFO )

       end if
   endif	
   do i = 1,n-1
	   etx(i)=sqrt(abs(etx(i)))
   end do
!   print*,' in alternate ',e(1:10)

	
end if	  

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
    implicit none

    call setup_localvectors
    call initialize_lanczos_vector

    return
		
end subroutine setup_for_lanczos
!*----------------------------------------------------------------
!   subroutine read_pivot_block
!
!   Purpose: This is just the modified readpivot subroutine from
!   bwfnlib.f90. The routine loads in a number of pivot vectors into
!   history equal to dimblock
!
!   the file of eigenvectors must be previously opened
!   and the header reader;
!   reads the # of kept vectors (nkeep),
!   prints out E J, T^2 of each vector
!   then asks which to use as a pivot
!
!   OUTPUT
!      v(n), n = 1,dimbasis = vector read from disk
!      dnorm = norm of pivot vector
!
!   CALLS
!     wnf_read_nkeep
!     wfn_readeigenvec
!     wfn_close_file
!
!*----------------------------------------------------------------
!
! adapted in 7.11.1 for using in option (h) to calculate hamiltonian matrix elements
!  if for_ham_me = T, then will have previously opened a file and figured out nkeep
!
subroutine read_pivot_block(for_ham_me)
   use localvectors
   use localblocks
   use precisions
   use basis
   use io
   use nodeinfo
   use bvectorlib_mod
   use mod_reorthog
   use wfn_mod
   implicit none

   logical :: for_ham_me
   real (kind = 8) :: dnorm
   integer nnkeep
   integer ikeep
   integer(4):: i,j,n, ierr
   real :: e,xj,xt2
   logical zeroflag
   integer :: aerr
   real(kind=8) :: scale
   integer :: ipar


   if(for_ham_me)then
	   dimblock = nkeep
   else
	   if(iproc==0)print*, " reading in pivot block"
	   
       call wfn_read_nkeep(oldwfnfile, nnkeep)  ! does BCAST
   
       if(nnkeep < dimblock)then
	      if(iproc==0)print*,' Insufficient pivot vectors, only ',nnkeep,', asked for ',dimblock
	      if(iproc==0)print*,' Block size must be <= ',nnkeep

#ifdef _MPI
	      call BMPI_ABORT(MPI_COMM_WORLD,101,ierr)
#endif
	      stop
	  end if
   end if

   do ikeep = 1, dimblock
!    call BMPI_BCAST(ikeep,1,0,icomm,ierr)
      vec1 = 0.d0
      call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1, ikeep, e, xj, xt2) ! KSM: updated
	  
	  if(for_ham_me)then
		  etx(ikeep)=e
		  j2tx(ikeep)=xj
		  t2tx(ikeep)=xt2
  		  call find_parity(ipar,1)
		  partx(ikeep)=ipar
		  
	  end if
      call br_grab_vec1()! push vec1 -> br_reg
! DON'T NORMALIZE--THIS COMES LATER	  
      !call br_normalize(scale)  
      call br_add2hist(1)
!write(logfile,*) "------------history matrix--------------"
!write(logfile,*) br_histbuf(:,:)
   enddo

   if(.not.for_ham_me)call wfn_close_file(oldwfnfile)
   
   !call exit
!print*, " reading in pivot block done"
   return
end subroutine read_pivot_block
!
!*----------------------------------------------------------------
!*----------------------------------------------------------------
!   subroutine read_pivot_block_chose
!
!   modified 7.9.4 from read_pivot_block
!
!   Purpose: This is just the modified readpivot subroutine from
!   bwfnlib.f90. The routine loads in a number of pivot vectors into
!   history equal to dimblock
!
!
!  
!
!   OUTPUT
!      v(n), n = 1,dimbasis = vector read from disk
!      dnorm = norm of pivot vector
!
!   CALLS
!     wnf_read_nkeep
!     wfn_readeigenvec
!     wfn_close_file
!
!*----------------------------------------------------------------
subroutine read_pivot_block_choose
   use localvectors
   use localblocks
   use precisions
   use basis
   use io
   use nodeinfo
   use bvectorlib_mod
   use mod_reorthog
   use wfn_mod
   use mod_balg
   implicit none

   real (kind = 8) :: dnorm
   integer nnkeep
   integer ikeep
   integer(4):: i,j,n, ierr
   integer(4):: pivstart,pivstop
   real :: e,xj,xt2
   logical zeroflag,founderr
   integer :: aerr
   real(kind=8) :: scale
   character :: listchar
   
   integer :: ilist
   
   if(iproc==0)print*, " reading in pivot block"
   allocate(pivot_choice(dimblock))
   allocate(dnormblock(dimblock))

   call wfn_read_nkeep(oldwfnfile, nnkeep)  ! does BCAST
      
   if(nnkeep < dimblock)then
	   if(iproc==0)then
		   print*,' Insufficient pivot vectors, only ',nnkeep,', asked for ',dimblock
		   print*,' (NEW)  Fill out rest with random vectors '
	   end if
	   
!	   if(iproc==0)print*,' Block size must be <= ',nnkeep
!	   call BMPI_ABORT(icomm,101,ierr)
!	   stop
   end if
   
   if(iproc==0) print*,' There are ', nnkeep,' wavefunctions '

   if(iproc==0)then
        print*,nkeep,' states '
        print*,' '
        print*,' STATE      E         J          <H > '
	endif
		
    do ikeep = 1,nnkeep
      i = ikeep
      ! new interface, we say which vector to read - it checks
      call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,e,xj,xt2) ! KSM: updated
      if(iproc==0)then
          write(6,101)i,e,xj,xt2

 101      format(i5,2x,4(1x,f11.4))
      end if
   end do ! ikeep
   founderr=.false.
   if(iproc==0)then  ! modified in 7.10.6
1044   print*,' '
	   print*,' Need to pick a set of states for the pivot block '
	   print*,' Choose either (c) a contiguous list or (s) list of specified states '
	   print*,' (For a contiguous list, just need first, last states; any missing will be random)'
	   print*,' (For specified list, you must list all the chosen states )'
	   
	   if(auto_input)then
		   read(autoinputfile,'(a)')listchar
		   select case (listchar)
		   case ('c','C')
		   print*,' Enter start, stop of list out of ',nnkeep,' states '
		   read(autoinputfile,*)pivstart,pivstop
		   print*,pivstart,pivstop

	       pivot_choice=0
		   do i =1,pivstop-pivstart+1
			 pivot_choice(i)=pivstart+i-1
		   end do
		   
		   case ('s','S')
		   print*,' Enter a list of ',dimblock,' states for the pivot block'
		   print*,' (Please enter in order; enter 0 to create random )'  ! modified in 7.9.11
		   if(dimblock > nnkeep)print*,' At least ',dimblock-nnkeep,' must be random '
  	       read(autoinputfile,*)pivot_choice(:)
		   
		   case default
		   
		       goto 1044
		   
	       end select
		   
		   
	   else
		   read(5,'(a)')listchar
		   select case (listchar)
		   
		   case ('c','C')
		   print*,' Enter start, stop of list out of ',nnkeep,' states '
		   read(5,*)pivstart,pivstop
		   
		   if(pivstart > nnkeep-1 .or. pivstop > nnkeep)then
			   print*,' You only have ',nnkeep,' states to choose from '
			   goto 1044
		   end if
	       pivot_choice=0
		   do i =1,pivstop-pivstart+1
			 pivot_choice(i)=pivstart+i-1
		   end do
		   write(autoinputfile,'(a1)')listchar
		   write(autoinputfile,*)pivstart,pivstop
		   
		   case ('s','S')
		   print*,' Enter a list of ',dimblock,' states for the pivot block'
		   print*,' (Please enter in order; enter 0 to create random )'  ! modified in 7.9.11
		   if(dimblock > nnkeep)print*,' At least ',dimblock-nnkeep,' must be random '
		   write(autoinputfile,'(a1)')listchar
  	       read(5,*)pivot_choice(:)
		   write(autoinputfile,*)pivot_choice
		   
		   case default
		   
		   print*,' that choice is not available '
		   goto 1044
		   
	       end select
		   
	   end if
	   
   end if   
      
#ifdef _MPI
   call BMPI_BCAST(pivot_choice,dimblock,0,MPI_COMM_WORLD,ierr)	   
   call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif   
   do ilist = 1,dimblock
	   ikeep= pivot_choice(ilist)
!   if(iproc==0)call BMPI_BCAST(ikeep,1,0,icomm,ierr)
!	   pivot_choice(ilist)=ikeep
	   
	   if(ikeep > nnkeep)then
           if(iproc==0)then
			   print *, "Selected vector is not in range (A):", nnkeep
			   print*,' Will replace with random vector '
		   end if
		   pivot_choice(ilist)=0
!		   stop 1
	   end if

   end do
   if(iproc==0)print*,pivot_choice,' chosen '

   
   do ilist = 1,dimblock
!    call BMPI_BCAST(ikeep,1,0,icomm,ierr)
      vec1 = 0.d0
	  ikeep = pivot_choice(ilist)
	  
	  if(ikeep < 1)then ! create a random vector
	  	call add_pivot_br_reg(0,ilist)	
	  	call br_normalize(dnorm)
	  	call br_orthogonalize(0)
	  	call br_normalize(dnorm)		
	    call br_add2hist(1)	 
	  else
          call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1, ikeep, e, xj, xt2) ! KSM: updated
          call br_grab_vec1()! push vec1 -> br_reg
!.... MUST ORTHONORMALIZE HERE AND MUST KEEP!	  
! 	  call br_normalize(dnorm)		
	  
      !call br_normalize(scale)
          call br_add2hist(1)
	  end if
		  
!write(logfile,*) "------------history matrix--------------"
!write(logfile,*) br_histbuf(:,:)
   enddo

   call wfn_close_file(oldwfnfile)
   
   !call exit
if(iproc==0)print*, " reading in pivot block done"
   return
end subroutine read_pivot_block_choose
!
!*----------------------------------------------------------------
!  Subroutine: initialize_pivot_block
!
!  Purpose:  Creates a pivot block of random orthonormal vectors within history.
!
!
!  Parameters:
!      dimblock(global IN)- size of block
!
!  Returns:
!*-------------------------------------------------------------------

subroutine initialize_pivot_block
use localblocks
use mod_reorthog
use mod_balg
implicit none
integer(kind=4) :: i
integer :: pivchoice
real(8) :: dnorm

pivchoice = 0

call br_set_histpos(0)
call br_reset_histbuf()
do i = 1,dimblock
	call add_pivot_br_reg(pivchoice,i)	
	call br_normalize(dnorm)
	call br_orthogonalize(0)
	call br_normalize(dnorm)		
    call br_add2hist(1)
enddo

end subroutine initialize_pivot_block
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
   use localblocks
   use sporbit,only:allsameW
   use W_info
   use sectors
!   use mod_balg
   implicit none 
   character(1) :: vchar
   real(kind=8) :: dnorm, dnorm0,da
   real :: rv
   integer startiter
   integer (kind = basis_prec) :: i
   integer (kind=4) :: iter
   logical smallflag
   integer :: pivotoption   ! different choices
   integer(8) :: countstates
!..... ADDED IN 7.8.7 to improve pivot....
   
   integer (kind=basis_prec) :: startvec,stopvec,npxsd,nnxsd
   integer :: isector,ncsector,icsector,csector
   integer :: ierr
   
   integer:: nseed,value(8)
   integer, allocatable:: seed(:)

   pivotoption = 0
   
   if(iproc==0)print*,' USING PIVOT OPTION ',pivotoption

   vchar = 'n'
   if ( strengthflag  .or. menu_char=='np' .or. menu_char=='pv'.or. greenflag) then
	   
    if(piv_flag .eq. 'y' .or. blockstrengthflag) then
	   call read_pivot_block_choose   ! now default
    else
      call readpivot   ! FIXED so it works with fragments

	  return
   endif


   else if ( startiter == 0 ) then

      dnorm = dsqrt(real(dimbasis,kind=8))
1111 continue

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

!..... ALTERNATING SIGN SUGGESTED BY WEO..............
       case (3)
	   
	  do i = v1s, v1e

          vec1(i) = real( 1.d0 * (-1)**i /dnorm, kind=lanc_prec)

       end do

!.................CHOOSING LOWEST W MAY MINIMIZE TRACE...
	   
       case (4)
	   
	   vec1 = 0.d0
	   if(allsameW)then
		   startvec=v1s
		   stopvec =v1e
	 	  do i = startvec,stopvec
	           vec1(i) = real( 1.d0 * (-1)**i /dnorm, kind=lanc_prec)
	       end do
	   else
		   countstates=0
		   do isector = 1,nsectors(1)
			   if(xsd(1)%sector(isector)%basisstart > v1e .or.  & 
			     xsd(1)%sector(isector)%basisend < v1s  )cycle
			   if(xsd(1)%sector(isector)%Wx /= minW(1))cycle
			   
			   do icsector = 1,xsd(1)%sector(isector)%ncsectors
				   csector = xsd(1)%sector(isector)%csector(icsector)
				   if(xsd(2)%sector(csector)%Wx /= minW(2))cycle
				   startvec = pstart( xsd(1)%sector(isector)%xsdstart) & 
				   + nstart( xsd(2)%sector(csector)%xsdstart)
				   stopvec =startvec+ & 
				    xsd(1)%sector(isector)%nxsd*xsd(2)%sector(csector)%nxsd-1
				   
				   startvec = max(startvec,v1s)
				   stopvec  = min(stopvec,v1e)
				   if(stopvec < startvec)cycle
				   
                   countstates=countstates+stopvec-startvec+1
!........ NOW LOOP OVER THESE SECTORS................				   
                   do i = startvec,stopvec					   
                      vec1(i) = real( 1.0 *(-1)**i /dnorm, kind=lanc_prec)
                   end do   
			   end do		   
		   end do
		   
#ifdef _MPI
		   call BMPI_ALLREDUCE(countstates,1,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
		   if(countstates<10)then
			   pivotoption = 3
			   goto 1111
		   end if
		   
		   	   
	   end if
	   
	   case (5,-5)
	   
!	   if(nproc==1)then
	     call random_seed(size=nseed)
	     allocate(seed(nseed))
		 if(pivotoption==-5 .and. nproc==1)then
	        print*,' Enter random seed '
	        read*,seed(1)
			
		else
			call date_and_time(VALUES=value)  ! kluge to set random seed to time
			seed(1)=value(8)
			if(nseed>1)seed(2)=value(8)
			
			
		end if
		if(iproc==0)print*,' Creating pivot with random seeds ',seed
	   
	     call random_seed(put=seed)
!	   end if
!	   print*,' seed ',nseed,seed
 	  do i = v1s, v1e

 !....... TO PREVENT ROUNDOFF ERRORS AND ALSO TO MINIMIZE
 !        ACCIDENTAL STOPS, HAVE DECREASING WEIGHT
           call random_number(rv) 
!if(i< 11)print*,i,rv

           vec1(i) = real( (0.5d0-rv) /dnorm, kind=lanc_prec)
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


!subroutine dynamic_vector_initializer
	! added 7.6.3, deleted 7.9.11
!====================================================================
!
! routine is mostly not used and should be phase out
! original idea was to go back and forth between vec1 and vec2
! where vchar='n' means normal:  H vec1 = vec2
!  and  vchar='r' means reverse: H vec2 = vec1
! but this idea doesn't lend itself well to MPI 
!
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
  use mod_balg
  use localblocks
  implicit none

  integer(4) :: iter,i_thick
!  real    :: alpha(niter), beta(niter)
!  real(kind=egv_prec)    :: e(niter), eiglvec(niter,niter)
  real (kind =8 )::dnorm0,da    ! magnitude for strengths

  real(4) :: xj,xt
  integer :: xpar
!  logical :: flyeig
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
  if(strengthflag .or. blockstrengthflag)then
      call find_lanczos_eigenvalues(iter,.true.,i_thick)
  else
    call find_lanczos_eigenvalues(iter,get_JT,i_thick)
  end if
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
  if(.not.skip_T2)then
      call setup4obsmaster('T')
  end if
!.................. ADDED DEC 2010 by CWJ .........................
! IF USING 3-BODY (OR POSSIBLY SKIPPING ZEROES)
! NEED TO RESET ALL THE JUMPS

  call reset_jumps4obs
  end if
!...................................................................

  if ( iproc == 0 ) print*,' '
  if ( strengthflag ) then
	  
	  if(.not.blockstrengthflag)then
          if ( iproc == 0 ) then

              write(resultfile,*)' Energy   Strength    '
              write(resultfile,*)' ______   ________  '

              write(6,*)' Energy   Strength  '
              write(6,*)' ______   ________  '

              do j = 1,iter
                  write(6,88)e(j),eiglvec(1,j)**2*dnorm0*dnorm0
                 if(writeout)write(resultfile,88)e(j),eiglvec(1,j)**2*dnorm0*dnorm0
              enddo
              write(resultfile,*)'_______   ________   '
			  write(resultfile,*)dnorm0**2,'= total strength, included above '
			  write(resultfile,*)' '
              write(6,*)'_______   ________  '
			  write(6,*)dnorm0**2,'= total strength, included above '
			  write(6,*)' '
              write(6,*)' '       
           end if
	 
	   else   ! BLOCK STRENGTH-- NEED TO CONVERT
!		   print*,' normalization ',orthonorm_pivot
	       do i = 1,dimblock
			   if(iproc==0)then
				   
				   dtmp  = 0.d0
				   
				   do j = 1,dimblock
					   
					   dtmp = dtmp + orthonorm_pivot(j,i)**2
				   end do
				   write(6,*)'  '
				   write(6,*)'################### '
				   write(6,*)'  '
				   write(6,*)' For pivot state ',pivot_choice(i),' total strength = ',dtmp
	               write(6,*)' Energy   Strength   '
	               write(6,*)' ______   ________  '		
				   write(resultfile,*)'###################'
				   write(resultfile,*)' For pivot state ',pivot_choice(i),' total strength = ',dtmp
	               write(resultfile,*)' Energy   Strength  '
	               write(resultfile,*)' ______   ________  '
				  	   
				   do j = 1,iter
				      dtmp = 0.d0
				      do k = 1,dimblock
					      dtmp = dtmp +   orthonorm_pivot(k,i)*eiglvec(k,j)
				       end do
				   
				       write(6,88)e(j),dtmp*dtmp !, dtmp*dtmp 
			           if(writeout)write(resultfile,88)e(j),dtmp*dtmp !,dtmp*dtmp
			   
		            end do ! j
				end if  ! iproc == 0
		   end do
	 
       endif 
!     return
  endif
88 format(3f12.6)

!------------------ APPLY RESOLVENT/GREEN FUNCTION AND WRITE OUT----------------------
!   (Added 7.9.10
if(greenflag)then
	eiglvec = 0.
	do j = 1,iter
		if(complex_green)then
 		   eiglvec(j,1)=gvec(j)*dnorm0
 		   eiglvec(j,2)=gvec_i(j)*dnorm0
		else
		   eiglvec(j,1)=gvec(j)*dnorm0
	    end if
	end do
    call br_transform_basis(eiglvec, nkeep, iter)
    call br_retrieve_hist(1)
    ! It turns out that we need the vector loaded into both vec1 and vec2
    ! of course, they are different slices on each node
    call br_restore_vec1()
    call wfn_writeeigenvec(wfnfile, frag1, vec1, 1, 0.0, 0.0, 0.0)
	
	if(complex_green)then
	    call br_retrieve_hist(2)
	    ! It turns out that we need the vector loaded into both vec1 and vec2
	    ! of course, they are different slices on each node
	    call br_restore_vec1()
	    call wfn_writeeigenvec(wfnfile, frag1, vec1, 2, 0.0, 0.0, 0.0)		
	end if
	
    if(iproc==0 .and. write_wfn)then
  	  print*,' '
  	  print*,' Note: output wfns not normalized  '
	  if(complex_green)then
		  print*,' '
		  print*,' Note: *** real part written to vector 1, imaginary to vector 2 *** '
		  
	  end if
    end if
	
	
	call wfn_close_file(wfnfile)
	
	if(dotflag)then
		if(iproc==0)print*,' Finally, read in a previously computed file to take dot product '
		call overlap(.true.)
	end if
	
	return
	
end if
!--------------------------------------------------------------

!---------------- WRITE OUT SINGLE PARTICLE STATES with occupation report --------
  if(spoccflag .and. .not.densityflag)  call pocc_write_orbits()


  if( strengthflag .and. .not.write_wfn  )then
	  if( iproc==0 )then
		  print*,' '
		  print*,' Note: not writing wfns to file, nor computing J, T '
		  print*,' '
		  write(logfile,*)' Note: not writing wfns to file, nor computing J, T '
		  write(resultfile,*)' Note: not writing wfns to file, nor computing J, T '
	  
	  end if
	  return
	  
  end if
  call pocc_write_table_header() !  E   Ex ...
  
  allocate(energy(nkeep), xjlist(nkeep),xtlist(nkeep),parity(nkeep), stat=aerr )
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
#ifdef _MPI
              call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
#endif
           end if
        end if
		
        call br_load_vec2_from_vec1()
        energy(i) = real(e(i), kind(energy(1)))
!------------ COMPUTE J2, T2.................................
        call clocker('aob','sta')
        xj2 = 0.0e0_obs_prec
        xt2 = 0.0e0_obs_prec
        twoobsflag = .true.  ! compute both J^2 and T^2
		
		if(skip_T2)twoobsflag=.false. ! special case
		
        call applyobsbundled(1)
        call clocker('aob','end')
        xj = real( -0.5 + sqrt(xj2 + 0.25), kind=4)
        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)

        xjlist(i) = xj
        xtlist(i) = xt
   	    !if(print_parity)
		call find_parity(xpar,1)
		parity(i)=xpar
		
!!----------------WRITE OUT WFN..............................................
!!  KSM: do this before computing occupations because data like xt2 will be overwritten
!!       not clear why we are saving xt2 instead of xt, but I'm not changing it.
        call clocker('wev','sta')

        if ( writeout .and. write_wfn) then
		 
   		   if(.not.strengthnorm .and. strengthflag .and. .not. blockstrengthflag)then ! added 7.9.7
 !  			  print*,eiglvec(1,i)*dnorm0,' testing '
   			  vec1=vec1*abs(eiglvec(1,i)*dnorm0)
   		   end if
           call wfn_writeeigenvec(wfnfile, frag1, vec1, i, real(e(i),kind=4), xj, real(xt2,kind=4))
        end if
        call clocker('wev','end')

!----------- COMPUTE SINGLE PARTICLE OCCUPATIONS---------- added in 7.3.7 
        call pocc_compute_spocc(i, .true.)  ! true for restore J/T setup

!----------------- WRITE OUT RESULTS.............................
           if ( iproc == 0 ) then
           select case (outchoice)

               case('d')
                 call pocc_write_ejt(i, e, xj, xt,xpar)

               case('b')
                 ! write(6,12)i,e(i), e(i) - e(1),xj,xt,(spoccx(1,i,iorb),iorb=1,numorb(1))
                 call pocc_write_ejt(i, e, xj, xt,xpar)
                 call pocc_write_occvec(6, spoccx, i, 1, "    p occ:")
                 call pocc_write_occvec(6, spoccx, i, 2, "    n occ:")
                 if ( writeout ) then
                    call pocc_write_occvec(resultfile, spoccx, i, 1, "    p occ:")
                    call pocc_write_occvec(resultfile, spoccx, i, 2, "    n_occ:")
                 end if
               case('p')
                 ! write(6,12)i,e(i), e(i) - e(1),xj,xt,(spoccx(1,i,iorb),iorb=1,numorb(1))
                 call pocc_write_ejt(i, e, xj, xt,xpar)
                 call pocc_write_occvec(6, spoccx, i, 1, "    p occ:")

                 if ( writeout ) then
                    !  write(resultfile,12)i,e(i),e(i) - e(1),xj,xt,(spoccx(1,i,iorb),iorb=1,numorb(1))
                    call pocc_write_occvec(resultfile, spoccx, i, 1, "    p occ:")
                 end if

               case('n')
                 ! write(6,13)i,e(i), e(i) - e(1),xj,xt,(spoccx(2,i,iorb),iorb=1,numorb(2))
                 call pocc_write_ejt(i, e, xj, xt,xpar)
                 call pocc_write_occvec(6, spoccx, i, 2, "    n occ:")

                 if ( writeout ) then
                    ! write(resultfile,13)i,e(i),e(i) - e(1),xj,xt,(spoccx(1,i,iorb),iorb=1,numorb(2))
                    call pocc_write_ejt(i, e, xj, xt,xpar)
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
#ifdef _MPI
        call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
#endif
     end if
     call clocker('eig','end')

!------------ COMPUTE J2, T2.................................
     call clocker('aob','sta')
     xj2 = 0.0e0_obs_prec
     xt2 = 0.0e0_obs_prec
     twoobsflag = .true.  ! compute both J^2 and T^2
	 
	 if(skip_T2)twoobsflag=.false. ! special case
	 
     call applyobsbundled(1)
     call clocker('aob','end')
     xj = real( -0.5 + sqrt(xj2 + 0.25), kind=4)
     xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
	 
	 !if(print_parity)
	 call find_parity(xpar,1)
	 parity(i)=xpar

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
		 
		 if(.not.strengthnorm .and. strengthflag .and. .not. blockstrengthflag)then ! added 7.9.7
			 print*,eiglvec(1,i)*dnorm0,' testing '
			 vec1=vec1*abs(eiglvec(1,i)*dnorm0)
		 end if
		 
        call wfn_writeeigenvec(wfnfile, frag1, vec1, i, real(e(i),kind=4), xj, real(xt2,kind=4))
     end if
     call clocker('wev','end')
  end do
  end if
  
  else    ! just print out energies for get_JT=.false.
	  if(iproc==0)then
          do i = 1,nkeep
             write(6,11)i,e(i), e(i) - e(1)
             if ( writeout ) write(resultfile,11)i,e(i),e(i) - e(1)
	      end do
      end if
	  
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

if(densityflag) then
	if(iproc==0)then
		write(denresultfile,*)' '
		write(occresultfile,*)' '
		
	    call pocc_write_orbits_alt(denresultfile)
	    call pocc_write_orbits_alt(occresultfile)
    endif
	 
 end if

  if(.not.allsamew .and. get_JT)call Wcounter    ! added in 7.2.5; computes and prints out truncation components

  call close_lanczosfile

!........... ADDED in 7.9.7: OPTION TO CARRY FORWARD NORMALIZATIONS....


  if(.not.strengthnorm .and. strengthflag .and. .not. blockstrengthflag .and. iproc==0 .and.write_wfn)then
	  print*,' '
	  print*,' Note: output wfns not normalized;  '
	  print*,' they carry forward input normalization and strength'
	  print*,'  (To change this, use menu option (sn))'
	  print*,' '
	  write(logfile,*)' '
	  write(logfile,*)' Note: output wfns not normalized;  '
	  write(logfile,*)' they carry forward input normalization and strength'
	  write(logfile,*)'  (To change this, use menu option (sn))'
	  write(logfile,*)' '
	  write(resultfile,*)' '
	  write(resultfile,*)' Note: output wfns not normalized;  '
	  write(resultfile,*)' they carry forward input normalization and strength'
	  write(resultfile,*)'  (To change this, use menu option (sn))'
	  write(resultfile,*)' '
	  
  end if
  if(strengthnorm .and. strengthflag .and. .not. blockstrengthflag .and. iproc==0 .and. write_wfn)then
	  print*,' '
	  print*,' Note: output wfns are normalized;  '
	  print*,'  (To change this, use menu option (su))'
	  print*,' '
	  write(logfile,*)' '
	  write(logfile,*)' Note: output wfns are normalized;  '
	  write(logfile,*)'  (To change this, use menu option (su))'
	  write(logfile,*)' '
	  write(resultfile,*)' '
	  write(resultfile,*)' Note: output wfns are normalized;  '
	  write(resultfile,*)'  (To change this, use menu option (su))'
	  write(resultfile,*)' '
	  
  end if
  
  if(densityflag)call density1b_output  !(e,eiglvec)

  if ( trdensout  ) then
     call clocker('trd','sta')
     call output_TRDENS(.false.)
     call clocker('trd','end')
  end if

  call pocc_cleanup()
  call wfn_close_file(wfnfile)
  

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
! CALLED BY:
!   lanczos_p
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
  use mod_balg,only:lanczos_mat
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
  integer(4):: i,j,k,l
  integer(4):: mymax  ! for targeted thick-restart
  integer(4):: iclosest
  real(kind=8) :: eclosest
  logical :: verbosetarget = .true.  ! used for debugging excited state thick-restart
  integer, parameter :: tmpfile = 58
  real (kind = lanc_prec) :: da
  integer(kind=basis_prec) :: kl
  integer(4) :: file_mode,file_info,ierror
  integer(4) :: inode
  real(8) :: lanczos_mat_tmp(niter,niter)

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
  
!----- ALTERNATE, MORE EXACT but more time-consuming METHOD ADDED IN 7.9.0
!      recompute the matrix elements by applying the Hamiltonian

if(lanczchar=='tx' )then
	
	if(txchoice==1)then
		lanczos_mat = 0.d0
		
	   do i = 1,nthick
		  lanczos_mat(i,i)=e(i)
		  lanczos_mat(i,nthick+1)=beta(iter-1)*eiglvec(iter-1,i)
		  lanczos_mat(nthick+1,i)=beta(iter-1)*eiglvec(iter-1,i)	
	   end do
	   
   else
	   print*,' That choice is not implemented'
	   stop 909

    end if	
end if
  
!------------ CHECK SIZE OF COUPLING -------------------------

  thickconverge = 0.d0
  if(lanczchar=='tx' .and. txchoice==2)then
      do i = 1,nkeep
  	     thickconverge = thickconverge + lanczos_mat(i,nthick+1)**2/real(nkeep,kind=8)
      end do	  
	  
  else
      do i = 1,nkeep
  	     thickconverge = thickconverge + beta(i)**2/real(nkeep,kind=8)
      end do
  end if
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
	 
	 if(lanczchar=='tx' .and. txchoice==2)then
         call br_transform_basis(eiglvectx, nthick, iter-1)  ! ? iter-1 or iter-2?
		 
	 else

        call br_transform_basis(eiglvec, nthick, iter-1)
		
	end if

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
#ifdef _MPI
        if(vchar == 'n') call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
        if(vchar == 'r') call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
#endif
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
#ifdef _MPI
print*,' SHOULD NOT BE HERE IN THICK_RESTART'
stop
!     call BMPI_FILE_CLOSE(tmpfile,ierror)
#endif
  end if

  return
end subroutine thick_restart_sub_p

!=============================================================================
!
!  thick-restart block lanczos:
!  based on thick_restart_sub_p 
!  -R.M.Z. SDSU 2022
!  after some iterations, take nthick lowest eigenvectors
!  and reintroduce them as lanczos vectors
!  then start up iterations again
! 
!  NB: in normal runs, up to this point have done iter = nthick+dimblock iterations
!
! nthick = # of vectors to keep for thick-restart
! iter = # of lanczos iterations in total (# block iter*dimblock)
! iter_thick = which iteration of thick_restart (starting at 0)
! maxiter = dimension of arrays (# block iter*dimblock) + dimblock
! e(:), eiglev(:,:) : eigenvectors
! vchar (no longer in significant use)
!
! ALSO USED but not passed in arguments: niter, dimension of truncated Hamiltonian to diagonalize
! NOT USED here: nkeep, the # of final eigenvectors
! CALLED BY:
!   lanczos_block
! CALLS:
!  find_lanczos_eigenvalues
!  br_transform_basis
!  br_pull_blreg_from_hist
!  br_push_blreg_to_hist
subroutine block_thick_restart_sub_p(nthick,iter,iter_thick,maxiter,vchar)

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
  use mod_balg
  use localblocks ! for BL
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
  integer(4):: i,j,k,l,m
  integer(4):: mymax  ! for targeted thick-restart
  integer(4):: iclosest
  real(kind=8) :: eclosest
  logical :: verbosetarget = .true.  ! used for debugging excited state thick-restart
  integer, parameter :: tmpfile = 58
  real (kind = lanc_prec) :: da
  integer(kind=basis_prec) :: kl
  integer(4) :: file_mode,file_info,ierror
  integer(4) :: inode
  real(8) :: lanczos_mat_tmp(niter,niter)

  real(kind = 8) :: dtmp

  integer :: ierr
  real(kind=8)  :: thickconverge
  real(kind = 8) :: Uk(dimblock,dimblock),beta_block_tmp(dimblock,dimblock)
  real(kind = 8) :: over_lap(dimblock,dimblock)

beta_block_tmp = 0.d0
Uk = 0.d0
over_lap = 0.d0


  
  if(iproc==0)print*,' In thick restart niter, nthick = ',niter,nthick

!-------------------SOLVE FOR TRUNCATED EIGENVALUES -----------

  if(iter_thick == 0)then
    nthicktmp = 0
  else
    nthicktmp = nthick
  endif
  call find_lanczos_eigenvalues(iter,.true.,nthicktmp) !due to hist with BL we can call this on iter

!--------------------CREATE NEW COEFICIENTS -------------------

		lanczos_mat = 0.d0
		! construct lancozs matrix for restart
	   do i = 1,nthick
		  lanczos_mat(i,i)=e(i) ! main diagonal of E(nthick)
	   end do
	        		
	   do i = 1,nthick/dimblock
        !get correct starting history position for each block
        k =(i-1)*dimblock + 1
        m=(nthick/dimblock)*dimblock + 1 ! the position of block after nthick
        beta_block_tmp = beta_block
        Uk = eiglvec(iter-dimblock+1:iter,k:k+dimblock-1)
        call DGEMM("N","N",dimblock,dimblock,dimblock,1.d0,beta_block_tmp,dimblock,Uk,dimblock,0.d0,over_lap,dimblock)


		  lanczos_mat(k:k+dimblock-1,m:m+dimblock-1)=transpose(over_lap)
		  lanczos_mat(m:m+dimblock-1,k:k+dimblock-1)=over_lap	
		!  lanczos_mat(k:k+dimblock-1,m:m+dimblock-1)=transpose(matmul(beta_block,eiglvec(iter-dimblock+1:iter,k:k+dimblock-1)))
		!  lanczos_mat(m:m+dimblock-1,k:k+dimblock-1)=matmul(beta_block,eiglvec(iter-dimblock+1:iter,k:k+dimblock-1))	
	   end do

  
!------------ CHECK SIZE OF COUPLING -------------------------
! something to try for comparison
  thickconverge = 0.d0
  !if(lanczchar=='tx' .and. txchoice==2)then
      do i = 1,nkeep
        m=(nthick/dimblock)*dimblock + 1 ! the position of block after nthick
  	     thickconverge = thickconverge + lanczos_mat(i,m+dimblock-1)**2/real(nkeep,kind=8)
      end do	  
	  
  !else
  !    do i = 1,nkeep
  !	     thickconverge = thickconverge + beta(i)**2/real(nkeep,kind=8)
   !   end do
  !end if
  if(iproc==0)print*,' Alternate thick-restart convergence : ',thickconverge

!.................. REWRITE TO FILE...... ! this is done in main bl loop...
!  if(iproc ==0)then
!     rewind(coef_file)
!     print*,' rewriting to file ',coef_file,nthick
!     do i = 1,nthick
!       write(coef_file,'(i10,2(2x,f20.15))')-i,alpha(i), beta(i)
!     end do
!  end if

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
	 
	 if(lanczchar=='tx' .and. txchoice==2)then
      if(iproc .eq. 0) then
         print *, "This feature is not installed with block Lanczos thick-restart"
         stop
     end if
         !call br_transform_basis(eiglvectx, nthick, iter-1)  ! ? iter-1 or iter-2?
		 
	 else

        call br_transform_basis(eiglvec, nthick, iter)
		
	end if

     ! init for next iteration. 
     ! load block_vec1 from the top of the history.

      call br_pull_blreg_from_hist(blk_s(iter/dimblock+1))
      call br_push_blreg_to_hist(blk_s(nthick/dimblock+1)) ! from br_reg to hist k_thick + 1
      !debug
      !print*, "pulling from hist pos", iter/dimblock+1,blk_s(iter/dimblock+1)
      !print *,"pushing breg block to pos", nthick/dimblock+1,blk_s(nthick/dimblock+1)


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
#ifdef _MPI		 
        if(vchar == 'n') call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
        if(vchar == 'r') call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
#endif
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
	  print*,' SHOULD NOT BE HERE '
	  stop
!     call BMPI_FILE_CLOSE(tmpfile,ierror)
  end if

  return
end subroutine block_thick_restart_sub_p





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
             if(iproc==0)print*,' wrong ws ',wp,wn, minwtot,maxwtot
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
#ifdef _MPI			  
              call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
#endif
           end if
        end if
     else
        ! new interface, we say which vector to read and it checks
        call wfn_readeigenvec(wfnfile,frag1, fcomm1_index, vec1,istate,ei,xj,xtt)
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
#ifdef _MPI
        call BMPI_ALLREDUCE(wfrac(:), SIZE(wfrac), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place reduce
#endif
     end if
     if(iproc==0)then
        write(6,'(i4,4x,10f7.2)')istate,(Wfrac(wlist(n))*100.,n=1,min(nWvals,10))
        if(writeout)write(resultfile,'(i4,4x,10f7.2)')istate,(Wfrac(wlist(n))*100.,n=1,min(10,nWvals))
		if(nWvals> 10)then
	        write(6,'(i4,4x,10f7.2)')istate,(Wfrac(wlist(n))*100.,n=11,nWvals)
	        if(writeout)write(resultfile,'(i4,4x,10f7.2)')istate,(Wfrac(wlist(n))*100.,n=11,nWvals)
			
			
		end if
     end if
  end do  ! istate

  return
end subroutine Wcounter

!======================================================================
!  added 7.9.10; computes 1/(Egf - H ) on pivot
!
!
subroutine resolvent(iter,ddot)
	
	implicit none
	integer, intent(in) :: iter
	real(kind=8) :: htmp(niter+2,niter+2),b(niter+2,1)
	integer :: i
	integer :: ipiv(niter+1)
	integer :: info
	real(8), intent(out) :: ddot
	real(8) :: norm1,norm2
	
	htmp = 0.d0

!........ SET UP MATRIX TO BE INVERTED....	
	do i = 1,iter
		htmp(i,i)=Egf-alpha(i)
		if(i==iter)cycle
		htmp(i+1,i)=-beta(i)
		htmp(i,i+1)=-beta(i)
	end do
	b = 0.d0
	b(1,1)=1.d0
	call dgesv(iter,1,htmp,niter+2,ipiv,b,niter+2,info)
	
	if(info/=0)then
		print*,' error in inverting E- H in subroutine resolvent '
		stop 919
	end if
	
	do i = 1,iter
		gvec(i)=b(i,1)
	end do
	ddot = 0.d0
	if(iter ==1)return
!........... CHECK FOR CONVERGENCE.............	
	norm1= 0.0
	norm2 = 0.0
	do i = 1,iter+2
		ddot = ddot + gvec(i)*gvec_old(i)
		norm1 = norm1+gvec(i)*gvec(i)
		norm2 = norm2+gvec_old(i)*gvec_old(i)
	end do
	
!	print*,' testing ',ddot,norm1,norm2
	ddot = abs(ddot) / sqrt(norm1*norm2)
	
!........... STORE CURRENT SOLUTION AS 'OLD' VECTOR.....	
	gvec_old = 0.d0
	do i = 1,iter+2
		gvec_old(i)=gvec(i)
	end do
	
	return
	
end subroutine resolvent

!======================================================================
!  added 7.9.10; computes 1/(Egf - H ) on pivot
!  complex version of resolvent, added 7.10.6
!
subroutine zresolvent(iter,ddot)
	
	implicit none
	integer, intent(in) :: iter
	complex(kind=8) :: htmp(niter+2,niter+2),b(niter+2,1)
	integer :: i
	integer :: ipiv(niter+1)
	integer :: info
	real(8), intent(out) :: ddot
	real(8) :: norm1,norm2
	
	htmp = 0.d0

!........ SET UP MATRIX TO BE INVERTED....	
	do i = 1,iter
		htmp(i,i)=cmplx(Egf-alpha(i),Egfi)
		if(i==iter)cycle
		htmp(i+1,i)=-cmplx(beta(i),0.d0)
		htmp(i,i+1)=-cmplx(beta(i),0.d0)
	end do
	b = (0.d0,0.d0)
	b(1,1)=(1.d0,0.d0)
	call zgesv(iter,1,htmp,niter+2,ipiv,b,niter+2,info)
	
	if(info/=0)then
		print*,' error in inverting E- H in subroutine resolvent '
		stop 919
	end if
	
	do i = 1,iter
		gvec(i)=real(b(i,1))
		gvec_i(i)=aimag(b(i,1))
	end do
	ddot = 0.d0
	if(iter ==1)return
!........... CHECK FOR CONVERGENCE.............	
	norm1= 0.0
	norm2 = 0.0
	do i = 1,iter+2
		ddot = ddot + gvec(i)*gvec_old(i)+gvec_i(i)*gvec_i_old(i)
		norm1 = norm1+gvec(i)*gvec(i)+ gvec_i(i)*gvec_i(i)
		norm2 = norm2+gvec_old(i)*gvec_old(i)+  gvec_i_old(i)*gvec_i_old(i)
	end do
	
!	print*,' testing ',ddot,norm1,norm2
	ddot = abs(ddot) / sqrt(norm1*norm2)
	
!........... STORE CURRENT SOLUTION AS 'OLD' VECTOR.....	
	gvec_old = 0.d0
	gvec_i_old = 0.d0
	do i = 1,iter+2
		gvec_old(i)=gvec(i)
		gvec_i_old(i)=gvec_i(i)
	end do
	
	return
	
end subroutine zresolvent

!======================================================================
end module lanczos_util

!/////////////////////////////////////////////////////////////////////

module convergence
	use lanczos_info
	use lanczos_util
	use precisions
	use flagger
	implicit none
    integer :: ncheck_ex        ! # of extra states to check for convergence
    integer :: thick_restart_factor , thick_restart_min

!...... CONVERGENCE IN LANCZOS........
  logical :: ask_converge_test = .false.        !  flag to allow question about convergence  
  integer :: converge_test_def = 0 
  integer :: converge_test                     
!....... THERE ARE SEVERAL OPTIONS FOR CONVERGENCE TEST.....
!        FOR LANCZOS BETWEEN TWO ITERATIONS:
!
!    CONVERGE_TEST = 0   test on average difference in energies
!                  = 1   test on max difference     in energies
!                  = 2   test on avg difference in wavefunctions
!                  = 3   test on max difference in wavefunctions
!    NOTE: 2,3 are good for observables, but could be problematic 
!    if there are degeneracies in the spectrum
!
   real,parameter :: converge_tol0_def =  1.0e-3
   real,parameter :: converge_tol1_def =  5.0e-4
   real,parameter :: converge_tol2_def =  1.0e-5 
   real,parameter :: converge_tol3_def =  2.0e-5
   real :: converge_tol	

!   real    :: ediff0,ediff1
   real(8)     :: vdiff2,vdiff3,vdot
	
contains
	
!======================================================================
!
!
!  UNITED CONVERGENCE TEST
!  ADDED in 7.8.8 by CWJ
!

subroutine convergence_test(iter,i_thick_restart,finished)
	use io,only:greenflag
	use nodeinfo
	use butil_mod
	use bmpi_mod
	use localblocks,only:dimblock,block_flag
	use io,only:complex_green
	use menu_choices
	implicit none
	integer,intent(IN) :: iter
	integer,intent(IN) :: i_thick_restart
	
	integer :: green_skip_iter = 5 ! check convergence on every 5th one
	integer :: green_skip_print = 10
!	real :: green_tol = 1.e-5
	logical   :: finished
	integer :: i,j,k,ierr
    real(kind=egv_prec), pointer ::  epoint(:)
	
	real(8) :: ddot
	
	
	if(greenflag .and. mod(iter,green_skip_iter)==0 .and. menu_char/='pv')then
		if(complex_green)then
			call zresolvent(iter,ddot)
			
		else
		    call resolvent(iter,ddot)
	    end if
		if(mod(iter,green_skip_print)==0 .and. iproc==0)then
			print*,iter,' Convergence of resolvent on pivot: ',ddot
			
		end if
		if(abs(ddot-1.0)< green_tol)then
			finished=.true.
			if(iproc==0)print*,' Within tolerance (green_tol) of ',green_tol
		end if
		
		if(iter==niter)then
			if(iproc==0)print*,' Did not reach desired tolerance (green_tol) of ',green_tol
			finished=.true.
		end if
#ifdef _MPI
        call BMPI_BCAST(finished,1,0,MPI_COMM_WORLD,ierr) 
#endif		
		if(finished)return
		
	end if
	
	if(lanczchar /= 'tx' .or. txchoice==1)then
		epoint =>e
	else
		epoint => etx
	end if

	 
!----------------- NORMAL CHECK -------------------------	 	 
     if ( iter == niter ) then
        finished = .true.
        
        if ( .not. fixed_iter .and. iproc == 0 ) then
           write(6,*)' Did not converge fully '
        end if
		return
     end if
     if ( .not. fixed_iter ) then
        if ( iter > ncheck  .and. (.not.skipeig .or. (skipeig .and. flyeig) .or. & 
		   (block_flag==1 .and. iter > dimblock))) then
!................ CHECK DIFFERENCES IN ENERGY..........................
           ediff0 = 0.0
           ediff1 = 0.0
!           do j = 1,ncheck
           do j = 1,ntest_converge  ! modified in 7.11.2
			   ! this modification allows for faster test of convergence
              ediff0 = ediff0 + abs(eold(j) - epoint(j) )
              ediff1 = bmax(ediff1, abs(  real( eold(j) -epoint(j),kind=4) ) )
           end do
!           ediff0 = ediff0/sqrt(float(nkeep+1))*(1+i_thick_restart)
           ediff0 = ediff0/sqrt(float(ntest_converge+1))*(1+i_thick_restart)

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
           if ( (flyeig .or. (block_flag==1 .and. iter > dimblock)  ).and. iproc == 0 ) then
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
           eold(j) = epoint(j)
           if(altconverge)then
                do k = 1,iter
                  eiglvecold(k,j) = eiglvec(k,j)
                end do
           end if
        end do
		
     end if
#ifdef _MPI
     call BMPI_BCAST(finished,1,0,MPI_COMM_WORLD,ierr) 
#endif
     return	
	
end subroutine convergence_test

end module convergence


