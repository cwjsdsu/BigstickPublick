!*----------------------------------------------------------------
!   Module: mod_balg -- Block Lanczos Subroutine Library
!   -Ryan M Zbikowski 2018
!
!   Purpose: this file contains all the main subroutines called within
!            blanczos_block.f90/blanczosmain.f90's block lanczos routine
!
!======================================================================
!   SUBROUTINES IN THIS FILE
!
!   initialize_pivot_block
!   setup_hist_index
!   br_retrieve_hist_in_place
!   br_add2hist_in_place
!   br_remove_prev_overlap_in_placev2
!   br_block_dot_in_place
!   grab_alpha_blockv2
!   compute_W
!   spectral_decomp
!   add_alpha_Tmat
!   add_beta_Tmat
!   read_pivot_block
!    * Debugging/printing routines at EOF *
!*----------------------------------------------------------------

module mod_balg
use mod_reorthog
use precisions
implicit none
integer, allocatable :: blk_s(:),blk_e(:)
real(kind=8), allocatable :: alpha_block(:,:),beta_block(:,:)
real(kind=8), allocatable ::  lanczos_mat(:,:), ev(:)

real(kind=lanc_prec), allocatable :: block_reg(:)   ! working 'register' of blocks
                                                    ! analog to br_reg in breothog.f90

integer, allocatable :: block_hsendcounts(:)   ! width for each hrank
integer, allocatable :: block_hdispls(:)       ! offset in FRAGMENT (so kind=4)for each hrank					
!--- ADDED IN 7.9.4 FOR BLOCK STRENGTH OPTION---

real(kind=8),allocatable :: orthonorm_pivot(:,:)								

contains



!*----------------------------------------------------------------
!
!  ORIGINAL VERSION-- REPLACED 
!  see routine in blanczoslib.f90
!
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

subroutine initialize_pivot_block_orig
use localblocks
use mod_reorthog
implicit none
integer(kind=4) :: i

call br_set_histpos(0)
call br_reset_histbuf()
do i = 1,dimblock
    call random_restart_with_history
    call br_add2hist(1)
enddo

end subroutine initialize_pivot_block_orig

! ---------------------------------------------
!
!  ADDED in 7.8.7: add a random vector for pivot
!
!  various pivot options
!  pivchoice = 0  random
!  pivchoice = 1  alternating pivot
!
subroutine add_pivot_br_reg(pivchoice,ivec)
	use mod_reorthog
	use localblocks
	use basis
    use sporbit,only:allsameW
    use W_info
    use sectors
	use nodeinfo
	implicit none
	integer(4) :: pivchoice,ivec
!	integer (kind=basis_prec):: dimsize
	real :: rv
	real :: mynorm
    integer(kind=basis_prec) :: i
	integer :: iphase
    integer (kind=basis_prec) :: startvec,stopvec,npxsd,nnxsd
    integer :: isector,ncsector,icsector,csector
	integer :: icount
	
	select case (pivchoice)
	
	case (0)
  	    mynorm = 1.0/sqrt(float(dimbasis))  ! rough normalization to simplify

	    do i = ostart,ostop
		   call random_number(rv)
	 	   br_reg(i)= (rv-0.5)*mynorm
	    end do
	
	case (1)
	    br_reg=0.d0
		mynorm = sqrt(float(dimblock)/float(dimbasis))  ! rough normalization to simplify
		
	    do i = ostart+ivec-1,ostop,dimblock
	 	   br_reg(i)= mynorm
	    end do
	
	
	case (2)
	    br_reg=0.d0
		mynorm = sqrt(float(dimblock)/float(dimbasis))  ! rough normalization to simplify
		iphase = 1
	    do i = ostart+ivec-1,ostop,dimblock
	 	   br_reg(i)= mynorm*iphase
		   iphase = -iphase
	    end do	
		
!............. WHEN USING W-BLOCKS, FILL ONLY SMALLEST........		
	case (3)
!	print*,' OSTART STOP ',ostart,ostop
    br_reg=0.d0
	mynorm = sqrt(float(dimblock)/float(dimbasis))  ! rough normalization to simplify
	iphase = 1
	   if(allsameW)then
   	    do i = ostart+ivec-1,ostop,dimblock
   	 	   br_reg(i)= mynorm*iphase
   		   iphase = -iphase
   	    end do	
	   else
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
				   
				   startvec = max(startvec,ostart)
				   stopvec  = min(stopvec,ostop)
				   if(stopvec < startvec)cycle
!				   print*, ' STOP START ',startvec,stopvec

!........ NOW LOOP OVER THESE SECTORS................				   
                   do i = startvec,stopvec
					   call random_number(rv)
					   
                      br_reg(i) = real( 1.0 *iphase*mynorm*rv, kind=lanc_prec)
					  iphase = -iphase
                   end do
				   
			   end do
				 			   
		   end do

	   end if
	   
   	case (44)
   	    br_reg=0.d0
!   		mynorm = sqrt(float(dimblock)/float(dimbasis))  ! rough normalization to simplify
   		iphase = 1
		
		br_reg(ivec)= 1.0
		
!   	    do i = ostart+ivec-1,ostop,dimblock
   	 	   !mynorm*iphase
!   		   iphase = -iphase
!   	    end do	

		
    end select
	
	return

end subroutine add_pivot_br_reg
!*----------------------------------------------------------------
!  Subroutine: setup_hist_index
!
!  Purpose: Creates two arrays to keep track of starting and
!           termination indexs of vblocks within history
!
!  Parameters:
!      ivec(IN) - dimblock* # of block iterations ie. the dimension of the lanczos matrix
!                 the routine then computes the start stop locations based on the spacing
!                 of dimblock.
!
!  Returns: blk_s (global) list of starting positions in hist.
!           blk_e (global) list of end positions in hist.
!
!*-------------------------------------------------------------------

subroutine setup_hist_index(ivec)
use localblocks
implicit none
integer :: i,ivec
do i=1,ivec
    blk_s(i)= 1 + (i - 1)*dimblock
    blk_e(i) = i*dimblock
enddo
end subroutine setup_hist_index

!*----------------------------------------------------------------
!  Subroutine: br_retrieve_hist_in_place
!
!  Purpose:  Loads a single vector from history into the rowspace block_1
!
!
!  Parameters:
!      ivec(IN) - row index for block_1
!      pos(IN) - vector position in history
!
!
!  Returns: N/A
!
!  MPI USES:  BMPI_GATHERV followed by MPI_BCAST
! cf br_restore_vec in breorthog.f90
!
!*-------------------------------------------------------------------

subroutine br_retrieve_hist_in_place(ivec,pos)
   use mod_reorthog
   use localblocks
   use precisions
   use localvectors
   use nodeinfo
   use bmpi_mod
   implicit none
   integer :: pos, ivec
   integer(kind=basis_prec) :: istate
   integer :: ierr

   call br_check_setup()
   call br_set_histpos(pos)
   if(pos < 1 .or. pos > br_histpos) then
      if(iproc == 0) print *, "br_retrieve_hist:  position ", pos, " is out of range (1:", br_histpos, ")"
      stop 23
   end if

! block_1(ivec,:)= br_histbuf(ostart:ostop,br_histpos)
!print*,' retrieving ',iproc,ivec,pos,ostart,ostop

   do istate = ostart,ostop
	   blockvec_1(ivec + (istate-1)*dimblock) = br_histbuf(istate,br_histpos)
   end do
   
!   call BMPI_REDUCE(blockvec_1, size(blockvec_1), MPI_SUM, 0, fcomm1, ierr) ! in place
#ifdef _MPI
   call BMPI_ALLREDUCE(blockvec_1, size(blockvec_1), MPI_SUM,  fcomm1, ierr) ! in place
#endif
!   do istate = v1s,v1e
!	   if(iproc==0)print*,istate,ivec,blockvec_1(ivec + (istate-1)*dimblock), ' retrieve '
!   end do
   return

end subroutine br_retrieve_hist_in_place

!*----------------------------------------------------------------
!  Subroutine: br_retrieve_hist_in_place
!
!  Purpose:  Loads a block of vectors from history into the rowspace block_1
!
!
!  Parameters:
!      ivec(IN) - row index for block_1
!      pos(IN) - vector position in history
!
!
!  Returns: N/A
!
!  MPI USES:  BMPI_GATHERV followed by MPI_BCAST
! cf br_restore_vec in breorthog.f90
!
!*-------------------------------------------------------------------

subroutine br_pull_block_from_hist(posin)
   use mod_reorthog
   use localblocks
   use precisions
   use localvectors
   use nodeinfo
   use bmpi_mod
   implicit none
   integer,intent(in) :: posin
   integer :: pos, ivec
   integer(kind=basis_prec) :: istate
   integer :: ierr

   pos=posin
   call br_check_setup()
   call br_set_histpos(pos)
   if(pos < 1 .or. pos > br_histpos) then
      if(iproc == 0) print *, "br_retrieve_hist:  position ", pos, " is out of range (1:", br_histpos, ")"
      stop 23
   end if

! block_1(ivec,:)= br_histbuf(ostart:ostop,br_histpos)

!if(iproc==0)print*,' start list ',br_ostoplist

!if(1+(ostart-1)*dimblock < 1+dimblock*(v1s-1) )print*,' problem start ',iproc, ostart,v1s
!if(ostop > v1e)print*,' problem end ',iproc, ostop,v1e
   do ivec = 1,dimblock
       do istate = ostart,ostop
	      blockvec_1(ivec + (istate-1)*dimblock) = br_histbuf(istate,br_histpos)
!		  if(istate==1 .and. ivec==1)print*,' here ',br_histbuf(istate,br_histpos),br_histpos
       end do
       pos = pos+1
	   call br_set_histpos(pos)
	   if(pos < 1 .or. pos > br_histpos) then
	      if(iproc == 0) print *, "br_retrieve_hist:  position ", pos, " is out of range (1:", br_histpos, ")"
	      stop 24
	   end if
	   
   end do
   
!   call BMPI_REDUCE(blockvec_1, size(blockvec_1), MPI_SUM, 0, fcomm1, ierr) ! in place

#ifdef _MPI
   call BMPI_ALLREDUCE(blockvec_1, size(blockvec_1), MPI_SUM,  fcomm1, ierr) ! in place
#endif
!   do ivec = 1,dimblock
!   do istate = v1s,v1e
!	   if(iproc==0)print*,istate,ivec,blockvec_1(ivec + (istate-1)*dimblock), ' retrieve '
!   end do
!end do
   return

end subroutine br_pull_block_from_hist

!*----------------------------------------------------------------
!  Subroutine: br_pull_block_reg_from_hist
!
!  Purpose:  Loads a block of vectors from history into 
!            the intermediary structure block_reg
!
!
!  Parameters:
!      posin - vector position in history
!
!  Returns: N/A
!
!*-------------------------------------------------------------------

subroutine br_pull_blreg_from_hist(posin)
   use mod_reorthog
   use localblocks
   use precisions
   use localvectors
   use nodeinfo
   use bmpi_mod
   implicit none
   integer,intent(in) :: posin
   integer :: pos, ivec
   integer(kind=basis_prec) :: istate
   integer :: ierr

   pos=posin
   call br_check_setup()
   call br_set_histpos(pos)
   if(pos < 1 .or. pos > br_histpos) then
      if(iproc == 0) print *, "br_retrieve_hist:  position ", pos, " is out of range (1:", br_histpos, ")"
      stop 23
   end if

   do ivec = 1,dimblock
       do istate = ostart,ostop
	      block_reg(ivec + (istate-1)*dimblock) = br_histbuf(istate,br_histpos)
       end do
       pos = pos+1
	   call br_set_histpos(pos)
	   if(pos < 1 .or. pos > br_histpos) then
	      if(iproc == 0) print *, "br_retrieve_hist:  position ", pos, " is out of range (1:", br_histpos, ")"
	      stop 24
	   end if
	   
   end do

   return

end subroutine br_pull_blreg_from_hist
!*----------------------------------------------------------------
!  Subroutine: br_push_block_reg_to_hist
!
!  Purpose:  Loads a block of vectors into history from 
!            the intermediary structure block_reg
!
!
!  Parameters:
!      posin - vector position in history
!
!  Returns: N/A
!
!*-------------------------------------------------------------------

subroutine br_push_blreg_to_hist(posin)
   use mod_reorthog
   use localblocks
   use precisions
   use localvectors
   use nodeinfo
   use bmpi_mod
   implicit none
   integer,intent(in) :: posin
   integer :: pos, ivec
   integer(kind=basis_prec) :: istate
   integer :: ierr

   pos=posin
   call br_check_setup()
   call br_set_histpos(pos)
   if(pos < 1 .or. pos > br_histpos) then
      if(iproc == 0) print *, "br_retrieve_hist:  position ", pos, " is out of range (1:", br_histpos, ")"
      stop 23
   end if

   do ivec = 1,dimblock
       do istate = ostart,ostop
	      br_histbuf(istate,br_histpos)=block_reg(ivec + (istate-1)*dimblock) 
       end do
       pos = pos+1
	   call br_set_histpos(pos)
	   if(pos < 1 .or. pos > br_histpos) then
	      if(iproc == 0) print *, "br_retrieve_hist:  position ", pos, " is out of range (1:", br_histpos, ")"
	      stop 24
	   end if
	   
   end do

   return

end subroutine br_push_blreg_to_hist
!*----------------------------------------------------------------
!  Subroutine: br_add2hist_in_place
!
!  Purpose: takes a vector from the rowspace of block_2 and adds it
!           history.
!
!
!  Parameters:
!      ivec(IN) - row index for block_2
!      pos(IN) - vector position in history
!
!
!  Returns: N/A
!
!
!*-------------------------------------------------------------------

subroutine br_add2hist_in_place(ivec,pos)
use mod_reorthog
use localblocks
use precisions
use localvectors
implicit none
integer :: pos, ivec
integer(kind=basis_prec) :: istate

call br_check_setup()
call br_set_histpos(pos-1)
br_histpos = br_histpos + 1 ! make room for new data
if(br_histpos > UBOUND(br_histbuf, 2)) then
if(iproc==0) print *, "Trying to save too many block Lanczos vectors.  Buffer depth=", UBOUND(br_histbuf, 2),pos,br_histpos
stop 2333
end if
!br_histbuf(ostart:ostop, br_histpos) = block_2(ivec,:)

do istate = ostart,ostop
	br_histbuf(istate, br_histpos)=blockvec_2(ivec + (istate-1)*dimblock)
	
end do 
return
end subroutine br_add2hist_in_place

!*----------------------------------------------------------------
!  Subroutine: br_remove_prev_overlap_in_placev2
!
!  Purpose: Performs the Generalized Projection
!           vec(jvec) = vec(jvec) - vec(ivec) *( vec(ivec).vec(jvec))
!           within history.
!
!  Called by: grab_alpha_blockv2
!
!  Parameters:
!      ivec(IN) - position of vector within history
!      jvec(IN) - position of 2nd vector within history
!      a(OUT) - removed overlap
!
!
!  Returns: a
!*-------------------------------------------------------------------
subroutine br_remove_prev_overlap_in_placev2(ivec,jvec,a)
   use bmpi_mod
   use nodeinfo
   implicit none
   integer,intent(in) :: ivec,jvec
   integer :: ierr
   real(kind=8) :: a
   real(kind=8) :: dp
   integer(kind=basis_prec) :: istate

   call br_check_setup()
   dp =0.d0
   dp = dot_product8(br_histbuf(ostart:ostop,ivec),br_histbuf(ostart:ostop,jvec))

   a=dp
#ifdef _MPI
   call BMPI_ALLREDUCE(dp, a, 1, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif   
!   br_reg = br_reg - a * br_histbuf(ostart:ostop, br_histpos)
! tried using pointer into br_histbuf,   Intel 12.0 fails.
!$omp parallel do &
!$omp    private(istate) &
!$omp    firstprivate(a, ostart, ostop) &
!$omp    shared(br_reg)


   do istate = ostart, ostop
      br_histbuf(istate,jvec) = br_histbuf(istate,jvec) - real(a * br_histbuf(istate,ivec), kind=lanc_prec)
   end do

!omp end parallel do
   return
end subroutine br_remove_prev_overlap_in_placev2

!*----------------------------------------------------------------
!  Subroutine: br_block_dot_in_place
!
!  Purpose: Performs dot product between two vectors within history
!
!  Called by: spectral_decomp
!
!  Parameters:
!      ivec(IN) - position of vector within history
!      jvec(IN) - position of 2nd vector within history
!      a(OUT) - overlap
!
!  Returns: a
!
! CALLS:  br_check_setup
!         dot_product8
!*-------------------------------------------------------------------

subroutine br_block_dot_in_place(ivec,jvec,a)
   use bmpi_mod
   implicit none
   integer,intent(in) :: ivec,jvec
   integer :: ierr
   real(kind=8) :: a
   real(kind=8) :: dp
   integer(kind=basis_prec) :: istate

   call br_check_setup()
   dp = dot_product8(br_histbuf(ostart:ostop,ivec),br_histbuf(ostart:ostop,jvec))
   a=dp
#ifdef _MPI
   call BMPI_ALLREDUCE(dp, a, 1, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   return
end subroutine br_block_dot_in_place

!*----------------------------------------------------------------
!  Subroutine: grab_alpha_blockv2
!
!  Purpose: Calls br_remove_prev_overlap_in_placev2 and uses the returned
!            overlap to construct alpha_block
!
!  Parameters:
!      pos(IN) - The current block iteration. The overlap is removed from
!                the next block *ahead(pos + 1)* in history. We will not need to
!                perform reorthogonalization on the current block after this.
!
!  Returns: N/A
!*-------------------------------------------------------------------

subroutine grab_alpha_blockv2(pos)
use localblocks
use mod_reorthog
integer :: pos, i, j, m, n
real(kind = 8) :: ov_l

m = 1

do i = blk_s(pos), blk_e(pos)
    n = 1
    do j = blk_s(pos+1), blk_e(pos+1)
    call br_remove_prev_overlap_in_placev2(i,j,ov_l)
    alpha_block(m,n) = ov_l
    n = n + 1
    enddo
m = m + 1
enddo

end subroutine grab_alpha_blockv2

!*----------------------------------------------------------------
!  Subroutine: compute_W
!
!  Purpose: enforces reorthogonalization on block farthest ahead in history
!           against the remaining blocks within history; ignoring the current.
!           This routine is called before the spectral decomposition.
!
!
!  Parameters:
!      pos(IN) - The current block iteration. The overlap is removed from
!                the next block *ahead (pos + 1)* in history. We will not need to
!                perform reorthogonalization on the current block after this.
!
!  Returns: N/A
!*-------------------------------------------------------------------
subroutine compute_W(pos)
use localblocks
use mod_reorthog
implicit none
integer:: i, pos

do i = blk_s(pos+1), blk_e(pos+1)

    call br_orthogonalize_in_place(i,1,blk_e(pos-1))

enddo

end subroutine compute_W

!*----------------------------------------------------------------
!  Subroutine: spectral_decomp
!
!  Purpose: This routine factorizatizes a block within history *ahead (pos + 1)*
!           of the current given block postion
!           to produce the next krylov subspace, add it to history and to produce a new beta block
!           H Vn = Vn-1 Bn-1 + Vn An + Vn+1 Bn
!           The process is a follows:
!           H Vn = W
!           An = Vn^T W
!           W - VnAn = W'   
!               <--- here is where this routine comes in
!           Overlap = W'^T W'  which is dimension dimblock x dimblock
!           Bn^T Bn = Overlap =  U lamda**2 U^T;      
!              U is a unitary matrix of eigenvectors of Overlap  note that U^T Overlap U = lambda**2 
!           Bn = lambda U^T
!           W’ = Vn+1Bn
!           W’ = Vn+1(LamdaU^T)
!           W’ULamda^-1= Vn+1
!
!
!  NOTE ADDED IN 7.9.4 for block-strength option
!  Need to orthonormlize initial block by this routine
! 
!  store the beta_block = lambda U^T 
!
!
!  Parameters:
!      pos(IN) - The current block iteration. The overlap matrix is computed from W' formed
!      after calling compute_W on the block farther ahead within history
!
!
!  Returns: N/A
!
!  CALLED BY:
!
! SUBROUTINES CALLED: br_block_dot_in_place
!      DSYEV
!      DGEMM
!      dot_product
!*-------------------------------------------------------------------
subroutine spectral_decomp(pos)
use localblocks
use mod_reorthog
use io
use basis
implicit none
integer :: pos, i, j, m, n, info, mrk, k

real(kind = 8) :: eigv(dimblock),over_lap(dimblock,dimblock),lambda(dimblock,dimblock),U(dimblock,dimblock)
real(kind = 8) :: work4(dimblock*(3+dimblock/2)),temp2(dimblock,dimblock),sap(dimblock,dimblock)
real(kind = 8) :: temp1(dimblock,dimblock)
real(kind = 8) :: temp3(dimblock), temp4(dimblock)
real(kind=8) :: d_out, dsqrtnorm,dsqrtnorminv

lambda = 0.d0
temp1 = 0.d0
temp2 = 0.d0
sap = 0.d0
U = 0.d0
over_lap = 0.d0

!print*, " Entering Spectral Decomp"

m = 1
! compute overlap matrix W'
do i = blk_s(pos+1), blk_e(pos+1)
    n = 1
    do j = blk_s(pos+1), blk_e(pos+1)
        call br_block_dot_in_place(i,j,d_out)
        over_lap(m,n) = d_out
        n = n + 1
    enddo
   m = m + 1
enddo

! store overlap for later
temp2 = over_lap

if(pos==0 .and. iproc==0)then
	print*,' norms of pivot vectors '
!	do i = 1,dimblock
		write(6,*)(over_lap(j,j),j=1,dimblock)
!	end do
    print*,' '
end if
!write(logfile,*) "------------overlap Matrix--------------"
!call printmat8(dimblock,dimblock, over_lap)

! compute eigenpairs of overlap matrix
call DSYEV( 'V','U', dimblock, over_lap,dimblock, eigv, work4, dimblock*(3+dimblock/2), info )
!write(logfile,*) "------------eig overlap Matrix--------------"
!call printmat8(dimblock,dimblock, over_lap)
!write(logfile,*) "------------eigs of overlap Matrix--------------"
!call printmat8(dimblock,1,eigv)

!form unitary matrix
! NOT NECESSARY--SHOULD BE NORMALIZED ALREADY
!do j = 1, dimblock

 !   U(:,j) = over_lap(:,j) * 1.d0/NORM2(over_lap(:,j))

!enddo

! compute lambda**2
!lambda = matmul(transpose(U),matmul(temp2,U))   ! DON'T USE MATMUL! SLOW!  use dgemm


! factor diagonal matrix
!  BUT IF DONE CORRECTLY LAMBDA SHOULD ALREADY BE DIAGONAL
!do i = 1,dimblock
!    do j = 1, dimblock
!        if( i .eq. j) then
!            lambda(i,j) = sqrt(lambda(i,j))
!       endif
!    enddo
!enddo

! beta_block = sqrt(diag(norm))*U ^T 

!write(logfile,*) "------------ lambda Matrix--------------"
!call printmat8(dimblock,dimblock, lambda)
! compute beta block
!call DGEMM("N","T",dimblock,dimblock,dimblock,1.d0,lambda,dimblock,U,dimblock,0.d0,sap,dimblock)

!if(pos==0 .and. iproc ==0)print*,' norm eigenvalues ', eigv 

do i = 1,dimblock
	dsqrtnorm = sqrt(eigv(i))
	dsqrtnorminv = 1.d0/dsqrtnorm
	do j = 1,dimblock
		beta_block(i,j)= dsqrtnorm*over_lap(j,i)
		sap(j,i) =over_lap(j,i)* dsqrtnorminv
	end do
end do


!beta_block = sap !matmul(lambda,transpose(U))

!compute lambda ^-1
!do i = 1, dimblock
 !   do j = 1, dimblock
  !      if(i .eq. j) then
   !         if(lambda(i,j) .ne. 0) then
    !        lambda(i,j) = 1.d0/(lambda(i,j))  ! BUT NEED TO WORRY ABOUT LAMBDA TINY
     !       endif
     !   endif
  ! enddo
!enddo

!write(logfile,*) "------------ inv lambda Matrix--------------"
!call printmat8(dimblock,dimblock, lambda)

!write(logfile,*) "------------ UUT Matrix--------------"
!call printmat8(dimblock,dimblock, matmul(U,transpose(U)))

!temp1 = matmul(U,lambda)
!sap =  0.d0

! sap = U*1/sqrt(diag(norm))
!call DGEMM("N","N",dimblock,dimblock,dimblock,1.d0,U,dimblock,lambda,dimblock,0.d0,sap,dimblock)
!write(logfile,*) "------------ temp1 Matrix--------------"
!call printmat(dimblock,dimblock, temp1)
!write(logfile,*) "------------ temp3 Matrix--------------"

!print*, "push Q to block_1"
!     mrk = 1
!do i = blk_s(pos+1), blk_e(pos+1)
!            call br_retrieve_hist_in_place(mrk,i)
!     mrk = mrk + 1
!enddo
!write(logfile,*) "block_1 in decomp--------------"
!call printmat(10,dimblock,transpose(block_1))


!block_2 = transpose(matmul(transpose(block_1),sap))

!write(logfile,*) "block_2 in decomp--------------"
!call printmat(10,dimblock,transpose(block_2))

!print*," add block_2 to history"
!     mrk = 1
!do i = blk_s(pos + 1), blk_e(pos + 1)
!            call br_add2hist_in_place(mrk,i)
!     mrk = mrk + 1
!enddo

!write(logfile,*) "block_2 in decomp--------------"


! perform Vn+1 = W'*U^T*inv(lambda)
do i = ostart, ostop
    m = 1
    temp3 = 0
    do k = blk_s(pos+1),blk_e(pos+1)
        temp3(m) = br_histbuf(i,k)
        m = m + 1
    enddo
            n = blk_s(pos+1)
            do j = 1, dimblock
                temp4 = 0
                d_out = 0.d0
                temp4 = sap(:,j)

                d_out = dot_product(temp3,temp4)
!                if( i .le. 10) then
!                    write(logfile,*) d_out
!                endif
                br_histbuf(i,n) = d_out
           n = n + 1
    enddo
enddo

end subroutine spectral_decomp

!*----------------------------------------------------------------
!  Subroutine: add_alpha_Tmat
!
!  Purpose: adds alpha matrix to lanczos matrix
!
!
!  Parameters:
!      pos(IN) - The current block iteration and the position within the lanczos matrix.
!
!
!  Returns: N/A
!*-------------------------------------------------------------------

subroutine add_alpha_Tmat(pos)
use mod_reorthog
use localblocks
implicit none
integer:: i, pos

i = (pos-1)*dimblock + 1

    lanczos_mat(i:i + dimblock - 1,i:i + dimblock - 1) = alpha_block
	
	return

end subroutine add_alpha_Tmat

!*----------------------------------------------------------------
!  Subroutine: add_beta_Tmat
!
!  Purpose: adds beta matrix to lanczos matrix
!
!
!  Parameters:
!      pos(IN) - The current block iteration and the position within the lanczos matrix.
!
!
!  Returns: N/A
!*-------------------------------------------------------------------
subroutine add_beta_Tmat(pos)
use localblocks
use mod_reorthog
implicit none
integer:: i, pos

 i = (pos-1)*dimblock + 1

lanczos_mat(i:i + dimblock - 1,i + dimblock :i + 2*dimblock - 1) = transpose(beta_block)
lanczos_mat(i + dimblock:i + 2*dimblock - 1,i:i + dimblock - 1) = beta_block

end subroutine add_beta_Tmat


!---------------------------ROUTINES USED FOR DEBUGGING----------------------------------
!Print pretty matrix single precision
subroutine printmat(m,n,a)
use io
implicit none

integer:: m,n,i,j
real(kind=4) a(m,n)

do i = 1,m

write(logfile,*)( a(i,j), j=1,n )
write(logfile,*) ""
enddo
return
end subroutine printmat

!------------------------------------------------------------------------------
!Print pretty matrix double precision
subroutine printmat8(m,n,a)
use io
implicit none

integer :: m,n,i,j
real(kind=8) a(m,n)

do i = 1,m

write(logfile,*)( a(i,j), j=1,n )
write(logfile,*) ""
enddo
return
end subroutine printmat8

! ---------------------------------------------------
! used to check orthogonality within history between two blocks
subroutine check_orthog(pos1,pos2)
use mod_reorthog
use localblocks
use basis
use io
implicit none
integer :: pos1, pos2, i, j, n, m
real(kind = 8) :: olp, temp1(dimblock,dimblock)

m = 1
do i = blk_s(pos2), blk_e(pos2)
    n = 1
    do j = blk_s(pos1), blk_e(pos1)
        olp = dot_product8(br_histbuf(ostart:ostop,i),br_histbuf(ostart:ostop,j))
        temp1(m,n) = olp
        n = n + 1
    enddo
    m = m + 1
enddo


!write(logfile,*) "------------ Check Q1/Q2 Orthogonality after decomp--------------"
!call printmat8(dimblock,dimblock,temp1)

end subroutine check_orthog

subroutine report_alpha_beta(pos)
use localblocks
use lanczos_info
implicit none
integer :: i,j,pos
real :: dherm
write(coef_file,*) "",pos
dherm = 0.0
do i = 1, dimblock
    write(coef_file,*)( alpha_block(i,j), j=1,dimblock),(beta_block(i,j),j=1,dimblock)
    write(coef_file,*) ""
    !write(coef_file,*)( beta_block(i,j), j=1,dimblock)
!	do j = i+1,dimblock
!		dherm = dherm + abs(alpha_block(i,j)-alpha_block(j,i))
!	end do
enddo
!print*,' hermiticity test ',dherm
return

end subroutine report_alpha_beta

!====================================================
!
! in part, analog to setup_mpi_hist in breorthog
!
subroutine setup_mpi_block_reg
	use mod_reorthog
	use localblocks
	use nodeinfo
	implicit none
	integer :: i
!	integer :: br_aerr
		
    allocate(block_reg(1+dimblock*(ostart-1):dimblock*ostop), stat=br_aerr) ! working register
    if(br_aerr /= 0) call memerror("setup_mpi_block_reg 1")


    allocate(block_hsendcounts(0:maxkey), stat=br_aerr) ! working register
    if(br_aerr /= 0) call memerror("setup_mpi_block_reg 2")
	
    allocate(block_hdispls(0:maxkey), stat=br_aerr) ! working register
    if(br_aerr /= 0) call memerror("setup_mpi_block_reg 3")
!... RESCALE FOR GATHERV,SCATTERV	
	do i = 0,maxkey
		block_hsendcounts(i)=br_hsendcounts(i)*dimblock
		block_hdispls(i)=br_hdispls(i)*dimblock  
		
	end do
	
!	if(iproc==0)print*,' displacements ',br_hdispls
!	if(iproc==0)print*,' send counts ',br_hsendcounts
	
	return
	
end subroutine setup_mpi_block_reg

!  takes a block FROM block_reg and puts into blockvec1
!
!  CALLED BY: 
!      lanczos_block
!  USES:  BMPI_GATHERV followed by MPI_BCAST
!
subroutine br_pull_block1_from_reg()
   use bmpi_mod
   use localblocks
   use nodeinfo
   implicit none
   integer :: ierr
   integer :: count

   call br_check_setup()
   ! if(iproc == 0) print *, "KSM:  in br_restore_vec1, after br_check_setup"
   ! Now have to send back and distribute across fragment members
   !
   ! print *, "KSM: iproc=", iproc, ", about to GATHERV, br_hcomm=", br_hcomm, "ostart=", ostart, ", ostop=", ostop
   ! print *, "KSM: iproc=", iproc, ", size(br_reg)=", size(br_reg), ", size(vec1)=", size(vec1)
   count = int((ostop - ostart+1)*dimblock , 4)  ! fragments are limited to 2B size
   if(nproc == 1) then
      blockvec_1 = block_reg
   else
#ifdef _MPI
      call BMPI_GATHERV(block_reg, count, blockvec_1, block_hsendcounts, block_hdispls, 0, br_hcomm, ierr)
#endif
   endif

   ! Now set data out to all members of the fragment
#ifdef _MPI
   call BMPI_BCAST(blockvec_1, size(blockvec_1) , 0, fcomm1, ierr)
#endif
end subroutine br_pull_block1_from_reg

!
! Assemble vec2 onto slices, that is, FROM vec2 ONTO vector br_reg
!
!  CALLED BY: lanczos_block
!
! USES: MPI_SCATTERV
!
subroutine br_push_block2_to_reg()
   use bmpi_mod
   use localblocks
   use nodeinfo
   implicit none
   integer :: ierr
   integer :: count

   call br_check_setup()
   ! Scatter from fragment source to slices
   count = int((ostop-ostart+1)*dimblock,4)
   if(nproc == 1) then
      block_reg = blockvec_2
   else
#ifdef _MPI
      call BMPI_SCATTERV(blockvec_2, block_hsendcounts, block_hdispls, block_reg, count, 0, br_hcomm, ierr)
#endif
   end if
end subroutine br_push_block2_to_reg


!--------------------------------------- + -------------------------------------
!
!  TRANSFERS blockvec1 -> blockvec2 ON DIAGONAL ROOT 
!  then use MPI_BCAST
! 
!  CALLED BY:
!      lanczos_output
!      exactdiag_p, exactdiag_MPI
!
subroutine br_load_block2_from_block1()
   use bmpi_mod
   use localblocks
   use fragments
   implicit none
   integer :: ierr

   if(nfragments == 1) then
      blockvec_2 = blockvec_1
   else
      ! we can transfer from vec2 to vec1 on the diagonal root nodes and then broadcast out
      ! within each fragment using fcomm1
      if(isfragroot) blockvec_2 = blockvec_1
#ifdef _MPI
      call BMPI_BCAST(blockvec_2, size(blockvec_2) , 0, fcomm2, ierr)
#endif
   end if
end subroutine br_load_block2_from_block1

end module mod_balg
