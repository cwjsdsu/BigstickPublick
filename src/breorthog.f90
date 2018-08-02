!====================================================================
!  Reorthogonalization for Lanczos algorithm
!  The Lanczos algorithm produces a sequence of normalized vectors that are
!  mathematically orthogonal, but due to numerical issues later vectors
!  slowly lose orthogonality.
!
!  The process is:
!  1) Each fragment cluster of nodes has a root node from which the fragment
!     is scattered to the history stack slices.
!  2) The new vector is made orthogonal to all vectors in the stack
!  3) The new vector is normalized
!  4) The new vector is added to the stack
!  5) The slices for a fragment are gathered back to the root of the 
!     fragment cluster
!  6) The fragment cluster root broadcasts the fragment vector to all nodes
!     in the cluster.  
!
!  Initial version:  Aug 9, 2014 - Ken McElvain
!
!  SOME IMPORTANT NOTES ON STORING VECTORS IN HISTORY UNDER MPI
!  (distributed memory)  added April 2018 in  7.8.4
!
!  BIGSTICK uses vectors in two operations: mat-vec multiplication
!  and storage and linear algebra on previous Lanczos vectors
!  
!  For mat-vec multiplication, we break up vectors into FRAGMENTS.
!  These are are LARGE as POSSIBLE because the factorization algorithm
!  for mat-vec multiplication is most efficient for large loops.
!  Fragments are defined in module fragments. Each MPI process
!  has a vec1 and vec2, defined in module localvectors, with a fixed 
!  fragment assigned to vec1 and a fixed fragment assigned to vec2.
!  (They generally are not the same fragments.) Multiple MPI processes
!  may have the same initial and final fragments. 
!
!  The start and stop of a fragment is defined by basestart(fragment)/basestop(fragment)
!  which is known by all processes. On a given MPI process the start/end 
!  for vec1 and vec2 are v1s,v1e and v2s,v2e, respectively.
!
!  For Lanczos vector storage/linear algebra, we use a different scheme.
!  These processes are given over to storing and acting on Lanczos vectors.
!  Each fragment is divided into slices.  Any given process 
!  is assigned a specific slice, and it stores the same slice for ALL 
!  Lanczos vectors. This slice is defined by the variables ostart and ostop.
!  Storing slices aids in efficient linear algebra. In essence, when 
!  carrying out dot products, one carries out the partial dot product within
!  a slice, possibly for all vectors simultaneously. Then only that partial 
!  dot product is broadcast and then combined to get the total dot product.
!
!  Each MPI process has a unique slice, defined by a unique ostart /ostop.
!  These slices are stored in the variable br_histbuf
!
!  The working slice of a vector is call br_reg.
!  The routines br_grab_vec1 and br_grab_vec2 take the appropriate slice of vec1/2
!  and put it onto br_reg, respectively
!  The routines br_restore_vec1 and br_restore_vec2 take the slice that is 
!  br_reg and put into vec1 /vec2, respectively
!
!  The routine br_add2hist takes br_reg and adds to br_histbuf
!  The routine br_retrieve_hist takes the slice from br_histbuf and puts into br_reg
!  
!  In otherwords we have 3 layers
!
!               br_histbuf (storage of slices of lanczos vectors)
!                  |                  ^
!                  V                  |
!               br_reg   (a specified slice)
!                  |                  ^
!                  V                  |
!                vec1     or  vec2
!
!  The routine setup_mpi_hist does the following:
!     -- compute ostart and ostop for ALL MPI processes
!     -- assign ostart and ostop for the current MPI process (iproc)
!  If you set verbose_setup=.true. in setup_mpi_hist, it will printout the slices.
!  These are ONLY nontrivial for MPI runs
!
!  Routine br_normalize normalizes br_reg across ALL slices, and returns its normalization
!  Routine br_orthogonalize is a master routine for reorthgonalization
!  Routine dot_product8 takes partial dot products of vectors
!  Routine br_reorthogonalizes uses dot_product8 to take partial dot products between br_reg and a vector in br_hist
!  and then subtracts it
!
!  
!===========================================================================
!  SUBROUTINES IN THIS FILE
!
!   br_check_setup : checks history has been set up
!
!   dot_product8   :  Takes dot product of two lanc_prec vectors with result
!                     in kind=8 precision 
!
!   setup_mpi_hist     :  Assigns slices of fragments to node and sets up communicators
!   setup_mpi_hist_comm:
!   br_grab_vec1
!   br_grab_vec2
!   br_remove_prev_overlap
!   br_normalize
!   br_orthogonalize
!   br_restore_vec1
!   br_restore_vec2
!
!===========================================================================
!
! 
!===========================================================================
module mod_reorthog
   use precisions
   use nodeinfo
   use basis
   use fragments
   use localvectors
   implicit none

   ! br_ostartlist, br_ostoplist, br_ofraglist arrays indexed by iproc
   integer(kind=basis_prec), allocatable :: br_ostartlist(:), br_ostoplist(:)
   integer, allocatable :: br_ofraglist(:)  ! fragment containing (ostart:ostop)
   ! local values for iproc
   integer(kind=basis_prec) :: ostart, ostop
   integer :: ofrag                      ! fragment stored by iproc
   integer br_rank_of_frag1_root            ! used for read/write to disk

   ! sources have to be diagnal nodes, meaning that they map fragment to same fragment
   integer, allocatable ::  br_osourcelist(:)  ! 0 means not a source, >0 means source of that fragment
   ! communicator for moving data from fragment to history
   integer :: br_hcomm
   ! indexed by hrank in br_hcomm  (length is MPI_COMM_SIZE)
   integer, allocatable :: br_hsendcounts(:)   ! width for each hrank
   integer, allocatable :: br_hdispls(:)       ! offset in FRAGMENT (so kind=4)for each hrank

   ! holds slices of previous vectors
   integer :: br_histpos   ! position of top of stack.  starts 0
   integer :: br_histmax   ! last position we can fill
   
   ! br_reg IS CURRENT WORKING VECTOR for reorthogonalization against old vectors
   real(kind=lanc_prec), allocatable :: br_reg(:)   ! working 'register'

   ! br_histbuf IS WHERE WE KEEP THE LANCZOS VECTORS as well as eigenvectors after reconstruction
   ! first index select slot in vector, second index indicates vector
   real(kind=lanc_prec), allocatable, target :: br_histbuf(:,:)  ! ( ostart:ostop, 1:br_histmax,)
   
   ! buffer for dot products of br_reg against br_histbuf
   real(kind=8), allocatable :: br_orthogdot(:), br_orthogdotf(:) ! (1:br_histmax)

   integer :: dbg_u
   logical :: br_setup_done = .false.
   
   integer :: br_aerr
contains


!---------------------------------------------+++------------------------------
!
! CHECKS THAT HISTORY HAS BEEN SET  UP
!
subroutine br_check_setup()
   implicit none

   if(br_setup_done) return
   print *, "Calling mpi/fragment based reorthogonalization routine without setup"
   stop 23
end subroutine br_check_setup
!---------------------------------------------+++------------------------------
!
! Compute dot product for slices of history and new slices
! Accumulation is done in kind=8 even when lanc_prec=4 to avoid
! loss of accuracy
!
!  CALLED BY
!     br_remove_prev_overlap
!     br_normalize
!     br_orthogonalize
!
real(kind=8) function dot_product8(x, y)
   implicit none
   real(kind=lanc_prec):: x(:), y(:)
   integer:: i
   real(kind=8) :: sum

   sum = 0.0
   ! Need to implement openmp loop
!$omp parallel do private(i), reduction(+ : sum)
   do i = 1, SIZE(x)
      sum = sum + real(x(i), kind=8) * real(y(i), kind=8)
   end do
!$omp end parallel do
   dot_product8 = sum
   return
end function dot_product8
!==========================================================
! Fragments have been set up
! storage slices have been assigned to processors.
! Each fragment has one diagional node identified for communications
! with the history stack.   This diagonal node will scatterv the
! Because we picked a diagonal node, a single communicator will work
! for both values of vchar.
!
!  CALLED BY: setup_mpi_hist
!  
subroutine setup_mpi_hist_comm()
   use bmpi_mod
   implicit none

   integer :: color, key, ierr, maxkey, pi, hsize
   integer :: hranklist(0:nproc-1)
   integer :: keylist(0:nproc-1)

   br_setup_done = .false.
   if(.not. useNewReorthog) then
      br_rank_of_frag1_root = 0
      isfragroot = (iproc == 0)
      return
   end if
   br_setup_done = .true.
   hranklist = -1
   keylist = -1

   color = br_ofraglist(iproc)  ! fragment this node is involved with
   ! find source of this fragment
   do pi = 0, nproc - 1
      if(br_osourcelist(pi) == color) then
         hranklist(0) = pi ! map hrank back to irank
         keylist(pi) = 0  ! this will be the hrank in the sub-comm br_hcomm
         exit
      end if
   end do
   key = 1
   do pi = 0, nproc-1
      if(br_ofraglist(pi) == color .and. keylist(pi) < 0) then
         hranklist(key) = pi
         keylist(pi) = key
         key = key + 1
      end if
   end do
   maxkey = key - 1
   key = keylist(iproc)  ! will be rank in br_hcomm

   call BMPI_COMM_SPLIT(icomm, color, key, br_hcomm, ierr)

   ! double check that we counted correctly
   call BMPI_COMM_SIZE(br_hcomm, hsize, ierr)
   if((maxkey+1) .ne. hsize) then
      print *, "comm with color=", color, ", maxkey=", maxkey, ", hsize=", hsize, "  mismatch!"
      stop 23
   end if

   allocate(br_hsendcounts(0:maxkey), stat=br_aerr) ! chunk size
   if(br_aerr /= 0) call memerror("setup_mpi_hist_comm 1")
   allocate(br_hdispls(0:maxkey), stat=br_aerr)     ! chunk offset in fragment
   if(br_aerr /= 0) call memerror("setup_mpi_hist_comm 2")
   do key = 0, maxkey
      pi = hranklist(key)
!        if(iproc == 0) print *, "key=", key, ", pi=", pi
      br_hsendcounts(key) = int(br_ostoplist(pi) - br_ostartlist(pi) + 1, 4)
      ! displacement is from the beginning of the fragment (color)
      ! Assuming that tvec starts at 1
      br_hdispls(key) = int(br_ostartlist(pi) - basestart(color),kind(br_hdispls(1)))   
   end do

   ! if the key (rank in created comm) is 0, then this node is the 
   ! root of the fragment.  Data is scatterv/gatherv'd from this
   ! node to the nodes storing the slices.
   ! We also the same node to broadcast the orthongonalized vector
   ! out to the members of the matvec fragment nodes 
   ! such nodes are rank 0 in br_hcomm,fcomm1 and fcomm2
   isfragroot = keylist(iproc) == 0
   ! Find the root of the first fragment.
   br_rank_of_frag1_root = -1
   do pi = 0, nproc-1
      if(nodal(pi)%ifragment == 1 .and. nodal(pi)%ffragment == 1) then
         br_rank_of_frag1_root = pi
         exit
      end if
   end do
   if(isfragroot .and. frag1 == 1) then
      if(br_rank_of_frag1_root /= iproc) then
         print *, "picked wrong root for read/write"
         stop 23
      end if
   end if
end subroutine setup_mpi_hist_comm
!=============================================================================
! Initialize Lanczos historical vectors data structures
! Issues:
!    (1) slices must be contained within fragments
!    (2) vchar applies
!    (3) Each fragment has a diagonal node mapping a fragment onto itself
!    (4) The diagional node will be the source for the scatterv
!    (5) First priority is that the source node is to store the first slice
!    (6) Other nodes producing the fragment come next
!    (7) Finish with unclaimed nodes
!
! Item (4) is important so that we can guarentee that communicators do not
! overlap.   We will make one set of communicators for pushing from vec1 and another
! set of communicators for pushing from vec2, but both will use the same source
!
! IMPORTANT: This routine defines variables ostart and ostop, which are the start and stop
!  of storing Lanczos vectors in history
!
!  CALLED BY: setup_localvectors
!
subroutine setup_mpi_hist(maxiter)
   use fragments  ! to get nodal
   use bmpi_mod
   implicit none

   integer :: maxiter
   integer :: ierr
   integer(kind=basis_prec) :: npos, epos
   integer :: pcnt, pi, fi, bestpi
   logical :: claimed(0:nproc)
   integer(kind=basis_prec) :: fs(nfragments)
   integer :: II
   integer ::  piinfragcnt(nfragments) ! number of processors assigned to fragment
   integer ::  csize(nfragments)       ! slice size for fragment
   integer ::  whichfrag(0:nproc)    ! assignment for processor
   integer :: bestfrag
   real :: g, bestg
   
   logical :: verbose_setup = .true.

   br_histmax = maxiter+2
   piinfragcnt = 0
   whichfrag = 0

!  if(iproc == 0) print *, "In setup_hist"
   allocate(br_ostartlist(0:nproc-1), stat=br_aerr)
   if(br_aerr /= 0) call memerror("setup_mpi_hist 1")
   allocate(br_ostoplist(0:nproc-1), stat=br_aerr)
   if(br_aerr /= 0) call memerror("setup_mpi_hist 2")
   allocate(br_ofraglist(0:nproc-1), stat=br_aerr)
   if(br_aerr /= 0) call memerror("setup_mpi_hist 3")
   allocate(br_osourcelist(0: nproc-1), stat=br_aerr)
   if(br_aerr /= 0) call memerror("setup_mpi_hist 4")
   ! initialize
   br_osourcelist = 0   ! not source of any fragment yet
   br_ofraglist = 0
   br_ostartlist = 0
   br_ostoplist = -1

   ! First assign one diagonal node for each fragment to be used as root
   do pi = 0, nproc-1
      fi = nodal(pi)%ifragment  ! KSM ifragment is frag1 in iproc, ffragment is frag2
	                            ! CWJ : nodal(:)%ifragment,fragment set in routine setnodaltribution
      if(fi == nodal(pi)%ffragment) then ! node is diagonal, use as source for fragment
         if(piinfragcnt(fi) == 0) then
            br_osourcelist(pi) = fi  ! This process sources fragment fi
            whichfrag(pi) = fi    ! slice for fragment pi
            piinfragcnt(fi) = 1
         end if
      end if
   end do
   ! Now assign unused MPI processes to the best fragments.   The
   ! best fragment is defined as the one that would have the largest slices
   ! given the current process count for that fragment.
   pi = 0
   do while(pi < nproc)
      bestfrag = 0
      bestg = 0
      do fi = 1, nfragments
         g = real(basestop(fi) - basestart(fi) + 1) / piinfragcnt(fi)
         if(g > bestg) then
            bestg = g
            bestfrag = fi
         end if
      end do
      if(bestfrag /= 0) then
         ! find unassigned fragment
         do while(pi < nproc .and.  whichfrag(pi) /= 0)
            pi = pi + 1
         end do
         if(pi < nproc) then
            ! got one, assign it and bump the count
            whichfrag(pi) = bestfrag
            piinfragcnt(bestfrag) = piinfragcnt(bestfrag) + 1
            pi = pi + 1
         end if
      end if
   end do

   ! Now every processor belongs to a fragment and we know the count
   fs = basestart
   ! comp slice size for each fragment
   csize =  int(((basestop - basestart + 1) + piinfragcnt - 1) / piinfragcnt,4)
   do pi =0,nproc-1
      fi = whichfrag(pi)
      if(fi == 0) then
         print *, "whichfrag not set"
         stop 23
      end if
      epos = fs(fi) + csize(fi) - 1  ! end pos of slice
      if(epos > basestop(fi)) epos = basestop(fi) ! trim if too large
      br_ofraglist(pi) = fi        ! fragment being stored
      br_ostartlist(pi) = fs(fi)   ! slice being stored
      br_ostoplist(pi) = epos
      fs(fi) = epos + 1         ! next position needing coverage
   end do

   call BMPI_BARRIER(icomm, ierr)

   ! get local info
   ostart = br_ostartlist(iproc)
   ostop = br_ostoplist(iproc)
   ofrag = br_ofraglist(iproc)

   ! storage for keeping previous iterations
   allocate(br_histbuf(ostart:ostop, 1:br_histmax), stat=br_aerr)
   if(br_aerr /= 0) call memerror("setup_mpi_hist 5")
   allocate(br_orthogdot(1:br_histmax), stat=br_aerr)
   if(br_aerr /= 0) call memerror("setup_mpi_hist 6")
   allocate(br_orthogdotf(1:br_histmax), stat=br_aerr)
   if(br_aerr /= 0) call memerror("setup_mpi_hist 7")
   allocate(br_reg(ostart:ostop), stat=br_aerr) ! working register
   if(br_aerr /= 0) call memerror("setup_mpi_hist 8")

   br_histpos = 0 ! last written position

   call setup_mpi_hist_comm()  ! init  br_hcomm, tables for scatterv
   
   if(verbose_setup .and. iproc==0 )then
	   open(unit=83, file='slices.bigstick',status='unknown')
	   do pi = 0,nproc-1
		   write(83,*)pi,br_ostartlist(pi),br_ostoplist(pi),br_ofraglist(pi)
		   
	   end do
	   close(83)
	   
	   
   end if
	   
   return
end subroutine setup_mpi_hist
!=============================================================================
!  sets top of history stack (location) to pos
!
! CALLED BY:
!   thick_restart_sub_p
!
subroutine br_set_histpos(pos)
   implicit none
   integer, intent(in) :: pos

   if(pos < 0 .or. pos > br_histmax) then
      print *, "setting hist pos out of range: ", pos
      stop 23
   end if
   br_histpos = pos;
end subroutine br_set_histpos

!
! load a history value into (vector) br_reg
!
subroutine br_retrieve_hist(pos)
   implicit none
   integer :: pos

   call br_check_setup()
   if(pos < 1 .or. pos > br_histpos) then
      if(iproc == 0) print *, "br_retrieve_hist:  position ", pos, " is out of range (1:", br_histpos, ")"
      stop 23
   end if
   br_reg = br_histbuf(ostart:ostop, pos)
end subroutine br_retrieve_hist

!
! Assemble vec2 onto slices, that is, FROM vec2 ONTO vector br_reg
!
!  CALLED BY: lanczos_p
!
subroutine br_grab_vec2()
   use bmpi_mod
   implicit none
   integer :: ierr
   integer :: count

   call br_check_setup()
   ! Scatter from fragment source to slices
   count = int((ostop-ostart)+1,4)
   if(nproc == 1) then
      br_reg = vec2
   else
      call BMPI_SCATTERV(vec2, br_hsendcounts, br_hdispls, br_reg, count, 0, br_hcomm, ierr)
   end if
end subroutine br_grab_vec2

! Assemble vec1 onto slices, that is, FROM vec1 ONTO vector br_reg
!
!  CALLED BY: lanczos_p
!
subroutine br_grab_vec1()
   use bmpi_mod
   implicit none
   integer :: ierr
   integer :: count

   call br_check_setup()
   ! Scatter from fragment root nodes to slices
   count = int((ostop-ostart)+1, 4)
   if(nproc == 1) then
      br_reg = vec1
   else
      call BMPI_SCATTERV(vec1, br_hsendcounts, br_hdispls,  br_reg, count, 0, br_hcomm, ierr)
   end if
end subroutine br_grab_vec1

!  takes a vector FROM br_reg and puts into vec2
!
subroutine br_restore_vec2()
   use bmpi_mod
   implicit none
   integer :: ierr
   integer :: count

   call br_check_setup()
   ! Now have to send back and distribute across fragment members
   !
   count = int(ostop - ostart + 1, 4)  ! fragments are limited to 2B size
   if(nproc == 1) then
      vec2 = br_reg
   else
      call BMPI_GATHERV(br_reg, count, vec2, br_hsendcounts, br_hdispls, 0, br_hcomm, ierr)
   end if

   ! Now set data out to all members of the fragment
   call BMPI_BCAST(vec2, size(vec2) ,  0, fcomm2, ierr)

end subroutine br_restore_vec2

!  takes a vector FROM br_reg and puts into vec1

subroutine br_restore_vec1()
   use bmpi_mod
   implicit none
   integer :: ierr
   integer :: count

   call br_check_setup()
   ! if(iproc == 0) print *, "KSM:  in br_restore_vec1, after br_check_setup"
   ! Now have to send back and distribute across fragment members
   !
   ! print *, "KSM: iproc=", iproc, ", about to GATHERV, br_hcomm=", br_hcomm, "ostart=", ostart, ", ostop=", ostop
   ! print *, "KSM: iproc=", iproc, ", size(br_reg)=", size(br_reg), ", size(vec1)=", size(vec1)
   count = int(ostop - ostart + 1, 4)  ! fragments are limited to 2B size
   if(nproc == 1) then
      vec1 = br_reg
   else
      call BMPI_GATHERV(br_reg, count, vec1, br_hsendcounts, br_hdispls, 0, br_hcomm, ierr)
   endif

   ! Now set data out to all members of the fragment
   call BMPI_BCAST(vec1, size(vec1) , 0, fcomm1, ierr)

   ! if(iproc == 0) print *, "KSM:  in br_restore_vec1, after BCAST"

end subroutine br_restore_vec1


!
! Remove overlap with the top of the history stack
! This holds  v_i when we are processing v_{i+1}
! we want v_{i+1} = v_{i+1} -  v_{i+1} . v_i
! v_{i+1} is still not normalized
!
!  CALLED BY: lanczos_p
!
!  CALLS: dot_product8
!  
subroutine br_remove_prev_overlap(a)
   use bmpi_mod
   implicit none
   integer :: ierr
   real(kind=8) :: a
   real(kind=8) :: dp
   integer(kind=basis_prec) :: pos

   call br_check_setup()
   dp = dot_product8(br_histbuf(ostart:ostop,br_histpos), br_reg)

   call BMPI_ALLREDUCE(dp, a, 1, MPI_SUM, icomm, ierr)
!   br_reg = br_reg - a * br_histbuf(ostart:ostop, br_histpos)
! tried using pointer into br_histbuf,   Intel 12.0 fails.   
!$omp parallel do &
!$omp    private(pos) &
!$omp    firstprivate(a, ostart, ostop) &
!$omp    shared(br_reg) 
   do pos = ostart, ostop
      br_reg(pos) = br_reg(pos) - real(a * br_histbuf(pos,br_histpos), kind=lanc_prec)
   end do
!$omp end parallel do
   return
end subroutine br_remove_prev_overlap

!
! Normalize br_reg and return scale factor
!  CALLED BY: lanczos_p
!  CALLS: dot_product8
!
subroutine br_normalize(a)
   use bmpi_mod
   implicit none
   integer :: ierr
   real(kind=8) :: a
   real(kind=8) :: dp, tdp
   integer(kind=basis_prec) :: i
   real(kind=8) :: m

   call br_check_setup()
   dp = dot_product8(br_reg, br_reg)
   call BMPI_ALLREDUCE(dp, tdp, 1, MPI_SUM, icomm, ierr)
   a = sqrt(tdp)
   m = 1.0d0 / a
   ! KSM TODO: write loop for OMP
   ! br_reg = br_reg * m
!$omp parallel do  private(i) firstprivate(m) shared(br_reg,ostart,ostop)
   do i = ostart, ostop
      br_reg(i) = real(REAL(br_reg(i), kind=8) * m, kind=lanc_prec)
   end do
!$omp end parallel do
end subroutine br_normalize

!
! Remove overlap with the history buffer
! off supports ending lopo at br_histpos-1 when we have
! already removed overlap with the previous vector to get \alpha
!
!  CALLED BY: lanczos_p
!
!  CALLS:
!      dot_product8
subroutine br_orthogonalize(off)
   use bmpi_mod
   implicit none
   integer :: off
   integer :: ierr
   integer :: hi, maxhi
   integer(kind=basis_prec) :: pos
   real(kind=8) :: dp, tdp, sum, m8
   real(kind=lanc_prec) :: m

   call br_check_setup()
   br_orthogdotf = 0.d0
   br_orthogdot = 0.d0
!   print*,' testing reorthogo ',br_histpos, ostart,ostop
   ! Compute local dot products
   maxhi = br_histpos - off
!   print*,' orthog before ',br_orthogdot(1:maxhi)
   
   if(maxhi > 2*ompNumThreads) then
!$omp parallel do &
!$omp    private(pos, hi, m8, sum) &
!$omp    firstprivate(ostart, ostop, maxhi) &
!$omp    shared(br_orthogdot, br_histbuf)
      do hi = 1, maxhi
         sum = 0.0
         do pos = ostart, ostop
            sum = sum + REAL(br_reg(pos), kind=8) &
                * REAL(br_histbuf(pos,hi), kind=8)
                  end do
         br_orthogdot(hi) = sum
      end do
!$omp end parallel do
   else
      ! old version, OpenMPI inside dot_product
      do hi = 1, maxhi
         br_orthogdot(hi) = dot_product8(br_histbuf(ostart:ostop, hi), br_reg)
      end do
   end if

   ! do all the ALLREDUCES at same time
!   print*,' orthog ',br_orthogdot(1:maxhi)

   call BMPI_ALLREDUCE(br_orthogdot, br_orthogdotf, maxhi, MPI_SUM, icomm, ierr)
!   print*,' orthogf ',br_orthogdotf(1:maxhi)

   ! Now sum overlaps before subtracting - better numerics to add a bunch of small overlaps before
   ! subtracting from large number
!  ORIG: for hi:   br_reg2 = br_reg2 + m * br_histbuf(ostart:ostop, hi)
!$omp parallel do private(pos, hi, sum)
   do pos = ostart, ostop ! pos on outside for locality
      sum = 0.d0
      do hi = 1, maxhi
         sum = sum + br_orthogdotf(hi) * br_histbuf(pos, hi)
      end do
      ! subtract totaled overlap
      br_reg(pos) = br_reg(pos) - REAL(sum, kind=lanc_prec)
   end do
!$omp end parallel do
end subroutine br_orthogonalize
!--------------------------------------- + -------------------------------------

!
! pushes br_reg onto history buffer
!
subroutine br_add2hist(iter)
   implicit none
   integer :: iter

   call br_check_setup()
   br_histpos = br_histpos + 1 ! make room for new data
!    thick restart breaks this rule.   iter hasn't been updated yet
!   if(iter /= br_histpos) then
!      if(iproc == 0) print *, "iter and br_histpos do not match!"
!      stop 23
!   end if
   if(br_histpos > UBOUND(br_histbuf, 2)) then
      if(iproc==0) print *, "Trying to save too many Lanczos vectors.  Buffer depth=", UBOUND(br_histbuf, 2)
      stop 23
   end if
   br_histbuf(ostart:ostop, br_histpos) = br_reg
end subroutine br_add2hist


subroutine br_load_vec1_from_vec2()
   use bmpi_mod
   implicit none
   integer :: ierr

   if(nfragments == 1) then
      vec1 = vec2
   else
      ! we can transfer from vec2 to vec1 on the diagonal root nodes and then broadcast out
      ! within each fragment using fcomm1
      if(isfragroot) vec1 = vec2
      call BMPI_BCAST(vec1, size(vec1) , 0, fcomm1, ierr)
   end if
end subroutine br_load_vec1_from_vec2

subroutine br_load_vec2_from_vec1()
   use bmpi_mod
   implicit none
   integer :: ierr

   if(nfragments == 1) then
      vec2 = vec1
   else
      ! we can transfer from vec2 to vec1 on the diagonal root nodes and then broadcast out
      ! within each fragment using fcomm1
      if(isfragroot) vec2 = vec1
      call BMPI_BCAST(vec2, size(vec2) , 0, fcomm2, ierr)
   end if
end subroutine br_load_vec2_from_vec1

subroutine br_open_u()
   implicit none

   dbg_u = 31
   open(UNIT=dbg_u, FILE="history.txt")
end subroutine br_open_u

subroutine br_close_u()
   close(dbg_u)
end subroutine br_close_u

! Testing routine
subroutine br_write_buf(u, buf)
   implicit none
   integer :: u  ! unit to write to
   real (kind = lanc_prec) :: buf(:)
   integer :: i, sz

   if(iproc == 0) then
      print *, "In br_write_buf"
      sz = size(buf)
      write(u, FMT="(A)") "VECTOR"
      do i = 1, sz
         if(i == sz) then
            write(u, FMT="(f9.6)") buf(i)
         else
            write(u, FMT="(f9.6,A)", advance='no') buf(i), ", "
            if( mod(i, 5) == 0 ) then
               write(u, fmt="(A)") ""
            end if
         end if
      end do
      flush(u)
   end if
end subroutine br_write_buf

! Testing routine
subroutine br_write(u)
   use bmpi_mod
   implicit none
   integer :: u  ! unit to write to
   integer :: ierr
   integer :: pi
   integer ::  sendcounts(0:nproc-1), displs(0:nproc-1)
   real(kind=lanc_prec) :: buf(dimbasis)
   integer :: count

   do pi = 0, nproc-1
      sendcounts(pi) = int(br_ostoplist(pi) - br_ostartlist(pi) + 1,4)
      displs(pi) = int(br_ostartlist(pi) - 1, kind(displs(0)))
   end do

   count = int(ostop - ostart + 1, 4) ! fragments are limited to 2B size
   if(nproc == 1) then
      buf = br_reg
   else
      call BMPI_GATHERV(br_reg, count, buf, sendcounts, displs, 0, icomm, ierr)
   end if

   ! We have now gathered the entire vector to iproc==0 so we can do the I/O
   call br_write_buf(u, buf)
end subroutine br_write

subroutine br_reset_reg()
   implicit none
   br_reg = 0
end subroutine br_reset_reg

subroutine br_reset_histbuf()
   implicit none
   br_histbuf = 0.0
end subroutine br_reset_histbuf

!
! We have a partial basis for the complete space stored in br_histbuf.  It has
! br_histpos vectors in it.   
!
! eiglvec is the set of eigenvectors of the projection of H on this space (the \alpha, \beta
! tri-diagonal matrix).   They have been sorted by energy.
! nkeep is the number of vectors we want to report at the end, so we only transform that
! many.

! We can now transform the history buffer to produce approximate
! eigenvectors in the complete space by taking a linear combination of histbuf(i, :) vectors
! weighted by the components of an eigenvector on the reduced space.
!
!  INPUT:
!   eiglvec : array of eigenvectors of truncated Lanczos tridiagonal
!   nkeep : # of eigenvectors to keep/transform
!   eign  : dimension of eiglvec
!
subroutine br_transform_basis(eiglvec, nkeep, eign)
   implicit none
   integer, intent(in) :: nkeep
   integer, intent(in) :: eign
   ! first index in eiglvec is vector member, second identifies the vector
   real(kind=egv_prec) :: eiglvec(:,:)
   real(kind=egv_prec) :: dtmp
   real(kind=lanc_prec) :: tmpvec(nkeep)
   integer :: pos, j
   integer(kind=basis_prec) :: i

   call br_check_setup()
   if(nkeep > br_histpos .or. eign > br_histpos) then
      if(iproc == 0) print *, "br_transform_basis: nkeep > br_histpos",nkeep,eign,br_histpos
      stop 23
   end if
   ! process bit by bit to reduce intermediate storage
   tmpvec = 0.0
   do i = ostart, ostop
      do pos = 1, nkeep  ! for each output vector we want
         dtmp = 0.d0
         do j = 1, eign
            dtmp = dtmp + real(br_histbuf(i, j),kind=8) * real(eiglvec(j, pos), kind=8)
         end do
         tmpvec(pos) = real(dtmp, kind=lanc_prec) ! dtmp is higher precision so sum is accurate
      end do
      ! now replace the column in the history buffer
      do pos = 1, nkeep
         br_histbuf(i, pos) = tmpvec(pos)
      end do
   end do
end subroutine br_transform_basis

subroutine dbg_br_reg_rng(s, e)
   implicit none
   integer(kind=basis_prec) :: i
   integer :: s, e

   do i = s, e
      print *, "br_reg(", i, ") = ", br_reg(i)
   end do
end subroutine dbg_br_reg_rng

subroutine dbg_br_reg()
   implicit none
   call dbg_br_reg_rng(1, 10)
end subroutine dbg_br_reg

!......... ADDED 7.6.8/7.7.1s by CWJ........
!    fetches the coefficients for a specific basis label , ibasis

subroutine br_fetch_coef(pos,nkeep,vamp)
	use precisions
	use localvectors
	use nodeinfo
	use bmpi_mod
	implicit none
    integer(kind=basis_prec) :: pos
	integer :: nkeep
	real(kind = lanc_prec) :: vamp(nkeep),tmpvamp(nkeep)
	integer :: i,ierr
	integer :: tag
	integer :: stat(MPI_STATUS_SIZE)
	integer :: fromproc
	
!	call br_check_setup()
	vamp(:)=0.0
	
	if(nproc==1)then
		do i = 1,nkeep
 		   vamp(i) = br_histbuf(pos, i)
		   tmpvamp(i)=vamp(i)
		end do
	else
		
!...... FIND OUT WHERE THE COEFFICIENTS ORIGINATE FROM.....		
		fromproc = -1
		do i = 0,nproc -1
		    if(pos >= br_ostartlist(i) .and. pos <= br_ostoplist(i))then
				if(fromproc >= 0)then
					print*,' AAGG some problem with fromproc '
					print*,i,fromproc
					stop
				end if
				fromproc = i
			end if
			
		end do
		
	    if(pos >= ostart .and. pos <= ostop .and. fromproc > 0)then
			if(iproc/=fromproc)print*,' MISMATCH IN fromproc ',fromproc,iproc

			   tag = iproc
			   call BMPI_SEND(tmpvamp,nkeep,0,tag,icomm,ierr)
 	    end if
		if(iproc==0)then
			if(fromproc /= 0)then
 			   tag = fromproc
 			   call BMPI_RECV(vamp,nkeep,fromproc,tag,icomm,stat,ierr)
	         end if
		end if
    end if


!	call BMPI_REDUCE(vamp,nkeep,MPI_SUM,0,icomm,ierr)
	return
	
end subroutine br_fetch_coef

! ---------------------------------------------
!
!  ADDED in 7.7.6: add a random vector for restart
!
subroutine add_random_br_reg(dimsize)
	implicit none
	integer (kind=basis_prec):: dimsize
	real :: rv
	real :: mynorm
    integer(kind=basis_prec) :: i
	
	mynorm = 1.0/sqrt(float(dimsize))  ! rough normalization to simplify
	
	do i = ostart,ostop
		call random_number(rv)
		br_reg(i)= (rv-0.5)*mynorm
	end do
	return

end subroutine add_random_br_reg
! ---------------------------------------------
!
! ADDED in 7.7.6 -- creates a random vector in br_reg
! orthogonal to all previous vectors in the history

subroutine random_restart_with_history
	use basis, only:dimbasis
	implicit none
	real(8) :: dnorm
	
	call add_random_br_reg(dimbasis)
	call br_normalize(dnorm)
	call br_orthogonalize(0)
	call br_normalize(dnorm)	
	
	return
end subroutine random_restart_with_history	

end module mod_reorthog
