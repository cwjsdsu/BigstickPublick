!  MPI master notes:
!    module nodeinfo :
!     nproc = # of MPI proceses
!     iproc = rank of each MPI process (from 0, which is master root node, to nproc-1)
!    (note: for simulations, nprocs = simulated # of MPI process, and iproc_dummy simulated rank)
!     icomm = communicator group for all processes
!
!    module localvectors
!      on each MPI process:
!      frag1, frag2 = assigned fragments for vec1, vec2; these are set in setup_localvectors
!      v1s, v1e are start, end for basis indices on vec1, so vec1(v1s:v1e)
!           note: same as basestart(frag1), basestop(frag1)
!      v2s, v2e are start, end for basis indices on vec2
! 
!      fcomm1, fcomm1 = communicators for reducing vec1 and vec2; set in BMPI_COMM_SPLIT (called by setup_localvectors)
!      isfragroot: logical flag to denote a "root" process assigned a fragment; set in setup_mpi_hist_comm
!
!=============================================================
!
!  for "new" parallelization, store vectors (or fragments of vectors) 
!  in this module; this reduces the # of information for subroutine calls
!
!  by making the vectors targets, we can use pointers to switch between them
!
!  when operating in parallel, can swap between vec1 and vec2
!  vecstore is used when storing the lanczos vectors all in RAM
!
!  With fragments, we only compute one direction
!
module localvectors
   use precisions
   implicit none
   ! useVec2Thread gives each thread an independent buffer to write to.
   logical :: wantUseVec2Thread = .false.   ! combine with numthreads > 1 to set useVec2Thread
   logical :: useHZSomp = .true.    ! flag for new OMP merged in 7.7.3 from HZS 
   double precision :: mpiStartTime
   
   logical :: useVec2Thread
   integer :: ompNumThreads   ! top level number of OpenMP threads

   real(kind=lanc_prec), allocatable, target :: vec1(:), vec2(:)
   ! KSM  Add thread local vec2 clones for performance testing
   !! real(kind=lanc_prec), allocatable, target :: vec2thread(:, :) ! (v2s:v2e, 0:num_threads-1)
   ! make new vec2thread flat to avoid bug in old fortran compilers
   integer(kind=basis_prec) :: vec2threadchunk, vec2threadchunkm1, vec2threadend
   real(kind=lanc_prec), allocatable, target :: vec2threadflat(:) ! (0 : num_threads * (v2e - v2s + 1))
   ! KSM  {v1s,v1e,v2s,v2e}={basestart(frag1), basestop(frag1), basestart(frag2), basestop(frag2)}
   integer(kind=basis_prec) :: v1s, v1e, v2s, v2e   ! slice start and end for vec1 and vec2
   integer, target :: frag1,frag2  ! which fragments each vector is
   ! note:   frag1 = nodal(iproc)%ifragment,  see setup_localvectors

   ! KSM - communicators for broadcast and reduction across
   ! one fragment.  We need one for each direction of vchar
   integer :: fcomm1, fcomm2  ! communicators to use for bcast/reduce on frag1/frag2
   logical :: isfragroot  ! node to scatterv/gatherv to slices, bcast back to fragment matvec nodes

!........... FOR DISTRIBUTED STORAGE AND REORTHOGONALIZATION....
!            the following arrays are used for controlling
!            distributed storage of Lanczos vectors
!            and their reorthgonalization
!
!  We do not store an entire Lanczos vector on a node
!  instead we store a "piece" (a fraction of a vector) 
!  and on a given node we store that same gene for all 
!  Lanczos vectors 
! 
   logical :: storelanczosincoreMPI     ! store lanczos vectors across MPI cores
   logical :: storelanczosincore1       ! store lanczos vectors on a sing
   integer   :: npiece   ! # of pieces of the vector
   integer(8):: Lpiece  ! size of a "piece" of the vector
   integer(8) :: Ldim, Lstart
   integer(8), allocatable :: Ldim_r(:), Lstart_r(:)
   integer(kind=8) :: buff_max
  
end module localvectors
!=============================================================

!
module nodeinfo
  implicit none
  integer(4) :: iproc,nproc,icomm

  integer(4) :: nprocs    ! used as a dummy
  integer(4) :: iproc_dummy   ! used for simulations

  logical :: distributeMPI =  .true.    
  logical :: simulateMPIdefault   =  .false.
  logical :: simulateMPI  

!...........................................................................
end module nodeinfo


module distrinfo2

  integer(8) :: totpnops, totppops, totnnops,totpppops,totppnops,totpnnops,totnnnops
  real(8)    :: time_pnop, time_ppop, time_nnop    ! approximate times for each operation
  real(8)    :: time_pppop, time_ppnop, time_pnnop, time_nnnop

end module distrinfo2

!-------------------------------------------------------------
!  module for information on distributing work for OpenMP

module distrOpenMP
 integer nthreads
 integer nsectorjumps
 logical :: openmp_ok = .false. ! .true.   
 real :: tottime 

 type opensectorjumps
    integer :: isector
    integer :: fsector
    character(2) :: jumptype   ! = pn, pp, nn
    integer :: sectorjump, csectorjump
    character(1) :: dir      ! = f or b
    integer :: nops
    real :: time
    integer :: thread
 end type opensectorjumps
 
 type (opensectorjumps), allocatable :: osj(:)

 type protonsectorlist
     real :: time
     integer :: nsectorjumps
     integer, allocatable :: jumplist(:)
 end type protonsectorlist

 type (protonsectorlist), allocatable :: pslist(:)

 type opensectorjumplist
    integer :: nsectorjumps
    logical, allocatable :: psector(:)  ! organize by final proton sector
    real :: time
    integer, allocatable :: jumplist(:)
 end type opensectorjumplist

 type (opensectorjumplist), allocatable :: osjlist(:)

end module distrOpenMP

!=========================================================
!   OPERATION_STATS
!      data on operations between sectors
!
!started 8/2011 by CWJ
!

module operation_stats
  implicit none

!.... DEFINE FUNDAMENTAL TYPE FOR # OF OPERATIONS
  type op_mat
     integer(8) :: nopPP, nopNN, nopPN, nopPPP, nopPPN, nopPNN, nopNNN,nopSPE
     integer(8) :: nopP, nopN   ! for densities
     real(8) :: optot
     integer :: nnodes
  end type op_mat

  type (op_mat), allocatable :: opstat(:,:)      ! count operations between sectors
  type (op_mat), allocatable :: opfragstat(:,:)   ! operations between fragments

!....... advanced: weights for operations
!        when they do not all take the same time

! add by Shan in 7.7.3ZZ
  real(4) :: opwtPPb, opwtPNb
  integer(8) :: allopPPb, allopPNb

  real(4) :: opwtPP, opwtNN, opwtPN, opwtPPP, opwtNNN, opwtPNN, opwtPPN,opwtSPE

  integer(8) :: allops, allopPP, allopPN, allopNN, allopPPP, allopPPN, allopPNN, allopNNN
  real(8) :: totalops
  real(8), allocatable :: frag2fragops(:,:)   ! alternate data to opfragstat
  logical :: usewts = .true.

!------- INFORMATION ON # OPERATIONS ON CURRENT MPI PROCESSES
!        previously in module tribution, moved in 7.6.7
  integer(8) :: nopsnode   ! # of operations per node

end module operation_stats
!
! MODULE fragments MOVED to file bfragments.f90
!
!=========================================================
!
!   NOTES: If we have up to Nlanc lanczos vectors, we must store each fragment Nlanc times.
!   This means that Nlanc * nfragments <= # of nodes.  
!   The ideal case would be that on each of those Nlanc nodes we are using that fragment
!   for multiplication.  Not clear if this will happen
!
!   some cases: suppose we have 30,000 nodes and basis dim = 500M. 
!   If we allow 50 lanczos vectors, then Nfragment < 30,000/50 = 600, and each fragment must 
!   be around 1 M, which is doable.
!   If we allow 300 lanczos vectors, then Nfragment < 30,000/300 = 100, and each fragment must be 
!   around 5 M, also doable. 
!
!   Now suppose we have 100,000 nodes and basis dim = 10 G.  
!   Nlanc = 50 means Nfrag < 100,000/50 = 2000, and each fragment roughly 10 G/ 2000 = 5 M, also doable.
!
!   The trick, however, is trying to get operations distributed so that this is reasonable. 
!   We want roughly the same number of operations to each fragment, which seems unlikely.  
!   What needs to be computed is the mismatch.
!
!   Alternately, for small systems where we can't store all the lanczos vectors in RAM, we can use
!   MPI-IO to write fragments to disk. 
!=========================================================
!
! information on distribution
! started 8/2011 by CWJ @ SDSU
!
!  DELETE IN 7.6.7 June 2016
!
!module tribution
!   implicit none
!  --- MOVED TO MODULE OPERATION_STATS
!   integer(8) :: nopsnode   ! # of operations per node
!----- MOVED TO MODULE FRAGMENTS:
!   type whatsonanode
!      integer :: ifragment,ffragment   ! which fragments are here
!      integer :: sfragment  ! which fragment is stored here; later might want to store more than one
!      real :: frac_ops   ! what fraction of possible ops are stored here      
!      logical :: ifirst ! first node in group servicing ifragment
!   end type whatsonanode
!   type (whatsonanode), allocatable :: nodal(:)

!end module tribution
!=========================================================
!
!  opbundles are new (9/2011) storage of information about operations  
!

module opbundles

  implicit none
  integer nopbundles
  
  integer, parameter :: MAXTHREADS =64   ! ADDED in 7.7.3 for new OpenMP

  type bund
     integer :: node
     integer :: isector, fsector    ! initial/final proton sectors
     integer :: insector,fnsector   ! neutron initial and final sectors
     integer :: ifragment, ffragment   
     
     character(3) :: optype
     character(1) :: hchar
     integer(kind=8) :: pxstart,pxend   !start, stop of either:
                                      ! proton jumps
                                      !  or proton SDs
     integer(kind=8) :: nxstart,nxend   ! start, stop of either:
                                      ! neutron jumps
                                      !  or neutron SDs
     integer(kind=8) :: nsortstart,nsortend  ! start, stop for sorting neutron jumps; 
                                             ! this solves a problem for the MPI distribution
     integer(kind=8) :: cstride               ! incrementing pSDs
     integer(kind=8), allocatable ::  startp_thread(:)
     integer(kind=8), allocatable ::  startn_thread(:)
     integer :: numpthreads,numnthreads
     integer(8) :: nops      ! # of operations in this bundle
     integer(8) :: min_nop  ! min # of operations; this is governed, for example, by 
                         ! # of conjugate SDs
     integer(8) :: nsetops   ! how many sets or block of operations there are
			! NOTE: min_nop x netsops = nops
     integer(8) :: njumps
!	 logical    :: annexed    ! added 7.8.0 by CWJ -- if annexed, then the work is added to a nearby bundle
!     integer(kind=8) :: threadStart(0:MAXTHREADS), threadStop(0:MAXTHREADS)   ! added in 7.7.3 by HZS
	 
  end type bund

  type (bund), allocatable,target :: opbundle(:)  ! storage of information for operations
  type (bund), allocatable,target :: draft_opbundle(:)  
             ! first round storage of information for operations

  integer, allocatable :: nopbundlesnode(:)  ! # of bundles on a given node
  integer, allocatable :: opbundlestart(:), opbundleend(:)  ! first, last bundles on a node
  integer, allocatable :: frag_draftopbundlestart(:,:), frag_draftopbundleend(:,:)  ! for given initial, final frags
                                                            ! which are first, last draft opbundles

end module opbundles
!=========================================================
!
!  information on limits on jumps
!  so only create jumps needed on each node
!    -- initiated 9/2012 by CWJ @ SDSU/TRIUMF
!
! FOR FLAGS GO TO MODULE FLAGGER IN bmodule_flags.f90
!
module jumplimits
  implicit none
  
!...... THESE ARRAYS ALLOCATED AND FILLED in routine setMPIjumpslimits in file bparallellib3.f90
!
  logical, allocatable :: makep1bjumps(:), maken1bjumps(:),makeppjumps(:), makennjumps(:)
  logical, allocatable :: makePPPjumps(:), makeNNNjumps(:)

  integer(8), allocatable :: startPPjumps(:), stopPPjumps(:), startp1bjumps(:), stopp1bjumps(:)
  integer(8), allocatable :: startNNjumps(:), stopNNjumps(:), startPPPjumps(:), stopPPPjumps(:)
  integer(8), allocatable ::  startNNNjumps(:), stopNNNjumps(:), startn1bjumps(:), stopn1bjumps(:)
  real :: maxjumpmemory  ! in Gb

!...... ADDED 7.4.2.....  REDIRECTION BY ELIMINATING INTRONS.....
!.......... new storage by MPI process, needed for modeling/instrumentation
  integer(8), allocatable,target :: PPjumplength(:), NNjumplength(:), P1Bjumplength(:),N1Bjumplength(:)
  integer(8), allocatable,target :: PPPjumplength(:), NNNjumplength(:)
!......... redirection arrays...............different for each MPI process...
!
!  for jump type X, if native (original) jump label i,  for some index n Xjumpcut(n)<= i < Xjumpcut(n+1)
!  new jump label iprime = i - Xjumpshift(n) 
!  new indices go from 1 to Xjumplength(iproc) (on MPI process iproc)
!  use Xendcut(:) to determine excluded jumps, i.e., for some (original) jumplabel i, 
!  if Xjumpcut(n) < Xendcut(n) < i < Xjumpcut(n+1), set iprime < 0 (which triggers skipping)
!
  integer(8), allocatable,target :: P1Bjumpcut(:),N1Bjumpcut(:)
  integer(8), allocatable,target :: PPjumpcut(:),NNjumpcut(:)
  integer(8), allocatable,target :: PPPjumpcut(:),NNNjumpcut(:)
  integer(8), allocatable,target  :: P1Bjumpshift(:),N1Bjumpshift(:),PPjumpshift(:),NNjumpshift(:)
  integer(8), allocatable,target :: PPPjumpshift(:),NNNjumpshift(:)
  integer(8), allocatable,target :: P1bendcut(:),N1Bendcut(:),PPendcut(:),NNendcut(:),PPPendcut(:),NNNendcut(:)
  
  logical :: compactjumpstorage
  logical :: compactjumpstorage_enabled = .true.
end module jumplimits
!=========================================================
module timing_parallel

implicit none

real(8)              :: timelast
real(8), allocatable :: time_Ham_MPI(:)

real(8)              :: timelastbundle,timelastprocop
real(8), allocatable :: time_bundle(:)
real(8), allocatable ::  time_procSPE(:)
real(8), allocatable ::  time_procPP(:)
real(8), allocatable ::  time_procNN(:)
real(8), allocatable ::  time_procPN(:)
real(8), allocatable ::  time_procPPP(:)
real(8), allocatable ::  time_procPPN(:)
real(8), allocatable ::  time_procPNN(:)
real(8), allocatable ::  time_procNNN(:)

real(8), allocatable ::  time_procPPb(:)
real(8), allocatable ::  time_procPNb(:)


end module timing_parallel

