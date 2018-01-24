!---------------------------------
!
!  this module contains most of the flags and parameters
!  to be set at compile time for running bigstick
!
!-------------------------------

module flagger
  implicit none

! Tracing options - usefull with crashes on edison/hopper
  logical :: wantnoisy0 = .false.
  logical :: noisy0      ! iproc==0 .and. wantnoisy0.   See main
!............ OUTPUT OPTIONS.....................  added in 7.3.6

  logical :: spoccflag     ! flag to compute and write out single-particle occupations
 
!........... STORAGE OF PN and PPN/PNN UNCOUPLED MATRIX ELEMENTS.....

  logical,parameter :: XYme_altstore_enabled = .true.
  logical           :: XYme_altstore
 
!...........HF PRECONDITIONING.......

  logical :: HF    = .false.       ! IF TRUE THEN COMPUTE HF POTENTIAL; mostly obsolete

!................. LANCZOS .........................

  logical :: checkeig  =.true.  ! if FALSE do not check intermediate energies
  integer :: itskipeig = 10    ! how often to print out intermediate energies during Lanczos  
  logical :: skipeig = .true.  ! if TRUE then only check eigenvalues every ITSKIPEIG times
  logical :: forcewritelanczos = .false.  ! writes Lanczos vectors to file even if storing lanczos vectors in core
                                         ! needed if one wants to restart
  integer :: fulllimit = 700  ! if dimension is less than this, suggest full diagonalization
  integer :: ncheck_ex_def = 5  ! default # of extra states to check for convergence
  integer :: ncheck_ex        ! # of extra states to check for convergence

!........... ONLY THICK-RESTART IF ITER > MAX( nkeep+thick_restart_min, thick_restart_factor*nkeep)
  integer :: thick_restart_factor_def = 2, thick_restart_min_def = 10
  integer :: thick_restart_factor , thick_restart_min
  integer :: nthick_add_def = 5   ! default # of extra states to add for thick-restart; must be >= ncheck_ex
  integer :: nthick_add   ! # of extra states to add for thick-restart; must be >= ncheck_ex
  
  real    :: restart_tol = 1.0e-1   ! tolerance of beta to force a random restart

!............. REORTHOGONALIZATION IN LANCZOS..................

  logical :: alpha_before_orthog = .true.    ! if T, then compute alpha before orthogonalization
                                             ! else compute alpha at end of orthgonalization
                                             ! In exact arithmetic this makes no difference
                                             ! In finite arithmetic can change alpha, beta coefficients
                                             ! Converged energies agree to better than 1 keV;
                                             ! converged wfns agree in overlap to 10^4 or better 

  logical :: orthog_sequential = .false.      ! only for MPI implementation when storing lanczos in RAM
                                             ! if T, then orthgonalize sequentially, as done when storing vectors on disk
                                             ! otherwise subtract simultaneously; this can lead to differences 
                                             ! in alpha, beta coefficients as above but same converged energies etc

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

!............ RESTRICTING STORAGE OF JUMPS.......
!  The "jumps" store the Hamiltonian in compact form
!  When parallelized via MPI, it is possible to store on each node
!  only the jumps needed  
! SEE ALSO MODULE jumplimits in bmodules_parallel.f90
!
!

   logical :: dosortjumps      = .true.   ! flag to turn off sorting of jumps for OpenMP
   logical :: opbundlesplitMPI = .true.     ! if true then allow splitting opbundles for MPI distribution
                                            ! otherwise, do not split
   logical :: restrictjumps = .true.         ! if .true. then only store relevant jumps on a node
   logical :: setjumpceiling = .true.         ! memory ceiling on storing jumps
!   real    :: maxjumpmemory_default  = 32.5    ! in Gb
   real    :: maxjumpmemory_default  = 16.0    ! in Gb  KSM - for DM runs on Edison
   
!............ STORAGE OF LANCZOS VECTORS IN CORE...........
!              FIRST FOR A SINGLE PROCESSOR
   logical :: enablelanczosincore1 = .true.
   real    :: maxlanczosstorage1 = 16.000    ! in Gb

!..............................................................
!              MULTIPLE MPI PROCESSORS
   logical,parameter :: enablelanczosincoreMPI = .true.  ! allow storing lanczos vectors in core if more than one processor
                                           ! otherwise they are written to disk
   logical,parameter :: asklancmemstore = .false. ! .true.
! following controls max storage of lanczos vectors on one MPI process
   integer(8),parameter :: maxpiecedimdefault = 4000000000_8   ! uses 1 Gb under normal defaults
!   integer(8) :: maxpiecedimdefault = 250000000   ! uses 1 Gb under normal defaults
   integer :: minpiece = 1000     ! suggested minimum size for pieces; if too small then inefficient

!.............. BREAKING LANCZOS VECTORS FOR HAMILTONIAN MAT-VEC MULTIPLY............

   logical,parameter :: break_vectors_enabled = .true.
   logical,parameter :: ask2break_vectors = .true.
   logical :: break_vectors  
   integer(8) :: maxfragmentsize_default = 200000000   
   integer(8) :: maxfragmentsize

!..............................OTHER.........................

   logical :: subsume_SPE = .true.   ! if TRUE, then subsume SPE into two-body matrix elements ADDED 7.8.0
   logical :: annex_bundles = .false. ! if TRUE, then seek to annex bordering bundles which share jumps ADDED 7.8.0
   
   logical  :: iffulluseph  = .true.    ! if TRUE, when a shell is filled, automatically turns to p-h conjugation

   logical :: pnsymorb =.false.   ! if TRUE, when reading in .pn.int file, assume proton and neutron orbits are the same
                                 !       in which case assume V_pn(ab,cd)=V_pn(ba,dc )* phase
								 !  only really used if you copy an isospin conserving force as a .pn.int file
								 ! if FALSE, must enter V_pn(ab,cd) and V_pn(ba,dc) separatedly  PREFERRED
   logical :: storeXXmesjumps= .false.  ! store XX/XXX matrix elements directly in jumps, rather than by index
                                    ! specifically in arrays p2b_me(), n2b_me(), p3b_me(), n3b_me()
   logical :: sumdiagonalXXmesjumps = .false.  ! if true, compute and store summed diagonal XX/XXX matrix elements
                                  ! this is important when many particles to loop over, > 10 
								  !  added in 7.6.0 ; may go away
   logical,parameter :: disableNoOMPloopswitch  = .true.  ! if T, disables switching order of loops when there is no OpenMP
   logical,parameter :: allow_alt_PNb_loops = .true. ! ADDED 7.8.0  if true, then can run PNb like Pnf
   integer,parameter :: pstridecut = 1000     ! (typically ~1000) in applying NN, if there is no OpenMP, 
                                    ! the order of the loops is switched. However,
                                    ! if there are a lot of neutrons, then the locality of proton SDs
                                    ! is degraded and one should not switch the order of the loops
                                    ! THIS IS HIGHLY COMPILER AND FLAG-DEPENDENT and determined empirically

end module flagger
