!==========================================================================2
! MODULES FOR BIGSTICK
!===========================================================================
!
!  MODULES:
!       ------- BASICS AND SINGLE-PARTICLE SPACE ----

!    system_parameters: many-body space defined by Z,N, Jz, parity...
!    verbosity : flags for writing out information, file numbers, etc
!    bitstuff : # of bits per word assumed
!    sporbit:  single particle ORBITs & quantum numbers 
!        key derived variable: orbqn
!    spstate:  single particle STATEs & quantum numbers
!        key derived variable: spsqn
!        [ NB: an orbit labeled by (nlj) while state by (nljm) ]
!      -------  HAIKUS ----------
!    W_info   : truncations on weighting W
!    spsections : partition of the s.p. space by W; foundation of templates
!    templates : partition of the many-body space by W; like configurations, 
!                  only not as fine-grained
!    builder : derived types used in constructing templates
!    haiku_info: information on haikus, mostly limits on W
!    haikus  : half-slater determinants used to construct...full slater determinants
!    blocks : organization of haikus
!     ---------- BASIS -------------
!    sdlimits: q# limits on SDs, used in constructing basis
!    basis : basic information on basis, including dim and start arrays
!    sectors: organization of basis
!    ----------- HAIKU HOPS ---------------
!    bugs   : derived types for hops
!        --- CONSTRUCTING JUMPS --------
!    descendents:
!    geneologies: 
!    chains
!       ----- JUMPS ------------
!    jumpdef
!    jumpNbody
!      ------- INTERACTION -----------------
!    pairdef
!    interaction
!    diagh
!      ------- APPLYING THE HAMILTONIAN
!    lanczos_info
!      ------- OUTPUT ---------------------
!===========================================================================
!
!------  Modules for MPI node information ----------------------------------


!-------------- OVERALL QUANTUM NUMBERS OF THE SYSTEMS----------------------
module system_parameters
  implicit none
  integer :: np(2) ! Number of particles
  integer :: jz
  integer :: iparity
  character(len=1) :: cparity
  integer :: Acore    ! a of the core
  
!...... INFORMATION FOR PARTICLE-HOLE CONJUGATION........
  logical :: phconj(2)
  integer :: npmax(2)
  integer :: npeff(2)  ! effective number of particles regarding particle-hole conjugation;
                       ! added in 7.5.6

end module system_parameters

!...........................................................................
module menu_choices
  implicit none
  character(len=3)  :: menu_char
  logical :: menu_dx_omathematica ! density in mathematica format
end module menu_choices

!------------------------------------------
! MODULE PROGRAM_INFO
!
! contains information about code
!

  module program_info
     implicit none

     character(6) :: version = '7.11.4'
     character(9) :: lastmodified = 'Aug 2025'

  end module program_info


!-------------------------------------------------------------------
!  MODULE VERBOSITY
!-------------------------------------------------------------------
!  Various flags for printing out information.
!  This allows one to easily print out or suppress information.
!
!  NOTE: verbose_X  means write out information to terminal
!        chatty_X  means write out information as being computed (rare)
!        print_X    means write to file
!        X_file     is unit number for file
!-----------------------------------------------------------------
module verbosity
  implicit none

  logical :: print4modelinfo       ! prints out information needed for modeling
  integer :: modelinfo_file    =  44
  logical :: verbose_winnow    = .false.
  logical :: print_groups      = .false.
  integer :: group_file        = 37
  logical :: print_hsps        = .false.
  integer :: hsps_file         = 37
  logical :: verbose_templates =  .false. 
  logical :: verbose_haikus    = .false.
  logical :: print_haikus      = .false.  
  integer :: haiku_file        =  37
  logical :: verbose_blocks    = .false.
  logical :: verbose_hops      = .false.
  logical :: print_hops        = .false.
  integer :: hop_file          =  38 
  logical :: verbose_sectors   = .false.
  logical :: print_sectors     = .false.
  integer :: sector_file       = 37  
  logical :: chatty_descent    = .false.
  logical :: chatty_links      = .false.

  logical :: print_genes       = .false.
  integer :: gene_file         = 38
  logical :: print_basis       = .false.
  integer :: basis_file        = 31

  integer :: xsd_file          = 32  ! added 4/2011 CWJ
  logical :: verbose_uncouple  = .true.
  
  logical :: verbose_jumps  = .true.

  logical :: print_matrix      = .false.
  integer :: matrix_file       = 55
  logical :: verbose_orthog    = .false.

  logical :: print_block_ops   = .false.
  logical :: print_node_load   = .false.

  logical :: print_sectorjumps = .false.    ! added 5/2010 CWJ
  integer :: sectorjump_file   = 39

  logical :: print_jumps = .false.   ! added 4/2011 CWJ
  integer :: jump_file = 46

  logical :: verbose_3body_readin = .false.
  logical :: verbose_3bodyjumps = .false.

  logical :: verbose_distrib = .false.  ! added 8/2011 CWJ
  logical :: verbose_bundles  = .true.  ! added 10/2011 CWJ

  logical :: verbose_fragments = .false.
  logical :: verbose_distro    =.false.
  logical :: verbose_distro_report = .false.
  

end module verbosity

!...........................................................................
module bitstuff
  implicit none
  integer max_bit_word         ! the max. # of usable bits in a
  parameter (max_bit_word=30)  ! given word to describe an SD
end module bitstuff

!---------------------------------------------------------------------------
!  declaring size of variables
!---------------------------------------------------------------------------
module precisions
  integer,parameter :: lanc_prec =    4 
  integer,parameter :: basis_prec   = 8
  integer,parameter :: obs_prec     = 8
  integer,parameter :: egv_prec     = 8
!  integer(4) :: MPI_lanc
end module precisions

!-----------------------------------------------------
!  INFORMATION ON TIME PER OPERATION
!  USEFUL FOR LOAD BALANCING
!  DEFAULTS CAN BE MODIFIED FOR A GIVEN MACHINE
!-------------------------------------------------------
module timeinfo

  implicit none
!.............. TIMING PER OPERATION (in ns)............
  integer :: timefile = 59

  logical :: modifytimewts = .true.   ! enables modification of timing
                                     ! by reading in timing from file "timebigstick.info"
                                     ! (if it exists; otherwise use defaults)
!   DEFAULT VALUES IN ABSENCE OF OTHER INFORMATION

  real(8),parameter :: tpoSPEdefault = 2.0
  real(8),parameter :: tpoPPdefault = 6.0
  real(8),parameter :: tpoNNdefault = 6.0
  real(8),parameter :: tpoPNdefault = 8.0
  real(8),parameter :: tpoPPPdefault = 6.0
  real(8),parameter :: tpoNNNdefault = 6.0
  real(8),parameter :: tpoPPNdefault = 8.0
  real(8),parameter :: tpoPNNdefault = 8.0
  
! added by shan for backward PN
  real(8),parameter :: tpoPNbdefault = 8.0
  real(8),parameter :: tpoPPbdefault = 8.0
  
! added by CWJ for backward PPN, PNN

real(8),parameter :: tpoPPNbdefault = 8.0
real(8),parameter :: tpoPNNbdefault = 8.0  

!  INFORMATION FROM PAST, CURRENT RUNS
  real(8) :: tpoSPEprior, tpoSPE
  real(8) :: tpoPPprior, tpoPP
  real(8) :: tpoNNprior, tpoNN
  real(8) :: tpoPNprior, tpoPN
  real(8) :: tpoPPPprior, tpoPPP
  real(8) :: tpoNNNprior, tpoNNN
  real(8) :: tpoPPNprior, tpoPPN
  real(8) :: tpoPNNprior, tpoPNN
  real(8) :: tpoPPNbprior, tpoPPNb
  real(8) :: tpoPNNbprior, tpoPNNb
  
  real(8) :: tpoPPbprior, tpoPPb
  real(8) :: tpoPNbprior, tpoPNb

  logical :: priortimeinfo

end module timeinfo
!...........................................................................
!
!  GATHERS MISCELLANEOUS INFORMATION FOR END-OF-RUN REPORT
!  
module reporter
   implicit none
   character(len=40) :: spfilename
   character(len=50) :: tbmefilename(50)    ! assumes no more than 50 input files
   real              :: tbmescale(50),spescale(50)      ! scaling of 2-body matrix elements,spes
   character(1)      :: tbmetype(50)        ! type of TBME
   integer           :: ntbmefiles          ! # of tbme files; needs to be initialized
   character(len=50) :: wfn_in_filename

end module reporter
!...........................................................................
module io
  use precisions
  implicit none
  
  character(500) :: outfile              ! filename for outputs
  logical       :: writeout
  integer       :: autoinputfile   = 8
  integer       :: resultfile      = 12
  integer       :: denresultfile   = 32
  integer       :: occresultfile   = 33
  integer       :: wfnfile         = 13
  integer       :: oldwfnfile      = 14
  integer       :: oldtrwfnfile    = 94
  integer       :: mdensfile       = 15 ! mathematica density file
  integer       :: logfile         = 36 ! human-readable log file
  integer       :: binlogfile      = 22 ! binary log file, used for speeding up creations
  integer       :: fragfile        = 46 ! file for fragment info
  integer       :: configfile      = 91  ! file for configuration occupation
  logical       :: auto_readin         ! flag to read in information from .wfn file
  logical       :: auto_input          ! flag to read from autoinput.bigstick file
  logical       :: ham_readin          ! flag to read in hamiltonian files
  logical       :: op_readin           ! flag to read in (nonscalar) one-body operator
  logical       :: strengthflag        ! flag to compute strength functions starting from a pivot
  logical       :: blockstrengthflag   ! flag to compute strength functions starting from a pivot block ADDED
  logical       :: strengthnorm        !  if TRUE then normalize output wfns for strength function ADDED 7.9.8
  logical      :: dotflag               ! if true, compute dot produce, else compute relative entropy adding 7.9.9
  logical      :: complex_green  ! if true, allow for complex energy for green's function; added 7.10.6
  
  logical       :: densityflag         ! to compute density matrices inline
  logical       :: trdensout           ! a flag to write to file suitable for Navratil's TRDENS code
  logical       :: greenflag           ! flag to apply resolve/Green's function to a pivot ADDED 7.8.8
  logical       :: restart_enabled = .true.  ! if set = .false., makes for easier restart
                                       ! else lanczos vectors automatically erased
  logical       :: write_wfn           ! flag to write wavefunctions to file (.wfn) (added by P.G.K)
  logical       :: get_JT               ! flag to enable calculation of J, T; added 7.6.8
  logical       :: baseASCII           ! write out .bas file in human-readable ASCII   added 7.9.1
  logical       :: uncoupled_dens    ! write/read in uncoupled one-body densities/operator to apply (added 7.10.5)
  
  logical       :: print_parity = .true. ! added 7.10.8; if TRUE then print parity -1 or 1 for spectrum
  logical       :: only_odd          ! added 7.10.8; option to compute densities only for odd-J levels; only for M > 0
                              ! this is needed to avoid vanishing CG coefficients
  !
  logical printouthamflag  ! flag for printing out nonzeroes of Hamiltonian
  logical printouttrans1flag  ! flag for printing out nonzeroes  of 1- body transition
  logical :: modeldensities  ! added 7.9.6, May2020
  logical :: stopafterdim    ! added 7.10.5 for computing dimensions 
  logical :: binden ! added 7.10.9, write out 1-b densities as binary

  integer(kind=basis_prec) :: dimbasischeck

  ! KSM:  Added to support writing large files to scratch directories
  !       on HPC machines
  character (len=64) :: nersc_host
  character (len=1024) :: base_scratch_dir  ! like /tmp or $SCRATCH
  character (len=1024) :: scratch_dir       ! add on resultfile
  
end module io

!--------------------------------------------------
!
! module for controlling hermiticity
!
!module hermit

!  logical hermitian

!end module hermit

!--------------------------------------------------------------------------
!  single-particle ORBIT information
!  key derived variable: orbqn
!---------------------------------------------------------------------------
module sporbit
   implicit none
   integer :: numorb(2)  ! # of s.p. orbits
   integer :: numorbmax	        ! max of numorb
   logical :: isoflag           ! in isospin formalism
   logical :: pnwtflag          ! flag to allow p,n to have different weights even if in isospin formalism
                                ! added in version 7.2.8 by CWJ
   logical :: spinless          ! for "spinless" fermions
   logical :: allsameparity     ! flag for all sp orbits same parity
   logical :: allsameW          ! flag for all sp orbits same W
!------------ CREATE A DEFINED TYPE-----------------------------------------
   type orb
      integer :: nr            ! radial quantum number
      integer :: j             ! 2 x j
      integer :: l             ! L
      integer :: par           ! parity = +/- 1
      integer :: w             ! excitation 
   end type orb
   type (orb),allocatable :: orbqn(:,:) ! orbqn(species,iorb)

   character*200 :: sps_path    ! environment variable for searching for directory containing single particle orbits
   integer :: length_sps_path
   integer :: maxorblabel       ! added 7.4.1; used for error traps in reading in matrix element files     
   
contains
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  function parmult
!
! "multiplies" parities where 1 = + and 2 = -
!
!  INPUT: par1,par2
!
      integer function parmult(par1,par2)

      implicit none
      integer par1,par2  ! input
!..........................................
      integer i1,i2,i3

      if(par1 < 1 .or. par2 > 2)then
        print*,' wrong parity ',i1,i2
        stop
      endif

      if(par2 < 1 .or. par2 > 2)then
        print*,' wrong parity (2)',i1,i2
        stop
      endif

      i1 = 3-2*par1
      i2 = 3-2*par2
      i3 = i1*i2
      parmult = (3-i3)/2

      return
      end function parmult
end module sporbit

!---------------------------------------------------------------------------
!  single-particle STATE information
!  key dervied variable: spsqn
!---------------------------------------------------------------------------
module spstate
  implicit none
  integer :: nsps(2) ! # of s.p. states
  integer :: nspsmax
!------------ CREATE A DEFINED TYPE-----------------------------------------
  type spst
     integer :: nr       ! radial quantum number
     integer :: j        ! 2 x j
     integer :: m        ! 2 x jz
     integer :: l        ! L
     integer :: w        ! excitation 
     integer :: par      ! parity
     integer :: orb      ! orbit label
     integer :: group    ! which group -- labels by jz, w, par
     integer :: tr       ! time-reversed state   (added Aug 2010 by CWJ @ SDSU)
  end type spst
  type (spst),allocatable :: spsqn(:,:)  ! spsqn(ispecies,istate)

!--------- INFO ON HAIKU S.P. STATES:---------------------------------------
!  haikus are divided up into "left" built from s.p. states with jz < 0
!                         and "right" built from s.p. states with jz >= 0
!------------------ CONVENTION for hspsqn (and elsewhere)-------------------
!                   if ispecies < 0, then "left" haiku  jz < 0
!                               >0        "right" haiku  jz >= 0
!...........................................................................
  integer :: nhsps0(-2:2),nhsps_eff(-2:2)  
  integer :: nhspsmax
  type (spst),allocatable :: hspsqn(:,:) ! hspsqn(ispecies,istate)

!---------------------- GROUPS -- group together hsps by par, jz, par-------
  integer :: ngroups(-2:2)

  type grouplist
     integer, pointer :: jz(:)
     integer, pointer :: w(:)
     integer, pointer :: par(:)
     integer, pointer :: start(:)
     integer, pointer :: fin(:)
     logical, pointer :: active(:)   
  end type grouplist
  
  type (grouplist) :: group(-2:2)
  logical :: groupwinnow = .true.  ! flag to winnow down groups
end module spstate

!---------------------------------------------------------------------------
!  Information on W (excitation weighting) 
! (equivalent to Nhw in REDSTICK, t in ANTOINE)
!---------------------------------------------------------------------------
module W_info
  implicit none
  integer :: maxW(2)      ! max excitation allowed for p, n
  integer :: minW(2)      ! min excitation allowed for p, n 
  integer :: maxWtot,minWtot     ! max,min total excitation
  integer :: wceiling(2)  ! ceiling to limit creating too many haikus
end module W_info

! NOTE: in 7.6.9 module spsections MOVED to file bsections.f90


!---------------------------------------------------------------------------
!  used to construct the basis
!---------------------------------------------------------------------------
!
! NOTE: in 7.6.9 module templates MOVED to file bbasis_template.f90
!
!...........................................................................
module builder
!
!  USEFUL (?) DERIVED TYPE  used in basislib1.f90 routines

  implicit none
  type alist
     integer,pointer :: loc(:)
  end type alist
  
  type blist
     type (alist), pointer :: list(:)
  end type blist
  
end module builder


!---------------------------------------------------------------------------
!  Information about haikus
!  a HAIKU is a slater determinant, constructed from s.p. states 
!  "left haikus" built from s.p. states with jz < 0
!  "right haikus" built from s.p. states with Jz >= 0
!---------------------------------------------------------------------------
module haiku_info
  implicit none

                          ! nhsps is set in subroutine limit_hspspace in bspstates.f90
  integer             :: nhsps(-2:2) !these needs to be made consistent
                                                   ! and so does this:
  integer             :: Nword(2)           ! # of words needed to construct a haiku

  integer             :: minNh(-2:2) ! min number of particles in any haiku
  integer             :: maxNh(-2:2) ! max number of particles in any haiku
  integer,allocatable :: minWh(:,:)                ! min w found in haiku w/ n particles
  integer,allocatable :: maxWh(:,:)                ! maxWh(ispecies, nh); ispecies = +/- 1, 2
                                                   ! if ispecies -1,-2  = "left" haikus  jz < 0
                                                   !  if ispecies +1,+2 = "right" haikus jz >= 0
end module haiku_info

!  MODULE haikus moved to file bbasis_haikus.f90 (as of 7.6.9)

!  MODULE blocks moved to file bbasis_blocks.f90 (as of 7.6.9)

!  MODULE slaterlimits moved to file bbasis_limits.f90 in 7.6.9

!  MODULE basis moved to file bbasislib.f90 in 7.6.9

!---------------------------------------------------------------------------
!  sectors of the basis are divided by (first) jz, then parity, finally W
!  sectors are constructed from blocks of haikus
!---------------------------------------------------------------------------
module sectors

  implicit none
  integer     :: nsectors(2)                ! # of sections
  type sectorinfo
     integer(8)  :: xsdstart                          ! where, in the list of X-species sds, does this start
     integer(8)  :: xsdend
     integer(8)  :: nxsd                              ! # of SDs in this sector
	 integer(8)  :: ncxsd
     integer  :: jzX
     integer  :: parX
     integer  :: wX
     integer  :: ncsectors                         ! # of opposite-species sections can couple to
     integer, pointer :: csector(:)                ! list of opposite-species sections
     integer  :: nhblocks                          ! # of pairs of haiku blocks used to construct this section
     integer, pointer :: rhblock(:),lhblock(:)     ! list of right- and left-haiku blocks 
     integer(8), pointer :: blockstart(:),blockend(:) ! where, for a given block, this starts and ends
     integer,pointer :: cut(:)                     ! ADDED 7.4.1 to be used to cut proton sectors for fragments
     integer         :: ncuts
     integer(kind=8) :: basisstart,basisend        ! ADDED 7.4.9: for protons, where basis starts and stops
     logical :: ibelong                            ! ADDED 7.5.9: denotes if a sector belongs on a given MPI process
                                                   ! especially when vectors are broken into fragments
  end type sectorinfo
  type mastersect
     type (sectorinfo), pointer :: sector(:)
  end type mastersect
  
  type (mastersect) :: xsd(2)
   
end module sectors

!.... MODULE bugs MOVED TO bhopslib.f90....(and renamed hoppy)

!.... MODULE nutple_info moved to bjumplib_ntuples.f90


!=============================================================
!  THE FOLLOWING MODULES ARE USED TO CONSTRUCT THE JUMPS
!==============================================================


!==============================================================

module geneologies
implicit none

type genlink
   type (genlink), pointer :: next
   integer, pointer :: gen(:)   ! stores list of block hops; if > 0 then right blockhop, if < 0 then left
   integer :: iblock,fblock
!-------------- NOTE: I MIGHT BE ABLE TO GET RID OF THESE THROUGH SOME REVISION
   integer :: irblock,ilblock   ! indices of initial, final blocks
   integer :: frblock,flblock   ! indices of initial, final blocks
   integer(8) :: id

end type genlink

type (genlink), target :: geneology
	
!$OMP threadprivate(geneology)

end module geneologies


!==============================================================
!
!  definition of derived types for sector jumps which overlay jumps
!
module jumpdef

implicit none


type conjump  ! information on conjugate jumps
   integer ncjmps
   integer, pointer :: cjump(:)
end type conjump

type jumper
   integer(8) :: njumps
   integer(8) :: nstart
!   integer(8) :: stride 
   logical diag
!............. ADDED IN 7.5.3....  FOR INTRON DELETION SHIFTS.........
   logical :: containsexons     ! if "exons" (active jumps) for current MPI proc are here
   integer(8) :: exstart,exstop  ! stop/stop of exons
   integer(8) :: exshift        ! shift of exons for storage
!............. ADDED in 7.9.0.... FOR SPLITTING SECTOR JUMPS.........
!           these are filled out in bjumplib_master.f90
!           and then used to reconstruct the sector jumps
!
   integer(4) :: nsplit   ! number of splits
   integer(8),allocatable :: splitstart(:)   ! new starts for splits    
   	
end type jumper

type jumpsect
   integer :: nsectjumps
   integer, pointer :: isector(:)
   integer, pointer :: fsector(:)
   type (jumper), pointer :: sjmp(:)
   type (conjump), pointer :: csjmp(:)   
end type jumpsect

end module jumpdef
!=====================================================

!
!  contains info on 1- and 2-body jumps
!
!  "JUMPS" are the fundamental constructs that define operations
!  a jump is always within a species and contains the following information:
!  -- index of initial and final (proton or neutron) Slater determinants
!  -- index/ices of operators that take us from initial to final Slater det
!     here either a 1-body (a^dagger a) or 2-body (a^dagger a^dagger a a) operator
!  -- phase information from anticommutation
!
! the master routines for creating jumps are found in bjumplib1.f90
!         master1bodyjumps  and master2bodyjumps
!  however jumps themselves are constructed (from "chains", an intermediate construct)
!  in subroutine weld in bjumplib3.f90
!

module jumpNbody
   use precisions
   use jumpdef
   implicit none
   
!---- MADE EXPLICIT IN 7.7.4 -- memory storage per jump
   integer(4) :: bytesper2Bjump = 25    ! 2 x 8 isd/fsd + 8 m.e. label + 1 phase
   integer(4) :: bytesper1Bjump = 21    ! 2 x 8 isd/fsd + 2 x 2 indexs + 1 phase 
      
   logical :: MPIjumps  = .true.
   logical :: skipzeros =  .false.  !  flag to skip jumps with "known" zero matrix elements 
                          ! this is straightforward for PP, NN jumps; not easy for others

!---------- SECTOR JUMPS define clumps of jumps between proton and neutron sectors
!           they are created in routines
!            set1bsectorjumps and set2bsectorjumps in bjumplib1.f90
!           and set3bsectorjumps in bjumplib5.f90
!
!           each sector jumps is associated with "conjugate" sector jumps
!           for the opposite species
!       "conjugate" sector jumps are computed in 
!           masterfindconjumps1b and masterfindconjumps2b in bjumplib3.f90
!           and masterfindconjumps3b  in bjumplib5.f90
!

   type (jumpsect), target :: x1bjump(2)
   type (jumpsect),target :: x2bjump(2)
!-------------------
   ! one-body operators stay within species and include kinetic energy + external
   ! potential.  The external potential could be a sum over a core.
   ! We also have components of H_{xy} which have the same form.
   integer(8) :: totn1bjumps
   ! n1b_isd:  map neutron jump id to initial proton slater determinant
   ! n1b_fsd:  map neutron jump id to final proton slater determinant
   integer(kind=basis_prec),allocatable, target :: n1b_isd(:),n1b_fsd(:)
   integer(kind=basis_prec),allocatable, target :: n1b_isd0(:),n1b_fsd0(:)   ! used only when printing out jumps

   ! n1b_cop:  map neutron jump id to to creation op single particle label   
   ! n1b_dop:  map neutron jump id to to destruction op single particle label
   ! n1b_phase maps proton jump id to phase to apply to matrix element
   integer(kind=2),allocatable, target :: n1b_cop(:),n1b_dop(:)   ! used for 1-body jumps
   integer(kind=1), allocatable, target :: n1b_phase(:)

   integer(8) :: totp1bjumps
   ! p1b_isd:  map proton jump id to initial proton slater determinant
   ! p1b_fsd:  map proton jump id to final proton slater determinant
   integer(kind=basis_prec),allocatable, target :: p1b_isd(:),p1b_fsd(:)
   integer(kind=basis_prec),allocatable, target :: p1b_isd0(:),p1b_fsd0(:)  ! used only when printing out jumps

   ! p1b_cop:  map proton jump id to to creation op single particle label   
   ! p1b_dop:  map proton jump id to to destruction op single particle label
   !   they are used together    a^+_c a_d for hamiltonian entry
   ! p1b_phase maps proton jump id to phase to apply to matrix element
   integer(kind=2),allocatable, target :: p1b_cop(:),p1b_dop(:)   ! used for 1-body jumps
   integer(kind=1), allocatable, target :: p1b_phase(:)

   ! two body operators
   ! H_{xx} = (1/4) \sum {V_{ijkl}^{(xx)} {a^+_i(x) a^+_j(x) a_l(x) a_k(x)}   paper eqn 15
   ! must be diagonal in quantum numbers because species y is not involved.
   integer(8):: totn2bjumps,totp2bjumps
   integer(kind=basis_prec),allocatable, target :: n2b_isd(:),n2b_fsd(:),p2b_isd(:),p2b_fsd(:)
   integer(kind=basis_prec),allocatable, target :: n2b_isd0(:),n2b_fsd0(:),p2b_isd0(:),p2b_fsd0(:) ! used only when printing out jumps

   integer,allocatable, target :: n2b_op(:),p2b_op(:)
   real(4), allocatable, target :: n2b_me(:), p2b_me(:)
   integer(kind=1), allocatable, target :: n2b_phase(:),p2b_phase(:)
   
   logical :: onebodyonly    ! flag to signal only one-body; added 7.6.8
   logical :: alt_PNb_loop   ! flag to allow for alternate loop for PPb in mat-vec mult added 7.8.0
   
end module jumpNbody
!==============================================

!
!  INFORMATION ON COUPLED two-body pairs
!
!
      module pairdefcoupled

!
!  information on COUPLED two-body pairs 
!   restricting on W (truncation #) for storage of coupled TBMEs
!

      type coupair_qn
        integer :: par
        integer :: indx
        integer :: ia,ib
      end type coupair_qn

      type coupairinfo
        type (coupair_qn), allocatable :: pairc(:)
        integer :: meref(2), mestart(2)  ! i parity
      end type coupairinfo

!--- ADDED 7.9.3--- ALLOWS FOR SORTING INTO BLOCKS WITH W----  
      type coupairinfo_adv
        type (coupair_qn), allocatable :: pairc(:)
		integer :: nw,wstart,wstop
		integer,allocatable :: wlist(:)
        integer :: meref(3), mestart(3)  ! i parity
      end type coupairinfo_adv

      end module pairdefcoupled

!==============================================================
!
!  INFORMATION ON COUPLED MATRIX ELEMENTS (i.e., single -particle energies,
!   couple 2-body matrix elements, etc); 
!
!  ENCODING for Hamiltonian matrix elements
!    Vxx(ab,cd):
!     a >= b, c >= d,  (ab) >= (cd)
!       pair1 = a*(a-1)/2 + b  -> XXcouplemap(pair1)
!      pair2 = c*(c-1)/2 + d   -> XXcouplemap(pair2)
!    pair1 >= pair2
!   iref = XXcouples(it)%meref(ipar) 
!   istart = XXcouples(it)%mestart(ipar) 
!   indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!
!   Vpn(ab,cd)
!   a >= c; if a==c, b >= d
!      pair1 = numorb(2)*(a-1) + b  -> PNcouplemap(pair1)
!     pair2 = numorb(2)*(c-1) + d -> PNcouplemap(pair2)
!     iref = PNcouples%meref(ipar) 
!     istart = PNcouples%mestart(ipar) 
!    indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart   
!
!  ENCODING for two-body density matrix elements  ADDED 7.9.1
!    Vxx(ab,cd):
!     a >= b, c >= d
!       pair1 = a*(a-1)/2 + b  -> XXcouplemap(pair1)
!      pair2 = c*(c-1)/2 + d   -> XXcouplemap(pair2)
!    pair1 >= pair2
!   iref = XXcouples(it)%meref(ipar) 
!   istart = XXcouples(it)%mestart(ipar) 
!   indx = (pair1-iref)*numXXpairs(ipar) + pair2-iref+istart    
!
!   Vpn(ab,cd)
!   a >= c; if a==c, b >= d
!      pair1 = numorb(2)*(a-1) + b  -> PNcouplemap(pair1)
!     pair2 = numorb(2)*(c-1) + d -> PNcouplemap(pair2)
!     iref = PNcouples%meref(ipar) 
!     istart = PNcouples%mestart(ipar) 
!    indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart  
!
!  REQUIRES MODULE pairdefcoupled

   module coupledmatrixelements
	   use pairdefcoupled
	   implicit none
!               restrict on total W

	    integer ncouplesXX(2)
	    type (coupairinfo) :: XXcouples(2)
	    integer, allocatable, target :: PPcouplemap(:), NNcouplemap(:)
	    integer ncouplesPN
	    type (coupairinfo) :: PNcouples
     integer, allocatable :: PNcouplemap(:)
	 
     logical dens2bflag     ! flag for which encoding to use--Hamiltonian or 2-body density ADDED 7.9.1
	 logical diag_den2b     ! flag for writing out only 2-body densities for initial=final, only Jt=0  ADDED 7.9.5
!	 logical iso_den2b    ! flag to write out diagonal 2-body densities with good isospin (may leave out options)
!...........................................................................
!-----------   Logical variables to control how interaction matrix elements
!-----------   are put into pp, nn, and pn components when using isospin
!-----------   formalism
	     logical pp_int,nn_int,pn_int,iso_int,isov_int,isot_int, coul_int, mixed_T

!--------- INFORMATION ABOUT FOLDERS TO FIND MATRIX ELEMENTS

	      character*200 :: int_path
	      integer :: length_int_path

		  logical :: emptyHam     ! no interaction read
		  real  :: hw   ! for c.m.
		  real :: ephshift   ! overall shift in energy from particle-hole tranformation

		  character*55 intfilename
		  real,allocatable, target :: pspe(:),nspe(:)
		  real,allocatable, target :: psppot(:,:),nsppot(:,:)  ! added in 7.7.9; off-diagonal single-particle 
		                                                       ! potentials
		  logical :: call_spe  ! flag whether or not to call SPE routines
   
		  integer :: nv2bmedim(0:2)  ! 0 = pn, 1 = pp, 2 = nn  
		  integer :: norb_allow(2)  ! # of allowed orbits
		  integer :: nv2bme_tot

!----------- coupled  TBMES -----------------

		  type vjs
		     integer  :: Jmin,Jmax
		     real,pointer :: v(:)
		  end type vjs
		  type (vjs), allocatable,target :: ppme(:),pnme(:),nnme(:)
!----------- COUPLE two-body density matrices
          !
		  type vjs2b
		     integer  :: Jabmin,Jabmax,Jcdmin,Jcdmax,Jmin,Jmax
		     real,pointer :: v(:,:,:)
		  end type vjs2b
		  type (vjs2b), allocatable,target :: pp2bden(:),pn2bden(:),nn2bden(:)	  
			  
   end module coupledmatrixelements

!==============================================================
!
!  INFORMATION ON UNCOUPLED two-body pairs
!  
! REQUIRES MODULE pairdef
!
      module pairdef

      implicit none
!
!  information on UNCOUPLED two-body pairs | ia, ib : JM >
!
      type pair_qn
        integer :: M
        integer :: par
        integer :: W
        integer :: indx
        integer :: ia,ib  ! important: if ia, ib > 0, then m is >= 0 
                          ! else if ia < 0, then ma < 0
        integer :: tr        ! index of time-reversed pair
        integer(1) :: trphase   ! phase from time-reversal
		integer :: wb   ! added 7.9.4 for ordering by W of destruction operator in a^+ b
      end type pair_qn

      type pairinfo
        type (pair_qn), pointer :: pair(:)
        integer(kind=8), pointer       :: meref(:,:), mestart(:,:)  ! indices are M, parity
		integer(kind=8),pointer :: block(:,:)  ! ADDED in 7.9.1, size of blocks
      end type pairinfo
      end module pairdef

!==============================================================
!  storage of uncoupled two-body matrix elements
!    key derived  variables :vtbme; hmatXX,hmatpp,pairXX,pairpn
module interaction

  use pairdef
  use pairdefcoupled
  use precisions
  implicit none

!--- # of uncoupled TBMEs ---------------------
!  integer(8) :: nmatpp,nmatnn  NOT USED ANY MORE
  integer(8) :: nmatpn
  integer(8) :: nmatXX(2)
  
!---------- uncoupled TBMEs ------------------

  real, allocatable, target :: hmatpp(:),hmatnn(:)
  real, allocatable, target :: obsmatpp(:),obsmatnn(:)
  real, allocatable,target :: hmatpn(:),obsmatpn(:)
  real, pointer      :: vmatpn(:)
!------ ADDED 7.10.4 for Hermitian conjugates  
  
  
!--------- for 2-body density matrices, need double precision because adding many small numbers
  real(kind=obs_prec), allocatable, target :: dmatpp(:),dmatnn(:),dmatpn(:)  
!------ ADDED 7.10.4 for Hermitian conjugates  
real(kind=obs_prec), allocatable, target :: dmatpphc(:),dmatnnhc(:),dmatpnhc(:)  

!------------ coding of the pairs -------------
!  take states i,j,k,l  from list in hspsqn
!
!  We encode a_i a_j, or a^+_i a^+_j, as follows:
!
!  assume i >= j
!
!  then  (i,j) =>  i(i-1)/2 +j 

  integer npairXX(2)
  type (pairinfo) :: XX2(2)  
  integer :: npairpn
  type (pairinfo) :: PN2
  integer :: nrhophX(2)   
  type (pairinfo) :: rhophX(2)  ! added 7.9.4  -- one-body particle-hole operator for protons, neutrons  
  
  integer(kind=8), allocatable :: cpnpair(:,:),dpnpair(:,:)
  integer(kind=8), allocatable,target :: mappairPP(:),mappairNN(:),mappairpn(:)
!.... ADDITIONS in 7.10.2 for use of time-reversal symmetry....  
  logical :: useTRnew = .false.
  logical :: useTRdef = .true.
  logical :: useTR = .false.  ! use time-reversal to restrict some matrix elements
			      ! that is, Vijkl = V-i -j -k -l
  integer(1), allocatable :: phasepnpair(:,:)

!C--------------- ENCODE particle-hole pairs
  
  type (pair_qn),allocatable,target :: php(:),phn(:)  !particle-hole
  
  integer nphp,nphn           ! # of particle-hole operators
  integer,allocatable,target :: mapphp(:),mapphn(:)
  
!------ SPECIALIZED FLAG-- FOR RANDOMIZED UNCOUPLED MATRIX ELEMENTS----
!  SET FALSE EXCEPT FOR SPECIAL CASES

  logical :: random_mscheme = .false.  


end module interaction
!===========================================================================
!
!  information for deformed single-particle spaces; only used in special cases
!
module interaction_deformed
	implicit none
	logical :: deformed = .false.
	real, allocatable, target :: speXm(:,:)

end module interaction_deformed
!===========================================================================
! module diagh moved to bdiagh.f90

!===========================================================================
! Observables, densities, applied operators, etc.
!===========================================================================
module obs
  use precisions
  implicit none
  real(kind=obs_prec) :: xj2,xt2
  real(kind=4),allocatable :: energy(:),xjlist(:), xtlist(:)
  real(kind=obs_prec) :: xj2tmp,xt2tmp ! for global sum /allreduce/ (added by P.G.K.)
  integer(kind=4),allocatable :: parity(:)
  logical :: twoobsflag
  logical :: skip_T2  =.false. ! new flag in 7.11.2, in order to eliminate computing T2 of state to reduce memory
                              ! should be set = .false. by default
  
end module obs
   
!=======================================================================
!module densities moved to bdensities.f90
!
!===================================================================
!
! coupled matrix elements for transition operators   ADDED 7.5.2
! FOR application to a wavefunction
! (corresponding uncoupled matrix elements found in module onebodypot)

module opmatrixelements
	implicit none
	logical :: pnoperators    ! flag to signal if operators are in pn format or isospin format
	  real, allocatable :: op1bod(:,:)  ! reduced matrix elements for one-body operartor; to be deleted
	
	real,allocatable,target :: pop1bod(:,:),nop1bod(:,:)  ! singly-reduced matrix elements for proton, neutron 1-body operators
	integer :: jop                         ! J of transition operator
	integer :: top                         ! (optional) T of transition operator
	
    logical :: subtract_enabled=.false.       ! for MKGK application of r^2 Y00
    logical :: subtract
	
end module opmatrixelements

!===================================================================
!
! application of one-body (uncoupled) potentials, including self-consistent mean-field for starting
! the original coupled matrix elements found in module opmatrixelements
!
      module onebodypot

      implicit none

      real, pointer :: pPOT(:,:),nPOT(:,:)

      real, allocatable, target :: ppot_obs(:,:),npot_obs(:,:)
      real, allocatable, target :: ppot_h(:,:),npot_h(:,:)

      end module onebodypot

!===========================================================================
! LANCZOS_INFO
! NOTE some info found in module flagger in bmodule_flags.f90
!
!===========================================================================
module lanczos_info
  use precisions
  implicit none
  
  integer(4):: niter  ! # of iteractions
  integer nkeep  ! # iterations to keep
  integer :: actual_iterations   ! counts how many actual iterations have taken place
                                 ! added by CWJ 8/2010
  
  logical :: reorthog = .true.   ! logical flag for reorthogonalization
  logical :: storevectorsenabled = .false.  ! if false, store Lanczos vectors in RAM if possible
  logical :: writetodisk

  integer(4) :: lvec_file = 30
  integer :: coef_file = 21  ! file for alpha, beta coefficients
  real (kind = lanc_prec),allocatable :: lvec(:,:)
  character(2) :: lanczchar
  integer :: kpow  ! exponent for polynomial filtered Lanczos; otherwise == 1 
!------------- CONVERGENCE CONTROL ------
  logical :: fixed_iter     ! whether or not to used fixed # of iterations
!--------------THICK RESTART ---------------------
  integer :: nthick_keep     ! # of vectors to keep when restarting
  logical :: thick_restart
  integer :: max_num_thick_restart  ! max # of times to thick-restart
  real(8) :: Etarget        ! target energy for finding excited states
  logical :: targetX    
!----------- CONSTRUCTING INITIAL PIVOT --------------  ADDED in 7.6.3 ----
  logical :: initializing_enabled=.false.
  integer :: initial_maxiter     ! max # of initial iterations
  logical :: initializing        ! flag to signal initializing, so don't do everything
  logical :: applypntrace         ! apply  average pn diagonal
  logical :: applyXXonly         ! only apply PP,NN, not PN
  logical :: diagonalsectorsonly ! only apply to "diagonal" sector jumps
  integer :: nsave               ! # of vectors to keep for starting, typically 1 but could be more

!---------- MISC (MOSTLY OBSOLETE)
  logical :: coef_only    ! used to get just coefficients
  integer :: time_loop_max = 6400  ! # of loops for timing

end module lanczos_info

!===========================================================================
! TIMING
!===========================================================================
module timing
  implicit none

  real(8) :: startall, endall
  real(8) :: startallalt,endallalt
  real(8) :: startbasis,endbasis
  real(8) :: startham,endham
  real(8) :: startlanczos,endlanczos
  real(8) :: startobs,endobs,startwriteobs,endwriteobs
  real(8) :: startdens,enddens
  real(8) :: starts1b,ends1b
  real(8) :: startp2b,endp2b
  real(8) :: startn2b,endn2b

  real(8) :: time1body
  real(8) :: time2body
  real(8) :: timereorthog
  
! Added by PGK ( extra timing ).............................................
  real(8) :: time_applobs,time_writeobs
  real(8) :: time_reduce
  real(8) :: time_gather
  real(8) :: time_distr
  real(8) :: time_writeeigvec
  real(8) :: time_read_orthog
  real(8) :: time_pn
  real(8) :: time_pp
  real(8) :: time_nn
  real(8) :: time_ppb,time_pnb  ! 'backwards' 
  real(8) :: time_pro
  real(8) :: time_dot
  real(8) :: time_swp
  real(8) :: time_spe
  real(8) :: time_rre     ! reduce (reorthog. only)
  real(8) :: time_ror     ! read in reorthog.
  real(8) :: time_sort1b     ! read in reorthog.
  real(8) :: time_sort2b     ! read in reorthog.
  real(8) :: time_hops    ! time to generate hops
  real(8) :: time_ppp, time_nnn, time_pnn, time_ppn
  real(8) :: time_ppnf,time_ppnb, time_pnnf,time_pnnb  ! added  7.9.2
  real(8) :: time_jmpcnt, time_jmpcnt_last   ! for counting jumps
  real(8) :: time_meset,  time_meset_last    ! set up for matrix elements
  real(8) :: time_munc,  time_munc_last    ! uncouple matrix elements
  real(8) :: time_desc,  time_desc_last     ! compute descendents setting up for jumps
  real(8) :: time_intron,  time_intron_last    ! map and delete introns
  real(8) :: time_pivot, time_pivot_last       ! time to set up pivot
  real(8) :: time_p1b,time_p1b_last, time_n1b,time_n1b_last  ! time for one-body densities

  real(8) :: time_hmult
  real(8) :: time_eig   ! time to get eigenvalues
  real(kind=8) :: time_eigvec
  real(kind=8) :: time_trdens
  
! Separate "timelast" variables are necessary to prevent "time-mixing"......
  real(8) :: timelast1    ! 1-body
  real(8) :: timelast2    ! 2-body
  real(8) :: timelast3    ! reorthog
  real(8) :: timelast4    ! applobs
  real(8) :: timelast5    ! reduce (applyh only)
  real(8) :: timelast6    ! distribution
  real(8) :: timelast7    ! writeeigvec
  real(8) :: timelast8    ! read_orthog
  real(8) :: timelast9    ! gather
  real(8) :: time_pn_last
  real(8) :: time_pp_last ! 2-body (pp)
  real(8) :: time_nn_last ! 2-body (nn)
  real(8) :: time_ppb_last  ! backwards pp
  real(8) :: time_pnb_last ! backwards pn
  real(8) :: time_ppp_last, time_ppn_last, time_pnn_last, time_nnn_last
  real(8) :: time_ppnf_last, time_ppnb_last,time_pnnf_last,time_pnnb_last ! added 7.9.2
  real(8) :: time_hmult_last
  real(8) :: time_pro_last! project Lanczos vector
  real(8) :: time_dot_last! dot product
  real(8) :: time_swp_last! "swap" v <---> w
  real(8) :: time_spe_last! find eigen-values
  real(8) :: time_rre_last! reduce (reorthog. only)
  real(8) :: time_ror_last! read in reorthog.
  real(8) :: time_eigvec_last   ! build eigenvectors
  real(8) :: time_trd_last   ! build trdens eigenvectors
  real(8) :: time_eig_last

! Temporary timing .........................................................
  real(8) :: t1start,t1end
  real(8) :: t8start,t8end
  real(8) :: t2,t2last
  real(8) :: t3,t3last
  real(8) :: t4,t4last
  real(8) :: t5,t5last
  real(8) :: t6,t6last
  real(8) :: t7,t7last

  real(8) :: jumptime,jumptimestart   ! used to check on timing across MPI processes
  real(8) ::time_hmult_wait,time_wait_last    ! USED ONLY FOR DEBUGGING, TIME WAITING TO REDUCE
!...........................................................................
end module timing

!===========================================================================
!
! for computing trace
!
module tracy
  implicit none
  integer(8) :: noperation, ioperation

  integer, allocatable :: indx1(:), indx2(:)
  real, allocatable :: me(:)

  type mex
     integer(kind=8) :: ncol
     real, allocatable :: v(:)
     integer(kind=8), allocatable :: col(:)
  end type mex

  type (mex), allocatable :: metrace(:)

  logical :: tracecount
  
!.... ADDED IN 7.6.3: COMPUTE THE PN TRACE SECTOR BY SECTOR
!           IN 7.6.6: TRACE (actually: average) OF ALL SECTOR BY SECTOR
   real, allocatable :: sectortrace(:)
!....... ADDED IN 7.9.0.... so one only has to compute the diagonal...   
   logical :: getvariance = .false.
   real(kind=8) :: tr1,tr2
   
   

end module tracy

!===================================================================
