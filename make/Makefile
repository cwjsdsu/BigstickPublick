############################################################################
#  NOTE: FOR LIST OF OPTIONS 
#  "make help"
#                                                                          #
#                     Make file for BIGSTICK                               #
#                                                                          # 
#   1) Correctly insert the names of your compilers into the               #
#      enviornment variables at the top. FXX abd SERIAXX_COMPILER          #
#      should be the same. MPIXX_COMPILER is the name of you f90 and f77   #
#      MPI compiler.                                                       #
#                                                                          #
#   2) For OpenMP complilation, specify the proper compiler and linker     #
#      directives in the variables OPENMP_COMPILE and OPENMP_LINK, e.g.,   #
#      for ifort, this would be -openmp.                                   #
#                                                                          #
#   The target builds are:                                                 #
#                                                                          #
#   1) serial: 'make serial' builds the serial version with executable     #
#      name bigstick.x                                                     #
#                                                                          #
#   2) openmp: 'make openmp' builds the openmp version with executable     #
#      name bigstick-openmp.x                                              #
#                                                                          #
#   3) mpi: 'make mpi' builds the mpi version with executable              #
#      name bigstick-mpi.x                                                 #
#                                                                          #
#   4) mpi+openmp: 'make openmp-mpi'  builds the mpi+openmp version with   #
#      executable name bigstick-mpi-omp.x                                   #
#                                                                          #
#   In general, a simple 'make' will not work to build the serial version  #
#   as the mpif.h file may not be present, which is copied from the        #
#   file mpif.h.template. Like dummy_MPI_libs.f90, this is a set of        #
#   "dummy" set up information and routines in order to compile without    #
#   mpi activated.                                                         #
#                                                                          #
#   The current version has the ability to select the precision for        #
#   Lanczos vectors. This is done in the module Lanczos_precision          #
#   located in bmodules.f90. BIGSTICK is delivered with single-precision   #
#   as the default. To change to double precision, change the comments     #
#   in the module.                                                         #
#                                                                          #
############################################################################
src_dir=  ../src
bin_dir=  ../bin
F90= ifort #gfortran                   #  default FORTRAN 90 compiler name
F77= ifort # gfortran                   #  default FORTRAN 77 compiler name
PGI_COMPILER=pgfortran                  # Portland Group compiler
SERIAL90_COMPILER=  ifort #gfortran     #  serial FORTRAN 90 compiler name
SERIAL77_COMPILER=  ifort #gfortran     #  serial FORTRAN 77 compiler name
GFORTRAN_COMPILER =  gfortran #-fbounds-check "  # GNU fortran compiler; flag corrects problem with MPI calls
G95_COMPILER = g95 # # alternate gfortran compiler
MPI90_COMPILER=  mpiifort     #  MPI FORTRAN 90 compiler name
MPI77_COMPILER=   mpiifort     #  MPI FORTRAN 77 compiler name
MPI90_COMPILER_SIERRA= mpif90 # mpiifort     #  MPI FORTRAN 90 compiler name for SIERRA
MPI77_COMPILER_SIERRA= mpif90 # mpiifort     #  MPI FORTRAN 77 compiler name for SIERRA
MPI90_COMPILER_GRIZZLY= mpiifort     #  MPI FORTRAN 90 compiler name for SIERRA
MPI77_COMPILER_GRIZZLY= mpiifort     #  MPI FORTRAN 77 compiler name for SIERRA
MPI90_COMPILER_VULCAN=  mpixlf77_r  # mpigfortran  MPI FORTRAN 90 compiler name for VULCAN
MPI77_COMPILER_VULCAN=  mpixlf77_r     # mpigfortran # MPI FORTRAN 77 compiler name for VULCAN
MPI90_COMPILER_MIRA=  mpixlf90_r  # mpigfortran  MPI FORTRAN 90 compiler name for VULCAN
MPI77_COMPILER_MIRA=  mpixlf77_r     # mpigfortran # MPI FORTRAN 77 compiler name for VULCAN
EDISONF90_COMPILER=ftn     #  MPI FORTRAN 90 compiler name for EDISON
EDISONF77_COMPILER=ftn     #  MPI FORTRAN 77 compiler name for EDISON
MPI90_COMPILER_EDISON= ftn # mpiifort     #  MPI FORTRAN 90 compiler name for EDISON
MPI77_COMPILER_EDISON= ftn # mpiifort     #  MPI FORTRAN 77 compiler name for EDISON
MPI90_COMPILER_PERLMUTTER= ftn # mpiifort     #  MPI FORTRAN 90 compiler name for EDISON
MPI77_COMPILER_PERLMUTTER= ftn # mpiifort     #  MPI FORTRAN 77 compiler name for EDISON
OPT_edison="-c -fast -no-ipo"
OPT_perlmutter="-c -cpp -D_MPI=1 -Ofast -no-pie -fallow-argument-mismatch"
# NOTE optimization flag -O3 does not seem to speedup
FGCFLAGS = "-c -O2  -cpp" # "-c -O2 -std=legacy -cpp" #  -fbounds-check " # GFORTRAN COMPILER DIRECTIVE
F77CFLAGS =  -c -O2 -cpp #-fbounds-check      #  Fortran 77 compiler directives
F90CFLAGS      =  -c -O2 -cpp #-fbounds-check    #  Fortran 90 compiler directives
MPICFLAGS =" -O2 -cpp -D_MPI=1 " 
FGCFLAGS3B =" -c -O2 -cpp -D_BD3=1 " # Fortran 90 compiler directive for including 3-body 
#INTELCFLAGS    =  -c -O2    # default Intel compiler flags
PGCFLAGS       ="-c -fast  " # PGI Fortran compiler directives
LFLAGS_all     =  #-fbounds-check       #  Linker directives
LFLAGS_ifort   =    # -Wl
LFLAGS_VULCAN   =          # -Wl
LFLAGS_MIRA  =          # -Wl
OPENMP_COMPILE=-qopenmp      #  OpenMP compiler directive
OPENMP_LEGACY_COMPILE=-openmp      #  Legacy Intel OpenMP compiler directive
OPENMP_LINK=-qopenmp         #  OpenMP linker directive
OPENMP_LEGACY_LINK=-openmp        #  OpenMP Intel linker directive
OPENMP_PGI=-mp              # OpenMP directive for Portland Group
OPENMP_G95=-openmp              # OpenMP directive for GNU compilers
OPENMP_VULCAN_COMPILE= -qsmp=omp      # -fopenmp # OpenMP compiler directive
OPENMP_VULCAN_LINK=-qsmp=omp         #  -fopenmp # OpenMP linker directive
OPENMP_MIRA_COMPILE=" -qsmp=omp -qnostrict -O3 -q64 -qtune=qp -qarch=qp -qsimd=auto "     # -fopenmp # OpenMP compiler directive
OPENMP_MIRA_LINK=-qsmp=omp         #  -fopenmp # OpenMP linker directive
OPENMP_EDISON_COMPILE=-openmp # OpenMP directive for Edison
OPENMP_CORI_COMPILE=-qopenmp # OpenMP directive for Edison
OPENMP_PERLMUTTER_COMPILE=-fopenmp # OpenMP directive for Edison
OPENMP_EDISON_LINK="-openmp -lifcoremt" # tmp directive to account for error in ftn wrapper
OPENMP_CORI_LINK="-qopenmp -lifcoremt" 
OPENMP_PERLMUTTER_LINK="-fopenmp " 
JAGFLAGS       =   -c -fast # compiler directives for Jaguar


XSTICK_OBJECTS =     \
bmodules_main.o      \
bmodules_flags.o     \
bmodules_parallel.o  \
bmodules_3body.o     \
bmodules_block.o \
bsections.o \
bmpi.o               \
butil.o              \
bbasis_templates.o \
bbasis_haikus.o \
bbasis_blocks.o \
bbasis_limits.o \
bbasislib.o  \
bfragments.o \
bparallel_opstat.o \
breorthog.o          \
bvectorlib.o        \
bwfnlib.o           \
bspstates.o \
bhopslib.o \
bdiagh.o \
btbme_adv.o \
bparallel_jumpstart.o \
bparallel_opbundles.o \
bparallel_util.o \
bparallel_report.o \
bparallel_main.o \
bjumplib_descent.o \
bjumplib_ntuples.o \
bsppot.o  \
btbme.o      \
btbmepn.o    \
bsector.o  \
bjumplib_chains.o \
bjumplib_weld.o \
bjumplib_master.o \
bdensities.o \
bcentroidlib.o \
bpartitionlib1.o \
bapplycentroid.o \
bapplyhlib1.o \
bapplyhlibOMP.o \
bapplyhlibOMP_block.o \
bapplyhblocklib.o \
bapplyobslib1.o \
bparticleocc.o \
bblock_algebra.o \
blanczoslib.o \
blanczosmain.o \
bsortjumpslib.o \
b3bme.o \
bpanashlib.o \
bigstick_main.o      \
bmenu.o              \
btimelib.o           \
bloglib.o            \
binput.o             \
binput_tbme.o        \
bpandyalib.o        \
bwlib.o              \
bXYbmelib.o          \
bdeformedlib1.o  \
bdeformedlib2.o       \
b3bme_input.o        \
b3bme_input_cJ.o     \
bexpect.o            \
bdensmathlib.o   \
bwritebinden.o \
bdenslib1.o          \
bdenslib2.o          \
bdenslib3.o        \
bdenslib4.o \
bdenslib5.o \
bdenslib6.o \
boutputlib.o         \
b3bmelib2.o          \
b3from2bme.o         \
bintron.o          \
btracer.o            \
lapackeig.o          \
dlapackeig.o         \
dlapackeig_adv.o \
dgesvpackage.o \
zlapack.o \
dlamc3.o   \
libra.o \
backtracestub.o                

%.o : $(src_dir)/%.f90
#%.o : %.f90
	$(F90) $(F90CFLAGS) $(MP_FLAGS) -c $< -o $@  

%.o : $(src_dir)/%.f
#%.o : %.f
	$(F77) $(F77CFLAGS) $(MP_FLAGS) -c $< -o $@  

default: $(XSTICK_OBJECTS) $(OPENMP_OBJECTS) $(MPI_OBJECTS) 
	$(F90) -o ${EX}.x $(MP_LINK) $(XSTICK_OBJECTS) $(OPENMP_OBJECTS) $(MPI_OBJECTS) $(LFLAGS) $(LIBS)

serial :
	$(MAKE) F90=$(SERIAL90_COMPILER) F77=$(SERIAL77_COMPILER) LFLAGS=$(LFLAGS_ifort) OPENMP_OBJECTS=dummy_OpenMP_libs.o  EX=$(bin_dir)/bigstick

pgserial :
	cp $(src_dir)/mpif.h.template $(src_dir)/mpif.h
	$(MAKE) F90=$(PGI_COMPILER) F77=$(PGI_COMPILER) F90CFLAGS=$(PGCFLAGS) F77CFLAGS=$(PGCFLAGS)  OPENMP_OBJECTS=dummy_OpenMP_libs.o MPI_OBJECTS=dummy_MPI_libs.o EX=$(bin_dir)/bigstick
	rm $(src_dir)/mpif.h

gfortran :
#	cp $(src_dir)/mpif.h.template $(src_dir)/mpif.h
#	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) MP_FLAGS=-fopenmp MP_LINK=-fopenmp MPI_OBJECTS=dummy_MPI_libs.o EX=gbigstick
	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) F90CFLAGS=$(FGCFLAGS) LFLAGS=$(LFLAGS_all) OPENMP_OBJECTS=dummy_OpenMP_libs.o EX=$(bin_dir)/bigstick

sierra-serial :
	cp $(src_dir)/mpif.h.template $(src_dir)/mpif.h
	$(MAKE) F90=$(SERIAL90_COMPILER) F77=$(SERIAL77_COMPILER) LFLAGS=$(LFLAGS_ifort) OPENMP_OBJECTS=dummy_OpenMP_libs.o MPI_OBJECTS=dummy_MPI_libs.o EX=$(bin_dir)/bigstick
	rm $(src_dir)/mpif.h

openmp-legacy :
	$(MAKE) F90=$(SERIAL90_COMPILER) F77=$(SERIAL77_COMPILER) LFLAGS=$(LFLAGS_ifort) MPI_OBJECTS=dummy_MPI_libs.o \
MP_FLAGS=$(OPENMP_LEGACY_COMPILE) MP_LINK=$(OPENMP_LEGACY_LINK) EX=$(bin_dir)/bigstick-openmp

openmp :
	$(MAKE) F90=$(SERIAL90_COMPILER) F77=$(SERIAL77_COMPILER) LFLAGS=$(LFLAGS_ifort)  \
MP_FLAGS=$(OPENMP_COMPILE) MP_LINK=$(OPENMP_LINK) EX=$(bin_dir)/bigstick-openmp


pgopenmp :
	$(MAKE) F90=$(PGI_COMPILER) F77=$(PGI_COMPILER) MPI_OBJECTS=dummy_MPI_libs.o \
MP_FLAGS=$(OPENMP_PGI) MP_LINK=$(OPENMP_PGI) EX=$(bin_dir)/bigstick-openmp

gfortran-openmp :
#	cp $(src_dir)/mpif.h.template $(src_dir)/mpif.h
#	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) MP_FLAGS=-fopenmp MP_LINK=-fopenmp MPI_OBJECTS=dummy_MPI_libs.o EX=gbigstick
	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) LFLAGS=$(LFLAGS_all) F90CFLAGS=$(FGCFLAGS)   MP_FLAGS=-fopenmp MP_LINK=-fopenmp EX=$(bin_dir)/bigstick-openmp
#	rm $(src_dir)/mpif.h

gfortran-openmp-3body :
#	cp $(src_dir)/mpif.h.template $(src_dir)/mpif.h
#	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) MP_FLAGS=-fopenmp MP_LINK=-fopenmp MPI_OBJECTS=dummy_MPI_libs.o EX=gbigstick
	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) LFLAGS=$(LFLAGS_all) F90CFLAGS=$(FGCFLAGS3B)   MP_FLAGS=-fopenmp MP_LINK=-fopenmp EX=$(bin_dir)/bigstick-openmp3b

sierra-openmp :
	cp $(src_dir)/mpif.h.template $(src_dir)/mpif.h
	$(MAKE) F90=$(SERIAL90_COMPILER) LFLAGS=$(LFLAGS_ifort) F77=$(SERIAL77_COMPILER) MPI_OBJECTS=dummy_MPI_libs.o \
MP_FLAGS=$(OPENMP_COMPILE) MP_LINK=$(OPENMP_LINK) EX=$(bin_dir)/bigstick-openmp
	rm $(src_dir)/mpif.h

mpi :
	$(MAKE) F90=$(MPI90_COMPILER) F77=$(MPI77_COMPILER) F90CFLAGS=$(MPICFLAGS) LFLAGS=$(LFLAGS_all) OPENMP_OBJECTS=dummy_OpenMP_libs.o  EX=$(bin_dir)/bigstick-mpi

sierra-mpi :
	@echo ------------------------
	@echo  Before beginning be sure to load correct tools
	@echo /usr/global/tools/dotkit/init.sh
	@echo use mvapich2-intel-1.7
	@echo -------------------------
	$(MAKE) F90=$(MPI90_COMPILER_SIERRA) F77=$(MPI77_COMPILER_SIERRA) LFLAGS=$(LFLAGS_ifort) OPENMP_OBJECTS=dummy_OpenMP_libs.o  EX=$(bin_dir)/bigstick-mpi
	@echo  ----------------------------
	@echo  Note when submitting MPI jobs 
	@echo  be sure to include these lines: 
	@echo /usr/global/tools/dotkit/init.sh
	@echo use mvapich2-intel-1.7
	@echo export LD_LIBRARY_PATH=:/opt/intel-13.0/compiler/lib/intel64/:/usr/local/tools/mvapich2-intel-1.7/lib/
	@echo 	--------------------

openmp-mpi :
	$(MAKE) F90=$(MPI90_COMPILER) F77=$(MPI77_COMPILER) F90CFLAGS=$(MPICFLAGS) LFLAGS=$(LFLAGS_all) MP_FLAGS=$(OPENMP_COMPILE) \
MP_LINK=$(OPENMP_LINK) EX=$(bin_dir)/bigstick-mpi-omp

openmp-mpi-legacy :
	$(MAKE) F90=$(MPI90_COMPILER) F77=$(MPI77_COMPILER)LFLAGS=$(LFLAGS_ifort) MP_FLAGS=$(OPENMP_LEGACY_COMPILE) \
MP_LINK=$(OPENMP_LEGACY_LINK) EX=$(bin_dir)/bigstick-mpi-omp

grizzly :
	@echo ------------------------
	@echo  Before beginning be sure to load correct module
	@echo  module load intel intel-mpi 
	$(MAKE) F90=$(MPI90_COMPILER_GRIZZLY) F77=$(MPI77_COMPILER_GRIZZLY)LFLAGS=$(LFLAGS_ifort) MP_FLAGS=$(OPENMP_CORI_COMPILE) \
MP_LINK=$(OPENMP_CORI_COMPILE) EX=$(bin_dir)/bigstick-mpi-omp
	@echo  ----------------------------
	@echo  Note when submitting MPI jobs 
	@echo  be sure to include these lines: 
	@echo  module load intel intel-mpi 
	@echo  To help affinity
	@echo  srun --ntasks-per-node 4 -c 9 with proper numbers

sierra-openmp-mpi :
	$(MAKE) F90=$(MPI90_COMPILER_SIERRA) F77=$(MPI77_COMPILER_SIERRA) LFLAGS=$(LFLAGS_ifort) MP_FLAGS=$(OPENMP_COMPILE) \
MP_LINK=$(OPENMP_LINK) EX=$(bin_dir)/bigstick-mpi-omp

perlmutter-openmp-mpi :
	$(MAKE) F90=$(MPI90_COMPILER_PERLMUTTER) F77=$(MPI77_COMPILER_PERLMUTTER) F90CFLAGS=$(OPT_perlmutter) F77CFLAGS=$(OPT_perlmutter) MP_FLAGS=$(OPENMP_PERLMUTTER_COMPILE) \
MP_LINK=$(OPENMP_PERLMUTTER_LINK) EX= EX=$(bin_dir)/bigstick-mpi-omp

jaguar :
	$(MAKE) F90=$(JAGF90_COMPILER) F77=$(JAGF77_COMPILER) \
OPENMP_OBJECTS=dummy_OpenMP_libs.o EX=$(bin_dir)/bigstick-jaguar

vulcan-openmp-mpi :
	@echo ------------------------
	@echo  Before beginning be sure to load correct tools
	@echo  use bggcc-4.7.2
	$(MAKE) F90=$(MPI90_COMPILER_VULCAN) F77=$(MPI77_COMPILER_VULCAN) LFLAGS=$(LFLAGS_VULCAN) MP_FLAGS=$(OPENMP_VULCAN_COMPILE) \
MP_LINK=$(OPENMP_VULCAN_LINK) EX=$(bin_dir)/bigstick-mpi-openmp-vulcan

mira-openmp-mpi :
	$(MAKE) F90=$(MPI90_COMPILER_MIRA) F77=$(MPI77_COMPILER_MIRA) LFLAGS=$(LFLAGS_MIRA) MP_FLAGS=$(OPENMP_MIRA_COMPILE) \
MP_LINK=$(OPENMP_MIRA_LINK) EX=$(bin_dir)/bigstick-mpi-openmp

clean : 
	rm *.o *.mod $(bin_dir)/*bigstick*.x 

options  :
	@echo To see a list of options, make help

help  :
	@echo
	@echo Here are some compile options:
	@echo Note default compiler is intel ifort compiler
	@echo ----------------------------------------------
	@echo IMPORTANT if you change compiler option, e.g., serial to openmp or intel to gfortran
	@echo always make clean and then recompile
	@echo ---------------------------------------------
	@echo make clean -- deletes intermediate files, allows for fresh compile
	@echo make serial -- serial, bigstick.x
	@echo make openmp -- OpenMP parallelism,  bigstick-openmp.x
	@echo make openmp-legacy -- OpenMP parallelism on older Intel compiler,  bigstick-openmp.x
	@echo make mpi    -- MPI parallelism,  bigstick-mpi.x
	@echo make openmp-mpi -- hybrid OpenMP+MPI,  bigstick-mpi-omp.x
	@echo make openmp-mpi-legacy -- hybrid OpenMP+MPI on older Intel compiler,  bigstick-mpi-omp.x
	@echo -----------------------------------------------
	@echo make gfortran -- serial with GNU compiler gfortran, bigstick.x
	@echo make gfortran-openmp -- OpenMP with GNU compiler gfortran, bigstick-openmp.x
	@echo -----------------------------------------------
	@echo make gfortran-openmp-3body -- OpenMP with gfortran, with 3-body forces bigstick-openmp3b.x STILL IN BETA
	@echo -----------------------------------------------
	@echo make perlmutter-openmp-mpi  -- hybrid OpenMP+MPI on Perlmutter,  bigstick-mpi-omp.x
	@echo make help-legacy: list old options

help-legacy  :
	@echo make sierra -- serial with intel compiler on Sierra, bigstick.x
	@echo make sierra-mpi -- MPI with intel compiler on Sierra, bigstick-mpi.x
	@echo make sierra-openmp -- OpenMP with intel compiler on Sierra, bigstick-omp.x
	@echo make sierra-openmp-mpi  -- hybrid OpenMP+MPI on Sierra,  bigstick-mpi-omp.x
	@echo -----------------------------------------------
	@echo -----------------------------------------------
	@echo make grizzly -- hybrid OpenMP+MPI on grizzly,  bigstick-mpi-omp.x
	@echo  Before compiling AND running be sure to load correct module
	@echo  module load intel intel-mpi
	@echo  To help affinity when running set
	@echo  srun --ntasks-per-node 4 -c 9 with proper numbers
	
	
