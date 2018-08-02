############################################################################
#                                                                          #
#                     Make file for BIGSTICK                               #
#                                                                          # 
#   Be sure to switch the COMPILER Variable to your compiler.              #
#                                                                          #
#   For MPI execution follow these steps:                                  #
#      (1) Remove dummy_MPI_libs.o from the list of objects.               #
#          This is a set of dummy subrotuines to compile BIGSTICK for      #
#          serial execution.                                               #
#      (2) Move, rename, or delete the mpif.h header file.                 #
#          This is a placeholder for MPI declrations, such as MPI_REAL,    #
#          etc.                                                            #
#                                                                          #
#   Failure to follow these steps will not implement MPI protocals as      #
#   these files will override the standard build.                          #
#                                                                          #
#   The current version has the ability to select the precision for        #
#   Lanczos vectors. This is done in the module Lanczos_precision          #
#   located in bmodules.f90. BIGSTICK is delivered with single-precision   #
#   as the default. To change to double precision, change the comments     #
#   in the module.                                                         #
#                                                                          #
############################################################################
F77CFLAGS   =  -c -O2  #-openmp # -CB
F90CFLAGS   =  -c -O2  #-openmp #  -CB
LFLAGS      = # -openmp 
F77COMPILER =   gfortran #  ifort # mpif77
F90COMPILER =   gfortran # mpif90
MPIF77COMPILER = mpif77
MPIF90COMPILER = mpif90
PROTROP         =  tropic1b
PROGENSTRENGTH     = genstrength
PROGTSTRENGTH     = gtstrength

TR1B_OBJECTS =          \
trans1bmodules.o        \
trans1b_main.o            \
trans1binput.o             \
trans1b_subs.o           \
libra.o

GEN_STRENGTH_OBJECTS = \
gentrans.v10.o  \
libra.o

GT_STRENGTH_OBJECTS = \
GTstrength.v9.o  \
libra.o

${PROTROP}.x : $(TR1B_OBJECTS)
	$(F90COMPILER) -o ${PROTROP}.x $(TR1B_OBJECTS) $(LFLAGS)

${PROGENSTRENGTH}.x : $(GEN_STRENGTH_OBJECTS)
	$(F90COMPILER) -o ${PROGENSTRENGTH}.x $(GEN_STRENGTH_OBJECTS) $(LFLAGS)


${PROGTSTRENGTH}.x : $(GT_STRENGTH_OBJECTS)
	$(F90COMPILER) -o ${PROGTSTRENGTH}.x $(GT_STRENGTH_OBJECTS) $(LFLAGS)

%.o : %.f90
	$(F90COMPILER) $(F90CFLAGS) $(<F)  

%.o : %.f
	$(F77COMPILER) $(F77CFLAGS) $(<F)  
clean : 
	rm *.o *.mod *.x 
