#####################################################################
## HORSES - a High-Order Spectral Element Solver
##
##      Test case makefile. This Makefile compiles a shared library,
##   libproblemfile.so with user specific procedures, such as initial
##   and boundary conditions or solution checks.
##
##	Usage:
##		-> cd $(CASE_DIR)/SETUP
##		-> make COMPILER=gfortran/ifort MODE=RELEASE/DEBUG COMM=SERIAL/PARALLEL
##		-> cd $(CASE_DIR)
##		-> ldd ./HORSES3D (to check whether the solver points
##				to the default library (found in
##				$HOME_DIR) or to the local library,
##				./SETUP/libproblemfile.so)
##
#####################################################################
HOME_DIR=@HORSES_HOME_DIR
########################################
# 1.- User defined parameters
########################################
PLATFORM?=LINUX
COMPILER?=gfortran
MODE?=RELEASE
COMM?=SEQUENTIAL
ENABLE_THREADS?=YES
MAKE=make

########################################
# 2.- User defined paths
########################################
INCLUDE=-I$(HOME_DIR)/include
LIBS =
MACROS =
-include $(HOME_DIR)/make.inc
-include make.inc

########################################
# 3.- User defined compilation flags 
########################################
FFLAGS=-fPIC
GNU_HPC_FLAGS= -cpp -fallow-argument-mismatch -std=legacy -ffree-line-length-0 -O3 -ftree-vectorize \
		   -ftree-vectorizer-verbose=0 -fbackslash -D_has_Quad -march=znver3 -funroll-loops
GNU_RELEASE_FLAGS= -cpp -ffree-line-length-0 -O3 -ftree-vectorize  -ftree-vectorizer-verbose=0 \
		   -fbackslash -D_has_Quad
GNU_DEBUG_FLAGS=  -ffree-line-length-0 -cpp -O0 -g -fcheck=all -Wno-unused-variable \
		  -fbacktrace -fbounds-check -finit-real=snan -Wall \
		  -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries \
		  -fbackslash -fcheck=pointer -D_has_Quad
INTEL_HPC_FLAGS= -cpp -Ofast -msse4.2 -axAVX,CORE-AVX2 -assume bscc -D_has_Quad  -fp-model=fast=2 -ipo
INTEL_RELEASE_FLAGS= -cpp -O3 -assume bscc -D_has_Quad -fp-model precise
INTEL_DEBUG_FLAGS= -cpp -O0 -g -warn all -fpscomp logicals -traceback \
		   -check bounds -check uninit -ftrapuv -debug all \
		   -ftrace=full -ftz -fp-model precise \
		   -fp-speculation=off -assume bscc -D_has_Quad
ARFLAGS=crv

TARGET=./libproblemfile

##########################################
## Select platform
##########################################
ifeq ($(PLATFORM),MACOSX)
   LIBSUFFIX=.dylib
   TARGET_NS=$(TARGET)_ns$(LIBSUFFIX)
   TARGET_NSSA=$(TARGET)_nssa$(LIBSUFFIX)
   TARGET_NSTPT=$(TARGET)_nstpt$(LIBSUFFIX)
   TARGET_iNS=$(TARGET)_ins$(LIBSUFFIX)
   TARGET_CH=$(TARGET)_ch$(LIBSUFFIX)
   TARGET_MU=$(TARGET)_mu$(LIBSUFFIX)
   TARGET_CAA=$(TARGET)_caa$(LIBSUFFIX)
   DYNLIB_FLAG_NS=-fPIC -shared -install_name @rpath/$(TARGET_NS) -o $(TARGET_NS)#
   DYNLIB_FLAG_NSSA=-fPIC -shared -install_name @rpath/$(TARGET_NSSA) -o $(TARGET_NSSA)
   DYNLIB_FLAG_NSTPT=-fPIC -shared -install_name @rpath/$(TARGET_NSTPT) -o $(TARGET_NSTPT)#
   DYNLIB_FLAG_iNS=-fPIC -shared -install_name @rpath/$(TARGET_iNS) -o $(TARGET_iNS)
   DYNLIB_FLAG_CH=-fPIC -shared -install_name @rpath/$(TARGET_CH) -o $(TARGET_CH)
   DYNLIB_FLAG_MU=-fPIC -shared -install_name @rpath/$(TARGET_MU) -o $(TARGET_MU)
   DYNLIB_FLAG_CAA=-fPIC -shared -install_name @rpath/$(TARGET_CAA) -o $(TARGET_CAA)
else ifeq ($(PLATFORM),LINUX)
   LIBSUFFIX=.so
   TARGET_NS=$(TARGET)_ns$(LIBSUFFIX)
   TARGET_NSSA=$(TARGET)_nssa$(LIBSUFFIX)
   TARGET_NSTPT=$(TARGET)_nstpt$(LIBSUFFIX)
   TARGET_iNS=$(TARGET)_ins$(LIBSUFFIX)
   TARGET_CH=$(TARGET)_ch$(LIBSUFFIX)
   TARGET_MU=$(TARGET)_mu$(LIBSUFFIX)
   TARGET_CAA=$(TARGET)_caa$(LIBSUFFIX)
   DYNLIB_FLAG_NS=-fPIC -shared -o $(TARGET_NS)
   DYNLIB_FLAG_NSSA=-fPIC -shared -o $(TARGET_NSSA)
   DYNLIB_FLAG_iNS=-fPIC -shared -o $(TARGET_iNS)
   DYNLIB_FLAG_CH=-fPIC -shared -o $(TARGET_CH)
   DYNLIB_FLAG_MU=-fPIC -shared -o $(TARGET_MU)
   DYNLIB_FLAG_CAA=-fPIC -shared -o $(TARGET_CAA)
endif


##########################################
## Select compiler 
##########################################
ifeq ($(COMPILER),gfortran)
  FC:=gfortran
  MOD_FLAG:=-J
  WHOLE_FILE:=--whole-file
else ifeq ($(COMPILER),gnu)
  FC:=gfortran
  MOD_FLAG:=-J
  WHOLE_FILE:=--whole-file
else ifeq ($(COMPILER),GNU)
  FC:=gfortran
  MOD_FLAG:=-J
  WHOLE_FILE:=--whole-file
endif

ifeq ($(COMPILER),ifort)
  FC:=ifort
  MOD_FLAG:=-module
  WHOLE_FILE:=
else ifeq ($(COMPILER),intel)
  FC:=ifort
  MOD_FLAG:=-module
  WHOLE_FILE:=
else ifeq ($(COMPILER),INTEL)
  FC:=ifort
  MOD_FLAG:=-module
  WHOLE_FILE:=
endif

##########################################
## Select mode
##########################################
ifeq ($(MODE),RELEASE)
     ifeq ($(FC),gfortran)
         FFLAGS+=$(GNU_RELEASE_FLAGS)
     else ifeq ($(FC),ifort)
         FFLAGS+=$(INTEL_RELEASE_FLAGS)
     endif
else ifeq ($(MODE),HPC)
	 MACROS+= -D_HPC_MODE
     ifeq ($(FC),gfortran)
         FFLAGS+=$(GNU_HPC_FLAGS)
     else ifeq ($(FC),ifort)
         FFLAGS+=$(INTEL_HPC_FLAGS)
     endif	 
else ifeq ($(MODE),DEBUG)
     ifeq ($(FC),gfortran)
         FFLAGS+=$(GNU_DEBUG_FLAGS)
     else ifeq ($(FC),ifort)
         FFLAGS+=$(INTEL_DEBUG_FLAGS)
     endif
endif

ifeq ($(COMM),PARALLEL)
    ifeq ($(FC),gfortran)
        FC=mpif90
	MACROS+= -D_HAS_MPI_
    else ifeq ($(FC),ifort)
        FC=mpiifort
	MACROS+= -D_HAS_MPI_
    endif
endif

ifeq ($(ENABLE_THREADS),YES)
    FFLAGS+=-fopenmp
endif

##############################################################
##	Get objects 
##############################################################
HORSES_LIBS_NS=-ltimeintegrator_ns -ljacobian_ns -ldiscretization_ns -lmonitors_ns -lparticles_ns -lmesh_ns -lphysics_ns -lphysicsns -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_NSSA=-ltimeintegrator_nssa -ljacobian_nssa -ldiscretization_nssa -lmonitors_nssa -lparticles_nssa -lmesh_nssa -lphysics_nssa -lphysicsnssa -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_NSTPT=-ltimeintegrator_nstpt -ljacobian_nstpt -ldiscretization_nstpt -lmonitors_nstpt -lparticles_nstpt -lmesh_nstpt -lphysics_nstpt -lphysicsnstpt -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_iNS=-ltimeintegrator_ins -ljacobian_ins -ldiscretization_ins -lmonitors_ins -lparticles_ins -lmesh_ins -lphysics_ins -lphysicsns -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_CH=-ltimeintegrator_ch -ljacobian_ch -ldiscretization_ch -lmonitors_ch -lparticles_ch -lmesh_ch -lphysics_ch -lphysicsch -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_MU=-ltimeintegrator_mu -ljacobian_mu -ldiscretization_mu -lmonitors_mu -lparticles_mu -lmesh_mu -lphysics_mu -lphysicsch -lphysicsmu -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_CAA=-ltimeintegrator_caa -ljacobian_caa -ldiscretization_caa -lmonitors_caa -lparticles_caa -lmesh_caa -lphysics_caa -lphysicscaa -lio -lspectral -lmpiutils -lfoundation -lftobject

.DEFAULT_GOAL := all

all: $(TARGET_NS) $(TARGET_NSSA) $(TARGET_NSTPT) $(TARGET_iNS) $(TARGET_CH) $(TARGET_MU) $(TARGET_CAA)

$(TARGET_NS): FORCE ./ProblemFile_NS.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_NS) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ns ./ProblemFile_NS.o -L$(HOME_DIR)/lib $(HORSES_LIBS_NS)

$(TARGET_NSSA): FORCE ./ProblemFile_NSSA.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_NSSA) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_nssa ./ProblemFile_NSSA.o -L$(HOME_DIR)/lib $(HORSES_LIBS_NSSA)

$(TARGET_NSTPT): FORCE ./ProblemFile_NSTPT.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_NSTPT) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_nstpt ./ProblemFile_NSTPT.o -L$(HOME_DIR)/lib $(HORSES_LIBS_NSTPT)

$(TARGET_iNS): FORCE ./ProblemFile_iNS.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_iNS) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ins ./ProblemFile_iNS.o -L$(HOME_DIR)/lib $(HORSES_LIBS_iNS)

$(TARGET_CH): FORCE ./ProblemFile_CH.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_CH) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ch ./ProblemFile_CH.o -L$(HOME_DIR)/lib $(HORSES_LIBS_CH)
	
$(TARGET_MU): FORCE ./ProblemFile_MU.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_MU) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_mu ./ProblemFile_MU.o -L$(HOME_DIR)/lib $(HORSES_LIBS_MU)

$(TARGET_CAA): FORCE ./ProblemFile_CAA.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_CAA) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_caa ./ProblemFile_CAA.o -L$(HOME_DIR)/lib $(HORSES_LIBS_CAA)


# ---------- Executable macro -------------------
./ProblemFile_NS.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ns -c ./ProblemFile.f90 -o $@

./ProblemFile_NSSA.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES -DSPALARTALMARAS $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_nssa -c ./ProblemFile.f90 -o $@

./ProblemFile_NSTPT.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES -DTRANSPORT $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_nstpt -c ./ProblemFile.f90 -o $@

./ProblemFile_iNS.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DINCNS $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ins -c ./ProblemFile.f90 -o $@

./ProblemFile_CH.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DCAHNHILLIARD $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ch -c ./ProblemFile.f90 -o $@

./ProblemFile_MU.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DMULTIPHASE -DCAHNHILLIARD $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_mu -c ./ProblemFile.f90 -o $@

./ProblemFile_CAA.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DACOUSTIC $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_caa -c ./ProblemFile.f90 -o $@


clean: FORCE
	rm -f ./ProblemFile_NS.o ./ProblemFile_NSSA.o ./ProblemFile_NSTPT.o ./ProblemFile_iNS.o ./ProblemFile_CH.o ./ProblemFile_MU.o ./ProblemFile_CAA.o ./*.mod
	rm -f ./*__genmod.f90
	
allclean: FORCE clean
	rm -f ./libproblemfile.so
	rm -f ./libproblemfile.dylib
	rm -f ./libproblemfile.dSYM
	rm -f ./libproblemfile_ns.so
	rm -f ./libproblemfile_ns.dylib
	rm -f ./libproblemfile_ns.dSYM
	rm -f ./libproblemfile_nssa.so
	rm -f ./libproblemfile_nssa.dylib
	rm -f ./libproblemfile_nssa.dSYM
    rm -f ./libproblemfile_nstpt.so
	rm -f ./libproblemfile_nstpt.dylib
	rm -f ./libproblemfile_nstpt.dSYM
	rm -f ./libproblemfile_ins.so
	rm -f ./libproblemfile_ins.dylib
	rm -f ./libproblemfile_ins.dSYM
	rm -f ./libproblemfile_ch.so
	rm -f ./libproblemfile_ch.dylib
	rm -f ./libproblemfile_ch.dSYM
	rm -f ./libproblemfile_mu.so
	rm -f ./libproblemfile_mu.dylib
	rm -f ./libproblemfile_mu.dSYM
	rm -f ./libproblemfile_caa.so
	rm -f ./libproblemfile_caa.dylib
	rm -f ./libproblemfile_caa.dSYM

FORCE:
