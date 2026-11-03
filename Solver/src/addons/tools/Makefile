#####################################################################
## HORSES - a High-Order Spectral Element Solver
##
##	Main code Makefile
##
##
#####################################################################

########################################
# 1.- User defined parameters
########################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIN?=ranlib
MOD_FLAG?=-J
PLATFORM?=LINUX
COMPILER?=gfortran
MODE?=RELEASE
COMM?=PARALLEL
MAKE=make
AR=ar
RANLIB=ranlib
PROG=./bin/horses.tools

OBJS=	\
	FWHPostProc \
	SurfaceClass \
	FWHPreSurface \
	FWHTools \
	LocalRefinementClass \
	readSpecM \
	readHDF5 \
	readGMSH \
	LocalRefinementTool \
	LocalIBMRefinementTool \
	ConverStats \
	StatsMeshInterpolation \
	LambVectorMeshInterpolation \
	MeshInterpolationTool \
	main \

OBJS_BUILD=$(foreach obj,$(OBJS),./build/$(obj).o)


##########################################
##	Include dependencies 		##
##########################################
-include ./make.deps


##########################################################
##		COMPILATION RULES 			##
##########################################################

.DEFAULT_GOAL := $(PROG)

./build/%.o: ./%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(MACROS) $(EXTLIB_INC) -I./include $(INCLUDE) $(INCLUDE)_ns $(MACROS) -c $< -o $@ $(MOD_FLAG) ./include

$(PROG): header mkdirs $(OBJS_BUILD)
	@echo
	@echo
	@echo "---------------------------------"
	@echo "| Linking " $@ "|"
	@echo "---------------------------------"
	$(FC) $(FFLAGS) $(EXTLIB_INC) -o $(PROG) $(OBJS_BUILD) $(HORSES_LIBS_EXEC) $(LIBS)


install: $(PROG)
	cp -p $(PROG) $(INSTALL_DIR)/bin
	cp -p ./include/*.mod $(INSTALL_DIR)/include
	@echo
	@echo

##########################################
##	Extra auxiliary rules 		##
##########################################
clean: FORCE
	rm -f ./include/*.mod ./lib/*.a ./build/*.o
	rm -f $(PROG)

mkdirs: FORCE
	@mkdir -p ./build
	@mkdir -p ./include
	@mkdir -p ./bin
	@mkdir -p ./lib

header: FORCE
	@echo
	@echo "========================================"
	@echo ">>   Building horses.tools Solver  "
	@echo "========================================"
	@echo
	@echo

## Dummy target to make sure a rule is executed
FORCE: