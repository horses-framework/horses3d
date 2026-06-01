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
PROG=./bin/horses3d.nstpt

OBJS_HORSES=	SpatialDiscretization \
				main  

OBJS_HORSES_BUILD=$(foreach obj,$(OBJS_HORSES),./build/$(obj).o)


##########################################
##	Include dependencies 		##
##########################################
-include ./make.deps

./build/main.o: FORCE


##########################################################
##		COMPILATION RULES 			##
##########################################################

.DEFAULT_GOAL :=  $(PROG) 

./build/%.o: ./%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -I./include $(INCLUDE) $(INCLUDE)_nstpt $(MACROS) $(EXTLIB_INC) -c $< -o $@ $(MOD_FLAG) ./include

./build_nstpt/%.o: ./**/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -I./include $(INCLUDE) $(INCLUDE)_nstpt $(MACROS) $(EXTLIB_INC) -c $< -o $@ $(MOD_FLAG) ./include_nstpt

$(PROG): header mkdirs $(OBJS_HORSES_BUILD)
	@echo
	@echo
	@echo "---------------------------------"
	@echo "| Linking " $@ "|"
	@echo "---------------------------------"
	$(FC) $(FFLAGS) $(INCLUDE) $(INCLUDE)_nstpt -o $(PROG) $(OBJS_HORSES_BUILD) $(HORSES_LIBS_EXEC) $(LIBS)

install: $(PROG)
	cp -p $(PROG) $(INSTALL_DIR)/bin
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
	@echo "======================================================="
	@echo ">>   Building HORSES3D Navier-Stokes Transport Solver  "
	@echo "======================================================="
	@echo
	@echo

## Dummy target to make sure a rule is executed
FORCE: