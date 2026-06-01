#################################################
#	Foundation Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J




LIBOBJS = ./build/FluidData_NS.o \
		  ./build/PhysicsStorage_NS.o \
		  ./build/VariableConversion_NS.o \
		  ./build/ArtificialViscosity.o \
		  ./build/Physics_NS.o \
		  ./build/RiemannSolvers_NS.o \
		  ./build/ManufacturedSolutions.o 

LIB = physicsnstpt

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:= lib$(LIB).a


./build/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(MACROS) $(EXTLIB_INC) -I./include $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include

lib$(LIB).a: header mkdirs $(LIBOBJS)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS)
	@$(RANLIB) $@

install : lib$(LIB).a
	cp -p lib$(LIB).a $(INSTALL_DIR)/lib/
	cp -p ./include/*.mod $(INSTALL_DIR)/include/
	@echo
	@echo

header: FORCE
	@echo
	@echo "=============================================="
	@echo ">>   Building Navier-Stokes Transport physics "
	@echo "=============================================="
	@echo
	@echo
mkdirs: FORCE
	@mkdir -p ./include
	@mkdir -p ./build

clean: FORCE
	$(RM) lib$(LIB).a 
	$(RM) ./build/*.o ./include/*.mod
#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE: