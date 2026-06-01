#################################################
#	Foundation Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J

LIBOBJS_NS = 	./build_ns/FluidData.o \
		./build_ns/PhysicsStorage.o \
		./build_ns/VariableConversion.o \
		./build_ns/Physics.o \
		./build_ns/GenericBoundaryConditionClass.o \
		./build_ns/InflowBC.o \
		./build_ns/OutflowBC.o \
		./build_ns/NoSlipWallBC.o \
		./build_ns/FreeSlipWallBC.o \
		./build_ns/PeriodicBC.o \
		./build_ns/UserDefinedBC.o \
		./build_ns/WallFunctionDefinitions.o \
		./build_ns/WallFunctionBC.o \
		./build_ns/BoundaryConditions.o \
		./build_ns/LESModels.o

LIBOBJS_NSSA = 	./build_nssa/FluidData.o \
		./build_nssa/PhysicsStorage.o \
		./build_nssa/VariableConversion.o \
		./build_nssa/Physics.o \
		./build_nssa/GenericBoundaryConditionClass.o \
		./build_nssa/InflowBC.o \
		./build_nssa/OutflowBC.o \
		./build_nssa/NoSlipWallBC.o \
		./build_nssa/FreeSlipWallBC.o \
		./build_nssa/PeriodicBC.o \
		./build_nssa/UserDefinedBC.o \
		./build_nssa/WallFunctionDefinitions.o \
		./build_nssa/WallFunctionBC.o \
		./build_nssa/BoundaryConditions.o 

LIBOBJS_NSTPT = 	./build_nstpt/FluidData.o \
		./build_nstpt/PhysicsStorage.o \
		./build_nstpt/VariableConversion.o \
		./build_nstpt/Physics.o \
		./build_nstpt/GenericBoundaryConditionClass.o \
		./build_nstpt/InflowBC.o \
		./build_nstpt/OutflowBC.o \
		./build_nstpt/NoSlipWallBC.o \
		./build_nstpt/FreeSlipWallBC.o \
		./build_nstpt/PeriodicBC.o \
		./build_nstpt/UserDefinedBC.o \
		./build_nstpt/WallFunctionDefinitions.o \
		./build_nstpt/WallFunctionBC.o \
		./build_nstpt/BoundaryConditions.o \
		./build_nstpt/LESModels.o

LIBOBJS_iNS = 	./build_ins/FluidData.o \
		./build_ins/PhysicsStorage.o \
		./build_ins/VariableConversion.o \
		./build_ins/Physics.o \
		./build_ins/GenericBoundaryConditionClass.o \
		./build_ins/InflowBC.o \
		./build_ins/OutflowBC.o \
		./build_ins/NoSlipWallBC.o \
		./build_ins/FreeSlipWallBC.o \
		./build_ins/PeriodicBC.o \
		./build_ins/UserDefinedBC.o \
		./build_ins/BoundaryConditions.o \
		./build_ins/LESModels.o

LIBOBJS_CH = 	./build_ch/FluidData.o \
		./build_ch/PhysicsStorage.o \
		./build_ch/VariableConversion.o \
		./build_ch/Physics.o \
		./build_ch/GenericBoundaryConditionClass.o \
		./build_ch/InflowBC.o \
		./build_ch/OutflowBC.o \
		./build_ch/NoSlipWallBC.o \
		./build_ch/FreeSlipWallBC.o \
		./build_ch/PeriodicBC.o \
		./build_ch/UserDefinedBC.o \
		./build_ch/BoundaryConditions.o 

LIBOBJS_MU = 	./build_mu/FluidData.o \
		./build_mu/PhysicsStorage.o \
		./build_mu/VariableConversion.o \
		./build_mu/Physics.o \
		./build_mu/GenericBoundaryConditionClass.o \
		./build_mu/InflowBC.o \
		./build_mu/OutflowBC.o \
		./build_mu/NoSlipWallBC.o \
		./build_mu/FreeSlipWallBC.o \
		./build_mu/PeriodicBC.o \
		./build_mu/UserDefinedBC.o \
		./build_mu/BoundaryConditions.o \
		./build_mu/LESModels.o

LIBOBJS_CAA = 	./build_caa/FluidData.o \
		./build_caa/PhysicsStorage.o \
		./build_caa/VariableConversion.o \
		./build_caa/Physics.o \
		./build_caa/GenericBoundaryConditionClass.o \
		./build_caa/InflowBC.o \
		./build_caa/OutflowBC.o \
		./build_caa/NoSlipWallBC.o \
		./build_caa/FreeSlipWallBC.o \
		./build_caa/PeriodicBC.o \
		./build_caa/UserDefinedBC.o \
		./build_caa/BoundaryConditions.o 

LIB = physics

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL := lib$(LIB).a

all: lib$(LIB)_ns.a lib$(LIB)_nssa.a lib$(LIB)_nstpt.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a lib$(LIB)_caa.a

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_nssa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DSPALARTALMARAS -I./include_nssa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nssa -c $< -o $@ $(MOD_FLAG) ./include_nssa

./build_nstpt/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DTRANSPORT -I./include_nstpt $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nstpt -c $< -o $@ $(MOD_FLAG) ./include_nstpt

./build_ins/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_mu/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DMULTIPHASE -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $< -o $@ $(MOD_FLAG) ./include_mu

./build_caa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DACOUSTIC -I./include_caa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_caa -c $< -o $@ $(MOD_FLAG) ./include_caa

lib$(LIB)_ns.a: header mkdirs $(LIBOBJS_NS)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NS)
	@$(RANLIB) $@

lib$(LIB)_nssa.a: header mkdirs $(LIBOBJS_NSSA)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSSA)
	@$(RANLIB) $@

lib$(LIB)_nstpt.a: header mkdirs $(LIBOBJS_NSTPT)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSTPT)
	@$(RANLIB) $@

lib$(LIB)_ins.a: header mkdirs $(LIBOBJS_iNS)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_iNS)
	@$(RANLIB) $@

lib$(LIB)_ch.a: header mkdirs $(LIBOBJS_CH)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_CH)
	@$(RANLIB) $@

lib$(LIB)_mu.a: header mkdirs $(LIBOBJS_MU)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_MU)
	@$(RANLIB) $@

lib$(LIB)_caa.a: header mkdirs $(LIBOBJS_CAA)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_CAA) 
	@$(RANLIB) $@

install: all
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_nssa.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_nstpt.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ins.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_caa.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_nssa/*.mod $(INSTALL_DIR)/include_nssa/
	cp -p ./include_nstpt/*.mod $(INSTALL_DIR)/include_nstpt/
	cp -p ./include_ins/*.mod $(INSTALL_DIR)/include_ins/
	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	cp -p ./include_caa/*.mod $(INSTALL_DIR)/include_caa/
	@echo
	@echo

install-mu: lib$(LIB)_mu.a
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	@echo
	@echo

install-ns: lib$(LIB)_ns.a
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	@echo
	@echo

install-nssa: lib$(LIB)_nssa.a
	cp -p lib$(LIB)_nssa.a $(INSTALL_DIR)/lib/
	cp -p ./include_nssa/*.mod $(INSTALL_DIR)/include_nssa/
	@echo
	@echo

install-nstpt: lib$(LIB)_nstpt.a
	cp -p lib$(LIB)_nstpt.a $(INSTALL_DIR)/lib/
	cp -p ./include_nstpt/*.mod $(INSTALL_DIR)/include_nstpt/
	@echo
	@echo

install-caa: lib$(LIB)_caa.a
	cp -p lib$(LIB)_caa.a $(INSTALL_DIR)/lib/
	cp -p ./include_caa/*.mod $(INSTALL_DIR)/include_caa/
	@echo
	@echo

install-ins: lib$(LIB)_ins.a
	cp -p lib$(LIB)_ins.a $(INSTALL_DIR)/lib/
	cp -p ./include_ins/*.mod $(INSTALL_DIR)/include_ins/
	@echo
	@echo

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building common Physics    "
	@echo "================================"
	@echo
	@echo

mkdirs: FORCE
	@mkdir -p ./include_ns ./include_nssa ./include_nstpt ./include_ins ./include_ch ./include_mu ./include_caa
	@mkdir -p ./build_ns ./build_nssa ./build_nstpt ./build_ins ./build_ch ./build_mu ./build_caa

clean: FORCE
	$(RM) lib$(LIB)_ns.a lib$(LIB)_nssa.a lib$(LIB)_nstpt.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a lib$(LIB)_caa.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod
	$(RM) ./build_nssa/*.o ./include_nssa/*.mod
	$(RM) ./build_nstpt/*.o ./include_nstpt/*.mod
	$(RM) ./build_ins/*.o ./include_ins/*.mod
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_mu/*.o ./include_mu/*.mod
	$(RM) ./build_caa/*.o ./include_caa/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
