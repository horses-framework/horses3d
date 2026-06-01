#################################################
#	Spectral Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J

LIBOBJS_NS =    ./build_ns/TimeIntegratorDefinitions.o \
		./build_ns/AnalyticalJacobian.o \
		./build_ns/GenericLinSolverClass.o \
		./build_ns/PetscSolverClass.o \
		./build_ns/MatrixFreeSmootherClass.o \
		./build_ns/LinearMultigridSolverClass.o \
		./build_ns/MKLPardisoSolverClass.o \
		./build_ns/IterativeSolverClass.o \
		./build_ns/MatrixFreeGMRESClass.o \
		./build_ns/StaticCondensationSolverClass.o \
		./build_ns/LinearSolverClass.o \
		./build_ns/IMEXMethods.o \
		./build_ns/BDFTimeIntegrator.o \
		./build_ns/RosenbrockTimeIntegrator.o \
		./build_ns/MultigridTypes.o \
		./build_ns/TruncationErrorClass.o \
		./build_ns/ExplicitMethods.o \
		./build_ns/AnisFASMultigridClass.o \
		./build_ns/FASMultigridClass.o \
		./build_ns/MultiTauEstimationClass.o \
		./build_ns/AdaptiveTimeStepClass.o \
		./build_ns/pAdaptationClass.o \
		./build_ns/pAdaptationClassTE.o \
		./build_ns/pAdaptationClassRL.o \
		./build_ns/TimeIntegrator.o 

LIBOBJS_NSSA =    ./build_nssa/TimeIntegratorDefinitions.o \
		./build_nssa/AnalyticalJacobian.o \
		./build_nssa/GenericLinSolverClass.o \
		./build_nssa/PetscSolverClass.o \
		./build_nssa/MatrixFreeSmootherClass.o \
		./build_nssa/LinearMultigridSolverClass.o \
		./build_nssa/MKLPardisoSolverClass.o \
		./build_nssa/IterativeSolverClass.o \
		./build_nssa/MatrixFreeGMRESClass.o \
		./build_nssa/StaticCondensationSolverClass.o \
		./build_nssa/LinearSolverClass.o \
		./build_nssa/IMEXMethods.o \
		./build_nssa/BDFTimeIntegrator.o \
		./build_nssa/RosenbrockTimeIntegrator.o \
		./build_nssa/MultigridTypes.o \
		./build_nssa/TruncationErrorClass.o \
		./build_nssa/ExplicitMethods.o \
		./build_nssa/AnisFASMultigridClass.o \
		./build_nssa/FASMultigridClass.o \
		./build_nssa/MultiTauEstimationClass.o \
		./build_nssa/AdaptiveTimeStepClass.o \
		./build_nssa/pAdaptationClass.o \
		./build_nssa/pAdaptationClassTE.o \
		./build_nssa/pAdaptationClassRL.o \
		./build_nssa/TimeIntegrator.o 

LIBOBJS_NSTPT =    ./build_nstpt/TimeIntegratorDefinitions.o \
		./build_nstpt/AnalyticalJacobian.o \
		./build_nstpt/GenericLinSolverClass.o \
		./build_nstpt/PetscSolverClass.o \
		./build_nstpt/MatrixFreeSmootherClass.o \
		./build_nstpt/LinearMultigridSolverClass.o \
		./build_nstpt/MKLPardisoSolverClass.o \
		./build_nstpt/IterativeSolverClass.o \
		./build_nstpt/MatrixFreeGMRESClass.o \
		./build_nstpt/StaticCondensationSolverClass.o \
		./build_nstpt/LinearSolverClass.o \
		./build_nstpt/IMEXMethods.o \
		./build_nstpt/BDFTimeIntegrator.o \
		./build_nstpt/RosenbrockTimeIntegrator.o \
		./build_nstpt/MultigridTypes.o \
		./build_nstpt/TruncationErrorClass.o \
		./build_nstpt/ExplicitMethods.o \
		./build_nstpt/AnisFASMultigridClass.o \
		./build_nstpt/FASMultigridClass.o \
		./build_nstpt/MultiTauEstimationClass.o \
		./build_nstpt/AdaptiveTimeStepClass.o \
		./build_nstpt/pAdaptationClass.o \
		./build_nstpt/pAdaptationClassTE.o \
		./build_nstpt/pAdaptationClassRL.o \
		./build_nstpt/TimeIntegrator.o 

LIBOBJS_iNS =   ./build_ins/TimeIntegratorDefinitions.o \
		./build_ins/AnalyticalJacobian.o \
		./build_ins/GenericLinSolverClass.o \
		./build_ins/PetscSolverClass.o \
		./build_ins/MatrixFreeSmootherClass.o \
		./build_ins/LinearMultigridSolverClass.o \
		./build_ins/MKLPardisoSolverClass.o \
		./build_ins/IterativeSolverClass.o \
		./build_ins/MatrixFreeGMRESClass.o \
		./build_ins/StaticCondensationSolverClass.o \
		./build_ins/LinearSolverClass.o \
		./build_ins/IMEXMethods.o \
		./build_ins/BDFTimeIntegrator.o \
		./build_ins/RosenbrockTimeIntegrator.o \
		./build_ins/MultigridTypes.o \
		./build_ins/TruncationErrorClass.o \
		./build_ins/ExplicitMethods.o \
		./build_ins/AnisFASMultigridClass.o \
		./build_ins/FASMultigridClass.o \
		./build_ins/MultiTauEstimationClass.o \
		./build_ins/AdaptiveTimeStepClass.o \
		./build_ins/pAdaptationClass.o \
		./build_ins/pAdaptationClassTE.o \
		./build_ins/pAdaptationClassRL.o \
		./build_ins/TimeIntegrator.o 


LIBOBJS_CH =    ./build_ch/TimeIntegratorDefinitions.o \
		./build_ch/AnalyticalJacobian.o \
		./build_ch/GenericLinSolverClass.o \
		./build_ch/PetscSolverClass.o \
		./build_ch/MatrixFreeSmootherClass.o \
		./build_ch/LinearMultigridSolverClass.o \
		./build_ch/MKLPardisoSolverClass.o \
		./build_ch/IterativeSolverClass.o \
		./build_ch/MatrixFreeGMRESClass.o \
		./build_ch/StaticCondensationSolverClass.o \
		./build_ch/LinearSolverClass.o \
		./build_ch/IMEXMethods.o \
		./build_ch/BDFTimeIntegrator.o \
		./build_ch/RosenbrockTimeIntegrator.o \
		./build_ch/MultigridTypes.o \
		./build_ch/TruncationErrorClass.o \
		./build_ch/ExplicitMethods.o \
		./build_ch/AnisFASMultigridClass.o \
		./build_ch/FASMultigridClass.o \
		./build_ch/MultiTauEstimationClass.o \
		./build_ch/AdaptiveTimeStepClass.o \
		./build_ch/pAdaptationClass.o \
		./build_ch/pAdaptationClassTE.o \
		./build_ch/pAdaptationClassRL.o \
		./build_ch/TimeIntegrator.o 

LIBOBJS_MU =  ./build_mu/TimeIntegratorDefinitions.o \
		./build_mu/AnalyticalJacobian.o \
		./build_mu/GenericLinSolverClass.o \
		./build_mu/PetscSolverClass.o \
		./build_mu/MatrixFreeSmootherClass.o \
		./build_mu/LinearMultigridSolverClass.o \
		./build_mu/MKLPardisoSolverClass.o \
		./build_mu/IterativeSolverClass.o \
		./build_mu/MatrixFreeGMRESClass.o \
		./build_mu/StaticCondensationSolverClass.o \
		./build_mu/LinearSolverClass.o \
		./build_mu/IMEXMethods.o \
		./build_mu/BDFTimeIntegrator.o \
		./build_mu/RosenbrockTimeIntegrator.o \
		./build_mu/MultigridTypes.o \
		./build_mu/TruncationErrorClass.o \
		./build_mu/ExplicitMethods.o \
		./build_mu/AnisFASMultigridClass.o \
		./build_mu/FASMultigridClass.o \
		./build_mu/MultiTauEstimationClass.o \
		./build_mu/AdaptiveTimeStepClass.o \
		./build_mu/pAdaptationClass.o \
		./build_mu/pAdaptationClassTE.o \
		./build_mu/pAdaptationClassRL.o \
		./build_mu/TimeIntegrator.o 

LIBOBJS_CAA =  ./build_caa/TimeIntegratorDefinitions.o \
		./build_caa/AnalyticalJacobian.o \
		./build_caa/GenericLinSolverClass.o \
		./build_caa/PetscSolverClass.o \
		./build_caa/MatrixFreeSmootherClass.o \
		./build_caa/LinearMultigridSolverClass.o \
		./build_caa/MKLPardisoSolverClass.o \
		./build_caa/IterativeSolverClass.o \
		./build_caa/MatrixFreeGMRESClass.o \
		./build_caa/StaticCondensationSolverClass.o \
		./build_caa/LinearSolverClass.o \
		./build_caa/IMEXMethods.o \
		./build_caa/BDFTimeIntegrator.o \
		./build_caa/RosenbrockTimeIntegrator.o \
		./build_caa/MultigridTypes.o \
		./build_caa/TruncationErrorClass.o \
		./build_caa/ExplicitMethods.o \
		./build_caa/AnisFASMultigridClass.o \
		./build_caa/FASMultigridClass.o \
		./build_caa/MultiTauEstimationClass.o \
		./build_caa/AdaptiveTimeStepClass.o \
		./build_caa/pAdaptationClass.o \
		./build_caa/pAdaptationClassTE.o \
		./build_caa/pAdaptationClassRL.o \
		./build_caa/TimeIntegrator.o 

LIB = timeintegrator

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_nssa.a lib$(LIB)_nstpt.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a  lib$(LIB)_caa.a

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ns/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_nssa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DSPALARTALMARAS -I./include_nssa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nssa -c $< -o $@ $(MOD_FLAG) ./include_nssa

./build_nssa/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DSPALARTALMARAS -I./include_nssa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nssa -c $< -o $@ $(MOD_FLAG) ./include_nssa

./build_nstpt/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_nstpt $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nstpt -c $< -o $@ $(MOD_FLAG) ./include_nstpt

./build_nstpt/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_nstpt $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nstpt -c $< -o $@ $(MOD_FLAG) ./include_nstpt

./build_ins/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ins/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_ch/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_mu/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DMULTIPHASE -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $< -o $@ $(MOD_FLAG) ./include_mu

./build_mu/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DMULTIPHASE -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $<  -o $@ $(MOD_FLAG) ./include_mu

./build_caa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DACOUSTIC -I./include_caa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_caa -c $< -o $@ $(MOD_FLAG) ./include_caa

./build_caa/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DACOUSTIC -I./include_caa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_caa -c $<  -o $@ $(MOD_FLAG) ./include_caa

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
	@echo ">>   Building Time integrator   "
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
