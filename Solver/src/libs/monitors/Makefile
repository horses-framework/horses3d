#################################################
#	Monitors Makefile 		#
#################################################
INSTALL_DIR?=./
INCLUDE?=-I../../../include
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J


LIBOBJS_NS = 	./build_ns/MonitorDefinitions.o \
		./build_ns/Probe.o \
		./build_ns/StatisticsMonitor.o \
		./build_ns/SurfaceIntegrals.o \
		./build_ns/VolumeIntegrals.o \
		./build_ns/ResidualsMonitor.o \
		./build_ns/SurfaceMonitor.o \
		./build_ns/VolumeMonitor.o \
		./build_ns/LoadBalancingMonitor.o \
		./build_ns/Monitors.o \
		./build_ns/SamplingOperator.o \
		./build_ns/SurfaceSampling.o \
		./build_ns/PlaneSampling.o \
		./build_ns/SpatialMeanNode.o \
		./build_ns/Samplings.o \
		./build_ns/FWHDefinitions.o \
		./build_ns/FWHObseverClass.o \
		./build_ns/FWHGeneralClass.o

LIBOBJS_NSSA = 	./build_nssa/MonitorDefinitions.o \
		./build_nssa/Probe.o \
		./build_nssa/StatisticsMonitor.o \
		./build_nssa/SurfaceIntegrals.o \
		./build_nssa/VolumeIntegrals.o \
		./build_nssa/ResidualsMonitor.o \
		./build_nssa/SurfaceMonitor.o \
		./build_nssa/VolumeMonitor.o \
		./build_nssa/LoadBalancingMonitor.o \
		./build_nssa/Monitors.o \
		./build_nssa/SamplingOperator.o \
		./build_nssa/SurfaceSampling.o \
		./build_nssa/PlaneSampling.o \
		./build_nssa/SpatialMeanNode.o \
		./build_nssa/Samplings.o 

LIBOBJS_NSTPT = 	./build_nstpt/MonitorDefinitions.o \
		./build_nstpt/Probe.o \
		./build_nstpt/StatisticsMonitor.o \
		./build_nstpt/SurfaceIntegrals.o \
		./build_nstpt/VolumeIntegrals.o \
		./build_nstpt/ResidualsMonitor.o \
		./build_nstpt/SurfaceMonitor.o \
		./build_nstpt/VolumeMonitor.o \
		./build_nstpt/LoadBalancingMonitor.o \
		./build_nstpt/Monitors.o \
		./build_nstpt/SamplingOperator.o \
		./build_nstpt/SurfaceSampling.o \
		./build_nstpt/PlaneSampling.o \
		./build_nstpt/SpatialMeanNode.o \
		./build_nstpt/Samplings.o \
		./build_nstpt/FWHDefinitions.o \
		./build_nstpt/FWHObseverClass.o \
		./build_nstpt/FWHGeneralClass.o

LIBOBJS_iNS = 	./build_ins/MonitorDefinitions.o \
		./build_ins/Probe.o \
		./build_ins/StatisticsMonitor.o \
		./build_ins/SurfaceIntegrals.o \
		./build_ins/VolumeIntegrals.o \
		./build_ins/ResidualsMonitor.o \
		./build_ins/SurfaceMonitor.o \
		./build_ins/VolumeMonitor.o \
		./build_ins/LoadBalancingMonitor.o \
		./build_ins/Monitors.o \
		./build_ins/SamplingOperator.o \
		./build_ins/SurfaceSampling.o \
		./build_ins/PlaneSampling.o \
		./build_ins/SpatialMeanNode.o \
		./build_ins/Samplings.o 

LIBOBJS_CH = 	./build_ch/MonitorDefinitions.o \
		./build_ch/Probe.o\
		./build_ch/StatisticsMonitor.o \
		./build_ch/SurfaceIntegrals.o \
		./build_ch/VolumeIntegrals.o \
		./build_ch/ResidualsMonitor.o \
		./build_ch/SurfaceMonitor.o \
		./build_ch/VolumeMonitor.o \
		./build_ch/LoadBalancingMonitor.o \
		./build_ch/Monitors.o \
		./build_ch/SamplingOperator.o \
		./build_ch/SurfaceSampling.o \
		./build_ch/PlaneSampling.o \
		./build_ch/SpatialMeanNode.o \
		./build_ch/Samplings.o 

LIBOBJS_MU = 	./build_mu/MonitorDefinitions.o \
		./build_mu/Probe.o\
		./build_mu/StatisticsMonitor.o \
		./build_mu/SurfaceIntegrals.o \
		./build_mu/VolumeIntegrals.o \
		./build_mu/ResidualsMonitor.o \
		./build_mu/SurfaceMonitor.o \
		./build_mu/VolumeMonitor.o \
		./build_mu/LoadBalancingMonitor.o \
		./build_mu/Monitors.o \
		./build_mu/SamplingOperator.o \
		./build_mu/SurfaceSampling.o \
		./build_mu/PlaneSampling.o \
		./build_mu/SpatialMeanNode.o \
		./build_mu/Samplings.o 

LIBOBJS_CAA = 	./build_caa/MonitorDefinitions.o \
		./build_caa/Probe.o\
		./build_caa/StatisticsMonitor.o \
		./build_caa/SurfaceIntegrals.o \
		./build_caa/VolumeIntegrals.o \
		./build_caa/ResidualsMonitor.o \
		./build_caa/SurfaceMonitor.o \
		./build_caa/VolumeMonitor.o \
		./build_caa/LoadBalancingMonitor.o \
		./build_caa/Monitors.o \
		./build_caa/SamplingOperator.o \
		./build_caa/SurfaceSampling.o \
		./build_caa/PlaneSampling.o \
		./build_caa/SpatialMeanNode.o \
		./build_caa/Samplings.o 

LIB = monitors

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_nssa.a lib$(LIB)_nstpt.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a  lib$(LIB)_caa.a

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDEO) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_nssa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DSPALARTALMARAS -I./include_nssa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDEO) $(INCLUDE)_nssa -c $< -o $@ $(MOD_FLAG) ./include_nssa

./build_nstpt/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_nstpt $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDEO) $(INCLUDE)_nstpt -c $< -o $@ $(MOD_FLAG) ./include_nstpt

./build_ins/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDEO) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDEO) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_mu/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DMULTIPHASE -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDEO) $(INCLUDE)_mu -c $< -o $@ $(MOD_FLAG) ./include_mu

./build_caa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DACOUSTIC -I./include_caa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDEO) $(INCLUDE)_caa -c $< -o $@ $(MOD_FLAG) ./include_caa

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
	@echo ">>   Building Monitors    "
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
