################################################
#	Foundation Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J
PLATFORM?=LINUX


LIB = problemfile

ifeq ($(PLATFORM),MACOSX)
   DYNLIB_FLAG_NS=-fPIC -shared -install_name @rpath/libproblemfile_ns.dylib -o ./libproblemfile_ns.dylib
   DYNLIB_FLAG_NSSA=-fPIC -shared -install_name @rpath/libproblemfile_nssa.dylib -o ./libproblemfile_nssa.dylib
   DYNLIB_FLAG_NSTPT=-fPIC -shared -install_name @rpath/libproblemfile_nstpt.dylib -o ./libproblemfile_nstpt.dylib
   DYNLIB_FLAG_iNS=-fPIC -shared -install_name @rpath/libproblemfile_ins.dylib -o ./libproblemfile_ins.dylib
   DYNLIB_FLAG_CH=-fPIC -shared -install_name @rpath/libproblemfile_ch.dylib -o ./libproblemfile_ch.dylib
   DYNLIB_FLAG_MU=-fPIC -shared -install_name @rpath/libproblemfile_mu.dylib -o ./libproblemfile_mu.dylib
   DYNLIB_FLAG_CAA=-fPIC -shared -install_name @rpath/libproblemfile_caa.dylib -o ./libproblemfile_caa.dylib
   LIBSUFFIX=.dylib
else ifeq ($(PLATFORM),LINUX)
   DYNLIB_FLAG_NS=-fPIC -shared -o ./libproblemfile_ns.so
   DYNLIB_FLAG_NSSA=-fPIC -shared -o ./libproblemfile_nssa.so
   DYNLIB_FLAG_NSTPT=-fPIC -shared -o ./libproblemfile_nstpt.so
   DYNLIB_FLAG_iNS=-fPIC -shared -o ./libproblemfile_ins.so
   DYNLIB_FLAG_CH=-fPIC -shared -o ./libproblemfile_ch.so
   DYNLIB_FLAG_MU=-fPIC -shared -o ./libproblemfile_mu.so
   DYNLIB_FLAG_CAA=-fPIC -shared -o ./libproblemfile_caa.so
   LIBSUFFIX=.so
endif

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns$(LIBSUFFIX) lib$(LIB)_nssa$(LIBSUFFIX) lib$(LIB)_nstpt$(LIBSUFFIX) lib$(LIB)_ins$(LIBSUFFIX) lib$(LIB)_ch$(LIBSUFFIX) lib$(LIB)_mu$(LIBSUFFIX) lib$(LIB)_caa$(LIBSUFFIX)

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_nssa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include_nssa

./build_nstpt/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include_nstpt

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_mu/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include_mu

./build_ins/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_caa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include_caa

./build_ns/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns $(MOD_FLAG)./include_ns -c $< -o $@

./build_nssa/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES -DSPALARTALMARAS $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nssa $(MOD_FLAG)./include_nssa -c $< -o $@

./build_nstpt/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES -DTRANSPORT $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nstpt $(MOD_FLAG)./include_nstpt -c $< -o $@

./build_ins/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DINCNS $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins $(MOD_FLAG)./include_ins -c $< -o $@

./build_ch/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DCAHNHILLIARD $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch $(MOD_FLAG)./include_ch -c $< -o $@

./build_mu/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DMULTIPHASE -DCAHNHILLIARD $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu $(MOD_FLAG)./include_mu -c $< -o $@

./build_caa/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DACOUSTIC $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_caa $(MOD_FLAG)./include_caa -c $< -o $@

lib$(LIB)_ns$(LIBSUFFIX): header mkdirs ./build_ns/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_NS) $(FFLAGS) $(MACROS) $(INCLUDE)_ns ./build_ns/ProblemFile.o $(LIBS)

lib$(LIB)_nssa$(LIBSUFFIX): header mkdirs ./build_nssa/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_NSSA) $(FFLAGS) $(MACROS) $(INCLUDE)_nssa ./build_nssa/ProblemFile.o $(LIBS)

lib$(LIB)_nstpt$(LIBSUFFIX): header mkdirs ./build_nstpt/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_NSTPT) $(FFLAGS) $(MACROS) $(INCLUDE)_nstpt ./build_nstpt/ProblemFile.o $(LIBS)

lib$(LIB)_ins$(LIBSUFFIX): header mkdirs ./build_ins/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_iNS) $(FFLAGS) $(MACROS) $(INCLUDE)_ins ./build_ins/ProblemFile.o $(LIBS)

lib$(LIB)_ch$(LIBSUFFIX): header mkdirs ./build_ch/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_CH) $(FFLAGS) $(MACROS) $(INCLUDE)_ch ./build_ch/ProblemFile.o $(LIBS)

lib$(LIB)_mu$(LIBSUFFIX): header mkdirs ./build_mu/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_MU) $(FFLAGS) $(MACROS) $(INCLUDE)_mu ./build_mu/ProblemFile.o $(LIBS)

lib$(LIB)_caa$(LIBSUFFIX): header mkdirs ./build_caa/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_CAA) $(FFLAGS) $(MACROS) $(INCLUDE)_caa ./build_caa/ProblemFile.o $(LIBS)

install: all
	cp -p lib$(LIB)_ns$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_nssa$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_nstpt$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ins$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_mu$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_caa$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_nssa/*.mod $(INSTALL_DIR)/include_nssa/
	cp -p ./include_nstpt/*.mod $(INSTALL_DIR)/include_nstpt/
	cp -p ./include_ins/*.mod $(INSTALL_DIR)/include_ins/
	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	cp -p ./include_caa/*.mod $(INSTALL_DIR)/include_caa/
	@echo
	@echo

install-mu: lib$(LIB)_mu$(LIBSUFFIX)
	cp -p lib$(LIB)_mu$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	@echo

install-ns: lib$(LIB)_ns$(LIBSUFFIX)
	cp -p lib$(LIB)_ns$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	@echo

install-nssa: lib$(LIB)_nssa$(LIBSUFFIX)
	cp -p lib$(LIB)_nssa$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p ./include_nssa/*.mod $(INSTALL_DIR)/include_nssa/
	@echo

install-nstpt: lib$(LIB)_nstpt$(LIBSUFFIX)
	cp -p lib$(LIB)_nstpt$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p ./include_nstpt/*.mod $(INSTALL_DIR)/include_nstpt/
	@echo

install-caa: lib$(LIB)_caa$(LIBSUFFIX)
	cp -p lib$(LIB)_caa$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p ./include_caa/*.mod $(INSTALL_DIR)/include_caa/
	@echo

install-ins: lib$(LIB)_ins$(LIBSUFFIX)
	cp -p lib$(LIB)_ins$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p ./include_ins/*.mod $(INSTALL_DIR)/include_ins/
	@echo

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building Problem File      "
	@echo "================================"
	@echo
	@echo

mkdirs: FORCE
	@mkdir -p ./include_ns ./include_nssa ./include_nstpt ./include_ins ./include_ch ./include_mu ./include_caa
	@mkdir -p ./build_ns ./build_nssa ./build_nstpt ./build_ins ./build_ch ./build_mu ./build_caa

clean: FORCE
	$(RM) $(LIBNAME)
	$(RM) ./build_ns/*.o ./build_nssa/*.o ./build_nstpt/*.o ./build_ins/*.o ./build_ch/*.o ./build_mu/*.o ./build_caa/*.o ./include_ns/*.mod ./include_nssa/*.mod ./include_nstpt/*.mod ./include_ins/*.mod ./include_ch/*.mod ./include_mu/*.mod ./include_caa/*.mod ./*.so

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
