#################################################
#	Jacobian Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J




LIBOBJS_NS = ./build_ns/JacobianDefinitions.o \
	        ./build_ns/ColorsClass.o \
			    ./build_ns/mkl_spblas.o \
					./build_ns/GenericMatrixClass.o \
			    ./build_ns/LinkedListMatrixClass.o \
			    ./build_ns/CSRMatrixClass.o \
			    ./build_ns/PETScMatrixClass.o \
			    ./build_ns/DenseBlockDiagonalMatrixClass.o \
			    ./build_ns/SparseBlockDiagonalMatrixClass.o \
			    ./build_ns/StaticCondensedMatrixClass.o \
			    ./build_ns/MatrixFreeClass.o \
			    ./build_ns/MatrixClass.o \
				  ./build_ns/JacobianComputerClass.o \
				./build_ns/NumericalJacobian.o

LIBOBJS_NSSA = ./build_nssa/JacobianDefinitions.o \
	            ./build_nssa/ColorsClass.o \
			    ./build_nssa/mkl_spblas.o \
			    ./build_nssa/GenericMatrixClass.o \
			    ./build_nssa/LinkedListMatrixClass.o \
			    ./build_nssa/CSRMatrixClass.o \
			    ./build_nssa/PETScMatrixClass.o \
			    ./build_nssa/DenseBlockDiagonalMatrixClass.o \
			    ./build_nssa/SparseBlockDiagonalMatrixClass.o \
			    ./build_nssa/StaticCondensedMatrixClass.o \
			    ./build_nssa/MatrixFreeClass.o \
			    ./build_nssa/MatrixClass.o \
				./build_nssa/JacobianComputerClass.o \
				./build_nssa/NumericalJacobian.o

LIBOBJS_NSTPT = ./build_nstpt/JacobianDefinitions.o \
	        ./build_nstpt/ColorsClass.o \
			    ./build_nstpt/mkl_spblas.o \
					./build_nstpt/GenericMatrixClass.o \
			    ./build_nstpt/LinkedListMatrixClass.o \
			    ./build_nstpt/CSRMatrixClass.o \
			    ./build_nstpt/PETScMatrixClass.o \
			    ./build_nstpt/DenseBlockDiagonalMatrixClass.o \
			    ./build_nstpt/SparseBlockDiagonalMatrixClass.o \
			    ./build_nstpt/StaticCondensedMatrixClass.o \
			    ./build_nstpt/MatrixFreeClass.o \
			    ./build_nstpt/MatrixClass.o \
				  ./build_nstpt/JacobianComputerClass.o \
				./build_nstpt/NumericalJacobian.o

LIBOBJS_iNS = ./build_ins/JacobianDefinitions.o \
         	./build_ins/ColorsClass.o \
			    ./build_ins/mkl_spblas.o \
					./build_ins/GenericMatrixClass.o \
			    ./build_ins/LinkedListMatrixClass.o \
			    ./build_ins/CSRMatrixClass.o \
			    ./build_ins/PETScMatrixClass.o \
			    ./build_ins/DenseBlockDiagonalMatrixClass.o \
			    ./build_ins/SparseBlockDiagonalMatrixClass.o \
			    ./build_ins/StaticCondensedMatrixClass.o \
			    ./build_ins/MatrixFreeClass.o \
			    ./build_ins/MatrixClass.o \
				  ./build_ins/JacobianComputerClass.o \
				./build_ins/NumericalJacobian.o



LIBOBJS_CH = ./build_ch/JacobianDefinitions.o \
        	./build_ch/ColorsClass.o \
			    ./build_ch/mkl_spblas.o \
				  ./build_ch/GenericMatrixClass.o \
			    ./build_ch/LinkedListMatrixClass.o \
			    ./build_ch/CSRMatrixClass.o \
			    ./build_ch/PETScMatrixClass.o \
			    ./build_ch/DenseBlockDiagonalMatrixClass.o \
			    ./build_ch/SparseBlockDiagonalMatrixClass.o \
			    ./build_ch/StaticCondensedMatrixClass.o \
			    ./build_ch/MatrixFreeClass.o \
			    ./build_ch/MatrixClass.o \
				  ./build_ch/JacobianComputerClass.o \
				./build_ch/NumericalJacobian.o

LIBOBJS_MU = ./build_mu/JacobianDefinitions.o \
        	./build_mu/ColorsClass.o \
			    ./build_mu/mkl_spblas.o \
				./build_mu/GenericMatrixClass.o \
			    ./build_mu/LinkedListMatrixClass.o \
			    ./build_mu/CSRMatrixClass.o \
			    ./build_mu/PETScMatrixClass.o \
			    ./build_mu/DenseBlockDiagonalMatrixClass.o \
			    ./build_mu/SparseBlockDiagonalMatrixClass.o \
			    ./build_mu/StaticCondensedMatrixClass.o \
			    ./build_mu/MatrixFreeClass.o \
			    ./build_mu/MatrixClass.o \
				./build_mu/JacobianComputerClass.o \
				./build_mu/NumericalJacobian.o

LIBOBJS_CAA = ./build_caa/JacobianDefinitions.o \
        	./build_caa/ColorsClass.o \
			    ./build_caa/mkl_spblas.o \
				./build_caa/GenericMatrixClass.o \
			    ./build_caa/LinkedListMatrixClass.o \
			    ./build_caa/CSRMatrixClass.o \
			    ./build_caa/PETScMatrixClass.o \
			    ./build_caa/DenseBlockDiagonalMatrixClass.o \
			    ./build_caa/SparseBlockDiagonalMatrixClass.o \
			    ./build_caa/StaticCondensedMatrixClass.o \
			    ./build_caa/MatrixFreeClass.o \
			    ./build_caa/MatrixClass.o \
				./build_caa/JacobianComputerClass.o \
				./build_caa/NumericalJacobian.o

LIB = jacobian

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_nssa.a lib$(LIB)_nstpt.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a  lib$(LIB)_caa.a

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
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_nstpt $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nstpt -c $< -o $@ $(MOD_FLAG) ./include_nstpt

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
	@echo ">>   Building Jacobian    "
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
