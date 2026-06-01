#################################################
#	Mesh Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp -ffree-line-length-0
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J

LIBOBJS_NS = 	./build_ns/ConnectivityClass.o \
		./build_ns/MeshTypes.o \
		./build_ns/WallDistance.o \
		./build_ns/StorageClass.o \
		./build_ns/FacePatchClass.o \
		./build_ns/TransfiniteMaps3D.o \
		./build_ns/MappedGeometry.o \
		./build_ns/NodeClass.o \
		./build_ns/HexElementConnectivityDefinitions.o \
		./build_ns/FaceClass.o \
		./build_ns/HexElementClass.o \
		./build_ns/ZoneClass.o \
		./build_ns/TessellationTypes.o \
		./build_ns/OrientedBoundingBox.o \
		./build_ns/MPI_IBMUtilities.o \
		./build_ns/KDClass.o \
		./build_ns/IBMClass.o \
		./build_ns/HexMesh.o \
		./build_ns/MeshPartitioning.o \
		./build_ns/SurfaceMesh.o \
		./build_ns/Read_HDF5Mesh_HOPR.o \
		./build_ns/Read_SpecMesh.o \
		./build_ns/Read_GMSH.o \
		./build_ns/ReadMeshFile.o

LIBOBJS_NSSA = 	./build_nssa/ConnectivityClass.o \
		./build_nssa/MeshTypes.o \
		./build_nssa/WallDistance.o \
		./build_nssa/StorageClass.o \
		./build_nssa/FacePatchClass.o \
		./build_nssa/TransfiniteMaps3D.o \
		./build_nssa/MappedGeometry.o \
		./build_nssa/NodeClass.o \
		./build_nssa/HexElementConnectivityDefinitions.o \
		./build_nssa/FaceClass.o \
		./build_nssa/HexElementClass.o \
		./build_nssa/ZoneClass.o \
		./build_nssa/TessellationTypes.o \
		./build_nssa/OrientedBoundingBox.o \
		./build_nssa/MPI_IBMUtilities.o \
		./build_nssa/KDClass.o \
		./build_nssa/IBMClass.o \
		./build_nssa/HexMesh.o \
		./build_nssa/MeshPartitioning.o \
		./build_nssa/SurfaceMesh.o \
		./build_nssa/Read_HDF5Mesh_HOPR.o \
		./build_nssa/Read_SpecMesh.o \
		./build_nssa/Read_GMSH.o \
		./build_nssa/ReadMeshFile.o

LIBOBJS_NSTPT = 	./build_nstpt/ConnectivityClass.o \
		./build_nstpt/MeshTypes.o \
		./build_nstpt/WallDistance.o \
		./build_nstpt/StorageClass.o \
		./build_nstpt/FacePatchClass.o \
		./build_nstpt/TransfiniteMaps3D.o \
		./build_nstpt/MappedGeometry.o \
		./build_nstpt/NodeClass.o \
		./build_nstpt/HexElementConnectivityDefinitions.o \
		./build_nstpt/FaceClass.o \
		./build_nstpt/HexElementClass.o \
		./build_nstpt/ZoneClass.o \
		./build_nstpt/TessellationTypes.o \
		./build_nstpt/OrientedBoundingBox.o \
		./build_nstpt/MPI_IBMUtilities.o \
		./build_nstpt/KDClass.o \
		./build_nstpt/IBMClass.o \
		./build_nstpt/HexMesh.o \
		./build_nstpt/MeshPartitioning.o \
		./build_nstpt/SurfaceMesh.o \
		./build_nstpt/Read_HDF5Mesh_HOPR.o \
		./build_nstpt/Read_SpecMesh.o \
		./build_nstpt/Read_GMSH.o \
		./build_nstpt/ReadMeshFile.o

LIBOBJS_iNS = 	./build_ins/ConnectivityClass.o \
		./build_ins/MeshTypes.o \
		./build_ins/WallDistance.o \
		./build_ins/StorageClass.o \
		./build_ins/FacePatchClass.o \
		./build_ins/TransfiniteMaps3D.o \
		./build_ins/MappedGeometry.o \
		./build_ins/NodeClass.o \
		./build_ins/HexElementConnectivityDefinitions.o \
		./build_ins/FaceClass.o \
		./build_ins/HexElementClass.o \
		./build_ins/ZoneClass.o \
		./build_ins/TessellationTypes.o \
		./build_ins/OrientedBoundingBox.o \
		./build_ins/MPI_IBMUtilities.o \
		./build_ins/KDClass.o \
		./build_ins/IBMClass.o \
		./build_ins/HexMesh.o \
		./build_ins/MeshPartitioning.o \
		./build_ins/SurfaceMesh.o \
		./build_ins/Read_HDF5Mesh_HOPR.o \
		./build_ins/Read_SpecMesh.o \
		./build_ins/Read_GMSH.o \
		./build_ins/ReadMeshFile.o


LIBOBJS_CH = 	./build_ch/ConnectivityClass.o \
		./build_ch/MeshTypes.o \
		./build_ch/WallDistance.o \
		./build_ch/StorageClass.o \
		./build_ch/FacePatchClass.o \
		./build_ch/TransfiniteMaps3D.o \
		./build_ch/MappedGeometry.o \
		./build_ch/NodeClass.o \
		./build_ch/HexElementConnectivityDefinitions.o \
		./build_ch/FaceClass.o \
		./build_ch/HexElementClass.o \
		./build_ch/ZoneClass.o \
		./build_ch/TessellationTypes.o \
		./build_ch/OrientedBoundingBox.o \
		./build_ch/MPI_IBMUtilities.o \
		./build_ch/KDClass.o \
		./build_ch/IBMClass.o \
		./build_ch/HexMesh.o \
		./build_ch/MeshPartitioning.o \
		./build_ch/SurfaceMesh.o \
		./build_ch/Read_HDF5Mesh_HOPR.o \
		./build_ch/Read_SpecMesh.o \
		./build_ch/Read_GMSH.o \
		./build_ch/ReadMeshFile.o

LIBOBJS_MU =  ./build_mu/ConnectivityClass.o \
		./build_mu/MeshTypes.o \
		./build_mu/WallDistance.o \
		./build_mu/StorageClass.o \
		./build_mu/FacePatchClass.o \
		./build_mu/TransfiniteMaps3D.o \
		./build_mu/MappedGeometry.o \
		./build_mu/NodeClass.o \
		./build_mu/HexElementConnectivityDefinitions.o \
		./build_mu/FaceClass.o \
		./build_mu/HexElementClass.o \
		./build_mu/ZoneClass.o \
		./build_mu/TessellationTypes.o \
		./build_mu/OrientedBoundingBox.o \
		./build_mu/MPI_IBMUtilities.o \
		./build_mu/KDClass.o \
		./build_mu/IBMClass.o \
		./build_mu/HexMesh.o \
		./build_mu/MeshPartitioning.o \
		./build_mu/SurfaceMesh.o \
		./build_mu/Read_HDF5Mesh_HOPR.o \
		./build_mu/Read_SpecMesh.o \
		./build_mu/Read_GMSH.o \
		./build_mu/ReadMeshFile.o

LIBOBJS_CAA =  ./build_caa/ConnectivityClass.o \
		./build_caa/MeshTypes.o \
		./build_caa/WallDistance.o \
		./build_caa/StorageClass.o \
		./build_caa/FacePatchClass.o \
		./build_caa/TransfiniteMaps3D.o \
		./build_caa/MappedGeometry.o \
		./build_caa/NodeClass.o \
		./build_caa/HexElementConnectivityDefinitions.o \
		./build_caa/FaceClass.o \
		./build_caa/HexElementClass.o \
		./build_caa/ZoneClass.o \
		./build_caa/TessellationTypes.o \
		./build_caa/OrientedBoundingBox.o \
		./build_caa/MPI_IBMUtilities.o \
		./build_caa/KDClass.o \
		./build_caa/IBMClass.o \
		./build_caa/HexMesh.o \
		./build_caa/MeshPartitioning.o \
		./build_caa/SurfaceMesh.o \
		./build_caa/Read_HDF5Mesh_HOPR.o \
		./build_caa/Read_SpecMesh.o \
		./build_caa/Read_GMSH.o \
		./build_caa/ReadMeshFile.o

LIB = mesh

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_nstpt.a lib$(LIB)_nssa.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a  lib$(LIB)_caa.a

./build_ns/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_nssa/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DNAVIERSTOKES -DSPALARTALMARAS -I./include_nssa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nssa -c $< -o $@ $(MOD_FLAG) ./include_nssa

./build_nstpt/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DNAVIERSTOKES -I./include_nstpt $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nstpt -c $< -o $@ $(MOD_FLAG) ./include_nstpt

./build_ins/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ch/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_mu/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DNAVIERSTOKES -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $< -o $@ $(MOD_FLAG) ./include_mu

./build_caa/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DACOUSTIC -I./include_caa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_caa -c $< -o $@ $(MOD_FLAG) ./include_caa

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

lib$(LIB)_ns.a: header mkdirs $(LIBOBJS_NS) ./build_ns/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NS) ./build_ns/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_nssa.a: header mkdirs $(LIBOBJS_NSSA) ./build_nssa/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSSA) ./build_nssa/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_nstpt.a: header mkdirs $(LIBOBJS_NSTPT) ./build_nstpt/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSTPT) ./build_nstpt/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_ins.a: header mkdirs $(LIBOBJS_iNS) ./build_ins/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_iNS) ./build_ins/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_ch.a: header mkdirs $(LIBOBJS_CH) ./build_ch/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_CH) ./build_ch/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_mu.a: header mkdirs $(LIBOBJS_MU) ./build_mu/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_MU) ./build_mu/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_caa.a: header mkdirs $(LIBOBJS_CAA) ./build_caa/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_CAA) ./build_caa/METISPartitioning.o
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
	@echo ">>   Building Mesh        "
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
