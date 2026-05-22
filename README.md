# HORSES3D High-Order (DG) Spectral Element Solver

MIT License

Copyright (c) 2021 NUMATH https://numath.dmae.upm.es


## Synopsis

**HORSES3D** is a multiphysics environment where the compressible Navier-Stokes equations, the incompressible Navier–Stokes equations, the Cahn–Hilliard equation and entropy–stable variants are solved. Arbitrary high–order, p–anisotropic discretisations are used, including static and dynamic p–adaptation methods (feature-based and truncation error-based). Explicit and implicit time-steppers for steady and time-marching solutions are available, including efficient multigrid and preconditioners. Numerical and analytical Jacobian computations with a coloring algorithm have been implemented. Multiphase flows are solved using a diffuse interface model: Navier–Stokes/Cahn–Hilliard. Turbulent models implemented include RANS: Spalart-Allmaras and LES: Smagorinsky, Wale, Vreman; including wall models. Immersed boundary methods can be used, to avoid creating body fitted meshes. Acoustic propagation can be computed using Ffowcs-Williams and Hawkings models.

HORSES3D supports curvilinear, hexahedral, conforming meshes in GMSH, HDF5 and SpecMesh/HOHQMesh format. A hybrid CPU-based parallelisation strategy (shared and distributed memory) with OpenMP and MPI is followed.


## External libraries

The following external routines/libraries can be used with **HORSES3D**, but are not necessary: METIS, MPI, HDF5, MKL, PETSc.


## Compilers and third-party software

**HORSES3D** is an object-oriented Fortran 2008 solver, that can be compiled using gcc and the Intel compiler, in Unix-based operating systems.

- We recommend using recent versions of such compilers (2019 or newer).

- Make is necessary (e.g., Gnu's version, which is included in most linux distributions).

- Supported meshes are in GMSH, HDF5 (HOPR) and SpecMesh/HOHQMesh format.

- Post processing can be performed in tecplot or paraview.


## Compiling & Running

1. Go to the Solver folder and configure the project

    ```shell
    cd Solver
    ./configure
    ```

2. Build the solvers using `make` (see manual in the `/doc` folder)

    ```shell
    make clean
    make all [options]
    ```

    with the desired *options* (defaults are bold):

    - PLATFORM=MACOSX/**LINUX**

    - MODE=DEBUG/**RELEASE**

    - COMPILER=ifort/**gfortran**

    - COMM=PARALLEL/**SEQUENTIAL**

    - ENABLE_THREADS=NO/**YES**

    - WITH_PETSC=YES/**NO**

    - WITH_METIS=YES/**NO**

    - WITH_HDF5=YES/**NO**

    - WITH_MKL=YES/**NO**

    For example:

    ```shell
    make all COMPILER=ifort COMM=PARALLEL
    ```

3. Run the solver for the parameter file *file.control* (see manual in the `/doc` folder and examples in `/test`)

    ```shell
    ./horses3d.ns file.control
    ```

4. Test cases for various physics are provided in the folder `/test`


## Additional libraries

- PETSc:

    ```shell
    export PETSC_DIR=path_to_petsc
    ```

- METIS:

    ```shell
    export METIS_DIR=path_to_metis
    ```

- HDF5:

    ```shell
    export HDF5_DIR=path_to_hdf5
    ```

## Documentation

Read the complete documentation [here](https://loganoz.github.io/horses3d/).

## Acknowledgments
The authors acknowledge the funding received by the Grant DeepCFD (Project No. PID2022-137899OB-I00) funded by MICIU/AEI/10.13039/501100011033 and by ERDF, EU and also the
funding from the European Union (ERC, Off-coustics, project number 101086075) and (ROSAS, project number 101138319). Views and opinions expressed are, however, those of the authors only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.