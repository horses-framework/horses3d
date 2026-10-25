title: User Manual
ordered_subpage: input-and-output-files.md
ordered_subpage: running-a-simulation.md
ordered_subpage: spatial-discretization.md
ordered_subpage: physics-related-keyword.md
ordered_subpage: implicit-solvers-with-newton-linearization.md
ordered_subpage: explicit-solvers.md
ordered-subpage: a-nonlinear-p-multigrid-FAS.md
ordered-subpage: a-p-adaptation-methods.md
ordered-subpage: b-immersed-boundary-method.md
ordered-subpage: sliding-mesh.md
ordered-subpage: monitors.md
ordered-subpage: s-advanced-user-setup.md
ordered-subpage: s-postprocessing.md
ordered-subpage: t-appendix.md
---

# Details for compiling the code

- Clone the git repository or copy the source code into a desired folder.
- Go to the folder Solver.
- Run configure script.
```bash
    $ ./configure
```
- Install using the Makefile:
```bash
    $ make all <<Oprions>>
```
with the desired options (bold are default):

- MODE=DEBUG/HPC/**RELEASE**
- COMPILER=ifort/**gfortran**
- COMM=PARALLEL/**SEQUENTIAL**
- ENABLE_THREADS=NO/**YES**
- WITH_PETSC=YES/**NO**
- WITH_METIS=YES/**NO**
- WITH_HDF5=YES/**NO**
- WITH_MKL=YES/**NO**

For example:
```bash
$ make all COMPILER=ifort COMM=PARALLEL
```
- The **MODE=DEBUG/HPC/RELEASE** flag enables various compiler flags for different levels of code optimization. Furthermore, MODE=HPC disables residual file writing to improve performance. 

- The **ENABLE_THREADS=YES** flag enables shared memory simulations using OpenMP.

- The **COM=PARALLEL** flag enables distributed memory simulations using MPI.

- The **WITH_METIS=YES** flag activates METIS linking. To compile the code linking it with METIS (that is an option for creating the mesh partitions of MPI runs), it is needed that before compilation and running, an environment variable called METIS_DIR is found. This variable must contain the path to the METIS installation folder (it must have been compiled with the same compiler as HORSES3D).

- The **WITH_HDF5=YES** flag activates HDF5 linking. To compile the code linking it with HDF5 (necessary for reading HOPR meshes), it is needed that before compilation and running, an environment variable called HDF5_DIR is found. This variable must contain the path to the HDF5 installation folder (it must have been compiled with the same compiler as HORSES3D). In addition, the lib folder must be added to the environment variable LD_LIBRARY_PATH.

- The **WITH_PETSC=YES** flag activates PETSC linking. To compile the code linking it with PETSC, it is needed that before compilation and running, an environment variable called PETSC_DIR is found. This variable must contain the path to the PETSC installation folder.

- If you use *environment modules*, it is advised to use the HORSES3D module file:

```bash
$ export MODULEPATH=$HORSES_DIR/utils/modulefile:$MODULEPATH
```
where `\$HORSES\_DIR` is the installation directory.

- It is advised to run the `make clean` or `make allclean` command if some options of the compilation routine need to be changed and it has been compiled before.



