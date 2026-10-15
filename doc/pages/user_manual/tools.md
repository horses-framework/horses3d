title: Tools
---
[TOC]

# Horses tools

<!-- Advanced users can have additional control over a simulation without having to modify the source code and recompile the code. To do that, the user can provide a set of routines that are called in different stages of the simulation via the Problem file (*ProblemFile.f90*). A description of the routines of the Problem File can be found in the following section. -->

It is possible to exclusively **compile** the tools running, from the `Solver` folder:

```
make addons
```

## Stats mesh interpolation

Interpolate stats from one mesh to another. There is no need that meshes have the same geometry nor node type. The routine is completely general.

The algorithm creates a probe for each node of the output mesh `M_out`. This probe is searched in the input mesh `M_in`. Then, a high-order interpolation is computed.

**This tool should be run in sequential mode or in parallel with OpenMP. It does not work with MPI.**

### Keywords

The keywords located in the control file regarding the **input** mesh (the mesh where the fields are defined) are:

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `tool type` | Required to run this tool | `"stats mesh interpolation"`
| `stats solver` | Specifies which solver generated the stats | ns/ins/mu
| `stats geo mesh file name` | Input mesh geometry | `"path/to/mesh.msh"`
| `stats mesh file name` | Input mesh data in `hmesh` format | `"path/to/mesh.hmesh"`
| `stats qbase file name` | Stats file computed in the input mesh | `"path/to/file.stats.hsol"`
| `stats sound velocity squared file name` | In the case of Navier-Stokes simulations, stats file of the sound velocity computed in the input mesh | `"path/to/file.SoundVelocitySquared.stats.hsol"`
| `stats gradient sound velocity squared file name` | In the case of Navier-Stokes simulations, stats file of the gradient of the sound velocity computed in the input mesh | `"path/to/file.GradientSoundVelocitySquared.stats.hsol"`
| `stats lamb vector file name` | Stats file of the Lamb vector computed in the input mesh | `"path/to/file.Lamb.stats.hsol"`
<!-- | `lamb vector file name` | File of the Lamb vector computed in the input mesh | `"path/to/file.Lamb.hsol"` -->


The keywords located in the control file regarding the **output** mesh (the mesh where the fields should be interpolated into) are:

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `mesh file name` | Output mesh geometry | `"path/to/mesh.msh"`
| `discretization nodes` | Type of nodes of the output mesh | Gauss/Gauss-Lobatto
| `Polynomial order` | Polynomial order of the output mesh | 1

### Execution

To execute the `stats mesh interpolation` tool, it is required to set the `tool type = "stats mesh interpolation"` in the control file, and then

```
./horses.tools controlFile.control
```

## Lamb vector mesh interpolation

Interpolate Lamb vector data from one mesh to another. There is no need that meshes have the same geometry nor node type. The routine is completely general.

The algorithm creates a probe for each node of the output mesh `M_out`. This probe is searched in the input mesh `M_in`. Then, a high-order interpolation is computed.

The tool looks for all the `Lamb.hsol` files present in the directory. Then, checks if there is a mesh for each of the files. If there is only one mesh file, all the data is interpolated using that file as input mesh. If there are more than one mesh files (if p-adaptation have been performed), then there should be the same number of mesh files than Lamb vector files. The same iteration is used for both file sets. 

**This tool should be run in sequential mode or in parallel with OpenMP. It does not work with MPI.**

### Keywords

The keywords located in the control file regarding the **input** mesh (the mesh where the fields are defined) are:

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `tool type` | Required to run this tool | `"Lamb vector mesh interpolation"`
| `stats geo mesh file name` | Input mesh geometry | `"path/to/mesh.msh"`
| `stats mesh file name` | Input mesh data in `hmesh` format | `"path/to/mesh.hmesh"`
| `lamb vector dir` | Path to the directory with all the Lamb vector data | `"path/to/"`
| `lamb vector file name` | Basename of the Lamb vector data computed in the input mesh | `"file"`


The keywords located in the control file regarding the **output** mesh (the mesh where the fields should be interpolated into) are:

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `mesh file name` | Output mesh geometry | `"path/to/mesh.msh"`
| `discretization nodes` | Type of nodes of the output mesh | Gauss/Gauss-Lobatto
| `Polynomial order` | Polynomial order of the output mesh | 1

### Execution

To execute the `stats mesh interpolation` tool, it is required to set the `tool type = "stats mesh interpolation"` in the control file, and then

```
./horses.tools controlFile.control
```

## Mesh interpolation

Interpolate the data of `.hsol` files from one mesh (input mesh) onto another mesh (output mesh). The two meshes can be completely different: other elements, other polynomial orders (also different in each element and direction) and other discretization nodes (Gauss / Gauss-Lobatto). The interpolation is high order: each node of the output mesh is located in an element of the input mesh and the data is evaluated with the Lagrange polynomials of that element.

This tool gathers and generalizes the capabilities of the `stats mesh interpolation` and `Lamb vector mesh interpolation` tools, and of the `meshInterpolation` task of `horsesConverter`:

* **Any kind of file.** The layout of each file is detected from the file itself, so the same tool interpolates:
    * Solution files, with or without gradients, sensor, LES viscosities or source terms.
    * Statistics files of the NS, iNS and MU solvers, with or without gradients.
    * Sound velocity squared, gradient of the sound velocity squared and Lamb vector statistics files.
    * Lamb vector files.
    * `.hmesh` files: interpolating the node coordinates of the input mesh gives the `.hmesh` of the output mesh.
* **Several files in one run.** The position of the output nodes inside the input mesh is computed once and reused for all the files that share the same input mesh.
* **Time-dependent input meshes.** Optionally, each file can use its own input mesh (e.g. after p-adaptation).
* **Robust point search.** The search starts at the last element found and uses a background grid of bins. Nodes slightly outside an input element are accepted within a tolerance. Nodes outside the input mesh are projected onto the closest point of the closest input element and a warning is printed.
* The iteration, time and reference values of each input file are kept in the interpolated file.

**This tool should be run in sequential mode or in parallel with OpenMP (`OMP_NUM_THREADS`). It does not work with MPI.**

### Keywords

The keywords located in the control file regarding the **input** mesh (the mesh where the fields are defined) are:

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `tool type` | Required to run this tool | `"mesh interpolation"`
| `input hmesh file name` | Input mesh in `hmesh` format | `"path/to/mesh.hmesh"`
| `input hmesh per file` | *(Optional)* If `.true.`, each file `name_ITER.hsol` (or `name_ITER.*.hsol`) uses the input mesh `path/to/mesh_ITER.hmesh`, where `path/to/mesh.hmesh` is the value of `input hmesh file name`. Default: `.false.` | `.true.`/`.false.`

The **output** mesh (the mesh where the fields should be interpolated into) can be given in two ways:

**Option 1: `hmesh` file (recommended, see [Large meshes](#large-meshes-generate-the-output-hmesh-with-the-solver)).** The polynomial orders and the discretization nodes are those of the `hmesh` file.

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `output hmesh file name` | Output mesh in `hmesh` format | `"path/to/mesh.hmesh"`

**Option 2: mesh file**, with the same keywords used by the solver. The mesh is constructed with the mesh readers of the solver, so curved meshes are supported in all the formats. In this case the control file must also contain the `#define boundary` blocks of all the boundaries of the mesh, e.g. copied from the control file of the simulation.

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `mesh file name` | Output mesh geometry | `"path/to/mesh.msh"`
| `Polynomial order` | Polynomial order of the output mesh. `polynomial order i/j/k` and `polynomial order file` can also be used, as in the solver | 4
| `discretization nodes` | Type of nodes of the output mesh | Gauss/Gauss-Lobatto
| `reference length (m)` | *(Optional)* Reference length used to make the mesh dimensionless. It must be the same as in the simulation. Default: 1 | 1.0

The keywords regarding the **files** to interpolate are (they can be combined, and at least one file must be given):

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `input file name` | File to interpolate | `"path/to/file.hsol"`
| `input file name 2`, `input file name 3`, ... | *(Optional)* More files to interpolate (numbered consecutively) | `"path/to/file.stats.hsol"`
| `input file pattern` | *(Optional)* Shell pattern of files to interpolate | `"path/to/case_*.Lamb.hsol"`

Other **optional** keywords are:

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `interpolation field` | Kind of files expected. The run stops if a file is of another kind. Default: `auto` (any file) | auto/solution/stats/lamb vector
| `output directory` | Directory of the interpolated files (created if it does not exist). Default: `RESULTS` | `"path/to/dir"`
| `output file suffix` | Text added to the name of the interpolated files. Default: `.Interpolated` | `".Interpolated"`
| `outside point tolerance` | Tolerance, in reference coordinates, to accept nodes slightly outside an input element. Default: 0.1 | 0.1
| `stats solver` | Solver that wrote the statistics files. Only needed when the layout of a statistics file is ambiguous (MU statistics without gradients); the tool asks for it in that case | ns/ins/mu
| `interpolation variable blocks` | *(Advanced)* Number of variables of each array stored per element, when the layout of a file cannot be detected automatically | `[9,5,1]`

For compatibility, if the keywords above are not given, the keywords of the `stats mesh interpolation` and `Lamb vector mesh interpolation` tools are also accepted: `stats mesh file name` (input `hmesh`), `stats qbase file name`, `stats sound velocity squared file name`, `stats gradient sound velocity squared file name`, `stats lamb vector file name`, `lamb vector dir` and `lamb vector file name`.

### Output files

Each file is written to `output directory` with the suffix added before the extension:

| Input file | Interpolated file |
|---------|-------------|
| `case_0000000100.hsol` | `case_0000000100.Interpolated.hsol` |
| `case.stats.hsol` | `case.Interpolated.stats.hsol` |
| `case_0000000100.Lamb.hsol` | `case_0000000100.Lamb.Interpolated.hsol` |
| `case_0000000000.hmesh` | `case_0000000000.Interpolated.hmesh` |

The interpolated files have the same type and layout as the input files, so they can be used directly, e.g. to restart a simulation on the new mesh (`restart file name`) or with `horses2plt`.

For each input mesh, the tool reports how many output nodes were found inside the input mesh, how many were accepted within the tolerance and how many were projected onto the closest element. Projected nodes usually mean that the two meshes do not represent the same geometry in that region (e.g. a curved wall represented with straight-sided elements in one of the meshes).

### Example

Interpolate a solution, the statistics and a series of Lamb vector files onto a new mesh given as an `hmesh` file:

```
tool type               = mesh interpolation

input hmesh file name   = MESH/old_0000000000.hmesh
output hmesh file name  = MESH/new_0000000000.hmesh

input file name         = RESULTS/case_0000010000.hsol
input file name 2       = RESULTS/case.stats.hsol
input file name 3       = RESULTS/case.SoundVelocitySquared.stats.hsol
input file pattern      = RESULTS/case_*.Lamb.hsol

output directory        = RESULTS_NEW
```

### Large meshes: generate the output hmesh with the solver

When the output mesh is given as a mesh file (`mesh file name`), the tool builds it with the full mesh readers of the solver. This step is sequential and, for large meshes, it takes most of the execution time. As a reference, for a mesh of 2184 curved elements with polynomial order 4, building the output mesh from the `.msh` file took about 8 to 9 s, whereas reading the same number of elements from an `hmesh` file took about 0.05 s.

For large meshes it is therefore a good option to generate first the `hmesh` of the output mesh and then give it to the tool with `output hmesh file name`. The simplest way is to run the solver on the new mesh with **zero time steps**: it reads the mesh (also in parallel with MPI), writes the `hmesh` and stops. Use a copy of the control file of the simulation with the output mesh, polynomial order and discretization nodes that you want, and:

```
mesh file name        = "MESH/new.msh"
Polynomial order      = 4
discretization nodes  = Gauss-Lobatto
solution file name    = "RESULTS/new.hsol"
restart               = .false.
restart file name     = "RESULTS/new.hsol"
Number of time steps  = 0
```

```
mpiexec -n 64 ./horses3d.ns new_mesh.control
```

The solver writes the file `MESH/new_0000000000.hmesh` (the name of the solution file followed by the iteration). Then, in the control file of the tool:

```
output hmesh file name  = MESH/new_0000000000.hmesh
```

Besides being much faster, this has two more advantages:

* The `#define boundary` blocks and the physical parameters are not needed in the control file of the tool.
* The elements of the interpolated files are numbered exactly as the solver numbers them, so the files can be used to restart the simulation on the new mesh without any further check.

Alternatively, the `hmesh` of the output mesh can also be obtained with this tool itself: run it once with the mesh file option and add the input `hmesh` file to the list of files to interpolate. The interpolated coordinates (`*.Interpolated.hmesh`) are the `hmesh` of the output mesh and can be used in the following runs.

### Execution

To execute the `mesh interpolation` tool, it is required to set the `tool type = "mesh interpolation"` in the control file, and then

```
OMP_NUM_THREADS=8 ./horses.tools controlFile.control
```

### Test

The test case `Solver/test/AddOns/MeshInterpolation` checks the tool: a polynomial solution (with gradients) and a polynomial Lamb vector field are written by the solver on a Cartesian box with Gauss nodes and polynomial order 3, interpolated with this tool onto a different box with Gauss-Lobatto nodes and polynomial order 4, and recovered exactly after a restart of the solver on the new mesh.
