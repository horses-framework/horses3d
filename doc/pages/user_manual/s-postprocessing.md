title: Postprocessing 
---

[TOC]

For postprocessing the Simulation Results

## Visualization with Tecplot Format: *horses2plt*


`HORSES3D` provides a script for converting the native binary solution files (*.hsol) into tecplot ASCII format (*.tec), which can be visualized in Pareview or Tecplot. It can also export the solution to the more recent VTKHDF format; however, note that this feature does **not** export boundary information or mesh files. Usage:

```bash
	$ horses2plt SolutionFile.hsol MeshFile.hmesh <<Options>>
```

The options comprise following flags:

| Flag                | Description                                                                                                      | Default value |
|---------------------|------------------------------------------------------------------------------------------------------------------|---------------|
| --output-order=    | *INTEGER*: Output order nodes. The solution is interpolated into the desired number of points.                  | Not Present   |
| --output-basis=    | *CHARACTER*: Either *Homogeneous* (for equispaced nodes, or *Gauss*.)                                            | *Gauss*       |
| --output-mode=     | *CHARACTER*: Either *multizone* or *FE*. The option *multizone* generates a Tecplot zone for each element. The option *FE* generates only one Tecplot zone for the fluid and one for each boundary (if *--boundary-file* is defined). Each subcell is mapped as a linear finite element. This format is faster to read by Paraview and Tecplot. | *multizone*   |
| --output-variables=| *CHARACTER*: Output variables separated by commas. A complete description can be found in Section 2.             | Q             |
| --flow-eqs=| *CHARACTER*: Solver generating the solution file. Options are: ns, ins, ch, mu, caa.             | ns             |
| --dimensionless    | Specifies that the output quantities must be dimensionless.                                                     | Not Present   |
| --partition-file=  | *CHARACTER*: Specifies the path to the partition file (*.pmesh) to export the MPI ranks of the simulation.     | Not Present   |
| --mlrk-level-file= | *CHARACTER*: Specifies the path to the MLRK level file (*.lmesh), written once by multi-level RK runs, to export the level of each element (`mlrk_level`) and the level including neighbours (`mlrk_level_wn`). | Not Present   |
| --boundary-file=   | *CHARACTER*: Specifies the path to the boundary mesh file (*.bmesh) to export the surfaces as additional zones of the Tecplot file. | Not Present   |
| --output-type=     | *CHARACTER*: Specifies the type of output file: *tecplot* or *vtkhdf*.                                          | *tecplot*     |

* *Homogeneous* when *--output-order* is specified

Additionally, depending on the type of solution file, the user can specify additional options.

## Solution Files (*.hsol) 

For standard solution files, the user can specify which variables they want to be exported to the Tecplot file with the flag **-{}-output-variables=*.
The options are:

<div class="multicols" style="column-count: 5;">
  <ul>
    <li>\(Q\) (default)</li>
    <li>\(\rho\)</li>
    <li>\(u\)</li>
    <li>\(v\)</li>
    <li>\(w\)</li>
    <li>\(p\)</li>
    <li>\(T\)</li>
    <li>\(Mach\)</li>
    <li>\(s\)</li>
    <li>\(Vabs\)</li>
    <li>\(V\)</li>
    <li>\(Ht\)</li>
    <li>\(rhou\)</li>
    <li>\(rhov\)</li>
    <li>\(rhow\)</li>
    <li>\(rhoe\)</li>
    <li>\(c\)</li>
    <li>\(Nxi\)</li>
    <li>\(Neta\)</li>
    <li>\(Nzeta\)</li>
    <li>\(Nav\)</li>
    <li>\(N\)</li>
    <li>\(Ax\_Xi\)</li>
    <li>\(Ax\_Eta\)</li>
    <li>\(Ax\_Zeta\)</li>
    <li>\(ThreeAxes\)</li>
    <li>\(Axes\)</li>
    <li>\(lamb\)</li>
    <li>\(Qbase\)</li>
    <li>\(mpi\_rank\)</li>
    <li>\(eID\)</li>
    <li>\(gradV\)</li>
    <li>\(u\_x\)</li>
    <li>\(v\_x\)</li>
    <li>\(w\_x\)</li>
    <li>\(u\_y\)</li>
    <li>\(v\_y\)</li>
    <li>\(w\_y\)</li>
    <li>\(u\_z\)</li>
    <li>\(v\_z\)</li>
    <li>\(w\_z\)</li>
    <li>\(c\_x\)</li>
    <li>\(c\_y\)</li>
    <li>\(c\_z\)</li>
    <li>\(\omega\)</li>
    <li>\(\omega\_x\)</li>
    <li>\(\omega\_y\)</li>
    <li>\(\omega\_z\)</li>
    <li>\(\omega\_abs\)</li>
    <li>\(Qcrit\)</li>
  </ul>
</div>

When `--flow-eqs=mu`, only the whole set of variables \(Q\) is supported.

## Statistics Files (*.stats.hsol)
Statistics files can generate the standard variables as well as the following variables (being \(S_{ij}\) the components of the Reynolds Stress tensor):

<div class="multicols" style="column-count: 3;">
  <ul>
    <li>\(umean\)</li>
    <li>\(vmean\)</li>
    <li>\(wmean\)</li>
    <li>\(S_{xx}\)</li>
    <li>\(S_{yy}\)</li>
    <li>\(S_{zz}\)</li>
    <li>\(S_{xy}\)</li>
    <li>\(S_{xz}\)</li>
    <li>\(S_{yz}\)</li>
  </ul>
</div>



## Extract geometry
Under construction.

## Merge statistics tool

Tool to merge several statistics files. The usage is the following:

```bash
	$ horses.mergeStats *.hsol --initial-iteration=INTEGER --file-name=CHARACTER
```

Some remarks:

- Only usable with statistics files that are obtained with the "reset interval" keyword and/or with individual consecutive simulations.
- Only constant time-stepping is supported.
- If the hsol files have the gradients, the following flag must be used
```bash
	$ --has-gradients
```
- Dynamic p-adaptation is currently not supported.

## Mapping result to different mesh 
<a name="MaptomeshKey"></a>
`HORSES3D` addons, *horsesConverter*, has a capability to map result into different mesh file, with both have a consistent geometry. This is done by performing interpolation with the polynomial inside each element for each node point of the new mesh. The type of node quadrature will follow the quadrature defined in the .hmesh file with selected polynomial order in the control file. A control input file is required and must has name *horsesConverter.convert*. The template of control input file will be generated by default when executing *./horsesConverter* in a directory without the control file. Error message is given when at least one node point of the new mesh is not within any element of the old mesh. After completion, a new result file is generated and named *Result\_interpolation.hsol*. The required keywords in the control file are described in the table below. Command to execute:
```bash
	$ ./horsesConverter
```

| Keyword              | Description                                                        | Default value |
|----------------------|--------------------------------------------------------------------|---------------|
| Task                 | *meshInterpolation*                                               |               |
| Mesh Filename 1      | Location of the origin mesh (*.hmesh)                             |               |
| Boundary Filename 1  | Location of the origin boundary mesh (*.bmesh)                    |               |
| Result 1             | Location of the solution file with origin mesh (*.hsol)           |               |
| Mesh Filename 2      | Location of the target mesh (*.hmesh)                             |               |
| Boundary Filename 2  | Location of the target boundary mesh (*.bmesh)                    |               |
| Polynomial Order     | Polynomial order of the target mesh                                | (1, 1, 1)     |


## Generate OpenFOAM mesh 
<a name="GenerateOpenFOAMmeshKey"></a>
Another functionality of *horsesConverter* addons is to convert the mesh files, (\*.hmesh) and (\*.bmesh), into OpenFOAM format, the *polyMesh* folder. Each element is discretised into \(n_x \times n_y \times n_z\) cells distributed as Gauss-Lobatto nodes. The number of division of each element, (\(n_x\), \(n_y\), and \(n_z\)), is required in the control file, see [previous section](GenerateOpenFOAMmeshKey). After completion, a folder named `foamFiles` is generated. OpenFOAM mesh files, i.e. *points*, *faces*, *owner*, *neighbour*, and *boundary*, are located within `foamFiles/constant/polyMesh`. The required keywords in the control file are described in the table below Command to execute:

```bash
	$ ./horsesConverter
```

| Keyword              | Description                                                        | Default value |
|----------------------|--------------------------------------------------------------------|---------------|
| Task                 | *horsesMesh2OF*                                                   |               |
| Mesh Filename 1      | Location of the origin mesh (*.hmesh)                             |               |
| Boundary Filename 1  | Location of the origin boundary mesh (*.bmesh)                    |               |
| Polynomial Order     | Number of division of each element (\(n_x\), \(n_y\), and \(n_z\))      | (1, 1, 1)     |

NOTE: Before running the mesh in the OpenFOAM environment, the type of boundaries inside the boundary file needs to be adjusted according to the actual type (*patch*, *wall*, and *symmetry*).


## Generate HORSES3D solution file from OpenFOAM result
`HORSES3D` provides a capability to convert OpenFOAM result into `HORSES3D` solution file (\*.hsol). The mesh of the OpenFOAM result must be generated by converting HORSES3D mesh files, see the [previous section](GenerateOpenFOAMmeshKey). Beforehand, the OpenFOAM result must be converted into VTK format(*.vtk). This not only allows the result to be in the single file but also converts cell data into point data. In the OpenFOAM environment, the command for this conversion:  
```cpp
	$ foamToVTK -fields "(U p T rho)" -ascii -latestTime
```
The necessary file (.vtk) required for the control file input is inside VTK folder, see [previous section](GenerateOpenFOAMmeshKey) for the control file template. The `HORSES3D` solution file is named `Result\_OF.hsol`. The required keywords in the control file are described in the table below. Command to execute:

```bash
	$ ./horsesConverter
```

| Keyword                         | Description                                               | Default value |
|---------------------------------|-----------------------------------------------------------|---------------|
| Task                            | *OF2Horses*                                               |               |
| Mesh Filename 1                 | Location of the origin mesh (*.hmesh)                    |               |
| Boundary Filename 1             | Location of the origin boundary mesh (*.bmesh)           |               |
| Polynomial Order                | Polynomial order of the solution file (.hsol)             | (1, 1, 1)     |
| VTK file                        | Location of VTK file (.vtk)                               |               |
| Reynolds Number                 | Reynolds Number/m of the solution -- \(L_{ref}\)=1.0m       |               |
| Mach Number                     | Mach Number of the solution                               |               |
| Reference pressure (Pa)         | Reference Pressure                                        | 101325        |
| Reference temperature (K)       | Reference Temperature                                     | 288.889       |

