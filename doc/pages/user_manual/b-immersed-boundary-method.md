title: Immersed Boundary Method



The immersed boundary is activated during the simulation if the following lines are specified in the control file:

```markdown
#define IBM
   name                           = myIBM
   active                         = .true.
   penalization                   = 1.0d-6
   semi implicit                  = .false.
   number of objects              = 5 
   number of interpolation points = 15
   band region                    = .true.
   band region coefficient        = 1.3
   compute distance               = .true.
   clip axis                      =  1
   aab                            = .false.
   describe                       = .true.
   plot obb                       = .false.
   plot kdtree                    = .false.
   plot mask                      = .true.
   plot band points               = .false.
#end
```

A folder called 'IBM' must be created.

| Keyword                   | Description                                                                                                                                                                                | Default value |
|---------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| name                      | *CHARACTER*: Name assigned to immersed boundary method job.                                                                                                                                |               |
| active                    | *LOGICAL*: When .TRUE., the immersed boundary method is active.                                                                                                                            | .FALSE.       |
| penalization              | *REAL*: Specifies the value of the penalization term, \eta.                                                                                                                                | \(\Delta t\)    |
| theta energy source              | *REAL*: Specifies the value of the \(\theta\) parameter in the source term of the energy equation: \( - \lambda \rho ( \mathbf{u}_s \cdot (\mathbf{u} - \mathbf{u}_s) + \theta \lvert \lvert \mathbf{u} - \mathbf{u}_s \rvert \rvert^2 ) \). |  0.5    |
| semi implicit             | *LOGICAL*: The source term is treated in a semi-implicit manner.                                                                                                                           | .FALSE.       |
| number of objects        | *INTEGER*: Specifies the maximum number of objects inside a leaf of the KD-tree.                                                                                                            | 5             |
| number of interpolation points | *INTEGER*: Number of points used for the interpolation of the variables' values on the surface. It's needed for the computation of the forces.                                           | 15            |
| band region               | *LOGICAL*: If it's true, the band region is computed, otherwise it is not.                                                                                                                 | .FALSE.       |
| band region coefficient   | *INTEGER*: A region \(n\)-times the oriented bounding box is created (where \(n\) is the band region coefficient): all the points inside this region belong to the band region.            | 1.5           |
| compute distance          | *LOGICAL*: If it's true, the distance between the points in the band region and the STL file is computed, otherwise it is not. If the distance is not required, turn it off since it is an expansive operation. | .FALSE.       |
| clip axis                 | *INTEGER*: It's the axis along which the STL is cut. It is only needed if the forces are computed so that the integration of the variables is performed only on the portion of the STL surface lying inside the mesh. 1 corresponds to x-axis, 2 with y-axis, and 3 with z-axis. | 0 |
| aab                       | *LOGICAL*: The Axis Aligned Box is computed instead of the Oriented Bounding box. It is recommended when 'clip axis' \(\neq\) 0.                                                               | .FALSE.       |
| describe                  | *LOGICAL*: The immersed boundary parameters are printed on the screen.                                                                                                                    | .FALSE.       |
| plot obb                  | *LOGICAL*: The oriented-bounding box is plotted.                                                                                                                                           | .FALSE.       |
| plot kdtree               | *LOGICAL*: The kd-tree is plotted.                                                                                                                                                         | .FALSE.       |
| plot mask                 | *LOGICAL*: The degrees of freedom belonging to the mask are plotted.                                                                                                                       | .FALSE.       |
| plot band points          | *LOGICAL*: The band region's points are plotted.                                                                                                                                           | .FALSE.       |


## STL file

Immersed boundary requires, along with the mesh, a STL file. It must be put in the MESH folder with the mesh. The STL file name must be in lowercase character. In some programs, like AutoCAD, a STL file has always positive coordinates: the mesh should be built according to this consideration.
In the case of 2D simulations, the STL can be automatically cut by `horses3D` through the addition of the line *clip axis* (described in the previous section) so that only the STL portions inside the mesh are considered.

| Keyword              | Description                                                                                   | Default value       |
|----------------------|-----------------------------------------------------------------------------------------------|---------------------|
| number of stl =      | *INTEGER*: Number of stl files.                                                              | 1                   |
| stl file name =      | *CHARACTER*: The name of the STL file, without extension; it has to be inside the folder "MESH". | **Mandatory keyword** |
| stl file nameN =     | *CHARACTER*: The name of the \(\mathrm{N}^{th}\)-STL file (where N starts from 2), without extension; for the first STL just use "stl file name". It has to be inside the folder "MESH". | none                |


## Computing forces

In order to compute the forces on a body, the monitor should be defined as usual but the "Marker=" has to be equal to the name of the stl file on which the user wants to compute the forces. Given a STL file called "stlname", the monitor should be:

```markdown
#define surface monitor 1
   marker = stlname
   .
   .
   .
#end
```

The result of this operation is a `.tec` file inside the RESULTS folder. This file contains a scalar or a vector data projected on the body surface. 

## Moving bodies

If one or more of the stl files move, then the following lines must be added:

```markdown
#define stl motion 1
   stl name         = mySTL
   type             = rotation 
   angular velocity = 134.5d0
   motion axis      = 2
#end
\end{lstlisting}
```

| Keyword          | Description                                                                                        | Default value           |
|------------------|----------------------------------------------------------------------------------------------------|-------------------------|
| stl name         | *CHARACTER*: Name of the moving stl; it has to be equal to the name of one of the stl files.      | **Mandatory keyword** |
| type             | *CHARACTER*: Type of motion, it can be ROTATION or LINEAR.                                         | **Mandatory keyword** |
| angular velocity | *REAL*: Specifies the angular velocity. It must be in [Rad]/[s].                                    | **Mandatory keyword for rotation type** |
| velocity         | *REAL*: Specifies the translation velocity. It must be in [m]/[s].                                   | **Mandatory keyword for linear type** |
| motion axis      | *REAL*: Specifies the axis along which the rotation/translation occurs.                            | **Mandatory keyword** |
