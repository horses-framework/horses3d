title: Running a Simulation
---

[TOC]

## Control file (*.control) - Overview

The control file is the main file for running a simulation. A list of all the mandatory keywords for running a simulation and some basic optional keywords is presented in the [table](#runningkey) below. The specific keywords are listed in the other chapters.

<a name="runningkey"></a>

| Keyword | Description | Default value |
|---------|-------------|---------------|
| solution file name | *CHARACTER*: Path and name of the output file. The name of this file is used for naming other output files. | **Mandatory keyword** |
| simulation type | *CHARACTER*: Specifies if HORSES3D must perform a 'steady-state' or a 'time-accurate' simulation. | 'steady-state' |
| time integration | *CHARACTER*: Can be 'implicit', 'explicit', or 'FAS'. The latter uses the Full Algebraic Storage (FAS) multigrid scheme, which can have implicit or explicit smoothers. | 'explicit' |
| polynomial order | *INTEGER*: Polynomial order to be assigned uniformly to all the elements of the mesh. If the keyword *polynomial order file* is specified, the value of this keyword is overridden. | --* |
| polynomial order i <br> polynomial order j <br> polynomial order k | *INTEGER*: Polynomial order in the i, j, or k component for all the elements in the domain. If used, the three directions must be declared explicitly, unless you are using a polynomial order file. If the keyword *polynomial order file* is specified, the value of this keyword is overridden. | --* |
| polynomial order file | *CHARACTER*: Path to a file containing the polynomial order of each element in the domain. | --* |
| restart | *LOGICAL*: If .TRUE., initial conditions of simulation will be read from restart file specified using the keyword *restart file name*. | **Mandatory keyword** |
| cfl | *REAL*: A constant related with the *convective* Courant-Friedrichs-Lewy (CFL) condition that the program will use to compute the time step size. | --** |
| dcfl | *REAL*: A constant related with the *diffusive* Courant-Friedrichs-Lewy (DCFL) condition that the program will use to compute the time step size. | --** |
| dt | *REAL*: Constant time step size. | --** |
| final time | *REAL*: This keyword is mandatory for time-accurate solvers. | -- |
| mesh file name | *CHARACTER*: Name of the mesh file. The currently supported formats are *.mesh* (SpecMesh file format) and *.h5* (HOPR hdf5 file format). | **Mandatory keyword** |
| mesh inner curves | *LOGICAL*: Specifies if the mesh reader must suppose that the inner surfaces (faces connecting the elements of the mesh) are curved. This input variable only affects the hdf5 mesh reader. | .TRUE. |
| number of time steps | *INTEGER*: *Maximum* number of time steps that the program will compute. | **Mandatory keyword** |
| output interval | *INTEGER*: In steady-state, this keyword indicates the interval of time steps to display the residuals on screen. In time-accurate simulations, this keyword indicates how often a 3D output file must be stored. | **Mandatory keyword** |
| convergence tolerance | *REAL*: Residual convergence tolerance for steady-state cases. | **Mandatory keyword** |
| partitioning | *CHARACTER*: Specifies the method for partitioning the mesh in MPI simulations. Options are: 'metis' (the code must have been linked to METIS at compilation time, or 'SFC' (to use a space-filling curve method, no special compilation is needed for this option). | 'metis' |
| partitioning with weights | *LOGICAL*: Specifies if the method for partitioning the mesh in MPI simulations takes the local polynomial order as weights. Necessary for local polynomial refinement. | .TRUE. |
| manufactured solution | *CHARACTER*: Must have the value '2D' or '3D'. When this keyword is used, the program will add source terms for the conservative variables taken into account an exact analytic solution for each primitive variable j (\(\rho\), \(u\), \(v\), \(w\), $p$) of the form: \(j = j_C(1) + j_C(2) \sin(\pi j_C(5) x) + j_C(3) \sin(\pi j_C(6) y) + j_C(4) \sin(\pi j_C(7) z)\). Where \(j_C(i)\) are constants defined in the file *ManufacturedSolutions.f90*. Proper initial and boundary conditions must be imposed (see the test case). The mesh must be a unit cube. | -- |

\* One of these keywords must be specified.

\** For Euler simulations, the user must specify either the CFL number or the time-step size. For Navier-Stokes simulations, the user must specify the CFL and DCFL numbers *or* the time-step size.

---

## Boundary conditions

The boundary conditions are specified as blocks in the control file. The block starts with the keywords `#define' and ends with `\#end'. Inside the block, the options are specified as a pair of keywords and values, just like the normal body of the rest of the file.

Each boundary condition can be individually defined, or if multiple boundaries are set with the same definition, it could be done in the same block (with the name separated by a double underscore `$\_$' sign). The name of each boundary must match with the one specified in the mesh file.

The block in general can be seen below: 

```markdown
#define boundary myBoundary1__myBoundary2__myBoundary3
    type        = typeValue
    parameter 1 = value_1
    parameter 2 = value_2
# end
```

The table below shows the values for the type keyword, and the possible values for the parameters depend on the boundary condition.
```
| Keyword | Description | Default value |
|---------|-------------|---------------|
| type    | *CHARACTER*: Type of boundary condition to be applied. Options are: Inflow, Outflow, NoSlipWall, FreeSlipWall, Periodic, User-defined. | N/A |
```

By default, NoSlipWall is adiabatic. Isothermal wall can be activated with the following block:

```markdown
#define boundary myBoundary1__myBoundary2__myBoundary3
    type = NoSlipWall
    wall type (adiabatic/isothermal) = isothermal
    wall temperature = 2000.0d0  !Wall temperature in K
# end
```

It is also possible to set a moving wall with the keyword wall velocity = *value*.

For periodic boundary conditions, the second boundary that must be used as a complement must be specified by the keyword `coupled boundary`. These two boundaries must have the same node position in all directions but one. For mesh files generated by commercial software where this strict rule is not imposed, a comparison based on the minimum edge size of the face element can be used by a boolean parameter in the normal body of the control file (not in the block body), with the keyword `periodic relative tolerance`.

### Periodic revolution boundaries

Periodic revolution (rotary periodic) boundaries couple two boundaries that are images of each other by a rotation around an axis. Each pair is defined as a normal periodic boundary condition, and the periodicity by revolution is activated in the normal body of the control file:

```markdown
periodic enabled = .true.
rotation axis    = [0.0, 0.0, 1.0]
rotation center  = [0.0, 0.0, 0.0]

#define boundary sideA
   type             = Periodic
   coupled boundary = sideB
#end

#define boundary sideB
   type             = Periodic
   coupled boundary = sideA
#end
```

| Keyword          | Description | Default value |
|------------------|-------------|---------------|
| periodic enabled | *LOGICAL*: Match the periodic boundaries by a rotation around the rotation axis instead of a translation. | .false. |
| rotation axis    | *REAL ARRAY(3)*: Rotation axis. For periodic revolution it must be aligned with the x, y or z axis (otherwise periodic revolution is not activated). | [0.0, 0.0, 1.0] |
| rotation center  | *REAL ARRAY(3)*: A point of the rotation axis. | [0.0, 0.0, 0.0] |

Notes:

- The nodes of the two boundaries are matched by their radius and their axial position. The keyword `periodic relative tolerance` can also be used here.
- When `periodic enabled = .true.`, all the Periodic boundary pairs of the mesh are treated as periodic revolution boundaries. Translational and rotational periodicity cannot be combined in the same case.
- Periodic revolution boundaries are only available for the compressible Navier-Stokes solver.

<!-- 
%\emph{Juan's email (to be translated and adapted to the manual format as a complement):}
%
%Hola Gente,
%
%He tenido que hacer unas modificaciones bastante importantes en las BCs. Era la única parte del código que estaba “a la antigua” y no programada a objetos. Esto hacía que no fueran muy customizables, y por ejemplo las controlábamos con el número ese que siempre vale 0.0 jajaja. Ahora cada condición de contorno tiene los parámetros que necesitas y se pueden customizar. Lo malo es que ningún control file de los que tenéis van a seguir funcionando, pero os escribo los cambios para que sepáis adaptarlos, en cualquier caso, podéis pedirme ayuda y os cuento.
%
%Los cambios del código son:
%
%\begin{itemize}
%
%\item Las condiciones de contorno se definen igual que los monitores, con los $\#$define en la parte final del control file. Para definir una condicion de contorno se hace:
%\begin{lstlisting}
%        #define boundary name
%             type = Inflow/Outflow/NoSlipWall/FreeSlipWall/Periodic/User-defined
%             parametro1 = #valor
%             parametro2 = #valor
%        #end
%\end{lstlisting}
%\item Los parámetros1, … dependen de la condición de contorno que toque. Si no se especifica nada, pues está como estaba antes. Dos cambios importantes:
%
%          · He unificado las NoSlipWall (adiabatica e isoterma) en una sola. Por defecto es adiabática.
%          · En las periódicas es obligatorio ahora indicar a qué boundary se acopla (lo cual supone poco esfuerzo y reduce el tiempo de búsqueda al código)
%\begin{lstlisting}
%                 #define boundary name
%                       type = Periodic
%                      coupled boundary = nombredelboundaryalqueseacopla
%                 #end
%\end{lstlisting}
%\item Se pueden definir más de una condición de contorno del tirón, por ejemplo si boundary1, boundary2 y boundary3 son inflows se puede hacer:
%\begin{lstlisting}
%                 #define boundary boundary1__boundary2__boundary3
%                         type = Inflow
%                 #end
%\end{lstlisting}
%     es decir, separado con dos guiones bajos.
%
%\item Por pantalla, donde aparecía la info de las zones y tal, también aparece qué BC tiene y cuáles son los parámetros.
%
%\item La BC outflowspecifyP la llamo simplemente Outflow. Más que nada por que antes había algunos ficheros de control con la BC Outflow y no existía, pero por defecto se mandaba a Inflow. Para evitar problemas, pues Outflow.
%
%\item Los archivos de condición de contorno están en /physics/common en lugar de cada uno su archivo. Esto es por que al final son todas iguales y si se añade una nueva es más facil agregar un nuevo archivo que hacerlo individualmente en cada ecuación.
%
%\item Los bcTypeDictionary bcValueDictionary desaparecen. Las BC están en el module physics/common/BoundaryConditions.f90 como variable global, se llama BCs y dentro aloja todas las condiciones de contorno (una por zona, y en el mismo orden de las zonas).
%
%\end{itemize}
%
%Creo que eso es todo, lamento si os supone mucho cambio en vuestros ficheros de control que estéis corriendo a día de hoy, y si rompo algo que no reflejen los test. Pero estos cambios eran necesarios para darle más versatilidad (por ejemplo en multifase el inflow necesita bastante customización, para definir caudales de cada fase y cosas así). Además, creo que el enfoque OOP va en la dirección del resto del código.
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 -->

---

# Restarting a Case

| Keyword              | Description                                                                                                                                                                                                                             | Default value        |
|----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------|
| restart              | *LOGICAL*: If .TRUE., initial conditions of simulation will be read from restart file specified using the keyword *restart file name*.                                                                                                   | **Mandatory keyword** |
| restart file name    | *CHARACTER*: Name of the restart file to be written and, if keyword *restart* = .TRUE., also name of the restart file to be read for starting the simulation.                                                                         | **Mandatory keyword** |
| restart polorder     | *INTEGER*: Uniform polynomial order of the solution to restart from. This keyword is only needed when the restart solution is of a different order than the current case.                                                               | same as case's       |
| restart nodetype     | *CHARACTER*: Node type of the solution to restart from (Gauss or Gauss-Lobatto). This keyword is only needed when the restart node type is different than the current case.                                                          | same as case's       |
| restart polorder file| *CHARACTER*: File containing the polynomial orders of the solution to restart from. This keyword is only needed when the restart solution is of a different order than the current case.                                               | same as case's       |
| get discretization error of | *CHARACTER*: Path to solution file. This can be used to estimate the discretization error of a solution when restarting from a higher-order solution.                                                                                 | --                   |
| NS load from NSSA    | *LOGICAL*: Used only for NS simulations that are restarted from RANS SA.                                                                                                                                                                | .FALSE.              |



