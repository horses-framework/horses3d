title: Explicit Solvers
---


Explicit time integration schemes available in `HORSES3D`.
The main keywords to use it are shown in the table below

| Keyword             | Description                                                                                                                                                                                                                                                                                           | Default value |
|---------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| time integration    | *CHARACTER*: This is the main keyword to activate the multigrid solvers. The value of it should be set to 'FAS' for the Full Approximation Scheme (FAS) nonlinear multigrid  solvers and to 'AnisFAS' for anisotropic FAS schemes.                                                                 | 'explicit'    |
| simulation type     | *CHARACTER*: Specifies if HORSES3D must perform a ’steady-state’ or a ’time-accurate’. If 'time-accurate' the solver switches to BDF integration and uses FAS as a pseudo problem solver. Compatible only with 'FAS'.                                                                                  | 'steady-state' |
| explicit method     | *CHARACTER*: Select desired Runge-Kutta solver. Options are: 'Euler', 'RK3', 'RK5', 'RKOpt', 'SSPRK33', 'SSPRK43, 'MixedRK' and 'Multi Level RK3' (see the Multi-level RK3 section below)                                                                                                                                                                            | RK3           |
| rk order            | *INTEGER*: Order of Runge-Kutta method optimized for steady-state solver ('RKOpt'). Possible orders are from 2 to 7.                                                                                                                                                                                  | 2             |
| limit timestep      | *LOGICAL*: Activate the positivity limiter of Zhang and Shu (only for SSPRK methods).                                                                                                                                                                                                                | .false.       |
| limiter minimum     | *REAL*: Minimum value of density and pressure allowed by the limiter.                                                                                                                                                                                                                               | 1e-13         |
| compute time derivative after timestep | *LOGICAL*: Evaluate the time derivative again after each time step, so that the residuals and the monitors that use the time derivative or the gradients correspond to the new solution. Without it, they come from the last stage of the step. | .false. |

## Multi-level RK3

'Multi Level RK3' is a local time-stepping version of the low-storage 'RK3'. The elements are grouped into levels by their advective CFL number. Level 1 takes the global time step `dt`, and every RK3 stage of a level contains a complete RK3 step of the next level, so the elements of higher levels take smaller sub-steps. During the sub-steps of a level, only the elements of that level and higher (and their face neighbours) are evaluated.

| Keyword                  | Description                                                                                                                                                                                                 | Default value |
|--------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| explicit method          | *CHARACTER*: 'Multi Level RK3' activates the method.                                                                                                                                                        | RK3           |
| dt                       | *REAL*: Global time step, taken by level 1. A fixed time step is required: with a CFL-based time step the global step is set by the smallest elements and the method brings no gain.                        | --            |
| number of level          | *INTEGER*: Maximum number of levels.                                                                                                                                                                        | 3             |
| cfl cut-off              | *REAL*: CFL number that separates level 1 from level 2. The threshold of each following level is 2.5 times the previous one. The value is limited to the range [1e-4, 10].                                  | 0.5           |
| relevel iteration        | *INTEGER*: The levels are assigned at the first time step and then every 'relevel iteration' steps.                                                                                                         | 1000000000    |
| optimized partition      | *LOGICAL*: MPI only. Builds a level-aware partition, so that every MPI rank owns a share of the elements of each level. Works with both 'metis' and 'SFC' partitioning.                                     | .false.       |
| MLRK interface tolerance | *REAL*: Multiphase only. A node with a concentration between this tolerance and one minus it belongs to the interface. Interface elements are placed in the highest level.                                  | 1e-4          |
| MLRK user-defined levels | *LOGICAL*: Take the level of each element from `e % MLevel` as set by the user in `UserDefinedFinalSetup` (ProblemFile), instead of from the CFL number. The neighbour rule (at most one level difference) still applies. | .false.       |

How the levels are assigned:

- The advective CFL number of each element is computed with `dt` (viscous limits are not considered). An element is in level 1 if its CFL number is below 'cfl cut-off', in level 2 if it is below 2.5 times 'cfl cut-off', and so on. Elements above the last threshold go to the highest level, 'number of level'.
- Face neighbours of a higher-level element are raised so that neighbouring elements differ by at most one level.
- In multiphase simulations, the elements that contain the interface are placed in the highest level.
- The number of elements of each level is printed in the log under "Update Multi-Level Local Timestepping". If all the elements fall into level 1, the simulation stops: use a larger `dt` or a single-level method.

Notes:

- The evaluation restricted to a level is implemented in the Navier-Stokes (including Euler) and multiphase solvers. The other solvers accept the method, but evaluate the full time derivative at every sub-step, so there is no speed-up.
- As with 'RK3', the residuals and the monitors that use the time derivative or the gradients come from the last evaluation of each level. Set 'compute time derivative after timestep' to `.true.` to evaluate them after the step.

### Level file

At the first level assignment, the solver writes the level of every element to `./MESH/<solution file name>.lmesh`. The first line holds the number of elements, followed by one line per element, in global element order, with two integers:

- `MLevel`: the level the element is advanced with.
- `MLevelwN`: the level at which the element is needed, i.e. its own level or, for a neighbour of a higher-level element, that level. The level-aware partition uses this value.

The file can be added to the post-processing with `horses2plt` (see [Post-processing](s-postprocessing.html)), which exports the variables `mlrk_level` and `mlrk_level_wn`:
```
horses2plt MESH/mesh.hmesh RESULTS/solution.hsol --mlrk-level-file=MESH/solution.lmesh
```
