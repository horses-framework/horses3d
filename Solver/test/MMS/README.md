# MMS convergence tests

Method of manufactured solutions for the HORSES3D equation sets. A smooth periodic
solution is imposed with a source term, and the L2 error at the final time gives the
order of convergence under h-refinement.

| Case | Solver | Equations |
|---|---|---|
| `NS`    | `horses3d.ns`   | compressible Navier-Stokes with Sutherland's law |
| `Euler` | `horses3d.ns`   | compressible Euler |
| `NSSA`  | `horses3d.nssa` | compressible RANS with Spalart-Allmaras |
| `iNS`   | `horses3d.ins`  | incompressible Navier-Stokes, artificial compressibility |
| `MU`    | `horses3d.mu`   | Cahn-Hilliard + incompressible Navier-Stokes, two fluids |
| `CAA`   | `horses3d.caa`  | acoustic perturbation equations (APE-4) on a uniform base flow |

The Cahn-Hilliard solver alone (`horses3d.ch`) has no user source term, so it has no
case; its discretization is tested inside `MU`.

## Run on GitHub

Actions > MMS convergence > Run workflow. Choose the case (or all), the polynomial
orders, the meshes and, optionally, changes to the control file such as
`Reynolds number = 100; Riemann solver = Roe`. Each case builds its solver, writes the
table of errors and orders to the job summary, and uploads `errors.csv`,
`convergence.png` and the control files and logs of the runs.

GitHub shows the Run workflow button only for workflows on the default branch. On
another branch, once the workflow has run in the repository, use the GitHub CLI:

```
gh workflow run mms.yml --ref <branch> -f case=NS -f P="2 3 4" -f meshes="4 8 16" -f set="Reynolds number = 100"
```

## Run locally

Build the solver (`./configure && make ns` in `Solver`), then:

```
pip install sympy gmsh matplotlib
python3 run.py NS                                   # P = 2 3 on meshes 4 8
python3 run.py NS --P 2 3 4 5 --meshes 4 8 16       # full study
python3 run.py NS --set "Reynolds number = 100"     # change the control file
python3 run.py MU --check                           # compare with MU/reference.csv
python3 plot.py NS                                  # NS/convergence.png
```

`run.py` generates `<case>/SETUP/ProblemFile.f90`, builds it with the HORSES3D build
settings (`COMPILER`, `ENABLE_THREADS`, ... from the environment and `Solver/make.inc`),
and runs every mesh and P with `<case>/control_template.control`, filling in the mesh,
the solution file and the order. Everything else, including `dt` and the number of
steps, comes from the control file; edit it or override keys with `--set`. Results go
to `<case>/errors.csv`, and pairs already run with the same solver, problem file and
control file are skipped. `--check` fails if a run fails or, with the committed
control file, an error differs from `<case>/reference.csv` by more than 1e-3.

The meshes are periodic cubes [-1,1]^3 with N^3 elements, made by the cartesian-box
generator of [horses3d-meshes](https://github.com/horses-framework/horses3d-meshes)
(Gmsh Python API). `run.py` downloads the script at a fixed commit into `MESH/`, or uses
`$HORSES3D_MESHES/cartesian-box/` if that variable points to a clone. No mesh files
are kept in this repository.

## Files

* `generate.py`: manufactured solutions and equations of each case. The source term is
  the residual of the solution, emitted without symbolic simplification. The physical
  parameters are read from the solver at run time.
* `ProblemFile.f90.in`: problem file template. It sets the initial condition, adds the
  source term and writes `mms_l2_error.dat` with the total and per-variable L2 errors,
  normalised by the mesh volume (integral of 1).
* `<case>/control_template.control`: control file of each case.
* `<case>/reference.csv`: errors of the committed case on upstream master (gfortran,
  RELEASE). After changing a case on purpose, copy its `errors.csv` here.

## Results

Orders of the total error between meshes 8 and 16 (theory P+1), upstream master f751d86f5,
gfortran:

| P | NS | Euler | NSSA | iNS | MU | CAA |
|---|---|---|---|---|---|---|
| 2 | 2.65 | 2.71 | 2.39 | 2.17 | 2.00 | 2.71 |
| 3 | 3.90 | 3.76 | 3.30 | 3.35 | 3.04 | 4.04 |
| 4 | 4.87 | 4.37 | 4.86 | 4.41 | 4.06 | 4.92 |
| 5 | 5.69 | | | | | |

* The meshes up to 16 are not fully in the asymptotic range: between 16 and 32 elements
  per side NS gives 2.95, 3.89, 4.72 and 5.54.
* iNS and MU converge closer to P because the artificial compressibility multiplies the
  divergence error by c0^2: lowering c0^2 brings iNS at P = 3 from 3.35 to 3.83. In MU at
  P = 2 the momentum and the concentration converge more slowly than the pressure, which
  dominates the total. In NSSA at P = 2 the momentum converges at about P.
* NS runs at Re = 1e5, where the viscous terms are a small part of the source term;
  NSSA runs at Re = 100.

## Assumptions

* NS and NSSA use the default Sutherland constants, so do not set `sutherland temperature`
  or `sutherland reference temperature`. NSSA uses the default Spalart-Allmaras constants;
  with no walls the wall distance is huge and only the production cb1 |omega| rho theta
  remains, so the field keeps |omega| > 0. MU assumes per-fluid sound speeds, not
  `artificial sound speed square`. CAA assumes a uniform base flow and no Lamb vector.
  Gravity is not included.
* Keep the final time: with u ~ sin(t), as in NS and iNS, a short run leaves the flow
  nearly at rest, the source is dominated by du/dt and the measured order drops.
* The error sums are rank-local, so run without MPI (OpenMP is fine).
* The legacy MU solver is not bitwise reproducible with more than one OpenMP thread
  (differences around 1e-7); use `OMP_NUM_THREADS=1` to compare runs exactly.
