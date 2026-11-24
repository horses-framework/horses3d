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
case. `MU` solves Cahn-Hilliard with its own routines and free energy, and shares with
`horses3d.ch` only the BR1 discretization of the Laplacian.

## Step by step: using MMS in your development

Follow these steps to run a convergence study from a fresh checkout. The examples use the
NS case; for another case, replace `NS` with its name. Step 1 starts at the top of the
repository, and all other commands run from `Solver/test/MMS`.

### Before you begin

You need:

* Linux, because `run.py` builds the problem file as a shared library (`.so`).
* A Fortran compiler. The examples use gfortran, the default of the HORSES3D build.
* Python 3 with sympy, gmsh and matplotlib.
* Internet access the first time you run `run.py`, which downloads the mesh generator of
  horses3d-meshes.

To install the Python packages, run:

```bash
pip install sympy==1.14.0 gmsh==4.15.2 matplotlib   # the versions that the CI workflow uses
```

If pip isn't installed, or refuses to install packages for the system Python, as on
Ubuntu 24.04, do one of the following:

* On Ubuntu, install the system packages with
  `sudo apt install python3-sympy python3-gmsh python3-matplotlib`. These versions are
  older, but they work.
* Install venv with `sudo apt install python3-venv`, create a virtual environment with
  `python3 -m venv ~/mms && . ~/mms/bin/activate`, and then run the `pip` command above.

With pip, gmsh also needs the system libraries that the workflow installs:
`sudo apt install libglu1-mesa libgl1 libxcursor1 libxft2 libxinerama1 libfontconfig1`.

### 1. Build the solver

`run.py` runs the solver executable from `Solver/bin`, so you need to build the solver for
your case first. The make target is the last part of the executable name in the table at
the top of this page: `ns` for NS and Euler, and `nssa`, `ins`, `mu` or `caa` for the other
cases.

```bash
cd Solver
./configure      # once per checkout, and again after a pull that changes Makefile.in
make ns          # the defaults: gfortran, RELEASE mode, OpenMP threads and no MPI
```

Keep the default build without MPI (`COMM=SEQUENTIAL`). `run.py` starts the solver as a
single process, and the problem file sums the error only over the elements of its own
process, so the suite doesn't support MPI runs.

If you build with non-default options, pass the same options to `run.py`. It compiles the
problem file with the settings it finds in the environment (`COMPILER`, `MODE`,
`ENABLE_THREADS` and so on) and in `Solver/make.inc`. Make doesn't track compiler
options, so run `make clean` before you build with other options. For example, for a
debug build:

```bash
make clean
make ns MODE=DEBUG
cd test/MMS
MODE=DEBUG python3 run.py NS --P 3 --meshes 4 8
```

To go back to the default build, run `make clean` and `make ns` again in `Solver`.

### 2. Choose the manufactured solution

Each case imposes a known solution, called the manufactured solution, through a source
term. `generate.py` substitutes the solution into the equations that the solver
discretizes, and uses what is left over as the source term. The solution is then exact,
so the error of a run measures the error of the scheme.

The default solution makes every term of the equations active from the first time step.
For most changes to the solver, you can keep it. `FLOW` and `cells` in `generate.py`
define it, and it's periodic on the cube [-1,1]^3:

```
rho     = 1 + 0.2 sin(pi x) cos(pi y) sin(pi z) cos(t)
(u,v,w) = cos(t) [ (0.3 sin(pi z), 0.3 cos(pi z), 0) + 0.03 A + 0.15 B ]
p       = 1 + 0.1 sin(pi x) sin(pi y) sin(pi z) cos(t)

A = ( sin(pi x) cos(pi y) cos(pi z), cos(pi x) sin(pi y) cos(pi z), -2 cos(pi x) cos(pi y) sin(pi z) )
B = ( sin(pi x) cos(pi y) cos(pi z), cos(pi x) sin(pi y) cos(pi z),    cos(pi x) cos(pi y) sin(pi z) )
```

A is divergence-free and makes the vorticity vary. B is irrotational and compresses the
flow. The shear term (0.3 sin(pi z), 0.3 cos(pi z), 0) keeps |omega| > 0 everywhere,
which keeps the Spalart-Allmaras source term smooth (see the conditions below).

| Case | Fields | Where in `generate.py` |
|---|---|---|
| `NS`, `Euler`, `iNS` | rho, (u,v,w) and p above | `FLOW`, `cells` |
| `NSSA` | the same, and theta = 1 + 0.5 cos(pi x) cos(pi y) cos(pi z) cos(t) | `spalart_allmaras` |
| `MU` | c = 0.5 [1 + 0.8 cos(pi x) cos(pi y) cos(pi z) cos(t)], between 0.1 and 0.9, with (u,v,w) and p above; the density is the mixture c rho1 + (1 - c) rho2 of the two fluids in the control file | `multiphase` |
| `CAA` | the perturbations u' = cos(t) B and p' = sin(pi x) sin(pi y) sin(pi z) cos(t), and rho' = 0, on the uniform base flow in the control file | `acoustics` |

To use a different solution, edit these expressions in `generate.py`. The source term,
the initial condition and the error computation all follow from them. A new solution
must meet these conditions:

* It's smooth. The source term contains derivatives of the solution, up to the fourth
  derivative of c in MU. A jump would make the source term infinite, and the order of
  convergence would mean nothing. To model a shock, use a steep but smooth front, as in
  the example below. In NSSA, keep |omega| > 0 everywhere: with no walls, only the
  production term cb1 |omega| rho theta of Spalart-Allmaras remains, and it isn't smooth
  where omega = 0.
* It's periodic with period 2 in x, y and z, because the meshes are periodic cubes
  [-1,1]^3. Build it from sin(pi x), cos(pi x) and so on.
* It's physically valid: rho > 0 and p > 0 in the compressible cases, theta >= 0 in NSSA
  (the solver uses the theta >= 0 branch of Spalart-Allmaras), and c in [0, 1] in MU (the
  solver clips c to this range).
* It's stable with the time step of the control file. A faster flow might need a smaller
  `dt` (see step 5).
* It keeps every term active. A term that vanishes for your solution isn't tested.

#### Example: add a moving density front

This example replaces the density with a front that moves in x. The front is steep, like a
shock, but smooth. In `generate.py`, replace the first line of `FLOW`:

```python
FLOW = dict(rho=1 + 0.2*sp.sin(pi*x)*sp.cos(pi*y)*sp.sin(pi*z)*sp.cos(t),
```

with these two lines:

```python
delta = 0.3                    # width of the front: smaller is steeper
FLOW = dict(rho=1 + 0.2*sp.tanh(sp.sin(pi*(x - t))/delta),
```

The function tanh(sin(pi (x - t))/delta) is close to 1 or -1. It jumps smoothly wherever
sin(pi (x - t)) changes sign: up at x = t and down at x = t ± 1. The jumps move in x at
speed 1, and their width is proportional to delta. NS, Euler, NSSA and iNS all take their
density from `FLOW`.

Run the case without `--check`, because the reference errors are for the default
solution:

```bash
python3 run.py NS --P 3 --meshes 4 8 16
```

With delta = 0.3, the order of the total error at P = 3 is 0.6 between meshes 4 and 8,
and 2.4 between meshes 8 and 16. The total is mostly the energy error, though: the
density, which carries the front, converges at only about order 1 between meshes 8 and
16, because mesh 16 only starts to resolve the front. With delta = 0.1, the order of the
total between meshes 8 and 16 is 0.2, so a steeper front needs finer meshes.

Undo the change before you go on, or the `--check` runs of step 6 fail:

```bash
git checkout generate.py
```

### 3. Change the problem file (optional)

All cases use one problem file template, `ProblemFile.f90.in`. On every run,
`generate.py` inserts the solution (`MMS_State`) and the source term (`MMS_Source`) into
the template, and saves the result as `<case>/SETUP/ProblemFile.f90`. To change what the
problem file does, edit the template. If you edit `SETUP/ProblemFile.f90` instead, the
next run overwrites your change.

These routines in the template do the MMS work:

* `UserDefinedInitialCondition` sets the exact solution at t = 0.
* `UserDefinedSourceTermNS` adds the source term.
* `UserDefinedFinalize` computes the L2 error of each variable at the final time, prints
  it to the log and writes it to `mms_l2_error.dat`, which `run.py` reads.

Keep the arguments of these routines the same as in
`Solver/src/libs/problemfile/ProblemFile.f90`, the solver's own problem file.

#### Example: print the largest nodal error

To print the largest nodal error of each variable as well, edit `UserDefinedFinalize` in
`ProblemFile.f90.in`. Change the line marked `! changed`, and add the three lines marked
`! new`:

```fortran
            real(kind=RP) :: Qex(NCONS), wq, vol, sq(NCONS), emax(NCONS)    ! changed
            ...
            sq   = 0.0_RP
            emax = 0.0_RP                                                   ! new
            ...
                  sq  = sq + wq * (Qex - e % storage % Q(:,i,j,k))**2
                  emax = max(emax, abs(Qex - e % storage % Q(:,i,j,k)))     ! new
            ...
            print *, "MMS L2 error vars  : ", sqrt(sq / vol)
            print *, "MMS max error vars : ", emax                          ! new
```

Then run a case and open its log, `NS/CONTROL/N4_P3.log`:

```bash
python3 run.py NS --P 3 --meshes 4
```

The log has a line that starts with `MMS max error vars :`, followed by the largest error
of each variable. This change doesn't alter the L2 errors. To undo it, run
`git checkout ProblemFile.f90.in`.

### 4. Generate the problem file (optional)

`run.py` generates and compiles the problem file on every run, so you don't need to do it
yourself. To read the Fortran code that `generate.py` writes, for example to check the
source term of a new solution before you run it, run `generate.py` with the name of a
case (`NS`, `Euler`, `NSSA`, `iNS`, `MU` or `CAA`):

```bash
python3 generate.py NS          # writes NS/SETUP/ProblemFile.f90
```

### 5. Set up the run

The control file of each case, `<case>/control_template.control`, sets everything else:
the physical parameters, the scheme, `dt` and the number of time steps. `run.py` fills in
the mesh, the solution file and the polynomial order, which appear in the file as
`{MESH_FILE}`, `{SOLUTION_FILE}` and `{P}`.

Use a fixed time step. With `cfl` instead, the runs fail: NS stops at start-up unless
`dcfl` is also set, and with both the runs stop at another time, which `run.py` reports
as a failure.

To change a case permanently, edit its control file, and then update its reference (see
step 8). For example, `NS/control_template.control` contains these lines:

```
                         Reynolds number = 10.0
                    Discretization nodes = Gauss-Lobatto
                          Riemann solver = Matrix dissipation
                 Inviscid discretization = Split-form
                                      dt = 1.0e-4
                    Number of time steps = 200
```

To change keys for one run only, use `--set`:

```bash
python3 run.py NS --P 3 --meshes 4 8 --set "Reynolds number = 100"
python3 run.py NS --P 3 --meshes 4 8 --set "dt = 5.0e-5; Number of time steps = 400"
```

The first command changes the errors: on mesh 4, the L2 error is 4.0069e-04 instead of
4.2541e-04. The second halves the time step and keeps the final time. Its errors change by
less than 1e-9 relative, which shows that the time-integration error is negligible.

If you set a key that isn't in the control file, `run.py` adds it and prints a warning,
because HORSES3D ignores keys that it doesn't read. Check the spelling of the key,
including any units in its name, such as `Viscosity (Pa.s)` in iNS.

### 6. Run the case

`run.py` runs the case for each polynomial order that you give with `--P`, on each mesh of
N x N x N elements that you give with `--meshes`. It prints the L2 error of each run, and
then a table of errors and orders of convergence. If you leave out `--P` and `--meshes`,
`run.py` runs the default sweep in `sweeps.json`, P = 3 to 6 on meshes 4, 8 and 16, which
takes much longer.

`run.py` skips the runs that it has already done with the same build. To run them again,
use `--force`. `--keep` keeps the solution files only of the runs that `run.py` makes, so
combine it with `--force` for runs already done. Only one `run.py` at a time can work on
a case.

These are the most common commands:

```bash
python3 run.py NS --P 3 4 --meshes 4 8              # quick run, about 15 s on a laptop
python3 run.py NS --P 3 4 --meshes 4 8 --check      # the same, compared with NS/reference.csv
python3 run.py NS --check                           # default sweep, compared with NS/reference.csv
python3 run.py NS --P 3 --meshes 4 8 --force        # run again the runs already in NS/errors.csv
python3 run.py NS --P 3 --meshes 4 --keep --force   # run again, keeping the solution files in NS/RESULTS
python3 run.py MU --P 3 4 --meshes 4 8 --check      # another case (build it first: make mu)
python3 plot.py NS                                  # draw NS/convergence.png
python3 run.py --help                               # list every option
```

To run every case on GitHub instead, use the workflow (see Checking a change and Run on
GitHub). The commands below need the GitHub CLI, logged in with `gh auth login`, and your
branch pushed. Running a workflow needs write access to the repository: without it, push
the branch to your fork, enable Actions there, and add `-R <your-user>/horses3d` to the
commands.

```bash
gh workflow run mms.yml --ref <your-branch>                                   # every case, compared with the references
gh workflow run mms.yml --ref <your-branch> -f case=NS -f P="3 4" -f meshes="4 8 16"
gh workflow run mms.yml --ref <your-branch> -f case=NS -f set="Reynolds number = 100"
```

### 7. Read the results

After the runs, `run.py` prints a table for the case. For each P, the table lists every
mesh with:

* its number of degrees of freedom (NDOF)
* its total L2 error
* the order of convergence against the previous mesh, for the total error, for the
  variable with the lowest order (named in parentheses) and for each variable

Between meshes 8 and 16, the orders should be close to P or higher (see Results). With
`--check`, the last lines show how many runs `run.py` compared with the reference,
followed by `PASS` or `FAILED`.

For example, on a fresh checkout, `python3 run.py NS --P 3 4 --meshes 4 8 --check` prints
the following. If you ran the quick run first, the first four lines say instead that the
runs were skipped.

```
NS N=4 P=3: L2 = 4.2541e-04 (1 s)
NS N=8 P=3: L2 = 4.9913e-05 (3 s)
NS N=4 P=4: L2 = 5.0827e-05 (1 s)
NS N=8 P=4: L2 = 3.2747e-06 (6 s)

MMS NS: dt = 0.0001, t_final = 0.02

| P | N | NDOF | L2 total | order total | lowest order | order rho | order rhou | order rhov | order rhow | order rhoe |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 4 | 4096 | 4.2541e-04 |  |  |  |  |  |  |  |
| 3 | 8 | 32768 | 4.9913e-05 | 3.09 | 2.57 (rhou) | 3.04 | 2.57 | 3.06 | 3.03 | 3.13 |
| 4 | 4 | 8000 | 5.0827e-05 |  |  |  |  |  |  |  |
| 4 | 8 | 64000 | 3.2747e-06 | 3.96 | 3.43 (rhov) | 3.97 | 3.73 | 3.43 | 3.97 | 4.08 |

Compared 4 of 4 runs with reference.csv, relative tolerance 0.001
PASS
```

`run.py` also writes these files:

* `NS/errors.csv`: the results of every run
* `NS/CONTROL/`: the control file and the log of each run

`python3 plot.py NS` draws the results of the last runs in `NS/convergence.png`.

### 8. Update the reference

`<case>/reference.csv` contains the errors of the default sweep for the committed case. If
a deliberate change to the solution, the problem file or the control file changes the
errors of a case, `--check` fails until you write a new reference.

To update the reference, do the following:

1. Check that the orders of convergence are still right.
2. Write the new reference:

   ```bash
   python3 run.py NS --write-reference      # runs the default sweep, then writes NS/reference.csv
   ```

   `--write-reference` refuses `--set`, `--P` and `--meshes`: it runs the default sweep
   with the control file as it is.

3. Update the numbers in this README that changed, as described in Maintaining the suite.

To write the reference from a GitHub run instead, which is faster, see Maintaining the
suite.

## Checking a change

Run the workflow on the branch with the default inputs (see below). Every case then runs
P = 3 to 6 on meshes 4, 8 and 16, and the workflow fails if a run fails or if an error
differs from `<case>/reference.csv` by more than a relative 1e-3, which leaves room for
round-off from other compilers and versions. The summary of the report job of each case
gives the errors, the orders of convergence between consecutive meshes and the result of
the comparison. If a change to a case is intended, check that its orders are still those
of the Results section and regenerate its reference (see Maintaining the suite).

## Run on GitHub

Actions > MMS convergence > Run workflow. Choose the case (or all), the polynomial
orders and the meshes (empty for the defaults of each case in `sweeps.json`) and,
optionally, changes to the control file, which apply to every selected case (e.g.
`Reynolds number = 100; Riemann solver = Roe` with case NS: MU, iNS and CAA have no Roe
solver). Each case and order runs in its own job, in parallel; a job has a 6-hour limit,
which mesh 32 at high P may exceed. A final job per case merges them, writes the table of
errors and orders and the result of the comparison to the job summary, fails if a run
has no result, and uploads `errors.csv`, `convergence.png` and the control files, logs
and problem files of the runs. The default sweep, P = 3 to 6 on meshes 4, 8 and 16 for
every case, takes about 15 to 20 minutes for all cases.

GitHub shows the Run workflow button only for workflows on the default branch. On
another branch, once the workflow has run in the repository, use the GitHub CLI:

```
gh workflow run mms.yml --ref <branch> -f case=NS -f P="2 3 4" -f meshes="4 8 16" -f set="Reynolds number = 100"
```

## Run locally

The step-by-step guide above shows how to install the Python packages, build the solver
and run the cases; steps 5 and 6 cover the control file, `--set`, skipped runs and
`--force`. In more detail, `run.py` generates `<case>/SETUP/ProblemFile.f90`, builds it
with the HORSES3D build settings (`COMPILER`, `ENABLE_THREADS`, ... from the environment
and `Solver/make.inc`), and runs every mesh and P with `<case>/control_template.control`,
filling in the mesh, the solution file and the order. Results go to `<case>/errors.csv`.
`run.py` exits with an error if a run fails and, with `--check`, if an error differs by
more than a relative 1e-3 from the row of the same mesh and P in `<case>/reference.csv`,
or if a run of the default sweep has no row there. Runs with `--set`, and meshes or
orders outside the default sweep, are not compared; the last lines of the output say how
many runs were.

The meshes are periodic cubes [-1,1]^3 with N^3 elements, made by the cartesian-box
generator of [horses3d-meshes](https://github.com/horses-framework/horses3d-meshes)
(Gmsh Python API). `run.py` downloads the script at a fixed commit into `MESH/`, or uses
`$HORSES3D_MESHES_DIR/cartesian-box/` if that variable points to a clone, and runs it
with `$GMSH_PYTHON` if set, as in the CI workflows. The meshes are kept in `MESH/` and
reused: after changing the generator, delete `MESH/` and run with `--force`, since the
skipped pairs do not depend on the meshes. No mesh files are kept in this repository.

## Files

* `generate.py`: manufactured solutions and equations of each case. The solutions are at
  full amplitude from t = 0. Except in CAA, the density (concentration in MU) varies and
  the velocity is neither divergence-free nor irrotational, so that the terms of the
  equations are active from the first step; in CAA the density perturbation is 0, since
  APE-4 does not evolve it, and the velocity perturbation is irrotational. The source
  term is the residual of the solution, emitted without symbolic simplification. The
  physical parameters are read from the solver at run time.
* `ProblemFile.f90.in`: problem file template. It sets the initial condition, adds the
  source term and writes `mms_l2_error.dat` with the total and per-variable L2 errors,
  normalised by the mesh volume (integral of 1).
* `<case>/control_template.control`: control file of each case.
* `sweeps.json`: default orders and meshes of each case, used by `run.py` and the workflow.
* `<case>/reference.csv`: errors of the default sweep of the committed case, written by
  `run.py --write-reference` from a workflow run with the solver of upstream master
  6c750adef (gfortran 15.2, RELEASE). Maintaining the suite explains how to regenerate it.

## Results

Orders of convergence between meshes 8 and 16 in the reference runs (gfortran). Each cell
gives the order of the total error and, in brackets, the lowest order of any single
variable (CAA leaves out rho, which APE-4 does not evolve). The job summary of each run
gives the order of every variable between consecutive meshes, and the lowest of them.

| P | NS | Euler | NSSA | iNS | MU | CAA |
|---|---|---|---|---|---|---|
| 3 | 3.18 (3.08) | 3.68 (3.41) | 3.11 (3.07) | 3.82 (3.18) | 3.03 (2.47) | 3.16 (3.13) |
| 4 | 4.31 (3.92) | 4.59 (4.33) | 4.48 (4.09) | 4.83 (4.01) | 3.96 (2.52) | 4.34 (4.29) |
| 5 | 5.18 (5.04) | 5.67 (5.45) | 5.05 (5.00) | 5.68 (5.06) | 5.03 (4.95) | 5.32 (5.17) |
| 6 | 6.75 (6.16) | 6.91 (6.65) | 6.66 (6.12) | 6.98 (6.00) | 6.17 (6.14) | 6.43 (6.37) |

* These are the orders of short runs between two meshes, not asymptotic rates. The runs
  end at t = 0.02 (200 steps), or at t = 2e-4 in MU, whose explicit time step is limited
  by Cahn-Hilliard. The initial condition is exact at the nodes, where the error is
  measured, so the whole error is made during the run by the truncation error of the
  spatial operator at the nodes, which is the residual printed at iteration 0 of each log
  (maximum norm). Between meshes 8 and 16 its order is about P in CAA, between P - 0.7
  and P + 0.4 in Euler, between P - 1.6 and P + 0.2 in NS, NSSA and iNS (lowest at even
  P, from BR1), and between P - 2.7 and P + 0.1 in MU (lowest in c at even P). The
  dissipation of the scheme damps part of this error during the run, faster on finer
  meshes, so the orders change with the final time and with the meshes.
* Every variable converges at order P - 0.1 or higher between meshes 8 and 16, except in
  MU: the momentum at P = 3 and 4, and c at P = 4. In iNS and MU the total is mostly the
  pressure error (at least 90% of the squared total on mesh 16, 71% in iNS at P = 6), so
  the order of the total is the order of the pressure.
* One-off runs on finer meshes with the same final time give the orders below (CI does
  not run them). NS and iNS at P = 5 and 6 and NSSA at P = 6 used half the time step; MU
  at P = 5 and 6 used dt = 2e-7. Most orders rise. In NS and NSSA the odd P stay near P
  and the even P rise well above P, as expected from BR1. The iNS total falls at P = 5
  and 6; at P = 6 the density error, of order about P, grows to a third of the squared
  total on mesh 32. In MU, c stays below P at P = 3, and at P = 4 the momentum falls to
  3.6 between meshes 24 and 32. In Euler at P = 6, halving the step changes the mesh-32
  errors (about 4e-13) by at most 1.6e-4 relative.

| Case, meshes | P = 3 | P = 4 | P = 5 | P = 6 |
|---|---|---|---|---|
| NS, 16-32 | 3.25 (3.15) | 4.75 (4.61) | 5.11 (5.00) | 7.23 (6.88) |
| Euler, 16-32 | 3.62 (3.53) | 5.04 (4.58) | 6.05 (5.63) | 6.89 (6.11) |
| NSSA, 16-32 | 3.18 (3.11) | 4.82 (4.56) | 5.14 (5.06) | 7.28 (6.79) |
| iNS, 16-32 | 3.85 (3.14) | 5.00 (4.06) | 5.55 (5.06) | 6.83 (6.11) |
| MU, 16-24 | 3.05 (2.60) | 4.19 (4.18) | 5.08 (4.95) | 6.30 (6.24) |
| MU, 24-32 | 3.06 (2.73) | 4.45 (3.62) | | |
| CAA, 16-32 | 3.32 (3.15) | 4.62 (4.54) | 6.04 (5.89) | 6.91 (6.80) |

* Runs five times longer (t = 0.1) give these orders of the total error between meshes 8
  and 16 for P = 3 to 6: NS 3.51, 5.04, 5.12, 7.81; NSSA 3.34, 5.29, 5.25, 7.33; CAA 4.36,
  5.06, 6.02, 6.91. To run longer, raise `Number of time steps`.
* Halving the time step changes the errors at P = 6 on mesh 16 by at most 5.3e-5 (iNS),
  1.1e-6 (MU), 2.3e-7 (Euler), 7.4e-8 (CAA), 6.5e-8 (NSSA) and 2.6e-8 (NS), so they are
  errors of the spatial discretization.
* NS runs at Re = 10 and iNS with a viscosity of 0.01, so that the viscous terms are about
  a fifth of the momentum source; NSSA runs at Re = 100. Euler runs the defaults of
  `horses3d.ns` (Gauss nodes, standard volume integral, Roe flux), and NS the split form
  on Gauss-Lobatto nodes. The two MU fluids have different densities, viscosities and
  sound speeds, so that the mixture laws, which weight each property by the
  concentration, are exercised. CAA runs a base flow with rho0 = 1.2 and c0^2 = 1.5.
* Not tested: Riemann solvers other than matrix dissipation and Roe in the compressible
  cases, and other than the exact one in iNS and MU; two-point fluxes other than
  Chandrasekar (NS) and Skew-symmetric 2 (iNS); BR2 and IP; gradient variables other than
  the state; Gauss nodes with viscous terms; CAA fluxes other than Lax-Friedrichs; the
  rotational part of the APE-4 velocity flux, since the CAA velocity perturbation is
  irrotational; time integrators other than RK3; non-periodic boundary conditions; meshes
  other than uniform Cartesian ones (all elements are equal cubes with the same
  orientation, so the metric terms are constant and diagonal), non-conforming meshes, and
  orders that differ between elements or directions (mortars, p-adaptation); LES, SVV and
  shock capturing; the Spalart-Allmaras wall and destruction terms and its negative branch;
  the clipping of c in MU (c stays in [0.1, 0.9]); gravity; non-uniform CAA base flows and
  the Lamb source; `horses3d.ch` alone; MPI.

## Assumptions

* NS and NSSA use the default Sutherland constants, so do not set `sutherland temperature`
  or `sutherland reference temperature`. NSSA uses the default Spalart-Allmaras constants;
  with no walls the wall distance is huge and only the production cb1 |omega| rho theta
  remains, so the field keeps |omega| > 0. CAA assumes a uniform base flow and no Lamb
  vector. Gravity is not included.
* The time steps of the control files are stable up to P = 6 on mesh 16, the finest
  pair of the default sweeps, where twice the step is still stable in every case; higher
  orders or finer meshes may need a smaller `dt`, with more steps to keep the final time.
  In MU the explicit step is limited by Cahn-Hilliard: the limit scales like h^4, falls
  steeply with P and grows with the chemical characteristic time. With 20000 s (a small
  mobility), dt = 1e-6 is stable; at P = 6 on mesh 16, 2.5 times this step is stable and
  5 times diverges.
* The default sweeps start at P = 3: at P = 2 the truncation error of the MU momentum
  does not converge.
* The error sums are rank-local, so run without MPI (OpenMP is fine).
* The runs are bitwise reproducible, from run to run and with 1 or 4 OpenMP threads.

## Maintaining the suite

### How a run works

`generate.py <case>` writes `<case>/SETUP/ProblemFile.f90` from `ProblemFile.f90.in`,
with two routines generated from the manufactured solution: `MMS_State` gives the
solution and `MMS_Source` its residual at a point and time. The problem file sets the
initial condition with `MMS_State`, returns `MMS_Source` as the user source term, and at
the end of the run writes `mms_l2_error.dat` with the number of elements, P, the number
of degrees of freedom, the volume, the L2 errors and the final time. `run.py` builds it
into `<case>/SETUP/libproblemfile_<solver>.so`, which the solver loads from `./SETUP`,
runs the solver on each mesh and P in `<case>`, checks that file, and adds the errors to
`<case>/errors.csv`. The control file and log of each run are in `<case>/CONTROL/`.

The workflow (`.github/workflows/mms.yml`) has three jobs: `setup` reads the case list,
the orders and the meshes; `run` runs one case and P per job, with `--check` unless
turned off; `report` merges the `errors.csv` of the run jobs of a case with
`run.py --merge`, compares them with the reference if the check is on, and plots them.

### Changing a case or the defaults

The control files, `sweeps.json` and the manufactured solutions of `generate.py` define
the cases. After changing one of them on purpose:

1. Regenerate the reference from a default sweep with the committed control file. On
   GitHub, run the case with the check off, download the merged errors and write them:

   ```
   gh workflow run mms.yml --ref <branch> -f case=MU -f check=false
   gh run download <run id> -n mms-MU -D /tmp/mms-MU
   python3 run.py MU --merge /tmp/mms-MU/errors.csv --write-reference
   ```

   Locally, `python3 run.py MU --write-reference` runs the default sweep and writes
   `MU/reference.csv` from it. `--write-reference` refuses runs with `--set`, with other
   orders or meshes than those of `sweeps.json`, or with missing pairs.
2. Measure again the numbers of the Results section that change: the orders (job
   summary), the truncation orders (iteration 0 of the logs), the stability margin (twice
   `dt` with half the steps, at P = 6 on mesh 16) and the effect of halving `dt`, with
   `--set "dt = ...; Number of time steps = ..."`. The user manual page
   (`doc/pages/user_manual/s-verification.md`) quotes the range of the orders too.
3. Copy the NS and MU files to horses3d-gpu (see below).

### Adding a case

1. In `generate.py`, write a function that returns the variables of the solver, the
   residual of each equation, the Fortran names of the variables (`IRHO`, ...) and the
   solver parameters it uses, as Fortran expressions of the arguments of the problem-file
   routines. Add it to `CASES` with the solver and the names of the variables in the
   order in which the solver stores them, and add the solver to `MACROS` if it is new.
2. Write `<case>/control_template.control`. `run.py` fills `{MESH_FILE}`,
   `{SOLUTION_FILE}` and `{P}`. Give `dt` and `Number of time steps`, keep every key
   before the first `#define`, and make the six boundaries periodic: `Left`, `Right`,
   `Bottom`, `Top`, `Back` and `Front`, as the mesh generator names them.
3. Run it, check that twice the time step is stable, regenerate its reference, and add
   it to the `options` and to the `all` list of `.github/workflows/mms.yml`, to the
   tables of this file and of `doc/pages/user_manual/s-verification.md` and, if it needs
   other orders or meshes, to `sweeps.json`.

### Pitfalls

* The problem file is a separate library with its own copies of the solver's modules,
  which are never initialised: read the solver parameters from the arguments of the
  problem-file routines (`thermodynamics_`, `dimensionless_`, `refValues_`, `multiphase_`,
  `Qbase`), never from module variables.
* `ProblemFile.f90.in` copies the interfaces of the problem-file routines from
  `Solver/src/libs/problemfile/ProblemFile.f90`; change it when they change.
* HORSES3D ignores the keys it does not read, and some defaults change the scheme:
  without `Riemann solver`, iNS and MU use the central flux although the log says
  `Exact`, and without `Penalty parameter` MU has no BR1 interface penalty (the log
  shows `Penalty parameter: 0.000`).
* The vector-split (the default), exact-average and exact-jump Riemann solvers of CAA
  leave out the base-flow terms of the face flux
  ([#422](https://github.com/horses-framework/horses3d/issues/422)), so the CAA case uses
  `lax-friedrichs`.
* `--merge` takes the rows of the last build of each file and checks that their `dt` and
  final time are those of the control file with `--set`, but not that the files come from
  the same case or the same other keys.

### horses3d-gpu

The suite is also in horses3d-gpu, for NS and MU. `generate.py`, `run.py`, `plot.py`,
`ProblemFile.f90.in`, `sweeps.json`, `.gitignore` and the NS and MU control files and
references are the same files in both repositories: change them in both, and copy the
references from here. Only this README, `doc/pages/user_manual/s-verification.md` and
`mms.yml` (case lists, thread and library settings) differ.
