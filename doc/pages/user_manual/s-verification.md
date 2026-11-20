title: Verification with manufactured solutions
---

[TOC]

# MMS convergence tests

`Solver/test/MMS` is a convergence study by the method of manufactured solutions (MMS).
For each case, a smooth periodic solution is imposed with a source term in the problem
file, the solver runs on periodic cube meshes of N x N x N elements, and the L2 error at
the final time gives the order of convergence under mesh refinement.

| Case  | Solver        | Equations                                                       |
|-------|---------------|-----------------------------------------------------------------|
| NS    | horses3d.ns   | compressible Navier-Stokes with Sutherland's law                |
| Euler | horses3d.ns   | compressible Euler                                              |
| NSSA  | horses3d.nssa | compressible RANS with Spalart-Allmaras                         |
| iNS   | horses3d.ins  | incompressible Navier-Stokes, artificial compressibility        |
| MU    | horses3d.mu   | Cahn-Hilliard and incompressible Navier-Stokes, two fluids      |
| CAA   | horses3d.caa  | acoustic perturbation equations (APE-4) on a uniform base flow  |

By default every case runs P = 3 to 6 on meshes with N = 4, 8 and 16, and its errors are
compared with the reference errors of the case. Run it to check that a change to the code
keeps the errors and the orders of convergence of these cases.

## Running on GitHub

The workflow *MMS convergence* runs by hand: Actions > MMS convergence > Run workflow.
On a branch other than the default one, use the GitHub CLI:

```bash
gh workflow run mms.yml --ref <branch>
```

The inputs are the case (or all), the polynomial orders, the meshes, changes to the
control file (e.g. `Reynolds number = 100`) and whether to compare with the reference
errors, which is the default. Each case and order runs in its own job. A final job per
case writes the errors, the orders of convergence between consecutive meshes and the
result of the comparison to its summary, and uploads the errors, a convergence plot and
the control files and logs of the runs.

## Running locally

Build the solver of the case, e.g. `./configure && make ns` in `Solver`, and then, in
`Solver/test/MMS`:

```bash
pip install sympy==1.14.0 gmsh==4.15.2 matplotlib   # the versions of the workflow
python3 run.py NS --check                           # default sweep, compared with NS/reference.csv
python3 run.py NS --P 3 4 --meshes 4 8              # quick run
python3 plot.py NS                                  # NS/convergence.png
```

Where pip refuses to install packages for the system Python (e.g. Ubuntu 24.04), use a
virtual environment: `python3 -m venv ~/mms && . ~/mms/bin/activate`. The README lists the
system libraries that gmsh needs.

`run.py` generates and builds the problem file of the case, makes the meshes with Gmsh,
runs the solver for every mesh and order with the control file of the case and prints
the errors and orders. It needs Linux, since the problem file is built as a shared
library.

## Reading the results

* The comparison fails if an error differs by more than a relative 1e-3 from the
  reference error of the same mesh and order, or if a run of the default sweep has no
  reference error. A run that diverges, stops early or writes no errors fails too.
* The orders are those of short runs between two meshes, not asymptotic rates. Between
  meshes 8 and 16 the order of the total error is between P and P + 1 (3.96 in MU at
  P = 4), and every variable converges at order P - 0.1 or higher, except some in MU.
* After a change that alters the errors on purpose, check the orders and regenerate the
  reference errors of the case with `run.py --write-reference`.

The [README](https://github.com/horses-framework/horses3d/blob/master/Solver/test/MMS/README.md)
of `Solver/test/MMS` gives the orders of each case, the assumptions, what the cases do
not test, and how to change a case or add one.
