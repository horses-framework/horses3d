# MMS Convergence Verification Framework — HORSES3D

Automated Method of Manufactured Solutions (MMS) convergence study for one
solver per directory.  The three solver directories (NS, iNS, MU) share an
identical layout and workflow.

---

## Directory structure

Each solver folder (`NS/`, `iNS/`, `MU/`) is laid out identically.  Mesh
files are **not** stored per-solver; they live in a shared `TestMeshes/Cubes/`
directory two levels above the solver folder so all three solvers reuse the
same meshes.

```
<repo-root>/
├── bin/
│   └── horses3d.{ns,ins,mu}                  # compiled solver binaries
└── <parent>/
    ├── TestMeshes/
    │   └── Cubes/
    │       └── mesh<N>.h5                    # structured hex meshes (N^3 elements)  <- shared by all solvers
    └── <solver-group>/
        └── <solver>/                         # NS / iNS / MU
            ├── SETUP/
            │   ├── generate_problemfile_<solver>.py   # symbolic source-term generator  <- edit here
            │   └── Makefile                           # builds the problem-file shared library
            ├── CONTROL/
            │   ├── mms_N<N>_P<P>.control             # auto-generated per run  (do not edit)
            │   └── mms_N<N>_P<P>.log                 # solver log per run      (do not edit)
            ├── RESULTS/
            │   └── MMS_N<N>_P<P>.hsol                # solver output           (do not edit)
            ├── control_template.control               # solver settings  <- edit here
            ├── run_convergence_<solver>.py            # study config and orchestration  <- edit here
            ├── plot_convergence_<solver>.py           # plot styling (optional edit)
            └── errors.csv                             # accumulated results (auto-created)
```

The default paths in `run_convergence_<solver>.py` are:

| Constant | Default | Resolves to |
|---|---|---|
| `HORSES_BINARY` | `../../../bin/horses3d.<ext>` | `<repo-root>/bin/horses3d.<ext>` |
| `MESH_DIR`      | `./../../TestMeshes/Cubes`    | shared mesh directory two levels up |

Both are resolved relative to the script's own directory, so the runner can be
launched from any working directory.

---

## Quick start

### 1 — Generate the problem file

```bash
cd <solver>/SETUP
# Edit Section 1 and Section 2 of generate_problemfile_<solver>.py (see below)
python3 generate_problemfile_<solver>.py
make
```

This writes `ProblemFile.f90` with the MMS source terms and L2 error routine
hard-coded, then `make` builds it into a shared library
(e.g. `libproblemfile_ns.so`) that the existing HORSES3D binary loads at
runtime.  The HORSES3D executable itself is not relinked.

> **Always rebuild the problem-file library after changing the manufactured
> solution or physical parameters.**

### 2 — Prepare meshes

Meshes are **shared across all three solvers** (NS, iNS, MU) and live in
`TestMeshes/Cubes/` two levels above each solver folder.  Place structured
hex meshes there named `mesh<N>.h5`, where `N` is the number of elements per
side (e.g. `mesh3.h5`, `mesh4.h5`, `mesh5.h5`, `mesh6.h5`).

If your meshes are in a different location, update `MESH_DIR` in **Section 1**
of `run_convergence_<solver>.py`.  The path is resolved relative to the
script's own directory.

### 3 — Configure the run

Open `run_convergence_<solver>.py` and edit **Section 1** (paths) and
**Section 2** (study parameters):

| Parameter | Description |
|---|---|
| `HORSES_BINARY` | Path to compiled solver binary |
| `MESH_SIZES` | List of `N` values (elements per side) |
| `P_VALUES` | Polynomial orders to sweep |
| `DT` | Fixed time step (see note on time floor below) |
| `T_FINAL` | End time |

Additional solver settings (Riemann solver, viscosity model, etc.) live in
`control_template.control`.  Only these tokens are filled automatically at
runtime:

| Token | Replaced by |
|---|---|
| `{MESH_FILE}` | Path to `mesh<N>.h5` |
| `{SOLUTION_FILE}` | Path to output `.hsol` |
| `{P}` | Polynomial order |
| `{DT}` | Time step |
| `{N_STEPS}` | `max(1, ceil(T_FINAL / DT))` — ceiling-based to guarantee `T_FINAL` is reached |

### 4 — Run

```bash
python3 run_convergence_<solver>.py
python3 run_convergence_<solver>.py --dry-run   # preview without running
```

Results are appended to `errors.csv` one row per case.  The script is
**crash-safe**: if interrupted, re-running it skips already-completed
`(nelems, P)` pairs.

### 5 — Plot

```bash
python3 plot_convergence_<solver>.py
python3 plot_convergence_<solver>.py --csv path/to/errors.csv --outdir path/to/out
```

Outputs `convergence.pdf` and `convergence.png` with three panels:

- **Left** — L2 error vs sqrt(NDOF) (combined h/p view)
- **Centre** — h-convergence per p (log-log, best-fit slope vs theory p+1)
- **Right** — p-convergence per mesh (empirical reduction factors annotated)

---

## Expected convergence rates

For a smooth manufactured solution and a DG scheme of polynomial order p:

- **h-convergence** (fixed p, refine mesh): slope ~= p+1 in log(error) vs log(h)
- **p-convergence** (fixed mesh, increase p): exponential reduction

**Even-odd oscillations** in h-slopes are a known property of the DGSEM on
Gauss-Lobatto nodes (super-convergence alternates between p+1 and p+2 with the
parity of p).  This is not a defect — it may appear or disappear depending on
the symmetry of the chosen manufactured solution.

**Time-integration floor**: if the time step is too large, temporal error
dominates and curves flatten.  Faded markers in the plot flag cases where
consecutive p-errors differ by less than a factor of 2.  Fix: reduce `DT`
by 10x and re-run.

---

## Changing the manufactured solution

Open `generate_problemfile_<solver>.py` and edit **Section 2**.  Any smooth
SymPy expression that is periodic on `[-1,1]^3` is valid.  Constraints by
solver:

### Compressible NS

- **Pressure must be strictly positive everywhere** — add a constant offset.
  This avoids NaN in `T = p / (rho * R)`.
- **Density must be strictly positive** — same pattern.
- The velocity does **not** need to be divergence-free.  However, a
  non-solenoidal choice adds a large density source term that can make
  the problem stiff and require a smaller `DT`.

```python
# Variable-density example
rho = RHO0 * (1 + 0.1*sp.sin(pi*x)*sp.sin(pi*y)*sp.sin(pi*z)*sp.cos(t))
u   =  sp.sin(pi*x)*sp.cos(pi*y)*sp.cos(pi*z)*sp.sin(t)
v   =  sp.cos(pi*x)*sp.sin(pi*y)*sp.cos(pi*z)*sp.sin(t)
w   = -2*sp.cos(pi*x)*sp.cos(pi*y)*sp.sin(pi*z)*sp.sin(t)
p   =  1 + 0.1*sp.sin(pi*x)*sp.sin(pi*y)*sp.sin(pi*z)*sp.cos(t)
```

### Incompressible NS (iNS)

The velocity **must satisfy** `du/dx + dv/dy + dw/dz = 0` identically.
The easiest construction is to derive `u = curl(A)` from a vector potential `A`.
The default field satisfies this condition exactly.

```python
# Verify divergence-free: d/dx(sin*cos*cos) + d/dy(cos*sin*cos) + d/dz(-2cos*cos*sin)
# = pi*cos*cos*cos - pi*cos*cos*cos + 0 = 0  (check!)
u =  sp.sin(pi*x)*sp.cos(pi*y)*sp.cos(pi*z)*sp.sin(t)
v =  sp.cos(pi*x)*sp.sin(pi*y)*sp.cos(pi*z)*sp.sin(t)
w = -2*sp.cos(pi*x)*sp.cos(pi*y)*sp.sin(pi*z)*sp.sin(t)
```

### Multiphase (MU)

The phase field must satisfy `|phi| <= 1` pointwise.  Use an amplitude `|A| < 1`:

```python
phi =  0.5*sp.sin(pi*x)*sp.sin(pi*y)*sp.sin(pi*z)*sp.cos(t)
```

---

## Adding a body force (e.g. gravity)

To verify a module that adds gravity or another body force:

1. Locate the source-term derivation block in `generate_problemfile_<solver>.py`
   (look for the comment `Source term derivation`).

2. Add the body-force contribution to the appropriate symbolic equation.
   For NS gravity `g = (0, 0, -g_val)` acting on the momentum equations,
   the z-momentum source gains an extra `-rho * g_val`:

```python
# Inside the source-term derivation block, after the existing S expressions:
g_val = 9.81   # m/s^2 or non-dimensional equivalent
# Subtract from the z-momentum source (HORSES3D sign convention: S appears on RHS)
S_mom_z = S_mom_z - rho * g_val
```

3. Re-run the generator and rebuild the library (`make`).  The Fortran output
   will contain the correct modified source term automatically.  No other file
   needs changing.

---

## Verifying a new physics module

General procedure for any new term (rotating frame, immersed boundary,
turbulence model, surface tension, etc.):

1. Identify which equations and which terms are modified.
2. Add the corresponding SymPy expressions to the generator's source-term block.
3. Regenerate `ProblemFile.f90` and rebuild the library (`make`).
4. Run the sweep.

**Interpreting results:**

| Observation | Likely cause |
|---|---|
| All slopes match p+1 | New module correctly implemented |
| Slopes systematically low across all (N, p) | Bug in new module or source term mismatch |
| Only finest-mesh / highest-p points low (faded markers) | Time-integration floor — reduce DT |
| Slopes erratic, no clear trend | Check physical parameter consistency between generator and control file |

---

## Keeping physical parameters consistent

The parameters in `generate_problemfile_<solver>.py` (e.g. `RHO0`, `GAMMA`,
Reynolds number) **must match** those set in `control_template.control`.

A mismatch does **not** produce a runtime error, but the source terms will be
wrong and convergence will typically fail or give incorrect slopes.

---

## SymPy simplification options

In **Section 3** of the generator:

| `SIMPLIFY_STRATEGY` | Behaviour |
|---|---|
| `"auto"` (default) | Try `simplify()`, fall back to `trigsimp` after timeout |
| `"simplify"` | Always use `simplify()` — thorough but slow |
| `"trigsimp"` | Always use `trigsimp()` — fast; best for pure trig solutions |

The `SIMPLIFY_TIMEOUT` variable (default 60 s) controls the fallback.
If you see `simplify timed out` on stderr, this is normal for complex
expressions.  Set `SIMPLIFY_STRATEGY = "trigsimp"` to skip the timeout.

---

## File reference — what to edit

| File | Edit? | What to change |
|---|---|---|
| `SETUP/generate_problemfile_*.py` | **Yes** — Sections 1-2 | Physical params, manufactured solution, new physics terms |
| `run_convergence_*.py` | **Yes** — Sections 1-2 | Paths, mesh sizes, P values, DT, T_FINAL |
| `control_template.control` | **Yes** | Solver-level settings (Riemann solver, BCs, viscosity model, etc.) |
| `plot_convergence_*.py` | Optional — Section 1 | Colour scheme, reference slopes, floor threshold |
| `SETUP/Makefile` | Rarely | Compilation flags |
| `errors.csv` | **Never** | Auto-generated; delete to force a fresh sweep |
| `CONTROL/*.control` | **Never** | Auto-generated from template |
| `CONTROL/*.log` | **Never** | Auto-generated solver logs |

---

## Common issues

**`ERROR: template not found`**
Run the script from the directory containing `control_template.control`, or
update `TEMPLATE_FILE` in Section 1.

**`ERROR: could not read mms_l2_error.dat`**
The solver ran but `UserDefinedFinalize` did not write output.  Check the log
in `CONTROL/mms_N<N>_P<P>.log`.  Most common cause: the problem-file library
was built without the correct `#ifdef` flag (`NAVIERSTOKES`, `CAHNHILLIARD`,
etc.) — verify the flag in `SETUP/Makefile` matches the solver you are running.

**Slopes consistently below p+1 for all cases**
Reduce `DT` by 10x and re-run.  Also verify that physical parameters in the
generator match the control file.

**Negative pressure / NaN during NS run**
The manufactured pressure is not strictly positive.  Increase the constant
offset in `p = 1 + ...` or reduce the oscillation amplitude. 

**`simplify timed out`** (stderr during generator)
Normal for complex expressions.  Set `SIMPLIFY_STRATEGY = "trigsimp"`.

---

## Dependencies

| Tool | Purpose | Install |
|---|---|---|
| Python >= 3.9 | Runner and plotter | — |
| SymPy | Symbolic source-term derivation | `pip install sympy` |
| NumPy | Quadrature weight computation | `pip install numpy` |
| Matplotlib | Convergence plots | `pip install matplotlib` |
| HORSES3D binary | Solver execution | Build from source |
