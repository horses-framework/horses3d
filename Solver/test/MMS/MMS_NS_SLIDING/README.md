# MMS convergence test — sliding-mesh ALE fluxes

Method-of-Manufactured-Solutions convergence study for the **rotating-mesh (ALE)
fluxes**. The whole mesh rotates rigidly about the cylinder axis, so this isolates
the ALE flux discretisation from the non-conforming mortar/hanging-node machinery.

Rotation rate `omega` is swept (0, 1, 2, ... 256) at fixed `dt` and step count, so
p- and h-convergence can be compared across rotation regimes. `omega = 0` must
reproduce the non-rotating scheme exactly (reduction check).

## Layout

| File | Purpose |
|---|---|
| `SETUP/generate_problemfile_NS.py` | SymPy generator: manufactured solution -> source term + exact-state routines. **Edit the MS here**, then re-run to regenerate `ProblemFile.f90`. |
| `SETUP/ProblemFile.f90` | Generated. Provides the MS source term, the Dirichlet boundary state and the L2 error (written to `mms_l2_error.dat`). |
| `control_template_deep.control` | Control-file template (tokens `{MESH_FILE} {P} {DT} {N_STEPS} {SLIDING_BLOCK}`). |
| `run_deep_campaign.py` | Runs the (mesh, P, omega) grid in parallel, collects `errors_deep.csv`. |
| `plot_deep_grid.py` | Convergence plots from `errors_deep.csv`. |

## Running

```bash
./configure                       # from Solver/ : creates SETUP/Makefile + horses3d.ns link
cd Solver/test/MMS/MMS_NS_SLIDING/SETUP && make   # builds libproblemfile_ns.so
cd .. && python3 run_deep_campaign.py             # add --dry-run to preview
python3 plot_deep_grid.py
```

## Meshes

The study uses HOPR cylinder-annulus meshes with **several element layers along the
rotation axis** (so the solution's axial structure is resolved by interior faces
rather than imposed by the boundary conditions):

| name | elements | size | note |
|---|---|---|---|
| `CYLINDERDEEP_24_mesh.h5` | 24 x 8 x 4 | 4.6 MB | p-convergence + omega sweep |
| `CYLINDERDEEP_36_mesh.h5` | 36 x 12 x 6 | 16 MB | h-convergence |
| `CYLINDERDEEP_48_mesh.h5` | 48 x 16 x 8 | 37 MB | h-convergence |

> **:warning: LARGE BINARIES — REVIEW BEFORE MERGING UPSTREAM.**
> These three meshes total **~59 MB** and are committed here so the convergence study
> can be reproduced and checked on the fork. They are far larger than any other mesh
> in `Solver/test` (the next largest are a few MB). **Do not carry them into `main`
> as-is** — either drop them and regenerate with HOPR (parameters below), keep only
> `CYLINDERDEEP_24` (4.6 MB, enough for p-convergence and the omega sweep), or move
> them to Git LFS / an external mesh store.

The runner skips any mesh that is absent, so removing them degrades the study
gracefully rather than breaking it. HOPR generation parameters: `Mode=11`,
`WhichMapping=4` (full cylinder), `R_0=0.5`, `R_INF=1.0`, `nElems=(/N, N/3, K/)`,
`useCurveds=T`, `BoundaryOrder=6`, boundary names `SIDEA` (inner) / `SIDED` (outer)
/ `SIDEB`,`SIDEC` (axial ends).

Boundaries are **user-defined Dirichlet** (exact MS state), not periodic: the sliding
rebuild reconstructs the face topology and does not re-pair periodic faces.

## Notes on correctness of the measurement

- The error is the **volume-averaged L2 (RMS)** norm, `sqrt(integral(|e|^2)/volume)`.
  The volume is accumulated in the same quadrature loop. Normalising by the number
  of elements instead would change the normalisation between meshes and corrupt the
  h-convergence slope.
- The viscosity in the generated source term follows **Sutherland's law**, matching
  `get_laminar_mu_kappa` in the solver. A constant-viscosity source term would not
  be consistent with the solver's viscous flux.
- Node type in the generator (`NODE_TYPE`) must match `discretization nodes` in the
  control file — both are `Gauss` here — because it selects the quadrature weights
  used for the error integral.
