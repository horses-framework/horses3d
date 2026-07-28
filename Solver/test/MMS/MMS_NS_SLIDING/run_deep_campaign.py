"""
run_deep_campaign.py
====================
Unsteady-MMS convergence study for the sliding-mesh ALE fluxes on the MULTILAYER
CYLINDERDEEP meshes, whole mesh rotating about z (full rotation: radius > domain,
no hanging nodes / interface). Uses the python-generated MMS ProblemFile
(SETUP_MMS), which writes the L2 error to mms_l2_error.dat at the final step.

Every case is FULL ROTATION with angle = omega*dt (omega=0 => angle 0 => the
ALE fluxes reduce to the base scheme). Fixed dt / n_steps across all cases so
the convergence slopes in P and h reflect the (omega-stressed) spatial operator
with temporal error held below it.

Cases run in PARALLEL, each in its own isolated working directory (so the shared
mms_l2_error.dat can't clobber). Concurrency = NPROC / OMP_PER_CASE.

Output: errors_deep.csv  (mesh, P, omega, dt, n_steps, ndof, L2_error)
"""
import csv, os, shutil, subprocess, sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

HERE       = Path(__file__).resolve().parent
BINARY     = (HERE / "horses3d.ns").resolve()        # symlink created by ./configure
SETUP      = (HERE / "SETUP").resolve()               # holds libproblemfile_ns.so (build: make in SETUP)
MESH_DIR   = (HERE / "MESH").resolve()
TEMPLATE   = (HERE / "control_template_deep.control")
ERRORS_CSV = HERE / "errors_deep.csv"
WORK       = HERE / "WORK_deep"

# ── study parameters ─────────────────────────────────────────────────────────
DT       = 1.0e-5          # CFL-safe at omega=256 (calibrate in smoke test)
N_STEPS  = 40              # few steps; error measurable, temporal err negligible
OMEGAS   = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256]
MESHFILE = {"deep24": "MESH/CYLINDERDEEP_24_mesh.h5",
            "deep36": "MESH/CYLINDERDEEP_36_mesh.h5",
            "deep48": "MESH/CYLINDERDEEP_48_mesh.h5"}
OMP_PER_CASE = int(os.environ.get("SWEEP_OMP", "4"))
MAXPAR   = int(os.environ.get("SWEEP_MAXPAR",
               str(max(1, (os.cpu_count() or 4) // OMP_PER_CASE))))
# Optional comma-separated mesh-key filter (e.g. "deep24,deep36") so the heavy
# deep48 cases can be run in a separate low-concurrency phase to stay within RAM.
_MI = os.environ.get("SWEEP_MESHES", "").strip()
MESH_INCLUDE = set(_MI.split(",")) if _MI else None
_PI = os.environ.get("SWEEP_PS", "").strip()
P_INCLUDE = set(int(p) for p in _PI.split(",")) if _PI else None

# FULL grid: every (mesh, P, omega) combination — 3 x 3 x 10 = 90 cases.
# Resume-skip (errors_deep.csv) means already-run cases are not repeated.
CASES = [(m, P, om) for m in ("deep24", "deep36", "deep48")
         for P in (1, 2, 3, 4) for om in OMEGAS]

SLIDING = ("#define SlidingMesh\n  center = [0.0d0, 0.d0]\n  radius = 2.0d0\n"
           "  angle = {ANGLE}\n  rotation axis = z\n#end")


def case_name(mesh, P, om):
    return f"deep_{mesh}_P{P}_om{om:g}"


def load_done():
    d = set()
    if ERRORS_CSV.exists():
        for r in csv.DictReader(open(ERRORS_CSV)):
            d.add((r["mesh"], int(r["P"]), float(r["omega"])))
    return d


def run_one(mesh, P, om):
    name = case_name(mesh, P, om)
    wd = WORK / name
    if wd.exists():
        shutil.rmtree(wd)
    (wd / "RESULTS").mkdir(parents=True)
    os.symlink(SETUP, wd / "SETUP")
    # Private per-case MESH dir: symlink ONLY the needed .h5 so the solver's
    # generated .bmesh/.hmesh sidecars stay isolated. (A shared MESH symlink
    # makes concurrent same-mesh cases clobber the identical
    # CYLINDERDEEP_NN_mesh.bmesh -> intermittent "End of file" crash.)
    (wd / "MESH").mkdir()
    h5 = MESHFILE[mesh].split("/")[-1]
    os.symlink(os.path.realpath(MESH_DIR / h5), wd / "MESH" / h5)
    angle = om * DT
    block = SLIDING.format(ANGLE=f"{angle:.12e}".replace("e", "d"))
    ctl = (TEMPLATE.read_text()
           .replace("{MESH_FILE}", MESHFILE[mesh])
           .replace("{P}", str(P))
           .replace("{DT}", f"{DT:.12e}".replace("e", "d"))
           .replace("{N_STEPS}", str(N_STEPS))
           .replace("{SLIDING_BLOCK}", block))
    (wd / "case.control").write_text(ctl)
    env = dict(os.environ, OMP_NUM_THREADS=str(OMP_PER_CASE))
    with open(wd / "run.log", "w") as log:
        subprocess.run([str(BINARY), "case.control"], cwd=wd, stdout=log,
                       stderr=subprocess.STDOUT, env=env)
    ef = wd / "mms_l2_error.dat"
    if not ef.exists():
        return (mesh, P, om, None, "no mms_l2_error.dat (see run.log)")
    line = ef.read_text().split("\n")[0].split()
    # format: nelems, P, NDOF, L2_error, t_final
    try:
        ndof, l2 = line[2], line[3]
    except Exception:
        return (mesh, P, om, None, f"bad error line: {line}")
    return (mesh, P, om, (ndof, l2), None)


def main():
    dry = "--dry-run" in sys.argv
    WORK.mkdir(exist_ok=True)
    done = load_done()
    # Skip meshes that are not present: the CYLINDERDEEP_* HDF5 meshes are large
    # and are NOT committed. Generate them with HOPR (see README) and drop them in
    # Solver/test/TestMeshes/ (or this case's MESH/). Whatever is present runs.
    missing = [m for m, f in MESHFILE.items() if not (HERE / f).exists()]
    if missing:
        print(f"note: mesh(es) not found, skipping: {', '.join(sorted(missing))}")
    todo = [c for c in CASES if (c[0], c[1], float(c[2])) not in done
            and c[0] not in missing
            and (MESH_INCLUDE is None or c[0] in MESH_INCLUDE)
            and (P_INCLUDE is None or c[1] in P_INCLUDE)]
    if not todo and missing:
        print("nothing to run — no meshes available. See README for HOPR generation.")
    print(f"{len(todo)}/{len(CASES)} cases to run, {MAXPAR} in parallel "
          f"(OMP={OMP_PER_CASE}), {len(done)} already done")
    for c in todo:
        print("  queued:", case_name(*c))
    if dry:
        return
    new = not ERRORS_CSV.exists()
    fcsv = open(ERRORS_CSV, "a", newline="")
    w = csv.writer(fcsv)
    if new:
        w.writerow(["mesh", "P", "omega", "dt", "n_steps", "ndof", "L2_error"])
    with ThreadPoolExecutor(max_workers=MAXPAR) as ex:
        for mesh, P, om, res, err in ex.map(lambda c: run_one(*c), todo):
            nm = case_name(mesh, P, om)
            if res is None:
                print(f"  FAIL {nm}: {err}", flush=True)
                continue
            ndof, l2 = res
            w.writerow([mesh, P, om, DT, N_STEPS, ndof, l2]); fcsv.flush()
            print(f"  ok {nm}: L2={l2}", flush=True)
    fcsv.close()
    print("done.")


if __name__ == "__main__":
    main()
