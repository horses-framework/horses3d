#!/usr/bin/env python3
"""
MMS convergence test for HORSES3D.

    python3 run.py CASE [--P 2 3] [--meshes 4 8 16] [--set "Reynolds number = 100"] [--check]

Generates and builds <case>/SETUP/ProblemFile.f90, then runs the solver for
every mesh and P with <case>/control_template.control, where run.py fills the
mesh, the solution file and the polynomial order. --set changes or adds any
other key of the control file (repeat it, or separate pairs with ';').
The L2 errors go to <case>/errors.csv and the orders of convergence between
consecutive meshes are printed. Pairs already run with the same solver,
problem file and control file are skipped.

The meshes are periodic cubes [-1,1]^3 with N^3 elements, generated with the
cartesian-box script of horses-framework/horses3d-meshes (needs the Gmsh
Python API, pip install gmsh). Set HORSES3D_MESHES to a local clone of that
repository to use it instead of downloading the script.

With --check, it exits with an error if a run failed or, with the committed
control file, an error differs from <case>/reference.csv by more than 1e-3.
The problem file is built with the HORSES3D build settings taken from the
environment (e.g. COMPILER=gfortran ENABLE_THREADS=NO) and Solver/make.inc.
"""
import argparse
import csv
import hashlib
import math
import os
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import generate

HERE = Path(__file__).resolve().parent
MESH_DIR = HERE / "MESH"
MESH_GENERATOR = ("https://raw.githubusercontent.com/horses-framework/horses3d-meshes/"
                  "9fae9c071aeb5bb272fb2a40c835359e9f463015/cartesian-box/generate_cartesian_box_mesh.py")
RTOL = 1e-3
VOLUME = 8.0          # of the cube [-1,1]^3 made by mesh()


def mesh(N):
    out = MESH_DIR / f"cube_{N}.msh"
    if not out.exists():
        MESH_DIR.mkdir(exist_ok=True)
        if os.environ.get("HORSES3D_MESHES"):
            script = Path(os.environ["HORSES3D_MESHES"]) / "cartesian-box" / "generate_cartesian_box_mesh.py"
        else:
            script = MESH_DIR / "generate_cartesian_box_mesh.py"
            if not script.exists():
                urllib.request.urlretrieve(MESH_GENERATOR, script)
        subprocess.run([sys.executable, str(script), str(N), str(N), str(N), "--xlim", "-1", "1",
                        "--ylim", "-1", "1", "--zlim", "-1", "1", "-o", str(out)], check=True)
    return out


def build(case_dir, solver, home):
    setup = case_dir / "SETUP"
    generate.generate(case_dir.name)
    makefile = (home.parent / "utils" / "ProblemFile" / "Makefile").read_text()
    (setup / "Makefile").write_text(makefile.replace("@HORSES_HOME_DIR", str(home)))
    with open(setup / "make.log", "w") as log:
        if subprocess.run(["make", f"libproblemfile_{solver}.so"], cwd=setup, stdout=log, stderr=subprocess.STDOUT).returncode:
            sys.exit(f"Building the problem file failed, see {setup / 'make.log'}")


def control(case_dir, overrides):
    """Lines of the control file template with the --set changes, without indentation:
    the control file reader truncates lines at 132 characters."""
    lines = [l.strip() for l in (case_dir / "control_template.control").read_text().splitlines()]
    header = next(i for i, l in enumerate(lines) if l.startswith("#define"))   # keys come before the boundaries
    for key, value in overrides:
        line = f"{key} = {value}"
        hits = [i for i in range(header) if "=" in lines[i] and not lines[i].startswith("!")
                and lines[i].split("=", 1)[0].strip().lower() == key.lower()]
        for i in hits:
            lines[i] = line
        if not hits:
            lines.insert(header, line)
            header += 1
    return lines


def value(lines, key):
    return next(float(l.split("=", 1)[1].split("!")[0].replace("d", "e")) for l in lines
                if "=" in l and l.split("=", 1)[0].strip().lower() == key)


def run(case_dir, binary, lines, N, P, dt, t_final, keep):
    """Run one case and return [NDOF, volume, L2 total, L2 per variable...], or None if it failed."""
    name = f"N{N}_P{P}"
    for d in ("CONTROL", "RESULTS", "MESH"):        # MU writes the boundary mesh to ./MESH
        (case_dir / d).mkdir(exist_ok=True)
    text = "\n".join(lines) + "\n"
    for token, val in [("{MESH_FILE}", os.path.relpath(mesh(N), case_dir)),
                       ("{SOLUTION_FILE}", f"RESULTS/{name}.hsol"), ("{P}", str(P))]:
        text = text.replace(token, val)
    (case_dir / "CONTROL" / f"{name}.control").write_text(text)

    dat = case_dir / "mms_l2_error.dat"
    dat.unlink(missing_ok=True)
    with open(case_dir / "CONTROL" / f"{name}.log", "w") as log:
        rc = subprocess.run([str(binary), f"CONTROL/{name}.control"], cwd=case_dir, stdout=log, stderr=subprocess.STDOUT).returncode
    if not keep:
        for f in (case_dir / "RESULTS").glob(f"{name}*"):
            f.unlink()
    try:
        # nelems, P, NDOF, volume, total, each variable, time
        v = [float(s) for s in dat.read_text().split(",")]
    except (OSError, ValueError):
        return None
    if (rc != 0 or not all(map(math.isfinite, v)) or v[0] != N**3 or v[1] != P
            or abs(v[3] - VOLUME) > 1e-8 or abs(v[-1] - t_final) > 0.5 * dt):
        return None
    return [int(v[2])] + v[3:-1]


def orders(rows, P):
    """Rows of one P sorted by mesh, each with the orders against the previous mesh."""
    rows = sorted((r for r in rows if r["P"] == P), key=lambda r: r["N"])
    out = []
    for prev, r in zip([None] + rows, rows):
        eoc = {}
        if prev:
            for k in r:
                if k.startswith("L2_") and prev[k] > 0 and r[k] > 0:    # some variables are exact (CAA rho)
                    eoc[k] = math.log(prev[k] / r[k]) / math.log(r["N"] / prev["N"])
        out.append((r, eoc))
    return out


def report(case, rows, ps, title):
    variables = ["total"] + generate.CASES[case][2]
    lines = [f"MMS {case}: {title}", "",
             "| P | N | NDOF | L2 total | " + " | ".join(f"order {v}" for v in variables) + " |",
             "|---" * (4 + len(variables)) + "|"]
    for P in ps:
        for r, eoc in orders(rows, P):
            cells = [f"{eoc[f'L2_{v}']:.2f}" if f"L2_{v}" in eoc else "" for v in variables]
            lines.append(f"| {P} | {r['N']} | {r['NDOF']} | {r['L2_total']:.4e} | " + " | ".join(cells) + " |")
    text = "\n".join(lines) + "\n"
    print("\n" + text)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
            f.write(text + "\n")


def check(case_dir, rows, failed, overridden):
    ok = not failed
    for N, P in failed:
        print(f"FAIL: run N={N} P={P} failed")
    reference = case_dir / "reference.csv"
    if overridden or not reference.exists():
        print("No comparison with reference.csv: the control file was changed with --set" if overridden
              else "No reference.csv")
    else:
        ref = {(int(r["N"]), int(r["P"])): r for r in csv.DictReader(open(reference))}
        for r in rows:
            if (r["N"], r["P"]) not in ref:
                print(f"N={r['N']} P={r['P']}: no reference")
                continue
            diff = max(abs(r[k] - float(v)) / max(float(v), 1e-300) for k, v in ref[(r["N"], r["P"])].items()
                       if k.startswith("L2_"))
            if diff > RTOL:
                ok = False
                print(f"FAIL: N={r['N']} P={r['P']} differs from the reference by {diff:.1e}")
    print("PASS" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=generate.CASES)
    ap.add_argument("--P", type=int, nargs="+", default=[2, 3], help="polynomial orders (default 2 3)")
    ap.add_argument("--meshes", type=int, nargs="+", default=[4, 8], help="elements per side (default 4 8)")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="change a key of the control file; several pairs can be separated by ';'")
    ap.add_argument("--home", type=Path, default=HERE.parents[1], help="HORSES3D Solver directory")
    ap.add_argument("--check", action="store_true", help=f"fail on a failed run or a difference above {RTOL:g} with reference.csv")
    ap.add_argument("--keep", action="store_true", help="keep the solution files")
    args = ap.parse_args()

    overrides = [tuple(s.strip() for s in pair.split("=", 1)) for item in args.set for pair in item.split(";") if pair.strip()]
    if any(len(o) != 2 for o in overrides):
        sys.exit("--set expects KEY=VALUE")
    _, solver, variables = generate.CASES[args.case]
    case_dir = HERE / args.case
    home = args.home.resolve()
    binary = home / "bin" / f"horses3d.{solver}"
    if not binary.exists():
        sys.exit(f"{binary} not found: build HORSES3D first (make {solver})")
    build(case_dir, solver, home)
    h = hashlib.sha256()
    for f in (binary, case_dir / "SETUP" / "ProblemFile.f90", case_dir / "control_template.control"):
        h.update(f.read_bytes())
    h.update(repr(overrides).encode())
    build_id = h.hexdigest()[:12]
    lines = control(case_dir, overrides)
    dt = value(lines, "dt")
    t_final = dt * value(lines, "number of time steps")

    header = ["build", "N", "P", "dt", "t_final", "NDOF", "volume", "L2_total"] + [f"L2_{v}" for v in variables]
    errors_csv = case_dir / "errors.csv"
    rows = []
    if errors_csv.exists():
        with open(errors_csv) as f:
            reader = csv.DictReader(f)
            if reader.fieldnames == header:
                rows = [{k: v if k == "build" else int(v) if k in ("N", "P", "NDOF") else float(v) for k, v in r.items()}
                        for r in reader]

    failed = []
    for P in args.P:
        for N in args.meshes:
            if any((r["build"], r["N"], r["P"]) == (build_id, N, P) for r in rows):
                continue
            start = time.time()
            result = run(case_dir, binary, lines, N, P, dt, t_final, args.keep)
            if result is None:
                failed.append((N, P))
                print(f"{args.case} N={N} P={P}: FAILED, see {case_dir}/CONTROL/N{N}_P{P}.log", flush=True)
                continue
            rows.append(dict(zip(header, [build_id, N, P, dt, t_final] + result)))
            with open(errors_csv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=header)
                w.writeheader()
                w.writerows(rows)
            print(f"{args.case} N={N} P={P}: L2 = {rows[-1]['L2_total']:.4e} ({time.time() - start:.0f} s)", flush=True)

    rows = [r for r in rows if r["build"] == build_id and r["P"] in args.P and r["N"] in args.meshes]
    title = f"dt = {dt:g}, t_final = {t_final:g}"
    if overrides:
        title += ", " + "; ".join(f"{k} = {v}" for k, v in overrides)
    report(args.case, rows, args.P, title)
    if args.check and not check(case_dir, rows, failed, bool(overrides)):
        sys.exit(1)


if __name__ == "__main__":
    main()
