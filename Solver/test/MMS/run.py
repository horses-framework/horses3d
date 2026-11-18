#!/usr/bin/env python3
"""
MMS convergence test for HORSES3D.

    python3 run.py CASE [--P 3 4] [--meshes 4 8] [--set "Reynolds number = 100"] [--check]

Generates and builds <case>/SETUP/ProblemFile.f90, then runs the solver for
every mesh and P with <case>/control_template.control, where run.py fills the
mesh, the solution file and the polynomial order. Without --P or --meshes it
uses the orders and meshes of the case in sweeps.json. --set changes or adds any
other key of the control file (repeat it, or separate pairs with ';').
The L2 errors go to <case>/errors.csv and the orders of convergence between
consecutive meshes are printed. Pairs already run with the same solver,
problem file and control file are skipped.

The meshes are periodic cubes [-1,1]^3 with N^3 elements, generated with the
cartesian-box script of horses-framework/horses3d-meshes (needs the Gmsh
Python API, pip install gmsh). Set HORSES3D_MESHES_DIR to a local clone of
that repository to use it instead of downloading the script, and GMSH_PYTHON
to a Python interpreter with gmsh if it is not this one.

It exits with an error if a run fails and, with --check and the committed
control file, if an error differs from <case>/reference.csv by more than 1e-3.
The problem file is built with the HORSES3D build settings taken from the
environment (e.g. COMPILER=gfortran ENABLE_THREADS=NO) and Solver/make.inc.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import generate

HERE = Path(__file__).resolve().parent
SWEEPS = json.loads((HERE / "sweeps.json").read_text())    # default orders and meshes of each case
MESH_DIR = HERE / "MESH"
MESH_GENERATOR = ("https://raw.githubusercontent.com/horses-framework/horses3d-meshes/"
                  "9fae9c071aeb5bb272fb2a40c835359e9f463015/cartesian-box/generate_cartesian_box_mesh.py")
RTOL = 1e-3
VOLUME = 8.0          # of the cube [-1,1]^3 made by mesh()


def mesh(N):
    out = MESH_DIR / f"cube_{N}.msh"
    if not out.exists():
        MESH_DIR.mkdir(exist_ok=True)
        if os.environ.get("HORSES3D_MESHES_DIR"):
            script = Path(os.environ["HORSES3D_MESHES_DIR"]) / "cartesian-box" / "generate_cartesian_box_mesh.py"
        else:
            script = MESH_DIR / "generate_cartesian_box_mesh.py"
            if not script.exists():
                urllib.request.urlretrieve(MESH_GENERATOR, script)
        python = os.environ.get("GMSH_PYTHON") or sys.executable
        subprocess.run([python, str(script), str(N), str(N), str(N), "--xlim", "-1", "1",
                        "--ylim", "-1", "1", "--zlim", "-1", "1", "-o", str(out)], check=True)
    return out


def build(case_dir, solver, home):
    setup = case_dir / "SETUP"
    generate.generate(case_dir.name)
    makefile = (home.parent / "utils" / "ProblemFile" / "Makefile").read_text()
    (setup / "Makefile").write_text(makefile.replace("@HORSES_HOME_DIR", str(home)))
    log = setup / "make.log"
    with open(log, "w") as f:
        rc = subprocess.run(["make", f"libproblemfile_{solver}.so"], cwd=setup, stdout=f, stderr=subprocess.STDOUT).returncode
    if rc:
        print("".join(log.read_text().splitlines(keepends=True)[-40:]), file=sys.stderr)
        sys.exit(f"Building the problem file failed, see {log}")


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
    """Number given to a key of the control file, or None if the key is not there."""
    for l in lines:
        if "=" in l and not l.startswith("!") and l.split("=", 1)[0].strip().lower() == key:
            return float(l.split("=", 1)[1].split("!")[0].lower().replace("d", "e"))
    return None


def run(case_dir, binary, lines, N, P, dt, t_final, keep):
    """Run one case and return [NDOF, volume, L2 total, L2 per variable...], or why it failed."""
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
        return f"return code {rc}, no mms_l2_error.dat"
    why = [text for bad, text in [(rc != 0, f"return code {rc}"),
                                  (not all(map(math.isfinite, v)), "error not finite"),
                                  (v[0] != N**3, f"{v[0]:g} elements"), (v[1] != P, f"P = {v[1]:g}"),
                                  (abs(v[3] - VOLUME) > 1e-8, f"volume {v[3]:.10g}"),
                                  (abs(v[-1] - t_final) > 0.5 * dt, f"final time {v[-1]:g}, expected {t_final:g}")] if bad]
    return ", ".join(why) if why else [int(v[2])] + v[3:-1]


def read_rows(path, header, strict=False):
    """Rows of an errors.csv written with this header; none if it does not exist or has
    another one, or an error with strict."""
    if not path.exists():
        if strict:
            sys.exit(f"{path} not found")
        return []
    with open(path) as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != header:
            if strict:
                sys.exit(f"{path} is not an errors.csv of this case")
            return []
        return [{k: v if k == "build" else int(v) if k in ("N", "P", "NDOF") else float(v) for k, v in r.items()}
                for r in reader]


def write_rows(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)


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


def check(case_dir, rows, overridden):
    """Compare the errors with reference.csv; False if one differs by more than RTOL."""
    ok = True
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
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=[c for c in generate.CASES if (HERE / c).is_dir()])   # cases of this repository
    ap.add_argument("--P", type=int, nargs="+", help="polynomial orders (default: those of the case in sweeps.json)")
    ap.add_argument("--meshes", type=int, nargs="+", help="elements per side (default: those of the case in sweeps.json)")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="change a key of the control file; several pairs can be separated by ';'")
    ap.add_argument("--home", type=Path, default=HERE.parents[1], help="HORSES3D Solver directory")
    ap.add_argument("--check", action="store_true", help=f"also fail on a difference above {RTOL:g} with reference.csv")
    ap.add_argument("--keep", action="store_true", help="keep the solution files")
    ap.add_argument("--merge", nargs="+", type=Path, metavar="CSV",
                    help="do not run: merge these errors.csv files (e.g. of parallel jobs) into <case>/errors.csv and report")
    args = ap.parse_args()
    sweep = SWEEPS.get(args.case, SWEEPS["default"])
    args.P = args.P or sweep["P"]
    args.meshes = args.meshes or sweep["meshes"]

    overrides = [tuple(s.strip() for s in pair.split("=", 1)) for item in args.set for pair in item.split(";") if pair.strip()]
    if any(len(o) != 2 for o in overrides):
        sys.exit("--set expects KEY=VALUE")
    _, solver, variables = generate.CASES[args.case]
    case_dir = HERE / args.case
    lines = control(case_dir, overrides)
    dt, steps, final_time = (value(lines, k) for k in ("dt", "number of time steps", "final time"))
    if dt is None or steps is None:
        sys.exit("The control file must give dt and the number of time steps")
    t_final = dt * steps if final_time is None else min(dt * steps, final_time)    # the solver stops at the first
    title = f"dt = {dt:g}, t_final = {t_final:g}"
    if overrides:
        title += ", " + "; ".join(f"{k} = {v}" for k, v in overrides)
    header = ["build", "N", "P", "dt", "t_final", "NDOF", "volume", "L2_total"] + [f"L2_{v}" for v in variables]
    errors_csv = case_dir / "errors.csv"

    if args.merge:
        # One row per mesh and order, the last one if a pair was run more than once
        rows = list({(r["N"], r["P"]): dict(r, build="merged") for f in args.merge
                     for r in read_rows(f, header, strict=True)}.values())
        if len({(r["dt"], r["t_final"]) for r in rows}) > 1:
            sys.exit("The files to merge have different dt or final time")
        write_rows(errors_csv, header, rows)
        report(args.case, [r for r in rows if r["P"] in args.P and r["N"] in args.meshes], args.P, title)
        return

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
    rows = read_rows(errors_csv, header)

    failed = []
    for P in args.P:
        for N in args.meshes:
            if any((r["build"], r["N"], r["P"]) == (build_id, N, P) for r in rows):
                continue
            start = time.time()
            result = run(case_dir, binary, lines, N, P, dt, t_final, args.keep)
            if isinstance(result, str):
                failed.append((N, P))
                print(f"{args.case} N={N} P={P}: FAILED ({result}), see {case_dir}/CONTROL/N{N}_P{P}.log", flush=True)
                continue
            rows.append(dict(zip(header, [build_id, N, P, dt, t_final] + result)))
            write_rows(errors_csv, header, rows)
            print(f"{args.case} N={N} P={P}: L2 = {rows[-1]['L2_total']:.4e} ({time.time() - start:.0f} s)", flush=True)

    # The rows of this build go last, the ones plot.py draws
    rows = [r for r in rows if r["build"] != build_id] + [r for r in rows if r["build"] == build_id]
    write_rows(errors_csv, header, rows)
    rows = [r for r in rows if r["build"] == build_id and r["P"] in args.P and r["N"] in args.meshes]
    report(args.case, rows, args.P, title)
    for N, P in failed:
        print(f"FAIL: run N={N} P={P} failed")
    ok = not failed
    if args.check:
        ok = check(case_dir, rows, bool(overrides)) and ok
        print("PASS" if ok else "FAILED")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
