#!/usr/bin/env python3
"""
MMS convergence test for HORSES3D.

    python3 run.py CASE [--P 3 4] [--meshes 4 8] [--set "Reynolds number = 100"] [--check]

Generates and builds <case>/SETUP/ProblemFile.f90, then runs the solver for
every mesh and P with <case>/control_template.control, where run.py fills the
mesh, the solution file and the polynomial order. Without --P or --meshes it
uses the orders and meshes of the case in sweeps.json. --set changes or adds any
other key of the control file (repeat it, or separate pairs with ';'); a key that
is not in the template is added with a warning, since HORSES3D ignores the keys
it does not read. The L2 errors go to <case>/errors.csv and the orders of
convergence between consecutive meshes are printed. Pairs already run with the
same solver, problem file, build settings and control file are skipped, unless
--force is given.

The meshes are periodic cubes [-1,1]^3 with N^3 elements, generated with the
cartesian-box script of horses-framework/horses3d-meshes (needs the Gmsh
Python API, pip install gmsh). Set HORSES3D_MESHES_DIR to a local clone of
that repository to use it instead of downloading the script, and GMSH_PYTHON
to a Python interpreter with gmsh if it is not this one.

It exits with an error if a run fails and, with --check, if an error differs by
more than a relative 1e-3 from the row of the same mesh and P in
<case>/reference.csv, or if a run of the default sweep has no row there. Other
runs, and runs with --set, are not compared. --write-reference writes the errors
of the default sweep, run with the committed control file, to <case>/reference.csv.
The problem file is built with the HORSES3D build settings taken from the
environment (e.g. COMPILER=gfortran ENABLE_THREADS=NO) and Solver/make.inc.
Only one run.py at a time can work on a case, since its runs share files.
"""
import argparse
import csv
import fcntl
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
BUILD_SETTINGS = ("PLATFORM", "COMPILER", "MODE", "COMM", "ENABLE_THREADS", "EXTRA_FFLAGS")    # read by the problem-file Makefile
RTOL = 1e-3
VOLUME = 8.0          # of the cube [-1,1]^3 made by mesh()


def mesh(N):
    """The cube mesh with N^3 elements, made if it is not in MESH/. Every case uses MESH/,
    so making a mesh is locked, and a file appears there only when it is complete."""
    out = MESH_DIR / f"cube_{N}.msh"
    MESH_DIR.mkdir(exist_ok=True)
    with open(MESH_DIR / ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not out.exists():
            if os.environ.get("HORSES3D_MESHES_DIR"):
                script = Path(os.environ["HORSES3D_MESHES_DIR"]) / "cartesian-box" / "generate_cartesian_box_mesh.py"
            else:
                script = MESH_DIR / "generate_cartesian_box_mesh.py"
                if not script.exists():
                    urllib.request.urlretrieve(MESH_GENERATOR, MESH_DIR / "generator.tmp")
                    os.replace(MESH_DIR / "generator.tmp", script)
            tmp = MESH_DIR / f"cube_{N}.tmp.msh"        # gmsh takes the format from the extension
            python = os.environ.get("GMSH_PYTHON") or sys.executable
            subprocess.run([python, str(script), str(N), str(N), str(N), "--xlim", "-1", "1",
                            "--ylim", "-1", "1", "--zlim", "-1", "1", "-o", str(tmp)], check=True)
            os.replace(tmp, out)
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
    """Lines of the control file template with the --set changes, without indentation
    (the control file reader truncates lines at 132 characters), and the keys that
    were not in the template."""
    lines = [l.strip() for l in (case_dir / "control_template.control").read_text().splitlines()]
    header = next(i for i, l in enumerate(lines) if l.startswith("#define"))   # keys come before the boundaries
    added = []
    for key, value in overrides:
        line = f"{key} = {value}"
        hits = [i for i in range(header) if "=" in lines[i] and not lines[i].startswith("!")
                and lines[i].split("=", 1)[0].strip().lower() == key.lower()]
        for i in hits:
            lines[i] = line
        if not hits:
            lines.insert(header, line)
            header += 1
            added.append(key)
    return lines, added


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


def write_reference(path, header, rows):
    """reference.csv: N, P and the L2 errors to 11 significant digits, sorted by P and N."""
    cols = ["N", "P"] + [k for k in header if k.startswith("L2_")]
    with open(path, "w", newline="") as f:
        f.write(",".join(cols) + "\r\n")
        for r in sorted(rows, key=lambda r: (r["P"], r["N"])):
            f.write(",".join([str(r["N"]), str(r["P"])] + ["%.10e" % r[k] for k in cols[2:]]) + "\r\n")


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


def summary(text, markdown=None):
    """Print text and add it, or its markdown form, to the job summary of a GitHub workflow."""
    print(text)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
            f.write((markdown or text) + "\n")


def report(case, rows, ps, title):
    variables = generate.CASES[case][2]
    lines = [f"MMS {case}: {title}", "",
             "| P | N | NDOF | L2 total | order total | lowest order | " + " | ".join(f"order {v}" for v in variables) + " |",
             "|---" * (6 + len(variables)) + "|"]
    for P in ps:
        for r, eoc in orders(rows, P):
            each = {v: eoc[f"L2_{v}"] for v in variables if f"L2_{v}" in eoc}    # exact variables have no order
            low = min(each, key=each.get) if each else None
            cells = ([f"{eoc['L2_total']:.2f}" if "L2_total" in eoc else "", f"{each[low]:.2f} ({low})" if low else ""]
                     + [f"{each[v]:.2f}" if v in each else "" for v in variables])
            lines.append(f"| {P} | {r['N']} | {r['NDOF']} | {r['L2_total']:.4e} | " + " | ".join(cells) + " |")
    summary("\n" + "\n".join(lines) + "\n")


def check(case_dir, rows, overridden, sweep):
    """Compare the errors with the rows of the same mesh and P in reference.csv. A run of
    the default sweep without a row there fails, and other runs are not compared.
    Returns whether the check passed, and the lines to report."""
    reference = case_dir / "reference.csv"
    if overridden:
        return True, ["Not compared with reference.csv: the control file was changed with --set"]
    ref = {}
    if reference.exists():
        with open(reference) as f:
            ref = {(int(r["N"]), int(r["P"])): r for r in csv.DictReader(f)}
    if not ref:
        return False, [f"FAIL: {case_dir.name}/reference.csv is missing or empty"]
    out, compared, missing = [], 0, []
    for r in rows:
        if (r["N"], r["P"]) not in ref:
            if r["N"] in sweep["meshes"] and r["P"] in sweep["P"]:
                out.append(f"FAIL: N={r['N']} P={r['P']} is in the default sweep but not in reference.csv")
            else:
                missing.append(f"N={r['N']} P={r['P']}")
            continue
        compared += 1
        diff = max(abs(r[k] - float(v)) / max(float(v), 1e-300) for k, v in ref[(r["N"], r["P"])].items()
                   if k.startswith("L2_"))
        if diff > RTOL:
            out.append(f"FAIL: N={r['N']} P={r['P']} differs from the reference by {diff:.1e}")
    out.append(f"Compared {compared} of {len(rows)} runs with reference.csv, relative tolerance {RTOL:g}"
               + (f"; no reference row for {', '.join(missing)}" if missing else ""))
    return not any(l.startswith("FAIL") for l in out), out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=[c for c in generate.CASES if (HERE / c).is_dir()])   # cases of this repository
    ap.add_argument("--P", type=int, nargs="+", help="polynomial orders (default: those of the case in sweeps.json)")
    ap.add_argument("--meshes", type=int, nargs="+", help="elements per side (default: those of the case in sweeps.json)")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="change a key of the control file; several pairs can be separated by ';'")
    ap.add_argument("--home", type=Path, default=HERE.parents[1], help="HORSES3D Solver directory")
    ap.add_argument("--check", action="store_true", help=f"also fail on a difference above {RTOL:g} with reference.csv")
    ap.add_argument("--force", action="store_true", help="run again the pairs already in errors.csv for this build")
    ap.add_argument("--keep", action="store_true", help="keep the solution files")
    ap.add_argument("--merge", nargs="+", type=Path, metavar="CSV",
                    help="do not run: merge the rows of the last build of these errors.csv files (e.g. of parallel jobs) "
                         "into <case>/errors.csv and report")
    ap.add_argument("--write-reference", action="store_true",
                    help="write the errors to <case>/reference.csv (default sweep and committed control file only)")
    args = ap.parse_args()
    sweep = SWEEPS.get(args.case, SWEEPS["default"])
    args.P = args.P or sweep["P"]
    args.meshes = args.meshes or sweep["meshes"]

    overrides = [tuple(s.strip() for s in pair.split("=", 1)) for item in args.set for pair in item.split(";") if pair.strip()]
    if any(len(o) != 2 for o in overrides):
        sys.exit("--set expects KEY=VALUE")
    if args.write_reference and (overrides or set(args.P) != set(sweep["P"]) or set(args.meshes) != set(sweep["meshes"])):
        sys.exit("--write-reference writes the reference of the default sweep: run it without --set, --P and --meshes")
    _, solver, variables = generate.CASES[args.case]
    case_dir = HERE / args.case
    lock = open(case_dir / ".lock", "w")         # held until exit
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        sys.exit(f"Another run.py is working on {args.case}: its runs share the files of {case_dir}")
    lines, added = control(case_dir, overrides)
    for key in added:
        print(("::warning::" if os.environ.get("GITHUB_ACTIONS") else "WARNING: ")
              + f"'{key}' is not a key of {args.case}/control_template.control, so it was added; HORSES3D ignores"
              " the keys it does not read (check the spelling, and the units in the name, e.g. 'Viscosity (Pa.s)')")
    dt, steps, final_time = (value(lines, k) for k in ("dt", "number of time steps", "final time"))
    if dt is None or steps is None:
        sys.exit("The control file must give dt and the number of time steps")
    t_final = dt * steps if final_time is None else min(dt * steps, final_time)    # the solver stops at the first
    title = f"dt = {dt:g}, t_final = {t_final:g}"
    if overrides:
        title += ", " + "; ".join(f"{k} = {v}" + (" (new key)" if k in added else "") for k, v in overrides)
    header = ["build", "N", "P", "dt", "t_final", "NDOF", "volume", "L2_total"] + [f"L2_{v}" for v in variables]
    errors_csv = case_dir / "errors.csv"
    failed = {}

    if args.merge:
        # The rows of the last build of each file (run.py writes them last), one per mesh and
        # order, taken from the last file that has it
        merged = {}
        for f in args.merge:
            file_rows = read_rows(f, header, strict=True)
            for r in file_rows:
                if r["build"] == file_rows[-1]["build"]:
                    merged[(r["N"], r["P"])] = dict(r, build="merged")
        rows = list(merged.values())
        if any(abs(r["dt"] - dt) > 1e-9 * dt or abs(r["t_final"] - t_final) > 0.5 * dt for r in rows):
            sys.exit("The files to merge have another dt or final time than the control file: give the --set of their runs")
        write_rows(errors_csv, header, rows)
    else:
        home = args.home.resolve()
        binary = home / "bin" / f"horses3d.{solver}"
        if not binary.exists():
            sys.exit(f"{binary} not found: build HORSES3D first (make {solver})")
        build(case_dir, solver, home)
        h = hashlib.sha256()
        for f in (binary, case_dir / "SETUP" / "ProblemFile.f90", case_dir / "control_template.control", home / "make.inc"):
            if f.exists():
                h.update(f.read_bytes())
        h.update(repr(overrides).encode())
        h.update(repr([os.environ.get(v, "") for v in BUILD_SETTINGS]).encode())
        build_id = h.hexdigest()[:12]
        rows = read_rows(errors_csv, header)

        for P in args.P:
            for N in args.meshes:
                if any((r["build"], r["N"], r["P"]) == (build_id, N, P) for r in rows):
                    if not args.force:
                        print(f"{args.case} N={N} P={P}: already in errors.csv for this build, skipped (--force runs it again)")
                        continue
                    rows = [r for r in rows if (r["build"], r["N"], r["P"]) != (build_id, N, P)]
                start = time.time()
                result = run(case_dir, binary, lines, N, P, dt, t_final, args.keep)
                if isinstance(result, str):
                    failed[(N, P)] = result
                    print(f"{args.case} N={N} P={P}: FAILED ({result}), see {case_dir}/CONTROL/N{N}_P{P}.log", flush=True)
                    continue
                rows.append(dict(zip(header, [build_id, N, P, dt, t_final] + result)))
                write_rows(errors_csv, header, rows)
                print(f"{args.case} N={N} P={P}: L2 = {rows[-1]['L2_total']:.4e} ({time.time() - start:.0f} s)", flush=True)

        # The rows of this build go last, the ones plot.py draws
        rows = [r for r in rows if r["build"] != build_id] + [r for r in rows if r["build"] == build_id]
        write_rows(errors_csv, header, rows)
        rows = [r for r in rows if r["build"] == build_id]

    rows = [r for r in rows if r["P"] in args.P and r["N"] in args.meshes]
    report(args.case, rows, args.P, title)
    done = {(r["N"], r["P"]) for r in rows}
    missing = [(N, P) for P in args.P for N in args.meshes if (N, P) not in done]
    notes = [f"FAIL: N={N} P={P} " + (f"failed ({failed[(N, P)]})" if (N, P) in failed else "has no result (failed or not run)")
             for N, P in missing]
    ok = not missing
    if args.check:
        passed, text = check(case_dir, rows, bool(overrides), sweep)
        ok = passed and ok
        notes += text + ["PASS" if ok else "FAILED"]
    if args.write_reference:
        if missing:
            notes.append("reference.csv not written: some runs have no result")
        else:
            write_reference(case_dir / "reference.csv", header, rows)
            notes.append(f"Written {args.case}/reference.csv")
    if notes:
        summary("\n".join(notes), "".join(f"* {n}\n" for n in notes))
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
