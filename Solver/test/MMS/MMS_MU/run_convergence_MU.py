"""
run_convergence_MU.py
=====================
Runs the HORSES3D Multiphase MMS convergence study automatically.

USAGE
-----
1. Edit SECTION 1 — Paths (binary, meshes, results directory).
2. Edit SECTION 2 — Study parameters (mesh sizes, P values, dt, t_final).
3. Run:  python3 run_convergence_MU.py

The script will:
  - Generate one control file per (mesh, P) combination
  - Run the solver binary for each case sequentially
  - Read mms_l2_error.dat after each run
  - Append results to errors.csv
  - Skip cases already present in errors.csv (safe to restart after a crash)

Add --dry-run to print what would be run without executing anything:
  python3 run_convergence_MU.py --dry-run

OUTPUT
------
errors.csv — one row per case:
  nelems, P, NDOF, L2_error, t_final
"""

import os
import sys
import csv
import math
import subprocess
import argparse
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 1 — Paths
# ═══════════════════════════════════════════════════════════════════════════════

# Full path to the horses3d.MU binary
HORSES_BINARY = "../../../bin/horses3d.mu"

# Directory containing the mesh files (meshN.h5)
# Mesh files are expected to be named: mesh{N}.h5
MESH_DIR = "./../../TestMeshes/Cubes"

# Directory where solution files (.hsol) will be written
RESULTS_DIR = "./RESULTS"

# Control file template (relative to this script)
TEMPLATE_FILE = "control_template.control"

# Where each run's control file is written (overwritten each run)
CONTROL_DIR = "./CONTROL"

# Output CSV file
ERRORS_CSV = "errors.csv"

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 2 — Study parameters
# ═══════════════════════════════════════════════════════════════════════════════

# Mesh numbers, meshes are expected to be named meshN with changing N
MESH_SIZES = [4, 5, 6]

# Polynomial orders to test, cannot be higher than MAX_ORDER in the problem-file generator
# (MAX_ORDER = 8 by default)
P_VALUES = [2, 3, 4, 5, 6]

# Time step — kept fixed across all cases to avoid introducing
# variability in the errors. Choose conservative enough for the
# finest mesh and highest P in your study.
DT = 1.0e-7

# Final simulation time
T_FINAL = 1.0e-6

# ═══════════════════════════════════════════════════════════════════════════════
#  END OF USER CONFIGURATION 
# ═══════════════════════════════════════════════════════════════════════════════

CSV_HEADER = ["nelems", "P", "NDOF", "L2_error", "t_final"]


def n_steps(dt, t_final):
    """Compute number of time steps, ceiling-based to guarantee t_final is reached."""
    return max(1, math.ceil(t_final / dt))


def ndof(nelems, P):
    return nelems * (P + 1)**3


def load_completed(csv_path):
    """Return set of (nelems, P) tuples already in errors.csv."""
    completed = set()
    if not os.path.exists(csv_path):
        return completed
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                completed.add((int(row["nelems"]), int(row["P"])))
            except (KeyError, ValueError):
                pass
    return completed


def csv_t_finals(csv_path):
    """Return the set of distinct t_final values found in errors.csv (rounded
    to 12 significant digits to avoid float-noise false positives)."""
    seen = set()
    if not os.path.exists(csv_path):
        return seen
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                seen.add(float(f"{float(row['t_final']):.12g}"))
            except (KeyError, ValueError):
                pass
    return seen


def append_result(csv_path, row):
    """Append one result row to errors.csv, creating with header if needed."""
    write_header = not os.path.exists(csv_path)
    with open(csv_path, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CSV_HEADER)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def make_control_file(template_path, output_path, N, P, dt, t_final,
                      mesh_dir, results_dir):
    """Fill template tokens and write control file."""
    with open(template_path) as f:
        content = f.read()

    steps    = n_steps(dt, t_final)
    mesh_file  = os.path.join(mesh_dir,    f"mesh{N}.h5")
    sol_file   = os.path.join(results_dir, f"MMS_N{N}_P{P}.hsol")

    content = content.replace("{MESH_FILE}",     mesh_file)
    content = content.replace("{SOLUTION_FILE}", sol_file)
    content = content.replace("{P}",             str(P))
    content = content.replace("{DT}",            f"{dt:.6e}")
    content = content.replace("{N_STEPS}",       str(steps))

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(content)


def read_l2_error(dat_path):
    """
    Read mms_l2_error.dat written by UserDefinedFinalize.
    Format: nelems, P, NDOF, L2_error, t_final
    Returns a dict or None on failure.
    """
    if not os.path.exists(dat_path):
        return None
    with open(dat_path) as f:
        line = f.read().strip()
    parts = [p.strip() for p in line.split(',')]
    if len(parts) != 5:
        return None
    try:
        return {
            "nelems":   int(parts[0]),
            "P":        int(parts[1]),
            "NDOF":     int(parts[2]),
            "L2_error": float(parts[3]),
            "t_final":  float(parts[4]),
        }
    except ValueError:
        return None


def run_case(binary, control_file, run_dir, log_path):
    """
    Run one solver case.
    Returns (success, return_code).
    """
    cmd = [binary, control_file]
    with open(log_path, 'w') as log:
        result = subprocess.run(
            cmd,
            cwd=run_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    return result.returncode == 0, result.returncode


def build_case_list(mesh_sizes, p_values, completed):
    """Return list of (N, P) pairs not yet in completed."""
    cases = []
    for N in mesh_sizes:
        for P in p_values:
            if (N**3, P) not in completed:
                cases.append((N, P))
    return cases


def main():
    parser = argparse.ArgumentParser(description="MU MMS convergence runner")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print cases without running them")
    args = parser.parse_args()

    # ── Validate paths ──────────────────────────────────────────────────────
    script_dir = Path(__file__).resolve().parent

    # Resolve every path relative to the script's own directory so that all
    # file I/O, the subprocess cwd, and the paths written into the control
    # file are consistent regardless of where the user's shell cwd happens
    # to be when this script is launched.
    horses_binary = str((script_dir / HORSES_BINARY).resolve())
    mesh_dir      = str((script_dir / MESH_DIR).resolve())
    results_dir   = str((script_dir / RESULTS_DIR).resolve())
    control_dir   = str((script_dir / CONTROL_DIR).resolve())
    errors_csv    = str((script_dir / ERRORS_CSV).resolve())

    template_path = script_dir / TEMPLATE_FILE
    if not template_path.exists():
        print(f"ERROR: template not found: {template_path}")
        sys.exit(1)

    if not args.dry_run and not os.path.exists(horses_binary):
        print(f"ERROR: binary not found: {horses_binary}")
        sys.exit(1)

    # Validate mesh directory and per-N mesh files upfront so the user gets
    # a clear error before any case is launched.  Skipped under --dry-run so
    # users can preview the case list before meshes are in place.
    if not args.dry_run:
        if not os.path.isdir(mesh_dir):
            print(f"ERROR: mesh directory not found: {mesh_dir}")
            sys.exit(1)
        missing_meshes = [N for N in MESH_SIZES
                          if not os.path.exists(os.path.join(mesh_dir, f"mesh{N}.h5"))]
        if missing_meshes:
            print(f"ERROR: missing mesh files in {mesh_dir}:")
            for N in missing_meshes:
                print(f"  mesh{N}.h5")
            sys.exit(1)

    # ── Set up directories ──────────────────────────────────────────────────
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(control_dir, exist_ok=True)

    # ── Build run list ──────────────────────────────────────────────────────
    completed = load_completed(errors_csv)
    cases     = build_case_list(MESH_SIZES, P_VALUES, completed)

    # Warn if the existing CSV holds results from a different time-integration
    # setup.  Mixing such rows produces misleading convergence plots, but is
    # easy to do accidentally when iterating on study parameters.
    seen_tfinals = csv_t_finals(errors_csv)
    current_tfinal_rounded = float(f"{float(T_FINAL):.12g}")
    other_tfinals = seen_tfinals - {current_tfinal_rounded}
    if other_tfinals:
        print(f"WARNING: {errors_csv} contains rows with t_final values "
              f"{sorted(other_tfinals)} that differ from the current "
              f"T_FINAL = {T_FINAL}.  Mixing these rows in the same CSV "
              "yields misleading convergence comparisons.  Delete the CSV "
              "to force a fresh sweep.")
        print()

    total = len(MESH_SIZES) * len(P_VALUES)
    print("MMS Convergence Study — Multiphase (MU)")
    print(f"  Binary   : {horses_binary}")
    print(f"  Meshes   : {MESH_SIZES}")
    print(f"  P values : {P_VALUES}")
    print(f"  dt       : {DT:.2e}   t_final : {T_FINAL}")
    print(f"  Total cases    : {total}")
    print(f"  Already done   : {len(completed)}")
    print(f"  Remaining      : {len(cases)}")
    print()

    if args.dry_run:
        print("DRY RUN — cases that would be run:")
        for N, P in cases:
            steps = n_steps(DT, T_FINAL)
            print(f"  N={N}  P={P}  nelems={N**3}  "
                  f"NDOF={ndof(N**3, P)}  steps={steps}")
        return

    if not cases:
        print("All cases already complete. Nothing to run.")
        print(f"Results are in {errors_csv}")
        return

    # ── Run cases ────────────────────────────────────────────────────────────
    run_dir  = str(script_dir)
    dat_path = os.path.join(run_dir, "mms_l2_error.dat")

    for idx, (N, P) in enumerate(cases, 1):
        nelems = N**3
        steps  = n_steps(DT, T_FINAL)
        label  = f"N={N} P={P}"
        print(f"[{idx}/{len(cases)}] {label}  "
              f"(nelems={nelems}, NDOF={ndof(nelems,P)}, steps={steps})",
              end="  ", flush=True)

        # Write control file (absolute path so cwd of caller doesn't matter)
        ctrl_path = os.path.join(control_dir, f"mms_N{N}_P{P}.control")
        make_control_file(
            template_path = str(template_path),
            output_path   = ctrl_path,
            N=N, P=P, dt=DT, t_final=T_FINAL,
            mesh_dir    = mesh_dir,
            results_dir = results_dir,
        )

        # Remove stale dat file so a crashed run doesn't give a false result
        if os.path.exists(dat_path):
            os.remove(dat_path)

        # Run solver (absolute log path)
        log_path = os.path.join(control_dir, f"mms_N{N}_P{P}.log")
        success, retcode = run_case(horses_binary, ctrl_path, run_dir, log_path)

        if not success:
            print(f"FAILED (exit code {retcode}) — see {log_path}")
            continue

        # Read result
        result = read_l2_error(dat_path)
        if result is None:
            print(f"ERROR: could not read {dat_path}")
            continue

        # Append to CSV
        append_result(errors_csv, result)
        print(f"L2 = {result['L2_error']:.6e}")

    # ── Summary ──────────────────────────────────────────────────────────────
    print()
    completed_now = load_completed(errors_csv)
    print(f"Done. {len(completed_now)}/{total} cases in {errors_csv}")


if __name__ == "__main__":
    main()
