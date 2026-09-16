"""
plot_convergence_MU.py
======================
Reads errors.csv produced by run_convergence_MU.py and generates
convergence plots for the HORSES3D Multiphase MMS study.

USAGE
-----
    python3 plot_convergence_MU.py

Or with a specific CSV file:
    python3 plot_convergence_MU.py --csv path/to/errors.csv

OUTPUT
------
  convergence.pdf  —  vector figure
  convergence.png  —  raster version for quick viewing

PLOTS
-----
  Left   : L2 error vs sqrt(NDOF) — combined h- and p-convergence view
  Centre : L2 error vs h = 1/N (log-log) — h-convergence per P value
           Best-fit slope annotated. Expected slope = P+1 (theory).
           Note: default example exhibits even-odd behaviour
  Right  : L2 error vs polynomial order P — p-convergence per mesh
           Empirical rates annotated. Faded markers = time-floor.
"""

import csv
import argparse
import os
import sys
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 1 — Plot settings  (edit if needed)
# ═══════════════════════════════════════════════════════════════════════════════

FLOOR_RATIO_THRESHOLD = 2.0

REFERENCE_SLOPES = [
    (3, '#BBBBBB', '$O(h^{3})$'),
    (4, '#888888', '$O(h^{4})$'),
    (5, '#555555', '$O(h^{5})$'),
]

P_COLOURS = {2:'#E91E63', 3:'#2196F3', 4:'#4CAF50',
             5:'#FF9800', 6:'#9C27B0', 7:'#00BCD4', 8:'#795548'}
P_MARKERS = {2:'o', 3:'s', 4:'^', 5:'D', 6:'p', 7:'h', 8:'*'}

PDF_FILE = "convergence.pdf"
PNG_FILE = "convergence.png"

# ═══════════════════════════════════════════════════════════════════════════════
#  END OF USER CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

MESH_COLOURS = ["#F44336","#2196F3","#9C27B0","#E91E63",
                "#4CAF50","#FF9800","#00BCD4","#795548"]
MESH_MARKERS = ["v","o","p","s","^","D","h","*"]


def load_csv(csv_path):
    data   = defaultdict(dict)
    p_data = defaultdict(dict)
    t_finals = set()
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            nelems = int(row["nelems"])
            P      = int(row["P"])
            L2     = float(row["L2_error"])
            tf     = float(row["t_final"])
            N      = round(nelems ** (1/3))
            label  = f"{N}\u00d7{N}\u00d7{N}"
            data[label][P] = L2
            p_data[P][N]   = L2
            # Round to suppress float-repr noise before deduping
            t_finals.add(float(f"{tf:.12g}"))
    ordered = dict(sorted(data.items(),
                          key=lambda kv: int(kv[0].split('\u00d7')[0])))
    return ordered, p_data, sorted(t_finals)


def ndof_from_label(label, P):
    N = int(label.split('\u00d7')[0])
    return (N * (P + 1))**3


def detect_floor(errors_by_p):
    floor_ps = set()
    ps = sorted(errors_by_p.keys())
    for i in range(1, len(ps)):
        if errors_by_p[ps[i-1]] / errors_by_p[ps[i]] < FLOOR_RATIO_THRESHOLD:
            floor_ps.add(ps[i])
    return floor_ps


def empirical_rate_p(e1, e2, p1, p2):
    if e2 >= e1:
        return None
    return math.log10(e1 / e2) / math.log10((p2 + 1) / (p1 + 1))


def fit_slope(ns, errors):
    """Best-fit slope of log(error) vs log(h=1/N). Returns slope ~ P+1."""
    if len(ns) < 2:
        return None, None
    log_h = np.log10([1.0 / n for n in ns])
    log_e = np.log10(errors)
    coeffs = np.polyfit(log_h, log_e, 1)
    return coeffs[0], coeffs[1]   # slope, intercept


def make_plots(data, p_data, t_finals, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(20, 5.8))
    fig.patch.set_facecolor('#FAFAFA')
    fig.suptitle('Multiphase (MU) — MMS Convergence',
                 fontsize=14, fontweight='bold', y=1.01)

    for ax in axes:
        ax.set_facecolor('#FAFAFA')
        ax.grid(True, which='both', linestyle='--', linewidth=0.5,
                color='#CCCCCC', alpha=0.8)
        ax.set_yscale('log')
        ax.set_ylabel('$L^2$ error', fontsize=12)

    ax_dof, ax_h, ax_p = axes

    ax_dof.set_xlabel('$\\sqrt{\\mathrm{NDOF}}$', fontsize=12)
    ax_dof.set_xscale('log')
    ax_dof.set_title('Convergence vs DOF', fontsize=13, fontweight='bold')

    ax_h.set_xlabel('$h = 1/N$', fontsize=12)
    ax_h.set_xscale('log')
    ax_h.set_title('$h$-convergence per $p$  (slope $\\approx p+1$)',
                   fontsize=13, fontweight='bold')

    ax_p.set_xlabel('Polynomial order $p$', fontsize=12)
    ax_p.set_title('$p$-convergence per mesh', fontsize=13, fontweight='bold')
    ax_p.xaxis.set_major_locator(ticker.MultipleLocator(1))

    mesh_labels = list(data.keys())

    # ── Left and right panels (per mesh) ──────────────────────────────────────
    for mesh_idx, label in enumerate(mesh_labels):
        results  = data[label]
        col      = MESH_COLOURS[mesh_idx % len(MESH_COLOURS)]
        mk       = MESH_MARKERS[mesh_idx % len(MESH_MARKERS)]
        floor_ps = detect_floor(results)
        ps       = sorted(results.keys())
        errors   = [results[p] for p in ps]

        # Left
        xvals = [ndof_from_label(label, p)**0.5 for p in ps]
        ax_dof.plot(xvals, errors, '-', color=col, linewidth=1.8,
                    label=label, zorder=3)
        for xv, e, p in zip(xvals, errors, ps):
            fc    = '#EEEEEE' if p in floor_ps else 'white'
            alpha = 0.5 if p in floor_ps else 1.0
            ax_dof.plot(xv, e, mk, color=col, markersize=7,
                        markerfacecolor=fc, markeredgewidth=2,
                        alpha=alpha, zorder=4)
            ax_dof.annotate(f'$p={p}$', (xv, e),
                            textcoords='offset points', xytext=(5, 2),
                            fontsize=7.5, color=col)

        # Right
        ax_p.plot(ps, errors, '-', color=col, linewidth=1.8,
                  label=label, zorder=3)
        for p, e in zip(ps, errors):
            fc    = '#EEEEEE' if p in floor_ps else 'white'
            alpha = 0.5 if p in floor_ps else 1.0
            ax_p.plot(p, e, mk, color=col, markersize=7,
                      markerfacecolor=fc, markeredgewidth=2,
                      alpha=alpha, zorder=4)
        for i in range(len(ps) - 1):
            rate = empirical_rate_p(errors[i], errors[i+1], ps[i], ps[i+1])
            if rate is not None and rate > 0.5:
                xm = (ps[i] + ps[i+1]) / 2
                ym = (errors[i] * errors[i+1])**0.5
                ax_p.annotate(f'{rate:.1f}×', (xm, ym),
                              textcoords='offset points', xytext=(5, 0),
                              fontsize=7.5, color=col, alpha=0.85)

    # ── Centre panel: h-convergence per P ────────────────────────────────────
    for P in sorted(p_data.keys()):
        n_err  = p_data[P]
        ns     = sorted(n_err.keys())
        hs     = [1.0 / n for n in ns]
        errors = [n_err[n] for n in ns]
        col    = P_COLOURS.get(P, '#333333')
        mk     = P_MARKERS.get(P, 'o')

        # Data markers
        ax_h.plot(hs, errors, color=col, marker=mk, markersize=7,
                  markerfacecolor='white', markeredgewidth=2,
                  linewidth=0, zorder=4)

        # Best-fit line
        if len(ns) >= 2:
            slope, intercept = fit_slope(ns, errors)
            if slope is not None:
                log_h_arr = np.log10(hs)
                h_fit = np.logspace(log_h_arr.min() - 0.1,
                                    log_h_arr.max() + 0.1, 50)
                e_fit = 10**intercept * h_fit**slope
                ax_h.plot(h_fit, e_fit, '-', color=col, linewidth=1.8,
                          label=f'$p={P}$  '
                                f'fit={slope:.2f} / theory={P+1}',
                          zorder=3)
                # Annotate at midpoint of fit line
                mid = len(h_fit) // 2
                ax_h.annotate(f'$p={P}$\n{slope:.2f}',
                              (h_fit[mid], e_fit[mid]),
                              textcoords='offset points', xytext=(6, 2),
                              fontsize=7.5, color=col)

    # ── Reference slopes (left panel) ─────────────────────────────────────────
    first_label = mesh_labels[0]
    first_p     = sorted(data[first_label].keys())[0]
    ref_x = ndof_from_label(first_label, first_p)**0.5
    ref_e = data[first_label][first_p] * 3.5
    for slope, grey, lbl in REFERENCE_SLOPES:
        x2 = ref_x * 5
        e2 = ref_e * (x2 / ref_x)**(-slope)
        ax_dof.plot([ref_x, x2], [ref_e, e2], '--', color=grey,
                    linewidth=0.9, alpha=0.7, label=lbl)

    ax_dof.legend(fontsize=9, framealpha=0.9, loc='lower left', ncol=2)
    ax_h.legend(fontsize=8, framealpha=0.9, loc='lower right')
    ax_p.legend(fontsize=10, framealpha=0.9, loc='upper right')

    if len(t_finals) == 1:
        tf_str = f"t_final = {t_finals[0]:.3g}"
    elif len(t_finals) > 1:
        tf_str = ("t_final = MIXED " +
                  "(" + ", ".join(f"{tf:.3g}" for tf in t_finals) + ")")
    else:
        tf_str = "t_final = n/a"

    fig.text(0.5, 0.005,
             f'Faded markers: time-integration floor '
             f'(ratio < {FLOOR_RATIO_THRESHOLD}).  '
             f'{tf_str}',
             ha='center', fontsize=8, color='#777777', style='italic')

    plt.tight_layout(pad=2.0, rect=[0, 0.03, 1, 1])
    pdf_path = os.path.join(outdir, PDF_FILE)
    png_path = os.path.join(outdir, PNG_FILE)
    plt.savefig(pdf_path, bbox_inches='tight', dpi=150)
    plt.savefig(png_path, bbox_inches='tight', dpi=150)
    plt.close()
    return pdf_path, png_path


def main():
    parser = argparse.ArgumentParser(
        description="Plot MU MMS convergence from errors.csv")
    parser.add_argument("--csv", default=None,
                        help="Path to errors.csv (default: errors.csv next to this script)")
    parser.add_argument("--outdir", default=None,
                        help="Output directory for plots (default: same as CSV)")
    args = parser.parse_args()

    # When --csv is not given, default to errors.csv next to this script so
    # the plotter can be launched from any working directory.  An explicit
    # --csv argument is honoured as-is (relative paths resolve against cwd).
    script_dir = Path(__file__).resolve().parent
    if args.csv is None:
        csv_path = script_dir / "errors.csv"
    else:
        csv_path = Path(args.csv)

    if not csv_path.exists():
        print(f"ERROR: {csv_path} not found.")
        print("Run run_convergence_MU.py first to generate results.")
        sys.exit(1)

    outdir = args.outdir or str(csv_path.parent)
    os.makedirs(outdir, exist_ok=True)

    print(f"Reading {csv_path}...")
    data, p_data, t_finals = load_csv(str(csv_path))

    if not data:
        print("ERROR: no data found in CSV.")
        sys.exit(1)

    meshes = list(data.keys())
    print(f"  Meshes  : {meshes}")
    for label in meshes:
        ps = sorted(data[label].keys())
        print(f"  {label} : P = {ps}")
    print()

    pdf_path, png_path = make_plots(data, p_data, t_finals, outdir)
    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


if __name__ == "__main__":
    main()
