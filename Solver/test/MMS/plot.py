#!/usr/bin/env python3
"""
Plot <case>/errors.csv into <case>/convergence.png.

    python3 plot.py NS|iNS|MU

Left: error vs sqrt(NDOF) per mesh. Centre: error vs h = 1/N per P, with the
least-squares slope (theory P+1). Right: error vs P per mesh.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

P_COLOURS = {1: "#795548", 2: "#E91E63", 3: "#2196F3", 4: "#4CAF50", 5: "#FF9800", 6: "#9C27B0", 7: "#00BCD4"}
MESH_COLOURS = ["#F44336", "#2196F3", "#9C27B0", "#4CAF50", "#FF9800"]


def main(case):
    case_dir = Path(__file__).resolve().parent / case
    rows = list(csv.DictReader(open(case_dir / "errors.csv")))
    # Only the runs with the build, dt and final time of the last one
    key = lambda r: (r.get("build"), r["dt"], r["t_final"])
    rows = [r for r in rows if key(r) == key(rows[-1])]
    by_p, by_n = defaultdict(dict), defaultdict(dict)
    for r in rows:
        N, P, e = int(r["N"]), int(r["P"]), float(r["L2_total"])
        by_p[P][N] = e
        by_n[N][P] = e

    fig, (ax_dof, ax_h, ax_p) = plt.subplots(1, 3, figsize=(18, 5.5))
    fig.patch.set_facecolor("#FAFAFA")
    fig.suptitle(f"MMS convergence, {case}  (dt = {float(rows[-1]['dt']):g}, t_final = {float(rows[-1]['t_final']):g})",
                 fontweight="bold")
    for ax, xlabel, title in [(ax_dof, r"$\sqrt{\mathrm{NDOF}}$", "Error vs DOF"),
                              (ax_h, r"$h = 1/N$", r"$h$-convergence per $p$"),
                              (ax_p, r"$p$", r"$p$-convergence per mesh")]:
        ax.set_facecolor("#FAFAFA")
        ax.grid(True, which="both", linestyle="--", linewidth=0.5, color="#CCCCCC")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("$L^2$ error")
        ax.set_title(title, fontweight="bold")
    ax_dof.set_xscale("log")
    ax_h.set_xscale("log")

    for i, N in enumerate(sorted(by_n)):
        ps = sorted(by_n[N])
        e = [by_n[N][p] for p in ps]
        col = MESH_COLOURS[i % len(MESH_COLOURS)]
        ax_dof.plot([(N*(p + 1))**1.5 for p in ps], e, "o-", color=col, mfc="white", mew=2, label=f"{N}x{N}x{N}")
        ax_p.plot(ps, e, "o-", color=col, mfc="white", mew=2, label=f"{N}x{N}x{N}")

    for P in sorted(by_p):
        ns = sorted(by_p[P])
        h = np.array([1.0/n for n in ns])
        e = np.array([by_p[P][n] for n in ns])
        col = P_COLOURS.get(P, "#333333")
        ax_h.plot(h, e, "o", color=col, mfc="white", mew=2)
        if len(ns) > 1:
            slope, icpt = np.polyfit(np.log(h), np.log(e), 1)
            hf = np.geomspace(h.min()/1.2, h.max()*1.2, 20)
            ax_h.plot(hf, np.exp(icpt)*hf**slope, "-", color=col, label=f"$p={P}$: {slope:.2f} (theory {P + 1})")

    for ax in (ax_dof, ax_h, ax_p):
        ax.legend(fontsize=9)
    ax_p.xaxis.set_major_locator(MultipleLocator(1))
    plt.tight_layout()
    out = case_dir / "convergence.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Written {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} CASE")
    main(sys.argv[1])
