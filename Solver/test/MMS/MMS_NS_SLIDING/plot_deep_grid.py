#!/usr/bin/env python3
"""
plot_deep_grid.py — convergence figures for the deep-cylinder unsteady-MMS
rotation study (full mesh x P x omega grid in errors_deep.csv).

Produces, in PLOTS/:
  * ndof_om{W}.png         L2 vs NDOF for each omega:
                             - solid colored line joins same P (h-refinement)
                             - dotted gray line joins same mesh (p-refinement)
  * ndof_overlay.png       L2 vs NDOF, omega in {0,16,256} overlaid (color = omega)
  * ndof_grid.png          small-multiples of all omega
  * pconv_deep48.png       deep48 L2 vs P, one colored curve per omega
"""
import csv, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
OUT  = HERE / "PLOTS"; OUT.mkdir(exist_ok=True)

ELEMS  = {"deep24": 768, "deep36": 2592, "deep48": 6144}   # N x N/3 x layers
MESHES = ["deep24", "deep36", "deep48"]
PS     = [1, 2, 3, 4]
OMEGAS = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256]
PCOLOR = {1: "#CC79A7", 2: "#0072B2", 3: "#E69F00", 4: "#009E73"}   # colorblind-safe
MESHLBL = {"deep24": "24", "deep36": "36", "deep48": "48"}
# hue-family overlay: each omega gets a color FAMILY, shade = P (light P1 -> dark P4).
# Families are assigned to omegas in list order; 4 shades = P1,P2,P3,P4.
FAMILIES = [
    ("reds",    ["#fca5a5", "#f87171", "#dc2626", "#7f1d1d"]),
    ("blues",   ["#bfdbfe", "#60a5fa", "#2563eb", "#1e3a8a"]),
    ("greens",  ["#bbf7d0", "#4ade80", "#16a34a", "#14532d"]),
    ("oranges", ["#fed7aa", "#fb923c", "#ea580c", "#9a3412"]),
    ("purples", ["#e9d5ff", "#c084fc", "#9333ea", "#6b21a8"]),
    ("pinks",   ["#fbcfe8", "#f472b6", "#db2777", "#9d174d"]),
]

def ndof(mesh, P): return ELEMS[mesh] * (P + 1) ** 3

# ── load ────────────────────────────────────────────────────────────────────
D = {}
for r in csv.reader(open(HERE / "errors_deep.csv")):
    if not r or r[0] == "mesh":
        continue
    mesh, P, om = r[0], int(r[1]), int(float(r[2]))
    l2 = float(str(r[-1]).replace(",", "").strip())        # strip runner's stray comma
    D[(mesh, P, om)] = l2

def pts_sameP(om, P):    # h-refinement: fixed P, vary mesh
    xs = [(ndof(m, P), D[(m, P, om)]) for m in MESHES if (m, P, om) in D]
    return sorted(xs)

def pts_samemesh(om, m): # p-refinement: fixed mesh, vary P
    xs = [(ndof(m, P), D[(m, P, om)]) for P in PS if (m, P, om) in D]
    return sorted(xs)

NAZ = {"deep24": 24, "deep36": 36, "deep48": 48}   # azimuthal count; h ~ 1/N

def hpts(om, P):         # (N, error) over meshes present, sorted by N
    return sorted((NAZ[m], D[(m, P, om)]) for m in MESHES if (m, P, om) in D)

def conv_order(om, P):
    """Least-squares order of convergence: error ~ h^order, h ~ 1/N."""
    d = hpts(om, P)
    if len(d) < 2:
        return None
    N = np.array([a for a, _ in d], float)
    e = np.array([b for _, b in d], float)
    return -np.polyfit(np.log(N), np.log(e), 1)[0]

def order_pair(om, P, m1, m2):
    e1, e2 = D.get((m1, P, om)), D.get((m2, P, om))
    if e1 is None or e2 is None:
        return None
    return math.log(e1 / e2) / math.log(NAZ[m2] / NAZ[m1])

# ── per-omega NDOF plot ───────────────────────────────────────────────────────
def draw_ndof(ax, om, legend=True):
    # dotted gray: same mesh (p-refinement) — drawn first, behind
    for m in MESHES:
        xy = pts_samemesh(om, m)
        if len(xy) >= 2:
            xs, ys = zip(*xy)
            ax.plot(xs, ys, ls=":", color="0.6", lw=1.2, zorder=1)
            ax.annotate(MESHLBL[m], (xs[0], ys[0]), textcoords="offset points",
                        xytext=(-14, -2), fontsize=7.5, color="0.45")
    # solid colored: same P (h-refinement)
    for P in PS:
        xy = pts_sameP(om, P)
        if not xy:
            continue
        xs, ys = zip(*xy)
        ax.plot(xs, ys, "-o", color=PCOLOR[P], lw=2, ms=5.5,
                label=f"P{P}", zorder=3)
        o = conv_order(om, P)
        if o is not None:
            ax.annotate(f"$p$={o:.2f}", (xs[-1], ys[-1]),
                        textcoords="offset points", xytext=(7, -1),
                        fontsize=8.5, color=PCOLOR[P], zorder=4)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.grid(True, which="both", ls="-", alpha=0.15)
    ax.set_title(f"$\\omega$ = {om}", fontsize=11)
    if legend:
        h = ax.plot([], [], ls=":", color="0.6", lw=1.2, label="same mesh (p-ref.)")
        ax.legend(fontsize=8, frameon=False, loc="lower left")

for om in OMEGAS:
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    draw_ndof(ax, om)
    ax.set_xlabel("Degrees of freedom"); ax.set_ylabel("$L_2$ error")
    fig.tight_layout(); fig.savefig(OUT / f"ndof_om{om}.png", dpi=160)
    plt.close(fig)

# ── small-multiples grid of all omega ─────────────────────────────────────────
fig, axes = plt.subplots(2, 5, figsize=(20, 8), sharex=True, sharey=True)
for ax, om in zip(axes.flat, OMEGAS):
    draw_ndof(ax, om, legend=False)
for ax in axes[:, 0]:  ax.set_ylabel("$L_2$ error")
for ax in axes[-1, :]: ax.set_xlabel("Degrees of freedom")
axes.flat[0].legend(fontsize=8, frameon=False, loc="lower left")
fig.suptitle("L2 error vs NDOF — solid=same P (h-refinement), dotted=same mesh (p-refinement)", fontsize=13)
fig.tight_layout(); fig.savefig(OUT / "ndof_grid.png", dpi=140); plt.close(fig)

# ── overlay omega in {0,16,256} (color = omega) ───────────────────────────────
OV = [0, 16, 256]
ovcol = {0: "#4477AA", 16: "#EE6677", 256: "#228833"}
fig, ax = plt.subplots(figsize=(6.4, 5.2))
for om in OV:
    for P in PS:                         # solid P-lines (h-refinement fan), per omega
        xy = pts_sameP(om, P)
        if not xy: continue
        xs, ys = zip(*xy)
        ax.plot(xs, ys, "-o", color=ovcol[om], lw=1.8, ms=5, alpha=0.9,
                label=f"$\\omega$={om}" if P == 2 else None)
ax.set_xscale("log"); ax.set_yscale("log")
ax.grid(True, which="both", ls="-", alpha=0.15)
ax.set_xlabel("Degrees of freedom"); ax.set_ylabel("$L_2$ error")
ax.set_title("h-refinement fans (solid = fixed P) at $\\omega$ = 0, 16, 256")
ax.legend(fontsize=9, frameon=False)
fig.tight_layout(); fig.savefig(OUT / "ndof_overlay.png", dpi=160); plt.close(fig)

# ── hue-family overlay: color family = omega, shade = P ───────────────────────
def draw_hue_overlay(omegas, fname, dotted_alpha=0.3):
    assert len(omegas) <= len(FAMILIES), "add more FAMILIES for this many omegas"
    hue = {om: {P: FAMILIES[i][1][j] for j, P in enumerate(PS)}
           for i, om in enumerate(omegas)}
    wide = len(omegas) > 3
    fig, ax = plt.subplots(figsize=(8.6, 6.4) if wide else (7.2, 5.6))
    # faint gray dotted: same mesh (p-refinement), background structure
    for om in omegas:
        for m in MESHES:
            xy = pts_samemesh(om, m)
            if len(xy) >= 2:
                xs, ys = zip(*xy)
                ax.plot(xs, ys, ls=":", color="0.6", lw=0.8,
                        alpha=dotted_alpha, zorder=1)
    # solid hue lines: same P (h-refinement)
    handles = {}
    for om in omegas:
        for P in PS:
            xy = pts_sameP(om, P)
            if not xy:
                continue
            xs, ys = zip(*xy)
            (h,) = ax.plot(xs, ys, "-o", color=hue[om][P], lw=2, ms=5, zorder=3)
            handles[(om, P)] = h
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.grid(True, which="both", ls="-", alpha=0.15)
    ax.set_xlabel("Degrees of freedom"); ax.set_ylabel("$L_2$ error")
    ax.set_title("L2 vs NDOF — color family = $\\omega$, shade = $P$  "
                 "(solid: same P, dotted: same mesh)", fontsize=10.5)
    # legend below, columns = omega, rows = P. matplotlib fills columns
    # top-to-bottom, so order omega-outer / P-inner => each column is one omega.
    order = [(om, P) for om in omegas for P in PS if (om, P) in handles]
    ax.legend([handles[k] for k in order],
              [f"$\\omega$={k[0]}, P{k[1]}" for k in order],
              ncol=len(omegas), fontsize=7.5, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.10),
              columnspacing=1.0, handlelength=1.5)
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=160,
                                    bbox_inches="tight"); plt.close(fig)

draw_hue_overlay([0, 16, 256], "ndof_overlay_hue.png", dotted_alpha=0.35)
draw_hue_overlay([0, 4, 16, 64, 256], "ndof_overlay_hue5.png", dotted_alpha=0.15)

# ── color = P (family), lightness = omega ─────────────────────────────────────
PCMAP = {1: cm.Purples, 2: cm.Reds, 3: cm.Greens, 4: cm.Blues}   # P1..P4
def omega_shade(P, om):
    i = OMEGAS.index(om)
    return PCMAP[P](0.28 + 0.67 * i / (len(OMEGAS) - 1))          # light=low, dark=high
def grey_shade(om):
    i = OMEGAS.index(om)
    return cm.Greys(0.28 + 0.67 * i / (len(OMEGAS) - 1))

VOL = math.pi * (1.0**2 - 0.5**2) * 2.0                            # annulus r in[.5,1], z in[-1,1]
def h_of(m): return (VOL / ELEMS[m]) ** (1.0 / 3.0)               # avg element spacing

def _pw_legends(ax, omega_key=(0, 4, 16, 64, 256), mesh_marks=None,
                p_loc="lower left", omega_loc="upper right"):
    """P-color legend + omega-lightness key (+ optional mesh-marker legend)."""
    pl = [Line2D([], [], color=PCMAP[P](0.78), marker="o", lw=2.2, label=f"P{P}") for P in PS]
    l1 = ax.legend(handles=pl, title="order (color)", fontsize=9,
                   loc=p_loc, frameon=False)
    ax.add_artist(l1)
    ol = [Line2D([], [], color=grey_shade(o), marker="o", lw=2.2, label=f"{o}") for o in omega_key]
    l2 = ax.legend(handles=ol, title="$\\omega$ (lightness)", fontsize=8,
                   loc=omega_loc, frameon=False)
    ax.add_artist(l2)
    if mesh_marks:
        ml = [Line2D([], [], color="0.35", marker=mesh_marks[m], lw=0, ms=8,
                     label=f"N={NAZ[m]}") for m in MESHES]
        ax.legend(handles=ml, title="mesh (marker)", fontsize=8,
                  loc="lower right", frameon=False)

def draw_l2_vs_h(omegas, fname, tag="", invert=True):
    # invert=True -> x-axis is 1/h (paper convention: refinement increases rightward,
    # error decreases rightward, |slope| = order). invert=False -> x-axis is h.
    xv = {m: (1.0 / h_of(m) if invert else h_of(m)) for m in MESHES}
    fig, ax = plt.subplots(figsize=(8.4, 6.4))
    for P in PS:
        for om in omegas:
            pts = sorted((xv[m], D[(m, P, om)]) for m in MESHES if (m, P, om) in D)
            if len(pts) < 2:
                continue
            xs, ys = zip(*pts)
            ax.plot(xs, ys, "-o", color=omega_shade(P, om), lw=1.5, ms=4.5, zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([xv[m] for m in MESHES])
    ax.set_xticklabels([f"{xv[m]:.2f}\nN={NAZ[m]}" for m in MESHES])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())     # kill log minor labels
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.grid(True, which="both", ls="-", alpha=0.15)
    xlabel = "$1/h$  (inverse element spacing)" if invert else "average element spacing  $h$"
    ax.set_xlabel(xlabel); ax.set_ylabel("$L_2$ error")
    ax.set_title(f"$L_2$ vs {'1/h' if invert else 'mesh spacing'}{tag} — "
                 "color = P, lightness = $\\omega$  (|slope| = order)")
    # place legends in the two empty corners (mirrored between the h and 1/h layouts)
    if invert:  # coarse/high-error on LEFT
        _pw_legends(ax, p_loc="upper right", omega_loc="lower left")
    else:       # coarse/high-error on RIGHT
        _pw_legends(ax, p_loc="upper left", omega_loc="lower right")
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=160); plt.close(fig)

MARK = {"deep24": "o", "deep36": "s", "deep48": "X"}
def draw_ndof_markers(omegas, fname, tag=""):
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    # dotted gray: join a fixed (mesh, omega) across P (p-refinement). Higher P is
    # lower with a steeper h-slope, so this connector "falls" as P increases.
    for m in MESHES:
        for om in omegas:
            pts = sorted((ndof(m, P), D[(m, P, om)]) for P in PS if (m, P, om) in D)
            if len(pts) >= 2:
                xs, ys = zip(*pts)
                ax.plot(xs, ys, ls=":", color="0.55", lw=0.9, alpha=0.55, zorder=1)
    # solid colored: same (P, omega) across meshes (h-refinement)
    for P in PS:
        for om in omegas:
            pts = sorted((ndof(m, P), D[(m, P, om)], m) for m in MESHES if (m, P, om) in D)
            if len(pts) < 2:
                continue
            c = omega_shade(P, om)
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", color=c, lw=1.1, zorder=2)
            for x, y, m in pts:
                ax.scatter([x], [y], marker=MARK[m], color=c, s=34, zorder=3, edgecolors="none")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.grid(True, which="both", ls="-", alpha=0.15)
    ax.set_xlabel("Degrees of freedom"); ax.set_ylabel("$L_2$ error")
    ax.set_title(f"$L_2$ vs NDOF{tag} — color = P, lightness = $\\omega$, marker = mesh")
    _pw_legends(ax, mesh_marks=MARK)
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=160); plt.close(fig)

OM5 = [0, 4, 16, 64, 256]
draw_l2_vs_h(OMEGAS, "l2_vs_invh.png")                              # 1/h (paper convention)
draw_l2_vs_h(OM5,    "l2_vs_invh_5omega.png", tag=" (5 $\\omega$)")
draw_l2_vs_h(OMEGAS, "l2_vs_h.png", invert=False)                  # h (kept for reference)
draw_l2_vs_h(OM5,    "l2_vs_h_5omega.png", tag=" (5 $\\omega$)", invert=False)
draw_ndof_markers(OMEGAS, "ndof_markers.png")
draw_ndof_markers(OM5,    "ndof_markers_5omega.png", tag=" (5 $\\omega$)")

# ── single-omega versions (cleanest for reading convergence order) ────────────
def _pcol(P):
    return PCMAP[P](0.72)

def draw_l2_vs_h_single(om, fname, invert=True):
    xv = {m: (1.0 / h_of(m) if invert else h_of(m)) for m in MESHES}
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    for P in PS:
        pts = sorted((xv[m], D[(m, P, om)]) for m in MESHES if (m, P, om) in D)
        if len(pts) < 2:
            continue
        xs, ys = zip(*pts)
        o = conv_order(om, P)
        ax.plot(xs, ys, "-o", color=_pcol(P), lw=2, ms=6.5,
                label=f"P{P}   (order {o:.2f})")
        im = len(xs) // 2                              # slope label on the line
        ax.annotate(f"slope {o:.2f}", (xs[im], ys[im]), textcoords="offset points",
                    xytext=(6, 9), fontsize=9.5, color=_pcol(P), fontweight="bold", zorder=5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([xv[m] for m in MESHES])
    ax.set_xticklabels([f"{xv[m]:.2f}\nN={NAZ[m]}" for m in MESHES])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.grid(True, which="both", ls="-", alpha=0.15)
    ax.set_xlabel("$1/h$  (inverse element spacing)" if invert else "$h$")
    ax.set_ylabel("$L_2$ error")
    ax.set_title(f"$L_2$ vs {'1/h' if invert else 'h'} — $\\omega$ = {om}   (|slope| = order)")
    ax.legend(title="best-fit order = |slope|", fontsize=10, frameon=False,
              loc="upper right" if invert else "upper left")
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=160); plt.close(fig)

def draw_ndof_markers_single(om, fname):
    fig, ax = plt.subplots(figsize=(7.6, 5.8))
    for m in MESHES:                                  # dotted gray p-refinement
        pts = sorted((ndof(m, P), D[(m, P, om)]) for P in PS if (m, P, om) in D)
        if len(pts) >= 2:
            xs, ys = zip(*pts)
            ax.plot(xs, ys, ls=":", color="0.55", lw=1.0, alpha=0.6, zorder=1)
    for P in PS:                                      # solid P-lines h-refinement
        pts = sorted((ndof(m, P), D[(m, P, om)], m) for m in MESHES if (m, P, om) in D)
        if len(pts) < 2:
            continue
        c, o = _pcol(P), conv_order(om, P)
        xln = [p[0] for p in pts]; yln = [p[1] for p in pts]
        ax.plot(xln, yln, "-", color=c, lw=1.6, zorder=2, label=f"P{P}   (order {o:.2f})")
        for x, y, m in pts:
            ax.scatter([x], [y], marker=MARK[m], color=c, s=46, zorder=3, edgecolors="none")
        im = len(xln) // 2                             # slope label on the line
        ax.annotate(f"slope {o:.2f}", (xln[im], yln[im]), textcoords="offset points",
                    xytext=(7, 9), fontsize=9.5, color=c, fontweight="bold", zorder=5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.grid(True, which="both", ls="-", alpha=0.15)
    ax.set_xlabel("Degrees of freedom"); ax.set_ylabel("$L_2$ error")
    ax.set_title(f"$L_2$ vs NDOF — $\\omega$ = {om}   (solid: h-refine, dotted gray: p-refine)")
    l1 = ax.legend(title="best-fit order = |slope|", fontsize=10, frameon=False, loc="upper right")
    ax.add_artist(l1)
    ml = [Line2D([], [], color="0.35", marker=MARK[m], lw=0, ms=8, label=f"N={NAZ[m]}") for m in MESHES]
    ax.legend(handles=ml, title="mesh (marker)", fontsize=8, frameon=False, loc="lower left")
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=160); plt.close(fig)

for om in (0, 16, 256):
    draw_l2_vs_h_single(om, f"single_l2_vs_invh_om{om}.png")
    draw_ndof_markers_single(om, f"single_ndof_markers_om{om}.png")

# ── deep48: L2 vs P, one curve per omega ──────────────────────────────────────
fig, ax = plt.subplots(figsize=(6.4, 5.2))
cmap = plt.cm.viridis(np.linspace(0, 1, len(OMEGAS)))
for c, om in zip(cmap, OMEGAS):
    xy = [(P, D[("deep48", P, om)]) for P in PS if ("deep48", P, om) in D]
    if len(xy) < 1: continue
    xs, ys = zip(*xy)
    ax.plot(xs, ys, "-o", color=c, lw=1.8, ms=5.5, label=f"{om}")
ax.set_yscale("log")
ax.set_xticks(PS)
ax.grid(True, which="both", ls="-", alpha=0.15)
ax.set_xlabel("Polynomial order $P$"); ax.set_ylabel("$L_2$ error")
ax.set_title("deep48 — p-convergence at each rotation rate $\\omega$")
ax.legend(title="$\\omega$", fontsize=8, frameon=False, ncol=2)
fig.tight_layout(); fig.savefig(OUT / "pconv_deep48.png", dpi=160); plt.close(fig)

print("wrote figures to", OUT)
for f in sorted(OUT.glob("*.png")):
    print("  ", f.name)

# ── convergence-order table (h-refinement, error ~ h^order) ───────────────────
with open(HERE / "convergence_orders.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["omega", "P", "order_lsq", "order_24_36", "order_36_48"])
    for om in OMEGAS:
        for P in PS:
            w.writerow([om, P,
                        f"{conv_order(om, P):.3f}",
                        f"{order_pair(om, P, 'deep24', 'deep36'):.3f}",
                        f"{order_pair(om, P, 'deep36', 'deep48'):.3f}"])

print("\nOrder of convergence  (LSQ fit over deep24/36/48, error ~ h^order)")
print(f"{'omega':>6} |" + "".join(f" {'P'+str(P):>6}" for P in PS))
print("-" * (8 + 7 * len(PS)))
for om in OMEGAS:
    print(f"{om:>6} |" + "".join(
        f" {conv_order(om, P):>6.2f}" if conv_order(om, P) is not None else f" {'--':>6}"
        for P in PS))
print(f"{'P+1':>6} |" + "".join(f" {P+1:>6}" for P in PS) + "   <- design order")
print("\n(full pairwise breakdown in convergence_orders.csv)")
