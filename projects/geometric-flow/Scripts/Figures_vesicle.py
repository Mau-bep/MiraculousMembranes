#!/usr/bin/env python3
"""The phase space of the vesicle in the style of the manuscript (figstyle.py), written as png and eps to Results/Paper_figures/vesicle_phase_space_vNN.*.

usage:  python3 Figures_vesicle.py                 one column wide, the colour bar on top
        python3 Figures_vesicle.py --cols 1.5      another width (columns); --height_cm 8 another height
        python3 Figures_vesicle.py --cbar right    the colour bar at the side
        python3 Figures_vesicle.py --legs in       only the way along the KB list ("in"), only the way back ("out"), "lowest" (default: the lower Total_E of the two)

The data are the chains of Phase_one_vesicle (Mem_one_vesicle_up.sh: KB scanned at fixed KI, KA = 1; Coverage_data_*.csv of Collect_chain_results.sh), the map of Plots.py phase_one_vesicle_map:
the points are binned (24 bins log-spaced in x, 16 linear in y), the colour of a bin is the mean coverage (CoveredArea) of its points, bins without points are left empty.
x = a sqrt(sigma / K_B),   y = w / sigma = K_I / K_A    (a the bead radius; points with K_A = 0 are left out).
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figstyle                                    # noqa: E402  the style of the manuscript

RESULTS = os.path.join(os.path.dirname(HERE), "Results")
COLS = 1.0                                         # width in columns
BATCH = os.path.join(RESULTS, "Phase_one_vesicle")


def load_points(batch=BATCH, legs="lowest", column="CoveredArea"):
    """x = a sqrt(KA/KB), y = KI/KA and the coloured column of the Coverage_data_*.csv of a chain folder (a glob of csv's also goes). legs "lowest": where a chain has two numbers for the same
    scanned value (the way in and the way back) the one of lower Total_E is kept, the cold start is left out; a tuple of legs ("in", "out", "cold") keeps those rows as they are."""
    files = sorted(glob.glob(os.path.join(batch, "Coverage_data_*.csv") if os.path.isdir(batch) else batch))
    if not files:
        raise SystemExit("no Coverage_data_*.csv in %s" % batch)
    rows, best = [], {}
    for fn in files:
        with open(fn, newline="") as f:
            for r in csv.DictReader(f):
                if float(r["KA"]) <= 0:                                  # w / sigma = KI / KA is undefined at zero tension
                    continue
                if legs == "lowest":
                    if r["leg"] not in ("in", "out"):
                        continue
                    key = (fn, r["scan_var"], float(r["scan_value"]))
                    if key not in best or float(r["Total_E"]) < float(best[key]["Total_E"]):
                        best[key] = r
                elif r["leg"] in legs:
                    rows.append(r)
    rows = list(best.values()) if legs == "lowest" else rows
    KA, KB, KI, r0 = (np.array([float(r[k]) for r in rows]) for k in ("KA", "KB", "KI", "BeadRadius"))
    X, Y, C = r0 * np.sqrt(KA / KB), KI / KA, np.array([float(r[column]) for r in rows])
    good = np.isfinite(X) & np.isfinite(Y) & np.isfinite(C) & (X > 0)
    print("%d points from %d files" % (good.sum(), len(files)))
    return X[good], Y[good], C[good]


def fig_phase_space(cols=COLS, height_cm=None, cbar_location="top", legs="lowest", nx=24, ny=16, batch=BATCH):
    """Coverage over the (a sqrt(sigma/KB), w/sigma) plane: the mean coverage of the points in every bin (log-spaced in x), empty bins blank."""
    X, Y, C = load_points(batch, legs)
    x_edges = np.logspace(np.log10(X.min()), np.log10(X.max()), nx + 1)
    y_edges = np.linspace(Y.min(), Y.max(), ny + 1)
    sums, _, _ = np.histogram2d(X, Y, bins=[x_edges, y_edges], weights=C)
    counts, _, _ = np.histogram2d(X, Y, bins=[x_edges, y_edges])
    mean = np.full_like(sums, np.nan, dtype=float)
    mean[counts > 0] = sums[counts > 0] / counts[counts > 0]

    figstyle.use()
    # the colour bar on top takes about 1.2 cm of height: the figure is taller than 0.8 of the width, as the planar phase space
    fig, ax = figstyle.figure(cols=cols, aspect=0.8 if cbar_location != "top" else 0.93, height_cm=height_cm)
    mesh = ax.pcolormesh(x_edges, y_edges, np.ma.masked_invalid(mean.T), shading="flat", cmap="viridis", vmin=0.0, vmax=1.0, rasterized=True)
    cbar = fig.colorbar(mesh, ax=ax, location=cbar_location, pad=0.02, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
                        **({"aspect": 22} if cbar_location in ("left", "right") else {"aspect": 45, "shrink": 0.7}))
    cbar.set_label("Coverage")
    cbar.outline.set_linewidth(figstyle.AXES_LW)
    cbar.ax.tick_params(width=figstyle.AXES_LW, direction="out")
    ax.set_xscale("log")
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    ax.xaxis.set_major_locator(FixedLocator([1, 10]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: "$%g$" % v))
    ax.yaxis.set_major_locator(figstyle.mpl.ticker.MultipleLocator(1))
    ax.set_xlabel(figstyle.label(symbol=r"a\sqrt{\sigma/K_B}"))
    ax.set_ylabel(figstyle.label(symbol=r"w/\sigma"))
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cols", type=float, default=COLS, help="width in columns (default 1)")
    ap.add_argument("--height_cm", type=float, default=None, help="height in cm (default: 0.93 of the width with the bar on top, 0.8 at the side)")
    ap.add_argument("--cbar", choices=["top", "right"], default="top", help="where the colour bar goes (default top)")
    ap.add_argument("--legs", default="lowest", help="lowest (default), or a comma list of in, out, cold")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    legs = args.legs if args.legs == "lowest" else tuple(args.legs.split(","))
    figstyle.save_figure(fig_phase_space(cols=args.cols, height_cm=args.height_cm, cbar_location=args.cbar, legs=legs), "vesicle_phase_space")


if __name__ == "__main__":
    main()
