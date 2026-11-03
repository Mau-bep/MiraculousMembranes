#!/usr/bin/env python3
"""The two tube figures of the manuscript, each 2/3 of a column wide (5.7 cm) (figstyle.py: 9 pt Latin Modern, lines and markers over the minima), written as png and eps to
Results/Paper_figures/<name>_vNN.* (a new iteration each run, the earlier ones are kept).

usage:  python3 Figures_tube.py                 both figures
        python3 Figures_tube.py pulling         tube_pulling_relaxation: f / f0 against the length dL of the tube, one colour per surface tension
        python3 Figures_tube.py radius          tube_radius_vs_tension: the tube radius r against the surface tension sigma, one series per bending rigidity
        python3 Figures_tube.py radius_ratio    r / r0 against sigma (r0 = sqrt(KB / (4 sigma)) the theory): the points fall on 1
        python3 Figures_tube.py --cols 0.5      another width (columns); --height_cm 5 another height
"""
import argparse
import glob
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figstyle                                    # noqa: E402  the style of the manuscript
import Plot_tube_data as P                         # noqa: E402  the readers of the csv's (chain_summary.csv, Bending.csv)

RESULTS = os.path.join(os.path.dirname(HERE), "Results")
COLS = 2.0 / 3.0                                   # width of both figures, in columns (5.7 cm: 1/3 column, 2.8 cm, cannot hold the 3.2 cm legend entry "Surface tension sigma = 16")
SERIES_COLORS = ["#1b4fb8", "#9a66e0", "#e0559c"]  # blue, purple, pink (colour-blind safe, as Plot_tube_data FILE_COLORS)
RADIUS_COLORS = ["magenta", "#7a7a7a", "black"]    # KB = 1, 2, 4 (Plot_tube_data used magenta, pink, black; pink is too pale at this size, a grey keeps the three apart in print)


def fig_pulling(cols=COLS, height_cm=None, files=None, full_legend=True, size=1.3):
    """f / f0 against the length ΔL of the tube (x_tip of chain_summary.csv), the way out and the way back of the chain (a relaxation at every position), one colour per surface tension
    sigma = KA, f0 = 2 pi sqrt(KB KA) the force of the infinite tube (KB = 4), a dot per relaxation (a scatter plot, dots 1.3 pt = 0.46 mm across: Mau asked for them much smaller than the 1 mm minimum for this figure (the check warns)), the dashed line at f / f0 = 1 (the figure of Plot_tube_data --color same --theory --y F_over_f0 --x x_tip).
    full_legend  the entries read "Surface tension sigma = 16" (3.2 cm wide, about 0.6 of a column of axes); False: "sigma = 16"."""
    files = files or sorted(glob.glob(os.path.join(RESULTS, "Tube_material", "Chain_data_KA_*.csv")), key=lambda f: float(re.search(r"KA_([0-9.]+)", f).group(1)))
    figstyle.use()
    fig, ax = figstyle.figure(cols=cols, aspect=0.9, height_cm=height_cm)
    ax.axhline(1.0, color="black", ls="--", lw=figstyle.AXES_LW, zorder=1)
    for f, c in zip(files, SERIES_COLORS):
        sigma = float(re.search(r"KA_([0-9.]+)", f).group(1).rstrip("."))
        d = P.load_chain_csv(f)
        f0 = P.infinite_tube_force("KA_%g" % sigma, 4.0)
        for leg in ("out", "in"):                  # the way back first, the way out on top; the two lie on each other (no hysteresis)
            sel = d["leg"] == leg
            order = np.argsort(d["link"][sel])
            ax.plot(d["x_tip"][sel][order], (d["F"] / f0)[sel][order], "o", color=c, ms=size, mew=0, zorder=2,
                    label=((r"Surface tension $\sigma = %g$" if full_legend else r"$\sigma = %g$") % sigma) if leg == "in" else None)
    ax.set_xlabel(figstyle.label("Membrane elongation", r"\Delta L"))
    ax.set_ylabel(figstyle.label(symbol="f/f_0"))
    ax.set_xlim(0, 6.5)
    ax.xaxis.set_major_locator(figstyle.mpl.ticker.MultipleLocator(2))
    ax.set_ylim(0, 1.5)
    ax.yaxis.set_major_locator(figstyle.mpl.ticker.MultipleLocator(0.5))
    ax.legend(loc="lower right", handletextpad=0.1, markerscale=2.2)          # the legend dots stay about 1 mm so the colours can be told apart          # the empty corner: below the plateau, right of the rise
    return fig


def fig_radius(cols=COLS, height_cm=None, path=None):
    """The tube radius r against the surface tension sigma (= KA) of the barbell with beads of radius 1, one series per KB, log-log; the dashed lines are r0 = sqrt(KB / (4 sigma))."""
    path = path or os.path.join(RESULTS, "Tube_beads", "Bending.csv")
    rows = P.load_radius_csv(path)
    figstyle.use()
    fig, ax = figstyle.figure(cols=cols, aspect=1.0, height_cm=height_cm)
    for kb, c in zip(sorted({r["KB"] for r in rows}), RADIUS_COLORS):
        sel = sorted((r for r in rows if r["KB"] == kb), key=lambda r: r["KA"])
        xs = np.array([r["KA"] for r in sel])
        ax.plot(xs, np.sqrt(kb / (4.0 * xs)), "--", color=c, lw=figstyle.LW, zorder=1)
        ax.plot(xs, [r["R_tube"] for r in sel], "o", color=c, ms=figstyle.MARKER, label=r"$K_B = %g$" % kb, zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    ax.xaxis.set_major_locator(FixedLocator([2, 5, 10, 50]))
    ax.yaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5]))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_locator(NullLocator())
        axis.set_major_formatter(FuncFormatter(lambda v, _p: "$%g$" % v))
    ax.set_xlabel(figstyle.label("Surface tension", r"\sigma"))
    ax.set_ylabel(figstyle.label("Tube radius", "r"))
    ax.legend(loc="lower left")
    return fig


def fig_radius_ratio(cols=COLS, height_cm=None, path=None):
    """r / r0 against the surface tension sigma (= KA), r0 = sqrt(KB / (4 sigma)) the theoretical radius: the three series of fig_radius fall on 1 (dashed line)."""
    path = path or os.path.join(RESULTS, "Tube_beads", "Bending.csv")
    rows = P.load_radius_csv(path)
    figstyle.use()
    fig, ax = figstyle.figure(cols=cols, aspect=1.0, height_cm=height_cm)
    ax.axhline(1.0, color="black", ls="--", lw=figstyle.AXES_LW, zorder=1)
    for kb, c in zip(sorted({r["KB"] for r in rows}), RADIUS_COLORS):
        sel = sorted((r for r in rows if r["KB"] == kb), key=lambda r: r["KA"])
        xs = np.array([r["KA"] for r in sel])
        ax.plot(xs, np.array([r["R_tube"] for r in sel]) / np.sqrt(kb / (4.0 * xs)), "o", color=c, ms=figstyle.MARKER, label=r"$K_B = %g$" % kb, zorder=2)
    ax.set_xscale("log")
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    ax.xaxis.set_major_locator(FixedLocator([2, 5, 10, 50]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: "$%g$" % v))
    ax.set_xlabel(figstyle.label("Surface tension", r"\sigma"))
    ax.set_ylabel(figstyle.label("Tube radius", "r/r_0"))
    ax.set_ylim(0.98, 1.005)
    ax.yaxis.set_major_locator(figstyle.mpl.ticker.MultipleLocator(0.01))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: "$%.2f$" % v))
    ax.legend(loc="lower right")                          # the empty corner: the points of the three series lie within 1 percent of 1
    return fig


FIGURES = {"pulling": ("tube_pulling_relaxation", fig_pulling), "radius": ("tube_radius_vs_tension", fig_radius), "radius_ratio": ("tube_radius_ratio_vs_tension", fig_radius_ratio)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("which", nargs="*", choices=list(FIGURES) + [[]], default=[], help="the figures to make (default: both)")
    ap.add_argument("--cols", type=float, default=COLS, help="width in columns (default 2/3)")
    ap.add_argument("--short_legend", action="store_true", help="pulling figure: the legend entries read sigma = 16 only (the full Surface tension sigma = 16 is 3.2 cm wide)")
    ap.add_argument("--height_cm", type=float, default=None, help="height in cm (default: square)")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    for key in args.which or list(FIGURES):
        name, make = FIGURES[key]
        fig = make(cols=args.cols, height_cm=args.height_cm, **({"full_legend": not args.short_legend} if key == "pulling" else {}))
        figstyle.save_figure(fig, name)


if __name__ == "__main__":
    main()
