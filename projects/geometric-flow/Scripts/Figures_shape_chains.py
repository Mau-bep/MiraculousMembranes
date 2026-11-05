#!/usr/bin/env python3
"""Figs. 11 (c0 = 2.4) and 13 (c0 = 3) of Seifert, Berndl and Lipowsky in the style of the manuscript (figstyle.py), written as png and eps to Results/Paper_figures/<name>_vNN.* (a new
iteration each run, the earlier ones are kept): the bending energy F_b / 8 pi kappa of every branch of vesicle shapes against the reduced volume v, with the shapes on the curves.

The data are those of Analyze_shape_chains.py (the chains of Run_shape_chain.py in Results/Shape_chains/c0_<c0>/, the dashed unstable part of prolate-1 of Shape_axisym.py saddle in
Studies/shape_chains/fig11_data/saddles.csv and fig13_data/saddles.csv); this script only draws them under the rules of the journal: 9 pt Latin Modern, filled markers over 1 mm, lines
over 0.5 pt, no transparency (eps has none). Links that reached the step cap (open dots in Analyze_shape_chains.py) are drawn filled like the others (their status is in points.csv), the grey
diamonds are the points of the paper read off its plot, the squares are the states a branch falls to where it ends, the black diamonds exact values of the theory.

usage:  python3 Figures_shape_chains.py                  both figures: shape_branches_c0_2p4_vNN, shape_branches_c0_3p0_vNN
        python3 Figures_shape_chains.py 3.0              one of them
        python3 Figures_shape_chains.py --cols 2         another width (columns, default 1.5); --height_cm 9 another height
        python3 Figures_shape_chains.py --no_paper       without the grey diamonds of the paper; --no_snapshots without the shapes on the curves; --out_dir folder: write there
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figstyle                                    # noqa: E402  the style of the manuscript
import Analyze_shape_chains as A                   # noqa: E402  the readers of the chains, the snapshots and the points of the paper

REPO = os.path.dirname(HERE)
RESULTS = os.path.join(REPO, "Results")
STUDY = os.path.join(REPO, "Studies", "shape_chains")
COLS = 1.5                                         # width in columns (12.75 cm)
ASPECT = 0.70                                      # height / width, as the 6.6 x 4.6 in figures of Analyze_shape_chains.py
REF_WIDTH_IN = 6.6                                 # the width on which the snapshot sizes and shifts (inches) of Analyze_shape_chains.py were chosen: they scale with the width
SNAP_SCALE = 0.30                                  # inches per unit length of the shape in a snapshot, on that figure
DOT = figstyle.MARKER                              # the dots of the curves (4 pt = 1.4 mm)
SQUARE = 4.0                                       # the square of a branch end, pt
PAPER_MS = 5.0                                     # the grey diamonds of the paper, pt (behind the dots)
GREY = "0.35"
PAPER_GREY = "0.7"
ORDER = ["prolate1", "prolate2", "pear"]
ANCHOR_TEXT = {"sphere": "sphere", "L^pear": r"$L^{\mathrm{pear}}$"}

# per c0: the name of the figure, the dashed segments, where the legend goes (the empty corner) and the spacing of the y ticks
FIGS = {2.4: dict(name="shape_branches_c0_2p4", saddles="fig11_data/saddles.csv", legend="upper right", dy=0.05),
        3.0: dict(name="shape_branches_c0_3p0", saddles="fig13_data/saddles.csv", legend="upper left", dy=0.1)}


def fig_branches(c0, cols=COLS, height_cm=None, root=None, saddles=None, paper=True, snapshots=True):
    """F_b / 8 pi kappa against v for the branches prolate-1, prolate-2 and the pear-shaped one at c0 (2.4 or 3), a dot per relaxation of the chains, the dashed unstable part of prolate-1,
    the shapes of a few links next to the curves (their outline, the meridian cut)."""
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator
    cfg = FIGS[c0]
    root = root or os.path.join(RESULTS, "Shape_chains", "c0_%.1f" % c0)
    parts = A.load_parts(root)
    if not parts:
        raise SystemExit("no chain_summary.csv under %s" % root)
    chains = {b: [row for r in pr for row in r] for b, pr in parts.items()}
    saddles = saddles or os.path.join(STUDY, cfg["saddles"])
    sad = A.load_saddles(saddles) if os.path.exists(saddles) else []
    view = A.VIEW[c0]

    figstyle.use()
    fig, ax = figstyle.figure(cols=cols, aspect=ASPECT, height_cm=height_cm)
    ends = False
    for b in ORDER:
        for rows in parts.get(b, []):
            pts = A.curve(rows)
            if not pts:
                continue
            nu, F = [r["nu"] for r in pts], [r["F"] for r in pts]
            ax.plot(nu, F, "-", color=A.COLORS[b], lw=figstyle.LW, zorder=2)
            ax.plot(nu, F, "o", color=A.COLORS[b], ms=DOT, zorder=3)
            for r in rows:
                if r["kind"] == "jump" and r["F"] != "":                 # the state a branch falls to where it ends
                    ax.plot([r["nu"]], [r["F"]], "s", color=A.COLORS[b], ms=SQUARE, zorder=3)
                    ends = True
    if sad and "prolate1" in chains:
        nu_last = A.curve(chains["prolate1"])[0]["nu"]                      # the dashed line goes on from the lowest link of the 3D chain
        low = [r for r in sad if r["nu"] <= nu_last + 1e-6]
        ax.plot([r["nu"] for r in low], [r["F"] for r in low], "--", color=A.COLORS["prolate1"], lw=figstyle.LW, zorder=2)
    for label, nu, F in A.anchors(c0):
        ax.plot([nu], [F], "D", color="black", ms=figstyle.MARKER, zorder=4)
        right = nu > view["xlim"][0] + 0.85 * (view["xlim"][1] - view["xlim"][0])
        ax.annotate(ANCHOR_TEXT.get(label, label), (nu, F), textcoords="offset points", xytext=(-4, -4) if right else (4, -4), ha="right" if right else "left", va="top")
    if paper:
        for _label, nu, F in A.PAPER.get(round(c0, 6), []):
            ax.plot([nu], [F], "D", color=PAPER_GREY, ms=PAPER_MS, zorder=1)

    ax.set_xlim(*view["xlim"])
    ax.set_ylim(*view["ylim"])
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(cfg["dy"]))
    ax.set_xlabel(figstyle.label("Reduced volume", "v"))
    ax.set_ylabel(figstyle.label("Bending energy", r"F_\mathrm{b}/8\pi\kappa"))
    ax.set_title(r"$c_0 = %g$" % c0)
    handles = [Line2D([], [], color=A.COLORS[b], lw=figstyle.LW, marker="o", ms=DOT) for b in ORDER if b in chains]
    labels = [A.LABELS[b] for b in ORDER if b in chains]
    if ends:
        handles.append(Line2D([], [], color=GREY, ls="None", marker="s", ms=SQUARE))
        labels.append("branch end")
    if sad:
        handles.append(Line2D([], [], color=GREY, lw=figstyle.LW, ls="--"))
        labels.append("unstable (1D solver)")
    if paper:
        handles.append(Line2D([], [], color=PAPER_GREY, ls="None", marker="D", ms=PAPER_MS))
        labels.append("Seifert et al.")
    ax.legend(handles, labels, loc=cfg["legend"], handletextpad=0.4)

    if snapshots:
        fig.canvas.draw()                                                   # the layout is final: the snapshots are placed in figure coordinates
        fig.set_layout_engine("none")
        f = fig.get_size_inches()[0] / REF_WIDTH_IN
        for b in ORDER:
            pts = A.curve(chains.get(b, []))
            for nu, side, shift in A.SNAPSHOTS.get(round(c0, 6), {}).get(b, []):
                r = A.nearest(pts, nu)
                if r is not None:
                    A.add_snapshot(fig, ax, r, (r["nu"], r["F"]), side, shift * f, A.COLORS[b], SNAP_SCALE * f, gap_in=0.10 * f, outline_lw=figstyle.LW * 0.8, line_lw=figstyle.LW_MIN, line_alpha=1.0)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("c0", nargs="*", type=float, help="the figures to make: 2.4 and / or 3.0 (default: both)")
    ap.add_argument("--cols", type=float, default=COLS, help="width in columns (default %g)" % COLS)
    ap.add_argument("--height_cm", type=float, default=None, help="height in cm (default: %g of the width)" % ASPECT)
    ap.add_argument("--no_paper", action="store_true", help="without the grey diamonds (the points of the paper read off its plot)")
    ap.add_argument("--no_snapshots", action="store_true", help="without the shapes on the curves")
    ap.add_argument("--root", type=str, default=None, help="folder of the chains (one c0 only; default Results/Shape_chains/c0_<c0>)")
    ap.add_argument("--saddles", type=str, default=None, help="csv of the dashed segments (one c0 only; default Studies/shape_chains/fig11_data|fig13_data/saddles.csv)")
    ap.add_argument("--out_dir", type=str, default=None, help="folder for the files (default Results/Paper_figures)")
    args = ap.parse_args()
    bad = [c for c in args.c0 if c not in FIGS]
    if bad:
        raise SystemExit("c0 = %s: there are figures for c0 = %s" % (", ".join("%g" % c for c in bad), " and ".join("%g" % c for c in sorted(FIGS))))
    import matplotlib
    matplotlib.use("Agg")
    for c0 in args.c0 or sorted(FIGS):
        fig = fig_branches(c0, cols=args.cols, height_cm=args.height_cm, root=args.root, saddles=args.saddles, paper=not args.no_paper, snapshots=not args.no_snapshots)
        figstyle.save_figure(fig, FIGS[c0]["name"], folder=args.out_dir)


if __name__ == "__main__":
    main()
