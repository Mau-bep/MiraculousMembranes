#!/usr/bin/env python3
"""Scatter plots of the chain summaries of the tube chains (Run_tube_material_chain.py): any column against any column, the way out and the way back in different colours or in one.

usage:  python3 Plot_tube_data.py [options] [label=]CHAIN_SUMMARY.csv [[label=]CHAIN_SUMMARY.csv ...]
        python3 Plot_tube_data.py --x x --y F --color leg ../Results/Tube_material/Chain_data_KA_*.csv
        python3 Plot_tube_data.py --x x --y E_rel --color same --out E_vs_x.png a=KA16.csv b=KA32.csv
        python3 Plot_tube_data.py --list_columns

As a function: from Plot_tube_data import plot_chain_data; plot_chain_data(["Chain_data_KA_16.csv", "Chain_data_KA_32.csv"], x="x", y="F", color="leg", out="F.png")

The files are the chain_summary.csv of a chain (they can be renamed, nothing else of the chain folder is needed): one row per relaxation with the columns link, kind, leg, x, run,
rundir, step, seconds, nV, x_tip, d_bead, Area, E_total, E_bend, E_tension, E_bead, Fx, Fy, Fz, capped. "The correct links" are the rows of kind point: the positions of the list
(link 0, the start of the way out, and every position of the way out and of the way back). The waypoints between two positions and the half move of a retry (kind mid) are
relaxations on the way and are left out (--kinds point,waypoint,mid puts them back); when a link has more than one point row (a chain run again) the last one is kept. The leg is
"in" for the way out (0.2 -> 6.0) and "out" for the way back (6.0 -> 0.0), as in chain_summary.csv: the way to larger x is called in there. The label of a series is the file name
without Chain_data_ (Chain_data_KA_16.csv -> KA_16) unless you give label=file; a label with the KA is written sigma = 16 in the legend.

x and y are any column of the files or one of the derived columns (per series, from its first point, link 0):
    F          the pulling force, -Fx (the default y: positive, the plateau of the tube is +f0; Fx itself is the force on the bead as chain_summary.csv has it, negative)
    E_rel      E_total - E_total(first point)         Eb_rel, Et_rel, Ebead_rel: the same for E_bend, E_tension, E_bead
    A_rel      Area - Area(first point)
    tip_gap    x_tip - x, how far the membrane reaches beyond the centre of the bead
    F_over_f0  F / f0, f0 = 2 pi sqrt(KB KA) of the infinite tube (KA from the label ..._KA_16..., KB from --KB, default 4): every plateau is at 1
--color leg   the way out and the way back in two colours (the series of several files are told apart by the shape of the marker)
--color same  one colour per file, the way out and the way back alike (the shape of the marker still says the way back with --leg_marker): blue, purple, pink for the
              first three files (--colors changes them), validated for colour-blind readers; --color leg: the way out is blue, the way back pink
--theory      for a force (Fx, F) a dotted line at the force f0 = 2 pi sqrt(KB KA) of the infinite tube, KA read from the label (..._KA_16...), KB from --KB (default 4);
              for F_over_f0 one black dashed line at 1 (the f0 of every file)
--plain       without the LaTeX figure set-up of Plots.py (the default); --fontsize 22 the axis labels; no title unless --title
--y F,F_over_f0   several y at once: one figure each (the two versions of the force plot), with --out the name gets _<y> before the extension
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

import numpy as np

DERIVED = {
    "F": "the pulling force, -Fx",
    "E_rel": "E_total - E_total(first point of the series)",
    "Eb_rel": "E_bend - E_bend(first point)",
    "Et_rel": "E_tension - E_tension(first point)",
    "Ebead_rel": "E_bead - E_bead(first point)",
    "A_rel": "Area - Area(first point)",
    "tip_gap": "x_tip - x",
    "F_over_f0": "F / f0 of the infinite tube (KA from the label, KB from --KB)",
}
FILE_COLORS = ["#1b4fb8", "#9a66e0", "#e0559c"]      # blue, purple, pink: pass the CVD / normal-vision / contrast checks of the dataviz validator (also all pairs)
TEXT_COLUMNS = ("kind", "leg", "run", "rundir")
LEG_NAME = {"in": "way out", "out": "way back"}      # chain_summary.csv calls the way to larger x "in" and the way back "out"


def default_label(path):
    base = os.path.splitext(os.path.basename(path))[0]
    return re.sub(r"^Chain_data_", "", base)


def load_chain_csv(path, kinds=("point",)):
    """The rows of one chain_summary.csv as a dict of numpy arrays (numbers as floats, nan where blank), only the given kinds, one row per (kind, link) (the last wins), sorted by link,
    with the derived columns."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("%s has no rows" % path)
    if "x" not in rows[0] or "Fx" not in rows[0]:
        raise SystemExit("%s is not a chain_summary.csv of a tube chain (no x and Fx columns): the columns are %s" % (path, list(rows[0])))
    if "kind" in rows[0]:
        rows = [r for r in rows if r.get("kind") in kinds]
    keep = {}
    for r in rows:                                                    # a chain run again can repeat a link: the last row of the link wins
        keep[(r.get("kind", "point"), r.get("link"))] = r
    rows = sorted(keep.values(), key=lambda r: (float(r["link"]) if r.get("link") not in (None, "") else 0.0, 0 if r.get("kind") == "waypoint" else 1))
    if not rows:
        raise SystemExit("%s has no rows of the kinds %s" % (path, ",".join(kinds)))

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return float("nan")

    data = {}
    for key in rows[0]:
        data[key] = np.array([r[key] for r in rows], dtype=object) if key in TEXT_COLUMNS else np.array([num(r[key]) for r in rows])
    first = next((i for i, r in enumerate(rows) if r.get("kind", "point") == "point"), 0)       # the first point row: link 0, the start of the way out
    data["F"] = -data["Fx"]
    for name, col in (("E_rel", "E_total"), ("Eb_rel", "E_bend"), ("Et_rel", "E_tension"), ("Ebead_rel", "E_bead"), ("A_rel", "Area")):
        if col in data:
            data[name] = data[col] - data[col][first]
    if "x_tip" in data:
        data["tip_gap"] = data["x_tip"] - data["x"]
    return data


# the figure set-up of Plots.py (LaTeX text, Latin Modern), copied from its first lines
def load_pretty_figure_setup():
    import matplotlib as mpl
    # if mpl seems to work with the wrong latex executable, uncomment and adapt path according to your installation
    # os.environ["PATH"] = os.path.expanduser("~/texlive/bin/x86_64-linux") + ":" + os.environ["PATH"]
    
    _preamble_shared = R"""
        \usepackage{graphicx}
        \DeclareMathOperator{\arcsinh}{arcsinh}
        \DeclareMathOperator{\km}{k_\mathrm{m}}
        \DeclareMathOperator{\fbi}{f_\mathrm{bi}}
        \DeclareMathOperator{\eps}{\epsilon_{\mathrm{mc}}}
        \DeclareMathOperator{\epscrit}{\epsilon_{\mathrm{mc}}^*}
        \DeclareMathOperator{\uf}{u_{\mathrm{f}}}
        \DeclareMathOperator{\kBT}{k_\mathrm{B}T}
        """[
        1:
    ]

    def mpl_rcParams_avenir():
        rcParams = {}
        rcParams["font.family"] = "sans-serif"
        rcParams["font.cursive"] = ["Optima"]
        rcParams["text.usetex"] = True
        # rcParams['text.latex.unicode']= True
        rcParams["pgf.texsystem"] = "lualatex"
        rcParams["pgf.rcfonts"] = False
        rcParams["pgf.preamble"] = (
            R"""
        \usepackage[utf8x]{inputenc}
        \usepackage[T1]{fontenc}
        \usepackage{fontspec}
        \usepackage{amsmath}
        \setmainfont{Avenir}[Scale=.9]
        \renewcommand{\setmainfont}{}
        \renewcommand{\sffamily}{}
        """[
                1:
            ]
            + "\n"
            + _preamble_shared
        )
        return rcParams
    
    def rc_params_setup():
        mpl.rcParams["font.family"] = "serif"
        mpl.rcParams["text.usetex"] = True
        mpl.rcParams["figure.constrained_layout.use"] = True
        mpl.rcParams.update(mpl_rcParams_avenir())
        # mpl.rcParams["pgf.texsystem"] = "lualatex"
        # mpl.rcParams["text.latex.preamble"] = mpl.rcParams['pgf.preamble'] #R"\usepackage{amsmath}\usepackage{lmodern}"
        mpl.rcParams["text.latex.preamble"] = (
            R"""
        \usepackage{lmodern}
        \usepackage{amsmath}
        """
            + "\n"
            + _preamble_shared
        )

    rc_params_setup()
    print("Pretty figure set-up loaded.")


def infinite_tube_force(label, KB):
    """f0 = 2 pi sqrt(KB KA) of the infinite tube, KA read from the label (..._KA_16...); None when the label has no KA."""
    m = re.search(r"KA[_=]?([0-9.]+)", label)
    return 2.0 * math.pi * math.sqrt(float(KB) * float(m.group(1).rstrip("."))) if m else None


def legend_name(label):
    """The name of a series in the legend: a label with the KA (KA_16, KA=16) is shown as sigma = 16 (KA is the surface tension sigma), any other label as it is."""
    m = re.search(r"KA[_=]?([0-9.]+)", label)
    return r"$\sigma = %g$" % float(m.group(1).rstrip(".")) if m else tex_safe(label)


def tex_safe(s):
    """With LaTeX text a bare _ or % in a label (x_tip, E_total) is an error: escape them unless the string is maths."""
    import matplotlib as mpl
    if mpl.rcParams["text.usetex"] and "$" not in s:
        return s.replace("_", r"\_").replace("%", r"\%").replace("&", r"\&").replace("#", r"\#")
    return s


def resolve_out(out, files):
    """A figure name without a directory goes to the directory of the (first) csv, a name with a directory (also ./name.png) is used as it is."""
    if out and not os.path.dirname(out) and not out.startswith("."):
        return os.path.join(os.path.dirname(os.path.abspath(files[0])), out)
    return out


def plot_chain_data(files, x="x", y="F", color="leg", labels=None, kinds=("point",), out=None, show=False, size=14, alpha=0.85, lines=False, theory=False, KB=4.0,
                    xlim=None, ylim=None, title=None, leg_marker=False, legs=("in", "out"), ax=None, colors=None, pretty=True, fontsize=22):
    """Scatter plot of column y against column x of one or more chain_summary.csv files.

    files    paths (a path, or a list of them); labels a list of names (default from the file names)
    color    "leg": the way out and the way back in different colours (the files are told apart by the marker), "same": one colour per file for both legs
    kinds    the kinds of row that are points: ("point",) the positions of the list; add "waypoint" and "mid" for every relaxation
    legs     which legs to draw: ("in", "out") both, ("in",) only the way out, ("out",) only the way back
    out      file to write the figure to: a bare name (F.png) is written next to the first csv, a name with a directory as given; show=True opens a window; ax draws into an existing axes
    theory   a dotted line at the force of the infinite tube, 2 pi sqrt(KB KA), for y Fx or F (KA from the label, KB the given one); for y F_over_f0 (F / f0) one black dashed line at 1
    pretty   the figure set-up of Plots.py (LaTeX text); fontsize the size of the axis labels (the ticks and the legend are 4 points smaller); title None: no title
    colors   the colours of the files with color="same" (default blue, purple, pink); with color="leg" the way out is blue and the way back pink
    Returns (fig, ax, series) with series = {label: dict of arrays}."""
    import matplotlib
    if out and not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if pretty:
        load_pretty_figure_setup()

    if isinstance(files, (str, os.PathLike)):
        files = [files]
    out = resolve_out(out, list(files))
    labels = list(labels) if labels else [default_label(f) for f in files]
    if len(labels) != len(files):
        raise SystemExit("%d labels for %d files" % (len(labels), len(files)))
    if len(set(labels)) != len(labels):
        raise SystemExit("two series have the same label %s: give them labels, label=file" % labels)
    series = {lab: load_chain_csv(f, kinds) for lab, f in zip(labels, files)}
    for lab, d in series.items():
        f0 = infinite_tube_force(lab, KB)
        if f0 is not None:
            d["F_over_f0"] = d["F"] / f0
        elif "F_over_f0" in (x, y):
            raise SystemExit("F_over_f0 needs the KA in the label of every file (Chain_data_KA_16.csv -> KA_16, or KA_16=file.csv): %r has none" % lab)
    for lab, d in series.items():
        for col in (x, y):
            if col not in d or col in TEXT_COLUMNS:
                raise SystemExit("%r is not a numeric column of %s: use one of %s" % (col, lab, sorted(c for c in d if c not in TEXT_COLUMNS)))
    if color not in ("leg", "same"):
        raise SystemExit("--color is leg or same, not %r" % color)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8.2, 5.8))
    else:
        fig = ax.figure
    cycle = list(colors or FILE_COLORS) + plt.rcParams["axes.prop_cycle"].by_key()["color"]     # blue, purple, pink, then the matplotlib ones
    leg_color = {"in": cycle[0], "out": FILE_COLORS[2]}
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    for si, (lab, d) in enumerate(series.items()):
        labelled = False
        # with two colours the way back is drawn first and larger, the way out on top: where they coincide (no hysteresis) both stay visible
        for leg in (sorted(legs, key=lambda l: l != "out") if color == "leg" else legs):
            sel = d["leg"] == leg if "leg" in d else np.ones(len(d[x]), bool)
            if not sel.any():
                continue
            if color == "leg":
                c, mk = leg_color[leg], markers[si % len(markers)]
                name = "%s, %s" % (legend_name(lab), LEG_NAME[leg]) if len(series) > 1 else LEG_NAME[leg]
            else:
                c = cycle[si % len(cycle)]
                mk = ("o" if leg == "in" else "s") if leg_marker else "o"
                name = (legend_name(lab) + ", " + LEG_NAME[leg]) if leg_marker else (legend_name(lab) if not labelled else None)
                labelled = True
            ax.scatter(d[x][sel], d[y][sel], s=size * (2.4 if (color == "leg" and leg == "out" and "in" in legs) else 1.0), c=c, marker=mk, alpha=alpha, label=name, edgecolors="none")
            if lines:
                order = np.argsort(d["link"][sel]) if "link" in d else np.arange(sel.sum())
                ax.plot(d[x][sel][order], d[y][sel][order], "-", color=c, lw=0.6, alpha=0.5)
        if theory and y in ("Fx", "F"):
            f0 = infinite_tube_force(lab, KB)
            if f0 is not None:
                ax.axhline(-f0 if y == "Fx" else f0, color=cycle[si % len(cycle)] if color == "same" else "0.4", ls=":", lw=0.9)         # no legend entry: its colour is that of its file
    if theory and y == "F_over_f0":                       # f0 of every file is 1 here: one black dashed line, no legend entry
        ax.axhline(1.0, color="black", ls="--", lw=1.1, zorder=1)
    labels_axis = {"x_tip":r"$\Delta L$","x": "x of the bead", "Fx": r"$F_x$", "F": r"$F_x$", "F_over_f0": r"$f/f_0$", "E_rel": "E_total - E_total(first)", "A_rel": "area - area(first)",
                   "tip_gap": "x_tip - x"}
    ax.set_xlabel(tex_safe(labels_axis.get(x, x)), fontsize=fontsize)
    ax.set_ylabel(tex_safe(labels_axis.get(y, y)), fontsize=fontsize)
    ax.tick_params(labelsize=fontsize - 4)
    if xlim:
        ax.set_xlim(*xlim)
        if not ylim:                                  # the y range of the points inside the x range, not of the whole chain
            inside = [d[y][(d[x] >= xlim[0]) & (d[x] <= xlim[1]) & np.isfinite(d[y])] for d in series.values()]
            inside = np.concatenate([v for v in inside if len(v)]) if any(len(v) for v in inside) else np.array([])
            if len(inside):
                pad = 0.05 * (inside.max() - inside.min() or 1.0)
                ax.set_ylim(inside.min() - pad, inside.max() + pad)
    if ylim:
        ax.set_ylim(*ylim)
    ax.grid(alpha=0.25, lw=0.5)
    handles, names = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=fontsize - 4, markerscale=2.0)
    if title is not None:
        ax.set_title(tex_safe(title), fontsize=fontsize)
    if not plt.rcParams["figure.constrained_layout.use"]:
        fig.tight_layout()
    if out:
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        fig.savefig(out, dpi=150)
        print("wrote %s" % out)
    if show:
        plt.show()
    return fig, ax, series


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Tube radius of the barbell with beads of radius 1 (Create_subjob_barbell_beads.py, Mem_barbell_beads_up.sh, Analyze_barbell_beads.py -> Bending.csv, Tension.csv)
BENDING_ICONS = [(4, 4, "*"), (1, 50, "s"), (2, 11.5, "^")]      # (KB, sigma, marker) of the Bending runs marked on the plots: star (KB = 4, sigma = 4), square (KB = 1, sigma = 50), triangle (KB = 2, sigma = 11.5)
SERIES_COLORS = ["magenta", "pink", "black"]          # the colours of Plots.py: one per series (KB in the Bending batch, sigma = KA in the Tension batch)


def load_radius_csv(path):
    """The rows of a csv of Analyze_barbell_beads.py as a list of dicts of floats (the run column and the True/False column are dropped)."""
    with open(path, newline="") as f:
        return [{k: float(v) for k, v in r.items() if k != "run" and v not in ("True", "False", "")} for r in csv.DictReader(f)]


def plot_tube_radius(path, group, out=None, colors=None, loglog=True, fontsize=22, size=7, pretty=True, title=None, ax=None, show=False, legend_loc=None,
                     lines=False, slope_line=False, slope_color="skyblue", theory=True, highlight=None):
    """Tube radius against the scanned constant, one series per value of the fixed one, the theory r_0 = sqrt(KB / (4 sigma)) dashed in the colour of its series (KA is the surface tension sigma).

    path        the csv of one batch (Bending.csv: KB fixed per series, KA scanned; Tension.csv: KA fixed per series, KB scanned)
    group       "KB" (a series per KB, x = sigma: the Bending batch) or "KA" (a series per sigma, x = KB: the Tension batch)
    out         file to write (a bare name goes next to the csv); colors the colour of each series in the order of the fixed values (default magenta, pink, black)
    loglog      log axes (default: the theory is then a straight line of slope -1/2 against sigma, +1/2 against KB); False: linear
    legend_loc  where the legend goes inside the axes (default: the empty corner)
    lines       a thin solid line through the points of each series (default: points only, the dashed theory is then not hidden)
    theory      draw the dashed theory r_0 (and its legend entries); False: only the points (and the slope line)
    highlight   points drawn again with a distinct, larger marker: a list of (fixed value, scanned value, marker), e.g. (4, 4, "*"): the point of the series group == fixed value
                nearest to the scanned value (within 3 percent, else an error); the marker takes the colour of its series
    slope_line  an extra line (slope_color, light blue) with the slope of the power law on the log-log plot, -1/2 against sigma (Bending) or +1/2 against KB (Tension), above the data:
                only its slope counts, it is not the theory
    Returns (fig, ax)."""
    import matplotlib
    if out and not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if pretty:
        load_pretty_figure_setup()
    rows = load_radius_csv(path)
    other = "KA" if group == "KB" else "KB"
    legend_loc = legend_loc or ("lower left" if group == "KB" else "lower right")      # the empty corner: the Bending lines fall to the right, the Tension lines rise to the right
    fixed = sorted({r[group] for r in rows})
    colors = list(colors or SERIES_COLORS)
    if len(fixed) > len(colors):
        colors += plt.rcParams["axes.prop_cycle"].by_key()["color"]
    if ax is None:
        fig, ax = plt.subplots(figsize=(8.2, 5.8))
    else:
        fig = ax.figure
    for g, c in zip(fixed, colors):
        sel = sorted((r for r in rows if r[group] == g), key=lambda r: r[other])
        xs = np.array([r[other] for r in sel])
        name = r"\sigma = %g" % g if group == "KA" else r"K_B = %g" % g          # KA is the surface tension sigma
        if lines:
            ax.plot(xs, [r["R_tube"] for r in sel], "-", color=c, lw=1.0, zorder=1)
        ax.scatter(xs, [r["R_tube"] for r in sel], s=size ** 2, color=c, edgecolors=c, label="$%s$" % name, zorder=2)
        if theory:
            xt = np.geomspace(xs.min(), xs.max(), 200) if loglog else np.linspace(xs.min(), xs.max(), 200)
            KB, KA = (g, xt) if group == "KB" else (xt, g)
            ax.plot(xt, np.sqrt(KB / (4.0 * KA)), "--", color=c, lw=1.8, label=r"$r_0(%s)$" % name, zorder=3)
    for g, xv, mk in (highlight or []):
        if g not in fixed:
            raise SystemExit("highlight: no series with %s = %g in %s (have %s)" % (group, g, path, fixed))
        sel = [r for r in rows if r[group] == g]
        pt = min(sel, key=lambda r: abs(r[other] - xv))
        if abs(pt[other] - xv) > 0.03 * xv:
            raise SystemExit("highlight: no run near %s = %g with %s = %g in %s (nearest %g)" % (other, xv, group, g, path, pt[other]))
        if pt[other] != xv:
            print("highlight: %s = %g, %s = %g: the run is at %s = %g" % (group, g, other, xv, other, pt[other]))
        ax.scatter([pt[other]], [pt["R_tube"]], s=(1.9 * size) ** 2, marker=mk, color=colors[fixed.index(g)], edgecolors=colors[fixed.index(g)], zorder=4)
    if slope_line:
        # a reference line, 1.7 times above the theory of the top series: over x from 4 times the smallest sigma to the largest (Bending), from the smallest KB to a quarter of the largest (Tension)
        allx = np.array([r[other] for r in rows])
        slope = -0.5 if group == "KB" else 0.5
        xa, xb = (allx.min() * 4.0, allx.max()) if group == "KB" else (allx.min(), allx.max() / 4.0)
        top = max(fixed) if group == "KB" else min(fixed)
        r_top_a = math.sqrt(top / (4.0 * xa)) if group == "KB" else math.sqrt(xa / (4.0 * top))
        xg = np.geomspace(xa, xb, 50)
        ax.plot(xg, 1.7 * r_top_a * (xg / xa) ** slope, "-", color=slope_color, lw=2.8, zorder=1, label=r"slope $%s1/2$" % ("-" if slope < 0 else "+"))
    if loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")
        from matplotlib.ticker import FuncFormatter

        def tick_label(v, _pos):          # plain numbers (0.2, 3, 50) at 1, 2, 3, 5 of every decade instead of 2e-01
            m = v / 10.0 ** math.floor(math.log10(v) + 1e-9)
            return ("$%g$" % v) if round(m, 6) in (1.0, 2.0, 3.0, 5.0) else ""      # maths mode: the Latin Modern numbers of the other plots
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(FuncFormatter(tick_label))
            axis.set_minor_formatter(FuncFormatter(tick_label))
    ax.set_xlabel(r"$\sigma$" if other == "KA" else r"$K_B$", fontsize=fontsize)
    ax.set_ylabel(r"$r$", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize - 4, which="both")
    ax.grid(alpha=0.25, lw=0.5, which="both")
    # the legend: the simulation entries together in the first column, then the theory lines, then the slope line
    handles, names = ax.get_legend_handles_labels()
    sim_i = [i for i, n in enumerate(names) if n.startswith("$") and not n.startswith("$r")]
    rest_i = [i for i in range(len(names)) if i not in sim_i]
    hs, ns = [handles[i] for i in sim_i], [names[i] for i in sim_i]
    while len(hs) < len(rest_i) and len(rest_i) > len(sim_i):        # an empty entry keeps the simulation entries in the first column when the second one is longer
        hs.append(plt.Line2D([], [], lw=0))
        ns.append("")
    ax.legend(hs + [handles[i] for i in rest_i], ns + [names[i] for i in rest_i], fontsize=fontsize - 6, ncol=2 if theory else 1, loc=legend_loc, markerscale=0.9, columnspacing=1.2, handletextpad=0.4)
    if title is not None:
        ax.set_title(tex_safe(title), fontsize=fontsize)
    if not plt.rcParams["figure.constrained_layout.use"]:
        fig.tight_layout()
    out = resolve_out(out, [path]) if out else None
    if out:
        fig.savefig(out, dpi=150)
        print("wrote %s" % out)
    if show:
        plt.show()
    return fig, ax


def save_bending_icons(out_dir="../Results/Tube_beads/icons", size_px=400, bending_icons=None):
    """One transparent png per icon of BENDING_ICONS (the marker and colour of the highlighted runs of the Bending plots): star_black.png, square_magenta.png, triangle_pink.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = {"*": "star", "s": "square", "^": "triangle", "o": "circle", "v": "triangle_down", "D": "diamond"}
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for KB, sigma, mk in (bending_icons or BENDING_ICONS):
        colour = SERIES_COLORS[sorted({g for g, _, _ in (bending_icons or BENDING_ICONS)}).index(KB)]
        fig = plt.figure(figsize=(1, 1), dpi=size_px)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.scatter([0.5], [0.5], s=(0.8 * 72) ** 2, marker=mk, color=colour, edgecolors=colour)
        path = os.path.join(out_dir, "%s_%s.png" % (names.get(mk, mk), colour))
        fig.savefig(path, transparent=True)
        plt.close(fig)
        paths.append(path)
        print("wrote %s (KB = %g, sigma = %g)" % (path, KB, sigma))
    return paths


def plot_tube_radius_bending(path="../Results/Tube_beads/Bending.csv", out="Bending_radius.png", **kw):
    """Bending batch: one series per KB (magenta, pink, black), the radius against sigma, the theory dashed; three runs marked with a distinct icon (BENDING_ICONS)."""
    kw.setdefault("highlight", BENDING_ICONS)
    return plot_tube_radius(path, "KB", out=out, **kw)


def plot_tube_radius_tension(path="../Results/Tube_beads/Tension.csv", out="Tension_radius.png", **kw):
    """Tension batch: one series per sigma (magenta, pink, black), the radius against KB, the theory dashed."""
    return plot_tube_radius(path, "KA", out=out, **kw)


def plot_tube_radius_bending_lines(path="../Results/Tube_beads/Bending.csv", out="Bending_radius_lines.png", **kw):
    """The Bending plot with a line through the points of each series (no theory curves) and a light blue reference line of slope -1/2."""
    kw.setdefault("highlight", BENDING_ICONS)
    return plot_tube_radius(path, "KB", out=out, lines=True, slope_line=True, theory=False, **kw)


def plot_tube_radius_tension_lines(path="../Results/Tube_beads/Tension.csv", out="Tension_radius_lines.png", **kw):
    """The Tension plot with a line through the points of each series (no theory curves) and a light blue reference line of slope +1/2."""
    return plot_tube_radius(path, "KA", out=out, lines=True, slope_line=True, theory=False, **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", help="[label=]chain_summary.csv (wildcards are expanded; default: ../Results/Tube_material/Chain_data_KA_*.csv)")
    ap.add_argument("--x", default="x", help="column of the x axis (default x, the position of the bead)")
    ap.add_argument("--y", default="F", help="column of the y axis (default F = -Fx, the pulling force on the bead: positive; Fx is the force as the chain writes it, negative); "
                                             "F_over_f0 = F / f0 of the infinite tube; a comma list (F,F_over_f0) makes one figure each")
    ap.add_argument("--color", choices=["leg", "same"], default="leg", help="leg: the way out and the way back in different colours; same: one colour per file")
    ap.add_argument("--leg_marker", action="store_true", help="with --color same: circles for the way out and squares for the way back")
    ap.add_argument("--legs", default="in,out", help="which legs to draw: in,out (default), in (the way out only) or out (the way back only)")
    ap.add_argument("--kinds", default="point", help="kinds of row to draw (default point: the positions of the list; waypoint and mid are the relaxations in between)")
    ap.add_argument("--lines", action="store_true", help="also draw a thin line through the points of each leg, in the order of the chain")
    ap.add_argument("--theory", action="store_true", help="for a force: a dotted line at the force 2 pi sqrt(KB KA) of the infinite tube (KA from the label KA_16)")
    ap.add_argument("--KB", type=float, default=4.0, help="KB for --theory (default 4)")
    ap.add_argument("--colors", default=None, help="colours of the files with --color same, comma separated (default blue, purple, pink: %s)" % ",".join(FILE_COLORS))
    ap.add_argument("--plain", action="store_true", help="without the figure set-up of Plots.py (LaTeX text): matplotlib's own fonts")
    ap.add_argument("--fontsize", type=float, default=22.0, help="size of the axis labels (default 22; the ticks and the legend are 4 smaller)")
    ap.add_argument("--size", type=float, default=14.0, help="marker size")
    ap.add_argument("--alpha", type=float, default=0.85)
    ap.add_argument("--xlim", type=float, nargs=2, default=None)
    ap.add_argument("--ylim", type=float, nargs=2, default=None)
    ap.add_argument("--title", default=None, help="a title (default: none)")
    ap.add_argument("--out", default=None, help="write the figure here (png, pdf, ...); a bare name goes to the directory of the first csv, and without --out (and --show) the figure is "
                                              "tube_<y>_vs_<x>_<color>.png in the directory of the first csv")
    ap.add_argument("--show", action="store_true", help="open the figure in a window")
    ap.add_argument("--list_columns", action="store_true", help="print the columns of the files and the derived ones and stop")
    args = ap.parse_args()

    if args.list_columns and not args.files:
        print("columns of a chain_summary.csv: link kind leg x run rundir step seconds nV x_tip d_bead Area E_total E_bend E_tension E_bead Fx Fy Fz capped")
        for k, v in DERIVED.items():
            print("derived: %-10s %s" % (k, v))
        return
    items = args.files or sorted(glob.glob(os.path.join("..", "Results", "Tube_material", "Chain_data_KA_*.csv")))
    files, labels = [], []
    for item in items:
        head, eq, tail = item.partition("=")
        label, path = (head, tail) if eq and "/" not in head and not os.path.exists(item) else (None, item)
        for p in (sorted(glob.glob(path)) or [path]):
            if not os.path.exists(p):
                raise SystemExit("no file %s" % p)
            files.append(p)
            labels.append(label if label and len(glob.glob(path)) <= 1 else default_label(p))
    if not files:
        raise SystemExit("no chain_summary.csv files given (and none found in ../Results/Tube_material/)")
    if args.list_columns:
        d = load_chain_csv(files[0], tuple(args.kinds.split(",")))
        print(" ".join(sorted(d)))
        return
    ys = args.y.split(",")
    for y in ys:
        out = args.out if args.out or args.show else "tube_%s_vs_%s_%s.png" % (y, args.x, args.color)
        if args.out and len(ys) > 1:                      # one figure per y: the y goes before the extension of --out
            stem, ext = os.path.splitext(args.out)
            out = "%s_%s%s" % (stem, y, ext)
        plot_chain_data(files, x=args.x, y=y, color=args.color, labels=labels, kinds=tuple(args.kinds.split(",")), out=out, show=args.show, size=args.size, alpha=args.alpha,
                        lines=args.lines, theory=args.theory, KB=args.KB, xlim=args.xlim, ylim=args.ylim, title=args.title, leg_marker=args.leg_marker, legs=tuple(args.legs.split(",")),
                        colors=args.colors.split(",") if args.colors else None, pretty=not args.plain, fontsize=args.fontsize)


if __name__ == "__main__":
    # main()    # the scatter plots of the chain summaries of the tube chains (everything this script did before the barbell radius plots): uncomment to use the command line again
    plot_tube_radius_bending()
    plot_tube_radius_tension()
    plot_tube_radius_bending_lines()
    plot_tube_radius_tension_lines()
