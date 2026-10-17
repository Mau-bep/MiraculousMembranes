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
without Chain_data_ (Chain_data_KA_16.csv -> KA_16) unless you give label=file.

x and y are any column of the files or one of the derived columns (per series, from its first point, link 0):
    F          the pulling force, -Fx (the default y: positive, the plateau of the tube is +f0; Fx itself is the force on the bead as chain_summary.csv has it, negative)
    E_rel      E_total - E_total(first point)         Eb_rel, Et_rel, Ebead_rel: the same for E_bend, E_tension, E_bead
    A_rel      Area - Area(first point)
    tip_gap    x_tip - x, how far the membrane reaches beyond the centre of the bead
--color leg   the way out and the way back in two colours (the series of several files are told apart by the shape of the marker)
--color same  one colour per file, the way out and the way back alike (the shape of the marker still says the way back with --leg_marker)
--theory      for a force (Fx, F) a dotted line at the force f0 = 2 pi sqrt(KB KA) of the infinite tube, KA read from the label (..._KA_16...), KB from --KB (default 4)
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
}
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


def resolve_out(out, files):
    """A figure name without a directory goes to the directory of the (first) csv, a name with a directory (also ./name.png) is used as it is."""
    if out and not os.path.dirname(out) and not out.startswith("."):
        return os.path.join(os.path.dirname(os.path.abspath(files[0])), out)
    return out


def plot_chain_data(files, x="x", y="F", color="leg", labels=None, kinds=("point",), out=None, show=False, size=14, alpha=0.85, lines=False, theory=False, KB=4.0,
                    xlim=None, ylim=None, title=None, leg_marker=False, legs=("in", "out"), ax=None):
    """Scatter plot of column y against column x of one or more chain_summary.csv files.

    files    paths (a path, or a list of them); labels a list of names (default from the file names)
    color    "leg": the way out and the way back in different colours (the files are told apart by the marker), "same": one colour per file for both legs
    kinds    the kinds of row that are points: ("point",) the positions of the list; add "waypoint" and "mid" for every relaxation
    legs     which legs to draw: ("in", "out") both, ("in",) only the way out, ("out",) only the way back
    out      file to write the figure to: a bare name (F.png) is written next to the first csv, a name with a directory as given; show=True opens a window; ax draws into an existing axes
    theory   a dotted line at the force of the infinite tube, 2 pi sqrt(KB KA), for y Fx or F (KA from the label, KB the given one)
    Returns (fig, ax, series) with series = {label: dict of arrays}."""
    import matplotlib
    if out and not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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
        for col in (x, y):
            if col not in d or col in TEXT_COLUMNS:
                raise SystemExit("%r is not a numeric column of %s: use one of %s" % (col, lab, sorted(c for c in d if c not in TEXT_COLUMNS)))
    if color not in ("leg", "same"):
        raise SystemExit("--color is leg or same, not %r" % color)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.4, 5.2))
    else:
        fig = ax.figure
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    leg_color = {"in": "tab:blue", "out": "tab:red"}
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
                name = "%s, %s" % (lab, LEG_NAME[leg]) if len(series) > 1 else LEG_NAME[leg]
            else:
                c = cycle[si % len(cycle)]
                mk = ("o" if leg == "in" else "s") if leg_marker else "o"
                name = (lab + ", " + LEG_NAME[leg]) if leg_marker else (lab if not labelled else None)
                labelled = True
            ax.scatter(d[x][sel], d[y][sel], s=size * (2.4 if (color == "leg" and leg == "out" and "in" in legs) else 1.0), c=c, marker=mk, alpha=alpha, label=name, edgecolors="none")
            if lines:
                order = np.argsort(d["link"][sel]) if "link" in d else np.arange(sel.sum())
                ax.plot(d[x][sel][order], d[y][sel][order], "-", color=c, lw=0.6, alpha=0.5)
        if theory and y in ("Fx", "F"):
            m = re.search(r"KA[_=]?([0-9.]+)", lab)
            if m:
                f0 = 2.0 * math.pi * math.sqrt(float(KB) * float(m.group(1).rstrip(".")))
                ax.axhline(-f0 if y == "Fx" else f0, color=cycle[si % len(cycle)] if color == "same" else "0.4", ls=":", lw=0.9,
                           label="f0 of the infinite tube, %s" % lab if color == "same" else ("f0 of the infinite tube" if si == 0 else None))
    labels_axis = {"x": "x of the bead", "Fx": "force on the bead along x", "F": "pulling force  -Fx", "E_rel": "E_total - E_total(first)", "A_rel": "area - area(first)",
                   "tip_gap": "x_tip - x"}
    ax.set_xlabel(labels_axis.get(x, x))
    ax.set_ylabel(labels_axis.get(y, y))
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
        ax.legend(fontsize=8)
    ax.set_title(title if title is not None else "%s against %s (%s)" % (y, x, ", ".join(kinds)), fontsize=10)
    fig.tight_layout()
    if out:
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        fig.savefig(out, dpi=150)
        print("wrote %s" % out)
    if show:
        plt.show()
    return fig, ax, series


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", help="[label=]chain_summary.csv (wildcards are expanded; default: ../Results/Tube_material/Chain_data_KA_*.csv)")
    ap.add_argument("--x", default="x", help="column of the x axis (default x, the position of the bead)")
    ap.add_argument("--y", default="F", help="column of the y axis (default F = -Fx, the pulling force on the bead: positive; Fx is the force as the chain writes it, negative)")
    ap.add_argument("--color", choices=["leg", "same"], default="leg", help="leg: the way out and the way back in different colours; same: one colour per file")
    ap.add_argument("--leg_marker", action="store_true", help="with --color same: circles for the way out and squares for the way back")
    ap.add_argument("--legs", default="in,out", help="which legs to draw: in,out (default), in (the way out only) or out (the way back only)")
    ap.add_argument("--kinds", default="point", help="kinds of row to draw (default point: the positions of the list; waypoint and mid are the relaxations in between)")
    ap.add_argument("--lines", action="store_true", help="also draw a thin line through the points of each leg, in the order of the chain")
    ap.add_argument("--theory", action="store_true", help="for a force: a dotted line at the force 2 pi sqrt(KB KA) of the infinite tube (KA from the label KA_16)")
    ap.add_argument("--KB", type=float, default=4.0, help="KB for --theory (default 4)")
    ap.add_argument("--size", type=float, default=14.0, help="marker size")
    ap.add_argument("--alpha", type=float, default=0.85)
    ap.add_argument("--xlim", type=float, nargs=2, default=None)
    ap.add_argument("--ylim", type=float, nargs=2, default=None)
    ap.add_argument("--title", default=None)
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
    out = args.out if args.out or args.show else "tube_%s_vs_%s_%s.png" % (args.y, args.x, args.color)
    plot_chain_data(files, x=args.x, y=args.y, color=args.color, labels=labels, kinds=tuple(args.kinds.split(",")), out=out, show=args.show, size=args.size, alpha=args.alpha,
                    lines=args.lines, theory=args.theory, KB=args.KB, xlim=args.xlim, ylim=args.ylim, title=args.title, leg_marker=args.leg_marker, legs=tuple(args.legs.split(",")))


if __name__ == "__main__":
    main()
