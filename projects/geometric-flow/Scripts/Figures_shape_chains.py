#!/usr/bin/env python3
"""Figs. 11 (c0 = 2.4) and 13 (c0 = 3) of Seifert, Berndl and Lipowsky in the style of the manuscript (figstyle.py), written as png and eps to Results/Paper_figures/<name>_vNN.* (a new
iteration each run, the earlier ones are kept): the bending energy F_b / 8 pi kappa of every branch of vesicle shapes against the reduced volume v, with the shapes on the curves.

The data are those of Analyze_shape_chains.py (the chains of Run_shape_chain.py in Results/Shape_chains/c0_<c0>/, the links picked for the snapshots are its SNAPSHOTS); this script only draws
them under the rules of the journal: 9 pt Latin Modern, filled markers over 1 mm, lines over 0.5 pt, no transparency (eps has none). Links that reached the step cap (open dots in
Analyze_shape_chains.py) are drawn filled like the others (their status is in points.csv); the black diamonds are exact values of the theory (the sphere, the budding limit L^pear).

Each figure is its own file, 1/3 of a column wide (2.83 cm), to be joined afterwards. At that width no legend fits, so the colours are for the caption (blue prolate-1, green prolate-2, red
pear-shaped), the axes are labelled with the symbols alone (v and F_b / 8 pi kappa) and the ticks are few. The curves are the line through every link with a dot about every 2.4 dot
diameters (a dot per link would be a solid band). The shapes of the picked links sit where there is room, each with a thin line to its point, as large as
they can be without lying on a curve or on each other. From one column on the legend, the words of the axis labels, the names of the black diamonds and the dot per link come back, and the
shapes sit next to their points as in the pngs of Analyze_shape_chains.py.
Off unless asked for: --paper (the grey diamonds, the points of the paper read off its plot), --unstable (the dashed unstable part of prolate-1 of Shape_axisym.py saddle, in
Studies/shape_chains/fig11_data|fig13_data/saddles.csv), --ends (a filled square where a branch ends: the state it falls to).

usage:  python3 Figures_shape_chains.py                  both figures, 1/3 column each: shape_branches_c0_2p4_vNN and shape_branches_c0_3p0_vNN
        python3 Figures_shape_chains.py 3.0              one of them (2.4 or 3.0)
        python3 Figures_shape_chains.py --cols 1.5       another width (columns); --height_cm 5 another height
        python3 Figures_shape_chains.py --paper --unstable --ends --cols 1.5      everything, as the pngs of Analyze_shape_chains.py
        python3 Figures_shape_chains.py --no_snapshots   without the shapes; --out_dir folder  write there instead of Results/Paper_figures
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figstyle                                    # noqa: E402  the style of the manuscript
import Analyze_shape_chains as A                   # noqa: E402  the readers of the chains, the snapshots and the points of the paper

REPO = os.path.dirname(HERE)
RESULTS = os.path.join(REPO, "Results")
STUDY = os.path.join(REPO, "Studies", "shape_chains")
COLS = 1.0 / 3.0                                   # width of one figure, in columns (2.83 cm)
REF_WIDTH_IN = 6.6                                 # the width on which the snapshot sizes and shifts (inches) of Analyze_shape_chains.py were chosen: they scale with the width
SNAP_SCALE = 0.30                                  # inches per unit length of the shape in a snapshot, on that figure
SQUARE = 4.0                                       # the square of a branch end, pt
PAPER_MS = 5.0                                     # the grey diamonds of the paper, pt (behind the dots)
NARROW_DOT = 3.4                                   # pt (1.2 mm): the dots of a figure below one column
DOT_PITCH = 2.4                                    # distance between the dots of a narrow figure in dot diameters
NARROW_ASPECT = 2.0                                # height / width below one column (a plot 1.7 cm wide needs room for the shapes)
SNAP_SCALES = (0.22, 0.19, 0.17, 0.15, 0.13, 0.115, 0.10, 0.09)       # inches per unit length tried for the shapes of a narrow figure, the largest that fits is kept
GREY = "0.35"
PAPER_GREY = "0.7"
ORDER = ["prolate1", "prolate2", "pear"]
ANCHOR_TEXT = {"sphere": "sphere", "L^pear": r"$L^{\mathrm{pear}}$"}

# per c0: the name of the figure, the dashed segments, where the legend goes (the empty corner), the spacing of the y ticks (wide, narrow)
FIGS = {2.4: dict(name="shape_branches_c0_2p4", saddles="fig11_data/saddles.csv", legend="upper right", dy=(0.05, 0.1)),
        3.0: dict(name="shape_branches_c0_3p0", saddles="fig13_data/saddles.csv", legend="upper left", dy=(0.1, 0.2))}


# ---- the axes and the curves ------------------------------------------------------------------------------------------------------------------------------
def decorate(ax, c0, cols):
    """Limits, ticks, labels and title of the panel of c0, `cols` columns wide."""
    from matplotlib.ticker import MultipleLocator
    cfg, view = FIGS[c0], A.VIEW[c0]
    ax.set_xlim(*view["xlim"])
    ax.set_ylim(*view["ylim"])
    if cols < 1.0:
        ax.xaxis.set_major_locator(MultipleLocator(0.2))                    # three ticks at most fit under 1.7 cm of axes
        ax.yaxis.set_major_locator(MultipleLocator(cfg["dy"][1]))
        ax.set_xlabel(figstyle.label(symbol="v"))
        ax.set_ylabel(figstyle.label(symbol=r"F_\mathrm{b}/8\pi\kappa"))
    else:
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(cfg["dy"][0]))
        ax.set_xlabel(figstyle.label("Reduced volume", "v"))
        ax.set_ylabel(figstyle.label("Bending energy", r"F_\mathrm{b}/8\pi\kappa"))
    ax.set_title(r"$c_0 = %g$" % c0)


def thin(ax, nu, F, pitch_pt):
    """Indices of the points of a curve to put a dot on: the first, the last and as many between as keep the dots `pitch_pt` apart on the page (the layout must be final)."""
    xy = ax.transData.transform(np.column_stack([nu, F])) * 72.0 / ax.figure.dpi
    keep = [0]
    for i in range(1, len(xy) - 1):
        if np.hypot(*(xy[i] - xy[keep[-1]])) >= pitch_pt:
            keep.append(i)
    if len(xy) > 1:
        if np.hypot(*(xy[-1] - xy[keep[-1]])) < 0.5 * pitch_pt and len(keep) > 1:
            keep[-1] = len(xy) - 1
        else:
            keep.append(len(xy) - 1)
    return keep


def draw_branches(ax, c0, cols, root=None, saddles=None, paper=False, unstable=False, ends=False, labels=True, legend=None):
    """The curves of F_b / 8 pi kappa against v at c0 on ax, a panel `cols` columns wide whose layout is final: the line through every link and the dots (every link from one column on, spaced
    below), the black diamonds, the optional extras and the legend (from one column on). Returns (parts, chains): the rows of every chain folder of every branch, and of every branch."""
    from matplotlib.lines import Line2D
    cfg = FIGS[c0]
    narrow = cols < 1.0
    root = root or os.path.join(RESULTS, "Shape_chains", "c0_%.1f" % c0)
    parts = A.load_parts(root)
    if not parts:
        raise SystemExit("no chain_summary.csv under %s" % root)
    chains = {b: [row for r in pr for row in r] for b, pr in parts.items()}
    view = A.VIEW[c0]
    dot = NARROW_DOT if narrow else figstyle.MARKER
    legend = (not narrow) if legend is None else legend
    ends_drawn = False
    for b in ORDER:
        for rows in parts.get(b, []):
            pts = A.curve(rows)
            if not pts:
                continue
            nu, F = np.array([r["nu"] for r in pts]), np.array([r["F"] for r in pts])
            ax.plot(nu, F, "-", color=A.COLORS[b], lw=figstyle.LW, zorder=2)
            k = thin(ax, nu, F, DOT_PITCH * dot) if narrow else slice(None)
            ax.plot(nu[k], F[k], "o", color=A.COLORS[b], ms=dot, zorder=3)
            if ends:
                for r in rows:
                    if r["kind"] == "jump" and r["F"] != "":                 # the state a branch falls to where it ends
                        ax.plot([r["nu"]], [r["F"]], "s", color=A.COLORS[b], ms=SQUARE, zorder=3)
                        ends_drawn = True
    sad = []
    if unstable:
        path = saddles or os.path.join(STUDY, cfg["saddles"])
        sad = A.load_saddles(path) if os.path.exists(path) else []
        if sad and "prolate1" in chains:
            nu_last = A.curve(chains["prolate1"])[0]["nu"]                  # the dashed line goes on from the lowest link of the 3D chain
            low = [r for r in sad if r["nu"] <= nu_last + 1e-6]
            ax.plot([r["nu"] for r in low], [r["F"] for r in low], "--", color=A.COLORS["prolate1"], lw=figstyle.LW, zorder=2)
    for label, nu, F in A.anchors(c0):
        ax.plot([nu], [F], "D", color="black", ms=max(dot, 3.0), zorder=4, clip_on=False)      # the sphere sits on the right edge of the axes: whole
        if labels and not narrow:
            right = nu > view["xlim"][0] + 0.85 * (view["xlim"][1] - view["xlim"][0])
            ax.annotate(ANCHOR_TEXT.get(label, label), (nu, F), textcoords="offset points", xytext=(-4, -4) if right else (4, -4), ha="right" if right else "left", va="top")
    if paper:
        for _label, nu, F in A.PAPER.get(round(c0, 6), []):
            ax.plot([nu], [F], "D", color=PAPER_GREY, ms=PAPER_MS, zorder=1)
    if legend:
        handles = [Line2D([], [], color=A.COLORS[b], lw=figstyle.LW, marker="o", ms=dot) for b in ORDER if b in chains]
        names = [A.LABELS[b] for b in ORDER if b in chains]
        if ends_drawn:
            handles.append(Line2D([], [], color=GREY, ls="None", marker="s", ms=SQUARE))
            names.append("branch end")
        if sad:
            handles.append(Line2D([], [], color=GREY, lw=figstyle.LW, ls="--"))
            names.append("unstable (1D solver)")
        if paper:
            handles.append(Line2D([], [], color=PAPER_GREY, ls="None", marker="D", ms=PAPER_MS))
            names.append("Seifert et al.")
        ax.legend(handles, names, loc=cfg["legend"], handletextpad=0.4)
    return parts, chains


# ---- the shapes on the curves -------------------------------------------------------------------------------------------------------------------------------
def picked(c0, chains):
    """The links picked for the snapshots (Analyze_shape_chains.SNAPSHOTS): [(branch, row)]."""
    out = []
    for b in ORDER:
        pts = A.curve(chains.get(b, []))
        for nu, _side, _shift in A.SNAPSHOTS.get(round(c0, 6), {}).get(b, []):
            r = A.nearest(pts, nu)
            if r is not None:
                out.append((b, r))
    return out


def draw_snapshots_wide(fig, ax, c0, chains, cols, snap_scale=None):
    """The shapes of the picked links next to their points (side and shift of Analyze_shape_chains.SNAPSHOTS), sized for a panel `cols` columns wide; the layout must be final."""
    f = cols * figstyle.COLUMN_CM / 2.54 / REF_WIDTH_IN
    scale = snap_scale if snap_scale is not None else SNAP_SCALE * f
    shrink = scale / SNAP_SCALE                                             # the shifts and gaps (inches) follow the size of the shapes
    for b in ORDER:
        pts = A.curve(chains.get(b, []))
        for nu, side, shift in A.SNAPSHOTS.get(round(c0, 6), {}).get(b, []):
            r = A.nearest(pts, nu)
            if r is not None:
                A.add_snapshot(fig, ax, r, (r["nu"], r["F"]), side, shift * shrink, A.COLORS[b], scale, gap_in=0.10 * shrink, outline_lw=figstyle.LW * 0.8, line_lw=figstyle.LW_MIN,
                               line_alpha=1.0)


def _inches(ax, xy):
    """Data points as inches from the lower left corner of the axes."""
    fig = ax.figure
    bb = ax.get_window_extent()
    return np.asarray(ax.transData.transform(np.asarray(xy, float).reshape(-1, 2))) / fig.dpi - np.array([bb.x0, bb.y0]) / fig.dpi


def _obstacles(ax, parts, c0, dot_in, step=0.012):
    """Points along every curve (the dots and the line between them) and the black diamonds, in inches from the lower left corner of the axes: what a picture has to keep clear of."""
    pts = []
    for rows_list in parts.values():
        for rows in rows_list:
            c = A.curve(rows)
            if len(c) < 1:
                continue
            xy = _inches(ax, [(r["nu"], r["F"]) for r in c])
            for p, q in zip(xy[:-1], xy[1:]):
                n = max(int(np.hypot(*(q - p)) / step), 1)
                pts.append(p + (q - p) * (np.arange(n)[:, None] / n))
            pts.append(xy[-1:])
    for _label, nu, F in A.anchors(c0):
        pts.append(_inches(ax, [(nu, F)]))
    return np.vstack(pts)


def _leader(spec, c):
    """The line from the point of a picture to the nearest point of the picture (a box of size w x h centred on c): (px, py, qx, qy)."""
    px, py = spec["anchor"]
    return px, py, min(max(px, c[0] - spec["w"] / 2), c[0] + spec["w"] / 2), min(max(py, c[1] - spec["h"] / 2), c[1] + spec["h"] / 2)


def _crosses(a, b):
    """Do the segments a and b (px, py, qx, qy) cross?"""
    def orient(ax, ay, bx, by, cx, cy):
        return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    return orient(*a[:2], *a[2:], *b[:2]) * orient(*a[:2], *a[2:], *b[2:]) < 0 and orient(*b[:2], *b[2:], *a[:2]) * orient(*b[:2], *b[2:], *a[2:]) < 0


def _layout(specs, samples, W, H, pad, step=0.02, passes=4):
    """Centres (inches from the lower left corner of the axes) for the pictures `specs` (dicts with the size w, h and the anchor point): each one on the free spot that keeps it off the curve samples
    and the other pictures, and its line to the point short and off the curves and the other lines. Greedy, the largest first, then a few passes that move one picture at a time.
    Returns (centres, overlap, reaches, crossings): the overlap is 0 when nothing lies on anything, the reach the length of the line of every picture, crossings the pairs of lines that cross."""
    n = len(specs)
    centres = [None] * n
    sx, sy = samples[:, 0][None, :], samples[:, 1][None, :]

    def best(i, others):
        w, h = specs[i]["w"], specs[i]["h"]
        xs = np.arange(w / 2 + 0.01, W - w / 2 - 0.01 + 1e-9, step)
        ys = np.arange(h / 2 + 0.01, H - h / 2 - 0.01 + 1e-9, step)
        if len(xs) == 0 or len(ys) == 0:
            return None
        cx, cy = (a.ravel() for a in np.meshgrid(xs, ys))
        inside = (sx > (cx - w / 2 - pad)[:, None]) & (sx < (cx + w / 2 + pad)[:, None]) & (sy > (cy - h / 2 - pad)[:, None]) & (sy < (cy + h / 2 + pad)[:, None])
        cost = 20.0 * inside.sum(axis=1)
        for j in others:
            ox, oy = centres[j]
            cost = cost + 1000.0 * ((np.abs(cx - ox) < (w + specs[j]["w"]) / 2 + pad) & (np.abs(cy - oy) < (h + specs[j]["h"]) / 2 + pad))
        px, py = specs[i]["anchor"]
        qx, qy = np.clip(px, cx - w / 2, cx + w / 2), np.clip(py, cy - h / 2, cy + h / 2)
        cost = cost + 30.0 * np.hypot(qx - px, qy - py)
        dx, dy = qx - px, qy - py                                           # the line must not run over a curve (but its own point)
        t = np.clip(((sx - px) * dx[:, None] + (sy - py) * dy[:, None]) / (dx * dx + dy * dy + 1e-12)[:, None], 0.0, 1.0)
        near = np.hypot(sx - (px + t * dx[:, None]), sy - (py + t * dy[:, None])) < 0.012
        near &= (np.hypot(sx - px, sy - py) > 0.08)
        cost = cost + 4.0 * near.sum(axis=1)
        tt = np.linspace(0.1, 0.9, 9)[None, :]
        lx, ly = px + tt * dx[:, None], py + tt * dy[:, None]
        for j in others:                                                    # nor over another picture, nor a picture over the line of another
            ox, oy = centres[j]
            cost = cost + 300.0 * ((np.abs(lx - ox) < specs[j]["w"] / 2 + 0.01) & (np.abs(ly - oy) < specs[j]["h"] / 2 + 0.01)).any(axis=1)
            a, b = _leader(specs[j], centres[j])[:2], _leader(specs[j], centres[j])[2:]
            mx, my = a[0] + tt * (b[0] - a[0]), a[1] + tt * (b[1] - a[1])
            cost = cost + 300.0 * ((np.abs(mx - cx[:, None]) < w / 2 + 0.01) & (np.abs(my - cy[:, None]) < h / 2 + 0.01)).any(axis=1)
        for j in others:                                                    # nor over the line of another picture
            a, b = _leader(specs[j], centres[j])[:2], _leader(specs[j], centres[j])[2:]
            o1 = (qx - px) * (a[1] - py) - (qy - py) * (a[0] - px)
            o2 = (qx - px) * (b[1] - py) - (qy - py) * (b[0] - px)
            o3 = (b[0] - a[0]) * (py - a[1]) - (b[1] - a[1]) * (px - a[0])
            o4 = (b[0] - a[0]) * (qy - a[1]) - (b[1] - a[1]) * (qx - a[0])
            cost = cost + 300.0 * ((o1 * o2 < 0) & (o3 * o4 < 0))
        k = int(np.argmin(cost))
        return float(cx[k]), float(cy[k])

    order = sorted(range(n), key=lambda i: -specs[i]["w"] * specs[i]["h"])
    for i in order:
        centres[i] = best(i, [j for j in order if centres[j] is not None])
    for _ in range(passes):
        for i in order:
            centres[i] = best(i, [j for j in range(n) if j != i and centres[j] is not None])
    over, reaches, lines = 0.0, [], []                                      # what still lies on what in the final layout
    for i in range(n):
        if centres[i] is None:
            over += 1e6
            reaches.append(9.0)
            lines.append((0, 0, 0, 0))
            continue
        (cx, cy), w, h = centres[i], specs[i]["w"], specs[i]["h"]
        over += 20.0 * float(((np.abs(samples[:, 0] - cx) < w / 2 + pad) & (np.abs(samples[:, 1] - cy) < h / 2 + pad)).sum())
        for j in range(i + 1, n):
            if centres[j] is not None and abs(cx - centres[j][0]) < (w + specs[j]["w"]) / 2 + pad and abs(cy - centres[j][1]) < (h + specs[j]["h"]) / 2 + pad:
                over += 1000.0
        ln = _leader(specs[i], centres[i])
        lines.append(ln)
        reaches.append(float(np.hypot(ln[2] - ln[0], ln[3] - ln[1])))
    crossings = sum(1 for i in range(n) for j in range(i + 1, n) if _crosses(lines[i], lines[j]))
    tt = np.linspace(0.1, 0.9, 9)
    for i in range(n):                                                      # a line over a picture that is not its own
        for j in range(n):
            if i != j and centres[i] is not None and centres[j] is not None:
                lx, ly = lines[i][0] + tt * (lines[i][2] - lines[i][0]), lines[i][1] + tt * (lines[i][3] - lines[i][1])
                if ((np.abs(lx - centres[j][0]) < specs[j]["w"] / 2 + 0.01) & (np.abs(ly - centres[j][1]) < specs[j]["h"] / 2 + 0.01)).any():
                    crossings += 1
    return centres, over, reaches, crossings


def draw_snapshots_narrow(fig, ax, c0, parts, chains, snap_scale=None, max_reach=0.40):
    """The shapes of the picked links where there is room in a narrow panel, each with a thin line to its point: the largest size (SNAP_SCALES inches per unit length) at which none lies on a
    curve or on another, no two lines cross and no line is longer than max_reach inches. The layout of the figure must be final. Returns (scale, what is left of the problems)."""
    from matplotlib.lines import Line2D
    sel = picked(c0, chains)
    fw, fh = fig.get_size_inches()
    bb = ax.get_window_extent()
    W, H = bb.width / fig.dpi, bb.height / fig.dpi
    segs = [A.outline(r) for _b, r in sel]
    ext = [(2.0 * np.abs(s.reshape(-1, 2)[:, 0]).max(), s.reshape(-1, 2)[:, 1].min(), s.reshape(-1, 2)[:, 1].max()) for s in segs]
    anchors = [_inches(ax, [(r["nu"], r["F"])])[0] for _b, r in sel]
    samples = _obstacles(ax, parts, c0, NARROW_DOT / 72.0)
    pad = 0.5 * NARROW_DOT / 72.0 + 0.012
    tried = []
    for scale in ([snap_scale] if snap_scale else SNAP_SCALES):
        specs = [dict(w=e[0] * scale, h=(e[2] - e[1]) * scale, anchor=tuple(a)) for e, a in zip(ext, anchors)]
        centres, over, reaches, crossings = _layout(specs, samples, W, H, pad)
        tried.append((over >= 1.0, crossings, max(reaches) > max_reach, -scale, scale, specs, centres, over, crossings, max(reaches)))
        if over < 1.0 and crossings == 0 and max(reaches) <= max_reach:
            break
    else:
        tried.sort(key=lambda t: t[:4])                                     # none is clean: the one with the fewest problems, the largest of those
        scale, specs, centres, over, crossings = tried[0][4:9]
    x0, y0 = bb.x0 / fig.dpi, bb.y0 / fig.dpi
    for (b, r), sp, c, s, e in zip(sel, specs, centres, segs, ext):
        if c is None:
            continue
        cx, cy = (x0 + c[0]) / fw, (y0 + c[1]) / fh                         # the centre in figure fractions
        w, h = sp["w"] / fw, sp["h"] / fh
        box = fig.add_axes([cx - w / 2, cy - h / 2, w, h])
        box.set_xlim(-e[0] / 2.0, e[0] / 2.0)
        box.set_ylim(e[1], e[2])
        box.axis("off")
        box.patch.set_alpha(0.0)
        A.draw_outline(box, s, A.COLORS[b], lw=figstyle.LW * 0.8)
        px, py, qx, qy = _leader(sp, c)
        fig.add_artist(Line2D([(x0 + px) / fw, (x0 + qx) / fw], [(y0 + py) / fh, (y0 + qy) / fh], transform=fig.transFigure, color=A.COLORS[b], lw=figstyle.LW_MIN, alpha=1.0))
    return scale, (over, crossings)


# ---- the figures ------------------------------------------------------------------------------------------------------------------------------------------------
def _panels(fig, axes, c0s, cols, snapshots, snap_scale, **kw):
    """Decorate, fix the layout, draw the curves and the shapes of every panel."""
    for c0, ax in zip(c0s, axes):
        decorate(ax, c0, cols)
    fig.canvas.draw()                                                       # the layout is final: the curves are thinned and the shapes placed on the page
    fig.set_layout_engine("none")
    drawn = [(c0, ax) + draw_branches(ax, c0, cols, **kw) for c0, ax in zip(c0s, axes)]
    if (cols >= 1.0) if snapshots is None else snapshots:
        for c0, ax, parts, chains in drawn:
            if cols < 1.0:
                scale, (over, crossings) = draw_snapshots_narrow(fig, ax, c0, parts, chains, snap_scale)
                print("c0 = %g: shapes at %g in per unit length%s" % (c0, scale, "" if over < 1.0 and crossings == 0 else " (still %d lines cross, overlap %g)" % (crossings, over)))
            else:
                draw_snapshots_wide(fig, ax, c0, chains, cols, snap_scale)


def fig_branches(c0, cols=COLS, height_cm=None, aspect=None, snapshots=None, snap_scale=None, **kw):
    """One figure, `cols` columns wide: F_b / 8 pi kappa against v for the branches prolate-1, prolate-2 and the pear-shaped one at c0 (2.4 or 3). snapshots None: always. The other keywords go to
    draw_branches (root, saddles, paper, unstable, ends, labels, legend)."""
    figstyle.use()
    aspect = aspect if aspect is not None else (NARROW_ASPECT if cols < 1.0 else 0.70)
    fig, ax = figstyle.figure(cols=cols, aspect=aspect, height_cm=height_cm)
    _panels(fig, [ax], [c0], cols, True if snapshots is None else snapshots, snap_scale, **kw)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("which", nargs="*", help="the figures to make: 2.4 and / or 3.0 (default: both)")
    ap.add_argument("--cols", type=float, default=COLS, help="width in columns (default 1/3)")
    ap.add_argument("--height_cm", type=float, default=None, help="height in cm (default: %g times the width of one figure below one column, 0.7 times from one column on)" % NARROW_ASPECT)
    ap.add_argument("--aspect", type=float, default=None, help="height / width of one figure instead of --height_cm")
    ap.add_argument("--paper", action="store_true", help="the grey diamonds: the points of the paper read off its plot")
    ap.add_argument("--unstable", action="store_true", help="the dashed unstable symmetric part of prolate-1 (Shape_axisym.py saddle)")
    ap.add_argument("--ends", action="store_true", help="a filled square where a branch ends (the state it falls to)")
    ap.add_argument("--legend", dest="legend", action="store_true", default=None, help="the legend even when it does not fit (default: from one column on)")
    ap.add_argument("--no_legend", dest="legend", action="store_false", help="no legend")
    ap.add_argument("--no_snapshots", dest="snapshots", action="store_false", default=None, help="without the shapes on the curves")
    ap.add_argument("--snap_scale", type=float, default=None, help="inches per unit length of the shapes (default: wide figures %g times the width over %g in, narrow ones the largest that fits)" % (SNAP_SCALE, REF_WIDTH_IN))
    ap.add_argument("--root", type=str, default=None, help="folder of the chains (one c0 only; default Results/Shape_chains/c0_<c0>)")
    ap.add_argument("--saddles", type=str, default=None, help="csv of the dashed segments (one c0 only; default Studies/shape_chains/fig11_data|fig13_data/saddles.csv)")
    ap.add_argument("--out_dir", type=str, default=None, help="folder for the files (default Results/Paper_figures)")
    args = ap.parse_args()
    which = args.which or ["2.4", "3.0"]
    bad = [w for w in which if not w.replace(".", "", 1).isdigit() or float(w) not in FIGS]
    if bad:
        raise SystemExit("%s: the figures are 2.4 and 3.0" % ", ".join(bad))
    import matplotlib
    matplotlib.use("Agg")
    common = dict(height_cm=args.height_cm, aspect=args.aspect, snapshots=args.snapshots, snap_scale=args.snap_scale, paper=args.paper, unstable=args.unstable, ends=args.ends, legend=args.legend)
    for w in which:
        fig = fig_branches(float(w), cols=args.cols, root=args.root, saddles=args.saddles, **common)
        figstyle.save_figure(fig, FIGS[float(w)]["name"], folder=args.out_dir)


if __name__ == "__main__":
    main()
