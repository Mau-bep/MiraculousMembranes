#!/usr/bin/env python3
"""The bending energy of each BRANCH of vesicle shapes against the reduced volume at a fixed c0 (Seifert, Berndl and Lipowsky, Phys. Rev. A 44, 1182 (1991), Figs. 11 and 13), from the
chains of Run_shape_chain.py, with the outline of the shapes on the curves.

usage:  micromamba run -n mir_mem python3 Analyze_shape_chains.py --c0 2.4 --root ../Results/Shape_chains/c0_2.4 --out ../Studies/shape_chains/fig11
        python3 Analyze_shape_chains.py --c0 3 prolate1=../Results/A/prolate1 pear=../Results/A/pear --out fig13      (explicit folders, name=folder)

--root holds one folder per branch (prolate1, prolate2, pear), each with the chain_summary.csv of its chain. Writes <out>.png (and .svg), <out>_shapes.png (the outline of every link of
every branch in order of nu: the check that a chain stayed on its branch), <out>_table.csv (the points of the curves) and <out>_landmarks.md (the maxima, minima, D^pear and the ends of
the branches next to the paper's values).

The energy is F_b / (8 pi kappa) = Bending / (4 pi K_B) at the reduced volume nu that was REACHED (the constraints are penalties, so it is a little off the target; the energy of a penalty
minimum is the constrained one at the nu and c0 reached). Marked on the curves: filled dots are converged links, open dots links that reached the step cap, crosses the state a branch fell to
where it ended (kind "jump": a state of another branch, not a point of the curve). Dashed lines (--saddles, from Shape_axisym.py saddle) are the up-down symmetric part of prolate-1 below the last link of
its chain: there the symmetric shapes are saddles (unstable against becoming a pear), a minimiser leaves them, and the 1D solver follows them with the symmetry imposed, down to the minimal volume M^pro.
Black diamonds are exact values of the theory that need no figure: the sphere F = (1 - c0/2)^2 at nu = 1
and, for c0 > 2 sqrt 2, the budding limit L^pear (two spheres of different radii with 1/R1 + 1/R2 = C0, F = 2 - c0 s + c0^2/4 with s = R1 + R2). Grey open diamonds (--paper) are the points of the paper's figure read off the plot (about +-0.003 in F, +-0.005 in nu).
"""
import argparse
import csv
import math
import os
import sys

import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Run_shape_chain import meridian_outline, read_obj  # noqa: E402

BRANCHES = ["prolate1", "prolate2", "pear"]
LABELS = {"prolate1": "prolate-1", "prolate2": "prolate-2", "pear": "pear-shaped"}
SHORT = {"prolate1": "P1", "prolate2": "P2", "pear": "pear"}
COLORS = {"prolate1": "#1f77b4", "prolate2": "#2ca02c", "pear": "#d62728"}
NUMERIC = ("nu_target", "nu", "F", "E_total", "c0_eff", "H0", "steps", "seconds", "nV", "length", "rmax", "r_big", "r_small", "rneck", "neck_ratio", "lobe_ratio", "skew", "nlobes")

# axes of the paper's figures and the points read off them
VIEW = {2.4: dict(xlim=(0.6, 1.0), ylim=(0.0, 0.2)), 3.0: dict(xlim=(0.5, 1.0), ylim=(0.0, 0.5))}
PAPER = {
    2.4: [("prolate-1 maximum", 0.822, 0.124), ("D^pear", 0.749, 0.108), ("C1^pear", 0.718, 0.086), ("M^pro", 0.705, 0.068), ("prolate-2 at 0.6", 0.600, 0.125),
          ("prolate-2 minimum", 0.700, 0.037), ("C^pear", 0.705, 0.041), ("M^pear", 0.800, 0.200)],
    3.0: [("prolate-1 maximum", 0.890, 0.287), ("D^pear", 0.823, 0.270), ("C1^pear", 0.750, 0.210), ("prolate-1 end", 0.707, 0.100), ("prolate-2 at 0.5", 0.500, 0.090),
          ("prolate-2 minimum", 0.565, 0.033), ("prolate-2 maximum", 0.600, 0.070), ("M^pear", 0.890, 0.450)],
}
# snapshots: per c0 and branch, (nu of the link, side of the point: above, below, left or right, shift in inches along the side)
SNAPSHOTS = {
    2.4: {"prolate1": [(0.97, "above", 0.12), (0.90, "below", 0.0), (0.82, "below", 0.0), (0.745, "above", -0.30)],
          "prolate2": [(0.62, "above", 0.0), (0.66, "above", 0.0), (0.695, "below", -0.15)],
          "pear": [(0.716, "right", -0.10), (0.746, "right", -0.05), (0.786, "right", 0.05)]},
    3.0: {"prolate1": [(0.97, "above", 0.12), (0.90, "above", 0.0), (0.82, "above", 0.0), (0.76, "above", -0.15)],
          "prolate2": [(0.52, "above", 0.0), (0.66, "above", 0.0), (0.70, "above", 0.15)],
          "pear": [(0.775, "left", 0.0), (0.833, "right", -0.45), (0.872, "right", 0.0), (0.8887, "right", 0.0)]},
}


def anchors(c0):
    """Exact points of the theory (label, nu, F): the sphere and, for c0 > 2 sqrt 2, the budding limit L^pear (two spheres of different radii with 1/R1 + 1/R2 = C0 joined by an ideal neck)."""
    pts = [("sphere", 1.0, (1.0 - c0 / 2.0) ** 2)]
    if c0 > 2.0 * math.sqrt(2.0):
        p = (1.0 + math.sqrt(1.0 + c0 * c0)) / (c0 * c0)
        s = c0 * p
        pts.append(("L^pear", s ** 3 - 3.0 * p * s, 2.0 - c0 * s + c0 * c0 / 4.0))
    return pts


def parabola_extremum(nu, F, kind):
    """(nu, F) of the maximum or minimum of the points: the best point refined by the parabola through it and its two neighbours."""
    i = int(np.argmax(F) if kind == "max" else np.argmin(F))
    if 0 < i < len(nu) - 1:
        a, b, c = np.polyfit(nu[i - 1:i + 2], F[i - 1:i + 2], 2)
        if a != 0.0:
            x = -b / (2.0 * a)
            if nu[i - 1] <= x <= nu[i + 1]:
                return float(x), float(np.polyval([a, b, c], x))
    return float(nu[i]), float(F[i])


def crossing(nu1, F1, nu2, F2):
    """nu of the first crossing of two curves F1(nu), F2(nu) (linear interpolation) in their common range, or None."""
    lo, hi = max(nu1.min(), nu2.min()), min(nu1.max(), nu2.max())
    if hi <= lo:
        return None
    grid = np.linspace(lo, hi, 2000)
    d = np.interp(grid, nu1, F1) - np.interp(grid, nu2, F2)
    idx = np.where(np.sign(d[:-1]) * np.sign(d[1:]) < 0)[0]
    return float(grid[idx[0]]) if len(idx) else None


def landmarks(chains, c0):
    """Our values of the points the paper marks on its figures: (label, nu, F)."""
    out = []
    cur = {b: curve(r) for b, r in chains.items()}
    arr = {b: (np.array([r["nu"] for r in pts]), np.array([r["F"] for r in pts])) for b, pts in cur.items() if len(pts) >= 3}
    if "prolate1" in arr:
        nu, F = arr["prolate1"]
        out.append(("prolate-1 maximum",) + parabola_extremum(nu, F, "max"))
        out.append(("prolate-1 lowest v reached", float(nu[0]), float(F[0])))
    if "prolate2" in arr:
        nu, F = arr["prolate2"]
        # the first local minimum (at c0 = 3 the branch has a second, deeper end towards the neck closure)
        inner = [i for i in range(1, len(F) - 1) if F[i] < F[i - 1] and F[i] <= F[i + 1]]
        if inner:
            i = inner[0]
            lo, hi = max(i - 1, 0), min(i + 2, len(F))
            out.append(("prolate-2 minimum",) + parabola_extremum(nu[lo:hi], F[lo:hi], "min"))
        else:
            out.append(("prolate-2 minimum",) + parabola_extremum(nu, F, "min"))
        out.append(("prolate-2 lowest v", float(nu[0]), float(F[0])))
        out.append(("prolate-2 highest v reached", float(nu[-1]), float(F[-1])))
        if c0 > 2.8:
            # c0 = 3: the energy of the dumbbell part of the branch (two lobes, one neck) has a maximum near v = 0.6, before it falls to the fold at v = 0.58; below it the branch is the three-lobed one
            sel = (nu > 0.585) & (nu < 0.66)
            if sel.sum() >= 3:
                out.append(("prolate-2 maximum",) + parabola_extremum(nu[sel], F[sel], "max"))
    if "pear" in arr:
        nu, F = arr["pear"]
        out.append(("pear, lowest v", float(nu[0]), float(F[0])))
        out.append(("pear, highest v (M^pear)", float(nu[-1]), float(F[-1])))
        if "prolate1" in arr:
            x = crossing(nu, F, *arr["prolate1"])
            if x is not None:
                out.append(("D^pear (pear = prolate-1)", x, float(np.interp(x, *arr["prolate1"]))))
    return out


def load_chain(folder):
    path = os.path.join(folder, "chain_summary.csv")
    if not os.path.exists(path):
        return []
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            row = {k: (float(v) if (k in NUMERIC and v != "") else v) for k, v in r.items()}
            if row["mesh"] and not os.path.isabs(row["mesh"]):
                row["mesh"] = os.path.join(folder, row["mesh"])
            rows.append(row)
    return rows


def load_parts(root=None, explicit=()):
    """{branch: [rows of every chain folder of the branch]}: the folders of --root (prolate1, prolate2, pear; a branch may be walked by several chains over different ranges of v: folders prolate2,
    prolate2b, ...: their points make one curve each, drawn in the same colour) and the `branch=folder,folder` pairs of `explicit`."""
    folders = {}
    if root:
        for b in BRANCHES:
            found = [os.path.join(root, b + suffix) for suffix in ("", "b", "c") if os.path.isdir(os.path.join(root, b + suffix))]
            if found:
                folders[b] = found
    for item in explicit:
        name, _, path = item.partition("=")
        folders[name] = path.split(",")
    parts = {b: [r for r in (load_chain(p) for p in paths) if r] for b, paths in folders.items()}
    return {b: pr for b, pr in parts.items() if pr}


def curve(rows):
    """The accepted links of a branch (the seed and the points) sorted by the reduced volume that was reached. A link whose run 'ended' (the integrator quit before converging or reaching the step cap:
    the neck was closing and the shape stuck at a v that is not the target) is left out, also in old chains that accepted it."""
    pts = [r for r in rows if r["kind"] in ("seed", "point") and r["F"] != "" and r["status"] != "ended"]
    return sorted(pts, key=lambda r: r["nu"])


def load_saddles(path):
    """The up-down symmetric (unstable) part of prolate-1 from Shape_axisym.py saddle, sorted by decreasing nu: nu, F, the lowest eigenvalue against asymmetric deformations."""
    with open(path) as f:
        rows = [{"nu": float(r["nu"]), "F": float(r["F"]), "lam_asym": float(r["lowest_eigenvalue_asymmetric"])} for r in csv.DictReader(f)]
    return sorted(rows, key=lambda r: -r["nu"])


def saddle_landmarks(rows):
    """C1^pear (where the lowest eigenvalue against asymmetric deformations changes sign, linear interpolation: below it the symmetric shape is a saddle) and the end of the branch (M^pro)."""
    out = []
    for a, b in zip(rows[:-1], rows[1:]):
        if a["lam_asym"] > 0.0 >= b["lam_asym"]:
            t = a["lam_asym"] / (a["lam_asym"] - b["lam_asym"])
            out.append(("C1^pear (1D solver: the symmetric prolate-1 becomes a saddle)", a["nu"] + t * (b["nu"] - a["nu"]), a["F"] + t * (b["F"] - a["F"])))
            break
    out.append(("prolate-1 minimal volume (1D solver, unstable)", rows[-1]["nu"], rows[-1]["F"]))
    return out


def outline(row):
    V, F = read_obj(row["mesh"])
    return meridian_outline(V, F)


def draw_outline(ax, segs, color, lw=0.8):
    ax.add_collection(LineCollection(segs, colors=color, linewidths=lw, capstyle="round"))


def add_snapshot(fig, ax, row, xy, side, shift_in, color, scale_in, gap_in=0.10, outline_lw=0.8, line_lw=0.5, line_alpha=0.8):
    """An outline of the shape of `row` next to the data point xy, on the given side of the curve, with a thin line to the point. scale_in: inches per unit of length of the shape.
    line_alpha 1 for a format without transparency (eps)."""
    segs = outline(row)
    if len(segs) == 0:
        return None
    pts = segs.reshape(-1, 2)
    w_ext = 2.0 * np.abs(pts[:, 0]).max()
    s_lo, s_hi = pts[:, 1].min(), pts[:, 1].max()
    w_in, h_in = w_ext * scale_in, (s_hi - s_lo) * scale_in
    fw, fh = fig.get_size_inches()
    px, py = fig.transFigure.inverted().transform(ax.transData.transform(xy))
    if side in ("above", "below"):
        sign = 1.0 if side == "above" else -1.0
        cx, cy = px + shift_in / fw, py + sign * (gap_in + h_in / 2.0) / fh           # shift_in moves it along the curve (horizontally)
    else:
        sign = 1.0 if side == "right" else -1.0
        cx, cy = px + sign * (gap_in + w_in / 2.0) / fw, py + shift_in / fh           # shift_in moves it vertically
    box = fig.add_axes([cx - w_in / 2.0 / fw, cy - h_in / 2.0 / fh, w_in / fw, h_in / fh])
    box.set_xlim(-w_ext / 2.0, w_ext / 2.0)
    box.set_ylim(s_lo, s_hi)
    box.axis("off")
    box.patch.set_alpha(0.0)
    draw_outline(box, segs, color, lw=outline_lw)
    # the line from the point to the nearest point of the outline's box
    qx = min(max(px, cx - w_in / 2.0 / fw), cx + w_in / 2.0 / fw)
    qy = min(max(py, cy - h_in / 2.0 / fh), cy + h_in / 2.0 / fh)
    fig.add_artist(Line2D([px, qx], [py, qy], transform=fig.transFigure, color=color, lw=line_lw, alpha=line_alpha))
    return {"nu": xy[0], "F": xy[1], "side": side, "shift_in": shift_in, "scale_in": scale_in, "width_in": w_in, "height_in": h_in}


def nearest(curve_rows, nu, tol=0.012):
    best = min(curve_rows, key=lambda r: abs(r["nu"] - nu)) if curve_rows else None
    return best if best is not None and abs(best["nu"] - nu) <= tol else None


def export_data(args, parts, drawn):
    """The numbers behind the figure, for whoever redraws it (with renders instead of outlines, say): every point of every curve with its mesh, and the shapes drawn on the curves."""
    import shutil
    from Run_shape_chain import waist_radius
    os.makedirs(os.path.join(args.data_dir, "snapshots"), exist_ok=True)
    name = lambda b, r: "%s_nu%.4f_F%.5f" % (b, r["nu"], r["F"])
    cols = ["branch", "chain", "leg", "tag", "kind", "nu_target", "nu", "F", "c0_eff", "H0", "status", "steps", "nV", "length", "rmax", "neck_radius", "neck_over_rmax", "neck_ratio_60slabs", "lobe_ratio",
            "skew", "mesh_original", "mesh_exported"]
    with open(os.path.join(args.data_dir, "points.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for b, prs in parts.items():
            for k, rows in enumerate(prs):
                for r in curve(rows):
                    V, Fc = read_obj(r["mesh"])
                    neck, rmax, ratio = waist_radius(V, Fc)
                    dest = ""
                    if args.all_meshes:
                        os.makedirs(os.path.join(args.data_dir, "meshes", b), exist_ok=True)
                        dest = os.path.join("meshes", b, name(b, r) + ".obj")
                        shutil.copy(r["mesh"], os.path.join(args.data_dir, dest))
                    w.writerow([b, k, r["leg"], r["tag"], r["kind"], r["nu_target"], "%.6g" % r["nu"], "%.6g" % r["F"], "%.6g" % r["c0_eff"], "%.6g" % r["H0"], r["status"], int(r["steps"]), int(r["nV"]),
                                "%.5g" % r["length"], "%.5g" % rmax, "%.5g" % neck, "%.5g" % ratio, "%.4g" % r["neck_ratio"], "%.4g" % r["lobe_ratio"], "%.5g" % r["skew"], r["mesh"], dest])
    with open(os.path.join(args.data_dir, "snapshots.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["branch", "nu", "F", "c0_eff", "tag", "mesh_file", "anchor_nu", "anchor_F", "side", "shift_in", "inches_per_length", "width_in", "height_in"])
        for b, r, info in drawn:
            dest = os.path.join("snapshots", name(b, r) + ".obj")
            shutil.copy(r["mesh"], os.path.join(args.data_dir, dest))
            w.writerow([b, "%.6g" % r["nu"], "%.6g" % r["F"], "%.6g" % r["c0_eff"], r["tag"], dest, "%.6g" % info["nu"], "%.6g" % info["F"], info["side"], info["shift_in"], info["scale_in"],
                        "%.4g" % info["width_in"], "%.4g" % info["height_in"]])
    if args.saddles and os.path.abspath(args.saddles) != os.path.abspath(os.path.join(args.data_dir, "saddles.csv")):
        shutil.copy(args.saddles, os.path.join(args.data_dir, "saddles.csv"))
    print("wrote %s/points.csv, snapshots.csv%s and %d snapshot meshes%s" % (args.data_dir, ", saddles.csv" if args.saddles else "", len(drawn), " and every mesh" if args.all_meshes else ""))


def continuum_figure(args, parts):
    """The continuum energy at the nu and c0 of every point (Shape_axisym.py, started from the shape of the mesh so it stays on the branch): the curves of the 1D problem with the mesh points, and the
    relative difference of the mesh energy to it. The mesh energy is that of the 3D discretisation at the remesher settings of the chain; the difference is its error."""
    import matplotlib.pyplot as plt
    from Shape_axisym import from_mesh
    path = args.out + "_continuum.csv"
    done = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            done[(r["branch"], r["chain"], r["tag"])] = r
    rows_out = []
    for b, prs in parts.items():
        for k, rows in enumerate(prs):
            for r in curve(rows):
                key = (b, str(k), r["tag"])
                if key in done:
                    rows_out.append(done[key])
                    continue
                try:
                    Fc, psi, L, ok, res, nu = from_mesh(r["mesh"], r["c0_eff"], r["nu"])
                except Exception as exc:                                  # a shape the 1D contour cannot hold (an overhang), or a solver failure
                    print("continuum: %s %s failed: %r" % (b, r["tag"], exc))
                    continue
                rec = {"branch": b, "chain": str(k), "tag": r["tag"], "nu": "%.6g" % r["nu"], "c0_eff": "%.6g" % r["c0_eff"], "F_mesh": "%.6g" % r["F"], "F_continuum": "%.6g" % Fc,
                       "diff_percent": "%.4g" % ((r["F"] - Fc) / Fc * 100.0), "converged": str(bool(ok and res < 1e-8))}
                rows_out.append(rec)
                print("continuum %s %-20s nu %.4f  F mesh %.5f  continuum %.5f  %+.2f %%  %s" % (b, r["tag"], r["nu"], r["F"], Fc, float(rec["diff_percent"]), "" if rec["converged"] == "True" else "NOT CONVERGED"), flush=True)
    keys = ["branch", "chain", "tag", "nu", "c0_eff", "F_mesh", "F_continuum", "diff_percent", "converged"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows_out)
    good = [r for r in rows_out if r["converged"] == "True"]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    for b in parts:
        pts = sorted((r for r in good if r["branch"] == b), key=lambda r: float(r["nu"]))
        if not pts:
            continue
        col = COLORS.get(b, "k")
        nu = [float(r["nu"]) for r in pts]
        axs[0].plot(nu, [float(r["F_continuum"]) for r in pts], "-", color=col, lw=1.0, label=LABELS.get(b, b) + ": continuum")
        axs[0].plot(nu, [float(r["F_mesh"]) for r in pts], "o", color=col, ms=2.5)
        axs[1].plot(nu, [float(r["diff_percent"]) for r in pts], "o-", color=col, ms=2.5, lw=0.8, label=LABELS.get(b, b))
    axs[0].set_xlabel(r"$v$")
    axs[0].set_ylabel(r"$F_b\,/\,8\pi\kappa$")
    axs[0].legend(fontsize=7)
    axs[0].set_title(r"$c_0=%g$: continuum (line) and 3D mesh (dots)" % args.c0, fontsize=9)
    axs[1].axhline(0, color="0.7", lw=0.6)
    axs[1].set_xlabel(r"$v$")
    axs[1].set_ylabel("(F mesh - F continuum) / F continuum  [%]")
    axs[1].legend(fontsize=7)
    axs[1].set_title("error of the 3D energy at the chain's remesher settings", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out + "_continuum.png", dpi=170)
    plt.close(fig)
    print("wrote %s_continuum.png, %s_continuum.csv" % (args.out, args.out))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--c0", type=float, required=True)
    ap.add_argument("--root", type=str, default=None, help="folder with one subfolder per branch (prolate1, prolate2, pear)")
    ap.add_argument("--out", type=str, required=True, help="output prefix (<out>.png, <out>.svg, <out>_shapes.png, <out>_table.csv)")
    ap.add_argument("--paper", action="store_true", help="mark the points of the paper's figure (read off the plot)")
    ap.add_argument("--no_snapshots", action="store_true")
    ap.add_argument("--data_dir", type=str, default=None, help="write the numbers behind the figure here: points.csv (every point of every curve: branch, chain, nu, F, c0_eff, shape numbers, neck radius, mesh), "
                                                                "snapshots.csv (the shapes drawn on the curves with their place in the plot) and snapshots/*.obj, so the snapshots can be replaced by renders")
    ap.add_argument("--all_meshes", action="store_true", help="with --data_dir also copy the mesh of every point to <data_dir>/meshes/<branch>/")
    ap.add_argument("--saddles", type=str, default=None, help="csv of Shape_axisym.py saddle: the up-down symmetric part of prolate-1 below its chain, drawn dashed (also exported and in the landmarks)")
    ap.add_argument("--continuum", action="store_true", help="also solve the 1D axisymmetric problem (Shape_axisym.py) from the shape of every point: <out>_continuum.png and .csv (needs scipy; a few seconds a point)")
    ap.add_argument("--scale", type=float, default=0.30, help="inches per unit length of the snapshots (a sphere of this problem is 1.1 long)")
    ap.add_argument("--xlim", type=float, nargs=2, default=None)
    ap.add_argument("--ylim", type=float, nargs=2, default=None)
    ap.add_argument("--title", type=str, default=None)
    ap.add_argument("folders", nargs="*", help="branch=folder pairs instead of --root")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parts = load_parts(args.root, args.folders)
    chains = {b: [row for r in pr for row in r] for b, pr in parts.items()}
    if not chains:
        raise SystemExit("no chain_summary.csv found")
    saddles = load_saddles(args.saddles) if args.saddles else []
    view = VIEW.get(round(args.c0, 6), {})
    xlim = tuple(args.xlim) if args.xlim else view.get("xlim")
    ylim = tuple(args.ylim) if args.ylim else view.get("ylim")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.12, top=0.93)
    table = []
    for b, prs in parts.items():
        col = COLORS.get(b, "k")
        for k, rows in enumerate(prs):
            pts = curve(rows)
            if not pts:
                continue
            nu = np.array([r["nu"] for r in pts])
            F = np.array([r["F"] for r in pts])
            ax.plot(nu, F, "-", color=col, lw=1.2, zorder=2, label=LABELS.get(b, b) if k == 0 else None)
            capped = np.array([r["status"] != "converged" for r in pts])
            ax.plot(nu[~capped], F[~capped], "o", color=col, ms=3.2, zorder=3)
            ax.plot(nu[capped], F[capped], "o", mfc="white", mec=col, ms=3.6, mew=0.9, zorder=3)
            for r in rows:
                if r["kind"] == "jump" and r["F"] != "":
                    ax.plot([r["nu"]], [r["F"]], "x", color=col, ms=5, mew=1.0, zorder=3)
            for r in pts:
                table.append((b, r["leg"], r["nu_target"], r["nu"], r["F"], r["c0_eff"], r["status"], r["neck_ratio"], r["lobe_ratio"], r["skew"], r["steps"], r["tag"]))
    if any(r["kind"] == "jump" and r["F"] != "" for rows in chains.values() for r in rows):
        ax.plot([], [], "x", color="0.3", ms=5, mew=1.0, label="where a branch ends (the state it falls to)")
    if saddles and "prolate1" in chains:
        nu_last = curve(chains["prolate1"])[0]["nu"]                 # the dashed line goes on from the lowest link of the 3D chain
        low = [r for r in saddles if r["nu"] <= nu_last + 1e-6]
        ax.plot([r["nu"] for r in low], [r["F"] for r in low], "--", color=COLORS["prolate1"], lw=1.2, zorder=2)
        ax.plot([], [], "--", color="0.35", lw=1.2, label="unstable symmetric part (1D solver)")
    for label, nu, F in anchors(args.c0):
        ax.plot([nu], [F], "D", mfc="black", mec="black", ms=4, zorder=4)
        right = xlim is not None and nu > xlim[0] + 0.85 * (xlim[1] - xlim[0])
        ax.annotate(label, (nu, F), textcoords="offset points", xytext=(-5, -10) if right else (5, -10), ha="right" if right else "left", fontsize=6.5)
    if args.paper:
        for label, nu, F in PAPER.get(round(args.c0, 6), []):
            ax.plot([nu], [F], "D", mfc="none", mec="0.45", ms=4.5, mew=0.8, zorder=1)
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel(r"reduced volume $v$")
    ax.set_ylabel(r"$F_b\,/\,8\pi\kappa$")
    ax.set_title(args.title or (r"$c_0 = %g$" % args.c0))
    ax.legend(fontsize=7.5, loc="upper left", frameon=False)
    fig.canvas.draw()
    drawn = []
    if not args.no_snapshots:
        for b, rows in chains.items():
            pts = curve(rows)
            for nu, side, shift in SNAPSHOTS.get(round(args.c0, 6), {}).get(b, []):
                r = nearest(pts, nu)
                if r is not None:
                    info = add_snapshot(fig, ax, r, (r["nu"], r["F"]), side, shift, COLORS.get(b, "k"), args.scale)
                    if info is not None:
                        drawn.append((b, r, info))
    fig.savefig(args.out + ".png", dpi=200)
    fig.savefig(args.out + ".svg")
    plt.close(fig)

    if args.data_dir:
        export_data(args, parts, drawn)

    with open(args.out + "_table.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["branch", "leg", "nu_target", "nu", "F", "c0_eff", "status", "neck_ratio", "lobe_ratio", "skew", "steps", "tag"])
        w.writerows(table)

    if args.continuum:
        continuum_figure(args, parts)

    ours = landmarks(chains, args.c0)
    if saddles:
        ours += saddle_landmarks(saddles)
    paper = PAPER.get(round(args.c0, 6), [])
    with open(args.out + "_landmarks.md", "w") as f:
        f.write("c0 = %g: points of the curves (ours) next to the paper's figure read off the plot (+-0.003 in F, +-0.005 in v) and the exact points of the theory\n\n" % args.c0)
        f.write("| point | v (ours) | F (ours) | v (paper) | F (paper) |\n|---|---|---|---|---|\n")
        match = {"prolate-1 maximum": "prolate-1 maximum", "prolate-2 minimum": "prolate-2 minimum", "prolate-2 maximum": "prolate-2 maximum", "D^pear (pear = prolate-1)": "D^pear",
                 "pear, highest v (M^pear)": "M^pear", "C1^pear (1D solver: the symmetric prolate-1 becomes a saddle)": "C1^pear",
                 "prolate-1 minimal volume (1D solver, unstable)": "M^pro" if any(n == "M^pro" for (n, a, b) in paper) else "prolate-1 end"}
        for label, nu, F in ours:
            pl = next(((a, b) for (n, a, b) in paper if n == match.get(label)), ("", ""))
            f.write("| %s | %.4f | %.4f | %s | %s |\n" % (label, nu, F, pl[0] if pl[0] == "" else "%.3f" % pl[0], pl[1] if pl[1] == "" else "%.3f" % pl[1]))
        for label, nu, F in anchors(args.c0):
            branch = {"sphere": "prolate1", "L^pear": "pear"}.get(label)          # the branch that ends in that exact point
            pts = curve(chains[branch]) if branch in chains else []
            f.write("| exact: %s | %.4f | %.4f | | |\n" % (label, nu, F))
            if pts:
                near = min(pts, key=lambda r: abs(r["nu"] - nu))
                if abs(near["nu"] - nu) < 0.01:
                    f.write("| ... nearest point of the %s branch (%s) | %.4f | %.4f | | |\n" % (LABELS.get(branch, branch), near["tag"], near["nu"], near["F"]))
    print(open(args.out + "_landmarks.md").read())

    # every shape of every branch in order of nu
    n = max(len(curve(r)) for r in chains.values())
    cols = min(max(n, 1), 10)
    rows_fig = sum(int(math.ceil(len(curve(r)) / cols)) for r in chains.values())
    fig2, axes = plt.subplots(rows_fig, cols, figsize=(1.25 * cols, 1.7 * rows_fig), squeeze=False)
    for a in axes.ravel():
        a.axis("off")
    i = 0
    for b, rows in chains.items():
        pts = curve(rows)
        for k, r in enumerate(pts):
            a = axes[i + k // cols][k % cols]
            segs = outline(r)
            if len(segs) == 0:
                continue
            p = segs.reshape(-1, 2)
            a.set_xlim(-1.1, 1.1)
            a.set_ylim(-1.1, 1.1)
            a.set_aspect("equal")
            draw_outline(a, segs / max(np.abs(p).max(), 1e-9), COLORS.get(b, "k"), lw=0.7)
            a.set_title("%s v=%.3f\nF=%.4f" % (SHORT.get(b, b), r["nu"], r["F"]), fontsize=5.5, pad=1)
        i += int(math.ceil(len(pts) / cols))
    fig2.tight_layout(pad=0.3)
    fig2.savefig(args.out + "_shapes.png", dpi=170)
    plt.close(fig2)
    print("wrote %s.png, %s.svg, %s_shapes.png, %s_table.csv" % ((args.out,) * 4))


if __name__ == "__main__":
    main()
