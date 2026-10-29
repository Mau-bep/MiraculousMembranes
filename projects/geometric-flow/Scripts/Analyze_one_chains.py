#!/usr/bin/env python3
"""Per-chain analysis of a batch of one-bead chains (Create_subjob_one_vesicle.py / Run_one_vesicle_chain.py), e.g. the Phase_one_vesicle and Phase_one_flat batches of
Mem_one_vesicle_up.sh. The sister of Analyze_two_chains.py.

usage:  python3 Analyze_one_chains.py ../Results/Phase_one_vesicle                  # every chain folder (KI_0, KI_0.3333, ...) of the batch
        python3 Analyze_one_chains.py ../Results/Phase_one_vesicle --chains KI_1 KI_2 --plots coverage volume
        python3 Analyze_one_chains.py ../Results/Phase_one_flat --axes fig2         # the flat disks: w~ against sigma~

Needs numpy + matplotlib (micromamba env mir_mem: ~/micromamba/envs/mir_mem/bin/python); runs where the run folders are (the cluster) and writes a folder small enough to
download. A chain is a folder <batch>/<Tag>/ with chain_summary.csv and the link folders. The points of a chain are its links of the list, both legs (the way along the list,
"in", and the way back, "out"); not the waypoints (_w1) nor the half step of a retry (_mid), not the cold start (--keep_cold keeps it). Below --out (default <batch>/one_bead_analysis):

    phase_space/phase_space_<axes>.png     the whole phase space, every chain of the batch: where a chain has the way in AND the way back for the same value of the scanned
                                           parameter, the point with the LOWER total energy is the one drawn (the same pick as coverage_lowest_energy of Plots.py). fig5 (closed
                                           vesicles): w/sigma against a/lambda, log axis, as deserno_fig5. fig2 (flat disks): w~ against sigma~ as planar_phase_deserno.
                                           No theory curves; the axes fit the points. Coloured by z = 2 x coverage; circles = the way in had the lower energy, squares = the way back
    phase_space/phase_space_map.png        the same points binned in the fashion of contour() of Plots.py (24 x 16 bins, the mean coverage of a bin, no theory)
    phase_space/phase_space_lowest.csv     the points that were drawn
    coverage_vs_adhesion/<Tag>.png         coverage against the adhesion strength w~ = 4 W a^2/KB = 2 w a^2/kappa for one chain, both legs (circles in, squares back;
                                           a black ring marks the point of lower energy); coverage_vs_adhesion/all_chains.png has the lower-energy points of every chain
    volume_difference/<Tag>.png            V - V_target of the final states (V the volume of the mesh, V_target that of the Volume_constraint of the config) against the scanned
                                           parameter, and in percent of V_target; all_chains.png has every chain. Closed vesicles only (a flat disk has no volume)
    final_states_per_chain/<Tag>/<Tag>_Final_link_NNN.obj   the final mesh of every point of the chain, in the order of the chain, with Bead_0_data.txt and states.txt
                                           (frame -> link, leg, scanned value, run folder)
    viewer_runs/<Tag>/                     a run folder for main_visualize (Input_file.json, membrane_<frame>.obj as hard links, Bead_0_data.txt):
                                           main_visualize <out>/viewer_runs   (the slider goes through the final states; the folders of the menu are the chains)
    csv/<Tag>.csv, csv/all_points.csv      the numbers of every point (all_points.csv has a column "lowest": 1 for the points that were drawn in the phase space)

"Coverage" (--cov): union (default) = the fraction of directions from the bead centre hit by a face of the shell of the config (CoverageUnion of PostProcessing and of the
                    Coverage_data_*.csv; the fraction Deserno's z is built from);
                    covered = the cosine-weighted solid-angle fraction in the shell (CoveredArea; the quantity of contour());
                    cap = the plain solid-angle fraction in the shell (the geometric wrapped fraction of Analyze_one_vesicle.py);
                    adhesion = f = -E_bead/(4 pi W a^2) (what the adhesion energy of the bead sees).
The a/lambda of the closed-vesicle phase space is a sqrt(KA/(kappa_over_KB KB)): kappa_over_KB = 0.5 (kappa = KB/2) in the Fig. 5 scatter as deserno_fig5 and 1 in the
binned map as contour() (a sqrt(sigma/KB)); --kappa_over_KB sets both.
"""
import argparse
import csv
import json
import math
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import Analyze_one_vesicle as A1  # noqa: E402  the energies, the shell and the coverages of a run
import Analyze_two_vesicle as A2  # noqa: E402  mesh helpers
import Coverage_data_chain as CD  # noqa: E402  the points of a chain, its run folders and the ray coverage

BLUE, PINK = "#1b4fb8", "#e0559c"          # way in, way back
COV_KEY = {"union": "union", "covered": "covered", "cap": "cap", "adhesion": "f_ad"}
COV_NAME = {"union": "CoverageUnion", "covered": "CoveredArea", "cap": "cap_cfg", "adhesion": "f_ad"}
NAN = float("nan")


# ----------------------------------------------------------------------------------------------------------------------
# the points of a chain
# ----------------------------------------------------------------------------------------------------------------------
def analyze(chain_name, chain, row, args):
    run = CD.run_folder(chain, row["run"])
    if run is None:
        raise RuntimeError("no run folder in %s" % os.path.join(chain, row["run"]))
    q = A1.analyze_run(run, row, row["scan_var"])
    cfg = json.load(open(os.path.join(run, "Input_file.json")))
    energies = {e["Name"]: e["constants"] for e in cfg["Energies"]}
    bead = cfg["Beads"][0]
    a, KI, KB = bead["radius"], bead.get("inter_str", 0.0), q["KB"]
    KA = 0.0
    for name in ("Surface_tension", "Membrane_tension", "Excess_tension"):
        if name in energies:
            KA = energies[name][0]
    planar = bool(cfg.get("boundary"))
    V_target = NAN
    if "Volume_constraint" in energies and energies["Volume_constraint"][1] > 0:
        V_target = energies["Volume_constraint"][1]               # constants: [KV, V_target]; a negative target means the volume of the start mesh (not used by the chains)
    union = multi = NAN
    if args.cov == "union":
        V, F = A2.read_obj(A2.last_mesh(run)[0])
        pos = A2.read_bead_positions(run, 1)[0]
        union, multi = CD.ray_union(V, F, pos, a, q["shell_width"], q["shell_power"], -1.0 if bead.get("outside", 1) < 0 else 1.0)
    p = {"chain": chain_name, "run_name": row["run"], "run": run, "link": q["link"], "leg": q["leg"], "scan_var": q["scan_var"], "scan_value": q["scan_value"],
         "KA": KA, "KB": KB, "KI": KI, "a": a, "wt": 4.0 * KI * a * a / KB if KI > 0 else NAN, "w_over_sigma": KI / KA if KA > 0 else NAN, "sigma_t": 2.0 * KA * a * a / KB,
         "Total_E": q["E_total"] * q["kappa"], "E_total_kappa": q["E_total"], "E_bend_kappa": q["E_bend"], "E_adh_kappa": q["E_adh"],
         "union": union, "multilayer": multi, "covered": q["cov_cfg"], "cap": q["cap_cfg"], "f_ad": q["f_ad"],
         "Volume": NAN if planar else q["Volume"], "V_target": V_target, "Area": q["Area"], "nV": q["nV"], "step": q["step"], "finished": q["finished"], "planar": planar,
         "lowest": 0}
    p["cov"] = p[COV_KEY[args.cov]]
    p["dV"] = p["Volume"] - V_target
    p["dV_rel"] = p["dV"] / V_target if V_target > 0 else NAN
    return p


def collect(chain_name, chain, args):
    rows = CD.chain_points(chain)
    if rows is None:
        print("  a tube chain (chain_summary.csv with a column kind): no bead to wrap, skipped")
        return []
    points = []
    for row in rows:
        if row["leg"] == "cold" and not args.keep_cold:
            continue
        try:
            p = analyze(chain_name, chain, row, args)
        except Exception as exc:  # a run that is missing or unreadable: reported, the chain goes on
            print("  skipped %s: %s: %s" % (row["run"], type(exc).__name__, exc), file=sys.stderr)
            continue
        if args.only_finished and not p["finished"]:
            continue
        points.append(p)
    points.sort(key=lambda p: p["link"])
    return points


def pick_lowest(points):
    """Marks p["lowest"] = 1 on one point per (chain, scanned parameter, value): of the way in and the way back for the same value the one with the lower Total_E
    (a tie keeps the way in; the cold start, if kept, does not take part). Returns the marked points and the number of values that had both legs."""
    best, count = {}, {}
    for p in points:
        if p["leg"] not in ("in", "out"):
            continue
        key = (p["chain"], p["scan_var"], round(p["scan_value"], 12))
        count[key] = count.get(key, 0) + 1
        if key not in best or p["Total_E"] < best[key]["Total_E"]:
            best[key] = p
    for p in points:
        p["lowest"] = 0
    for p in best.values():
        p["lowest"] = 1
    both = [k for k, n in count.items() if n > 1]
    return list(best.values()), len(both), sum(best[k]["leg"] == "out" for k in both)


# ----------------------------------------------------------------------------------------------------------------------
# the plots
# ----------------------------------------------------------------------------------------------------------------------
def save(fig, args, plt, folder, filename, dpi=150):
    d = os.path.join(args.out, folder)
    os.makedirs(d, exist_ok=True)
    fig.savefig(os.path.join(d, filename), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_coverage(chain_name, points, args, plt):
    """Coverage against the adhesion strength w~ = 4 W a^2 / KB of one chain, both legs; a ring marks the point of lower energy."""
    pts = [p for p in points if p["wt"] == p["wt"]]
    if not pts:
        print("  no adhesion (W = 0): no coverage plot")
        return
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    for leg, marker, color, size, alpha in (("in", "o", BLUE, 34, 1.0), ("out", "s", PINK, 22, 0.9)):
        sel = sorted([p for p in pts if p["leg"] == leg], key=lambda p: p["wt"])
        if sel:
            ax.plot([p["wt"] for p in sel], [p["cov"] for p in sel], marker + "-", color=color, ms=size ** 0.5, lw=0.9, alpha=alpha, label="way " + ("in" if leg == "in" else "back"))
    low = [p for p in pts if p["lowest"]]
    if low and any(p["leg"] == "in" for p in pts) and any(p["leg"] == "out" for p in pts):
        ax.plot([p["wt"] for p in low], [p["cov"] for p in low], "o", mfc="none", mec="k", ms=9, mew=1.0, label="lower energy")
    ws = [p["wt"] for p in pts]
    if max(ws) / min(ws) > 20:
        ax.set_xscale("log")
    ax.set_xlabel(r"adhesion strength  $\tilde w = 4Wa^2/K_B = 2wa^2/\kappa$")
    ax.set_ylabel("coverage (%s)" % COV_NAME[args.cov])
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("%s   (scan of %s, sigma = %.4g)" % (chain_name, pts[0]["scan_var"], pts[0]["KA"]), fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25, lw=0.5)
    save(fig, args, plt, "coverage_vs_adhesion", chain_name + ".png")


def plot_volume(chain_name, points, args, plt):
    """V - V_target of the final states against the scanned parameter."""
    pts = [p for p in points if p["V_target"] == p["V_target"] and p["Volume"] == p["Volume"]]
    if not pts:
        print("  no volume constraint (a flat disk): no volume plot")
        return
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    for leg, marker, color, alpha in (("in", "o", BLUE, 1.0), ("out", "s", PINK, 0.9)):
        sel = sorted([p for p in pts if p["leg"] == leg], key=lambda p: p["scan_value"])
        if sel:
            ax.plot([p["scan_value"] for p in sel], [p["dV"] for p in sel], marker + "-", color=color, ms=4.5, lw=0.9, alpha=alpha, label="way " + ("in" if leg == "in" else "back"))
    ax.axhline(0, color="0.5", lw=0.6)
    xs = [p["scan_value"] for p in pts]
    if min(xs) > 0 and max(xs) / min(xs) > 20:
        ax.set_xscale("log")
    ax.set_xlabel(pts[0]["scan_var"])
    ax.set_ylabel(r"$V - V_{target}$")
    V0 = pts[0]["V_target"]
    sec = ax.secondary_yaxis("right", functions=(lambda y: 100.0 * y / V0, lambda y: y * V0 / 100.0))
    sec.set_ylabel(r"$(V - V_{target})/V_{target}$ [%%], $V_{target}=%.5g$" % V0)
    ax.set_title(chain_name, fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25, lw=0.5)
    save(fig, args, plt, "volume_difference", chain_name + ".png")


def chain_colors(chains, key):
    """One colour per chain from the parameter of the chain (the y of its phase space: w/sigma or sigma~)."""
    import matplotlib.pyplot as plt
    vals = {c: next((p[key] for p in pts if p[key] == p[key]), NAN) for c, pts in chains.items()}
    good = [v for v in vals.values() if v == v]
    norm = plt.Normalize(min(good), max(good)) if good and max(good) > min(good) else plt.Normalize(0, 1)
    return {c: (plt.cm.viridis(norm(v)) if v == v else "0.6") for c, v in vals.items()}, plt.cm.ScalarMappable(norm=norm, cmap="viridis")


def overview_plots(chains, args, plt, axes):
    key = "w_over_sigma" if axes == "fig5" else "sigma_t"
    label = r"$w/\sigma$" if axes == "fig5" else r"$\tilde\sigma = 2\sigma a^2/\kappa$"
    colors, sm = chain_colors(chains, key)
    if "coverage" in args.plots:
        fig, ax = plt.subplots(figsize=(6.8, 5.2))
        n = 0
        for c, pts in chains.items():
            sel = sorted([p for p in pts if p["lowest"] and p["wt"] == p["wt"]], key=lambda p: p["wt"])
            if sel:
                ax.plot([p["wt"] for p in sel], [p["cov"] for p in sel], "o-", color=colors[c], ms=3.5, lw=0.9)
                n += 1
        if n:
            ws = [p["wt"] for pts in chains.values() for p in pts if p["wt"] == p["wt"]]
            if max(ws) / min(ws) > 20:
                ax.set_xscale("log")
            fig.colorbar(sm, ax=ax).set_label(label)
            ax.set_xlabel(r"adhesion strength  $\tilde w = 4Wa^2/K_B$")
            ax.set_ylabel("coverage (%s), lower-energy points" % COV_NAME[args.cov])
            ax.set_ylim(-0.02, 1.02)
            ax.grid(alpha=0.25, lw=0.5)
            save(fig, args, plt, "coverage_vs_adhesion", "all_chains.png")
        else:
            plt.close(fig)
    if "volume" in args.plots:
        fig, ax = plt.subplots(figsize=(6.8, 5.2))
        n = 0
        for c, pts in chains.items():
            for leg, marker, alpha in (("in", "o", 1.0), ("out", "s", 0.6)):
                sel = sorted([p for p in pts if p["leg"] == leg and p["dV_rel"] == p["dV_rel"]], key=lambda p: p["scan_value"])
                if sel:
                    ax.plot([p["scan_value"] for p in sel], [100.0 * p["dV_rel"] for p in sel], marker + "-", color=colors[c], ms=3, lw=0.6, alpha=alpha)
                    n += 1
        if n:
            xs = [p["scan_value"] for pts in chains.values() for p in pts]
            if min(xs) > 0 and max(xs) / min(xs) > 20:
                ax.set_xscale("log")
            fig.colorbar(sm, ax=ax).set_label(label)
            ax.axhline(0, color="0.5", lw=0.6)
            ax.set_xlabel(next(iter(chains.values()))[0]["scan_var"])
            ax.set_ylabel(r"$(V - V_{target})/V_{target}$ [%]   (circles in, squares back)")
            ax.grid(alpha=0.25, lw=0.5)
            save(fig, args, plt, "volume_difference", "all_chains.png")
        else:
            plt.close(fig)


def phase_xy(p, axes, kappa_over_KB):
    """fig5: a/lambda = a sqrt(KA/(kappa_over_KB KB)) and w/sigma = KI/KA (nan without a tension); fig2: w~ = 4 KI a^2/KB and sigma~ = 2 KA a^2/KB."""
    if axes == "fig5":
        if not p["KA"] > 0:
            return NAN, NAN
        return p["a"] * math.sqrt(p["KA"] / (kappa_over_KB * p["KB"])), p["KI"] / p["KA"]
    return 4.0 * p["KI"] * p["a"] ** 2 / p["KB"], 2.0 * p["KA"] * p["a"] ** 2 / p["KB"]


def phase_space(lowest, args, plt, axes):
    """The scatter of the lower-energy points in the axes of deserno_fig5 (closed vesicles: w/sigma against a/lambda, log x) or planar_phase_deserno (flat disks: w~ against
    sigma~), without the theory curves. The axes fit the points unless --xlim / --ylim say."""
    k = args.kappa_over_KB if args.kappa_over_KB is not None else 0.5
    xy = np.array([phase_xy(p, axes, k) for p in lowest])
    z = 2.0 * np.array([p["cov"] for p in lowest])
    good = np.isfinite(xy).all(axis=1) & np.isfinite(z) & ((xy[:, 0] > 0) if axes == "fig5" else True)
    if not good.any():
        print("no point with a surface tension and an adhesion: no phase space plot")
        return
    fig, ax = plt.subplots(figsize=(7.0, 5.4))
    last = None
    for leg, marker, size, name in (("in", "o", 40, "way in had the lower energy"), ("out", "s", 36, "way back had the lower energy")):
        sel = np.array([good[i] and lowest[i]["leg"] == leg for i in range(len(lowest))])
        if sel.any():
            last = ax.scatter(xy[sel, 0], xy[sel, 1], c=z[sel], s=size, marker=marker, cmap="viridis", vmin=0.0, vmax=2.0, edgecolor="black", linewidth=0.3, zorder=3, label=name)
    if last is not None:
        fig.colorbar(last, ax=ax).set_label(r"$z = 2\times$ coverage (%s)" % COV_NAME[args.cov])
    X, Y = xy[good, 0], xy[good, 1]
    if axes == "fig5":
        ax.set_xscale("log")
        ax.set_xlim(X.min() / 1.25, X.max() * 1.25)
        ax.set_xlabel(r"$a/\lambda$")
        ax.set_ylabel(r"$w/\sigma$")
    else:
        ax.set_xlim(X.min() - 0.05 * (X.max() - X.min() or 1.0), X.max() + 0.05 * (X.max() - X.min() or 1.0))
        ax.set_xlabel(r"$\tilde w = 2 w a^2/\kappa$")
        ax.set_ylabel(r"$\tilde\sigma = \sigma a^2/\kappa$")
    ax.set_ylim(Y.min() - 0.05 * (Y.max() - Y.min() or 1.0), Y.max() + 0.05 * (Y.max() - Y.min() or 1.0))
    if args.xlim:
        ax.set_xlim(*args.xlim)
    if args.ylim:
        ax.set_ylim(*args.ylim)
    ax.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False)       # above the axes: not on top of the points
    save(fig, args, plt, "phase_space", "phase_space_%s.png" % axes)


def phase_map(lowest, args, plt, axes, nx=24, ny=16):
    """The lower-energy points binned as contour() of Plots.py does: nx x ny bins (log x for fig5), the mean coverage of a bin, empty bins left blank, no theory."""
    k = args.kappa_over_KB if args.kappa_over_KB is not None else (1.0 if axes == "fig5" else 0.5)
    xy = np.array([phase_xy(p, axes, k) for p in lowest])
    c = np.array([p["cov"] for p in lowest])
    good = np.isfinite(xy).all(axis=1) & np.isfinite(c) & ((xy[:, 0] > 0) if axes == "fig5" else True)
    if good.sum() < 2:
        return
    X, Y, C = xy[good, 0], xy[good, 1], c[good]
    x_edges = np.logspace(np.log10(X.min()), np.log10(X.max()), nx + 1) if axes == "fig5" else np.linspace(X.min(), X.max() if X.max() > X.min() else X.min() + 1, nx + 1)
    y_edges = np.linspace(Y.min(), Y.max() if Y.max() > Y.min() else Y.min() + 1, ny + 1)
    s, _, _ = np.histogram2d(X, Y, bins=[x_edges, y_edges], weights=C)
    n, _, _ = np.histogram2d(X, Y, bins=[x_edges, y_edges])
    avg = np.full_like(s, np.nan, dtype=float)
    avg[n > 0] = s[n > 0] / n[n > 0]
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    m = ax.pcolormesh(x_edges, y_edges, np.ma.masked_invalid(avg.T), shading="flat", cmap="viridis", vmin=0.0, vmax=1.0)
    if axes == "fig5":
        ax.set_xscale("log")
        ax.set_xlabel(r"$a\,\sqrt{ \sigma / \kappa_\text{B}}$" if k == 1.0 else r"$a/\lambda$")
        ax.set_ylabel(r"$w/\sigma$")
    else:
        ax.set_xlabel(r"$\tilde w = 2 w a^2/\kappa$")
        ax.set_ylabel(r"$\tilde\sigma = \sigma a^2/\kappa$")
    fig.colorbar(m, ticks=[0.1 * i for i in range(11)]).set_label("Coverage")
    save(fig, args, plt, "phase_space", "phase_space_map.png")


# ----------------------------------------------------------------------------------------------------------------------
# the final states and the csv's
# ----------------------------------------------------------------------------------------------------------------------
def export_states(chain_name, points, args):
    """The final meshes in the order of the chain, Bead_0_data.txt, and a run folder for main_visualize (the final states as the frames of the slider)."""
    d = os.path.join(args.out, "final_states_per_chain", chain_name)
    v = os.path.join(args.out, "viewer_runs", chain_name)
    os.makedirs(d, exist_ok=True)
    os.makedirs(v, exist_ok=True)
    for folder in (d, v):                         # a run again: the frames of the earlier one go
        for f in os.listdir(folder):
            if (f.startswith("membrane_") or "_Final_link_" in f) and f.endswith(".obj"):
                os.remove(os.path.join(folder, f))
    bead_rows = []
    with open(os.path.join(d, "states.txt"), "w") as st:
        st.write("# frame link leg scan_var scan_value Total_E lowest step run_folder\n")
        for k, p in enumerate(points):
            mesh, _ = A2.last_mesh(p["run"])
            obj = os.path.join(d, "%s_Final_link_%03d.obj" % (chain_name, p["link"]))
            shutil.copyfile(mesh, obj)
            frame = os.path.join(v, "membrane_%d.obj" % k)
            try:
                os.link(obj, frame)  # a hard link, the viewer only reads regular files and this costs no space
            except OSError:
                shutil.copyfile(obj, frame)
            rows = [l.split() for l in open(os.path.join(p["run"], "Bead_0_data.txt")) if l.strip() and not l.startswith("#")]
            bead_rows.append([str(k)] + rows[-1][1:])  # the last row: the position (and the force) of the final state
            st.write("%d %d %s %s %.8g %.8g %d %d %s\n" % (k, p["link"], p["leg"], p["scan_var"], p["scan_value"], p["Total_E"], p["lowest"], p["step"], os.path.abspath(p["run"])))
    for folder in (d, v):
        with open(os.path.join(folder, "Bead_0_data.txt"), "w") as f:
            f.write("####### one row per final state: frame (index in states.txt), x y z, force fx fy fz of the last step of the run; radius %g \n" % points[0]["a"])
            for r in bead_rows:
                f.write(" ".join(r) + " \n")
    # Input_file.json of the viewer run: the config of the first state, its mesh as the initial mesh (already at its final size: no rescale, no displacement)
    cfg = json.load(open(os.path.join(points[0]["run"], "Input_file.json")))
    cfg["init_file"] = os.path.abspath(os.path.join(v, "membrane_0.obj"))
    cfg["first_dir"] = os.path.abspath(os.path.dirname(v)) + "/"
    cfg["rescale"] = 1.0
    cfg.pop("Displacement", None)
    cfg.pop("first_dir_subfolder", None)
    cfg["Beads"][0]["Pos"] = [float(x) for x in bead_rows[0][1:4]]
    with open(os.path.join(v, "Input_file.json"), "w") as f:
        json.dump(cfg, f, indent=2)
    print("  exported %d final states: %s, %s" % (len(points), d, v))


FIELDS = ["chain", "run_name", "link", "leg", "lowest", "scan_var", "scan_value", "KA", "KB", "KI", "a", "wt", "w_over_sigma", "sigma_t", "Total_E", "E_total_kappa", "E_bend_kappa", "E_adh_kappa",
          "union", "multilayer", "covered", "cap", "f_ad", "Volume", "V_target", "dV", "dV_rel", "Area", "nV", "step", "finished", "run"]


def write_csv(path, points):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(FIELDS)
        for p in points:
            w.writerow([("%.8g" % p[k]) if isinstance(p[k], float) else p[k] for k in FIELDS])


# ----------------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("batch", nargs="?", default=os.path.join(HERE, "..", "Results", "Phase_one_vesicle"), help="the folder that holds the chain folders (default Results/Phase_one_vesicle)")
    ap.add_argument("--chains", nargs="+", default=None, help="only these chain folders (the Tags, e.g. KI_1 KI_2); the phase space then has only their points")
    ap.add_argument("--out", default=None, help="output folder (default <batch>/one_bead_analysis)")
    ap.add_argument("--plots", nargs="+", choices=["phase", "coverage", "volume", "export", "csv"], default=["phase", "coverage", "volume", "export", "csv"],
                    help="what to write (default all)")
    ap.add_argument("--cov", choices=list(COV_KEY), default="union", help="coverage definition, see above (default union; the others skip the ray coverage, which is the slow part)")
    ap.add_argument("--axes", choices=["auto", "fig5", "fig2"], default="auto", help="phase space axes: fig5 = w/sigma against a/lambda (closed vesicles), fig2 = w~ against sigma~ (flat disks); auto: by the batch")
    ap.add_argument("--kappa_over_KB", type=float, default=None, help="kappa/KB of a/lambda (default 0.5 in the Fig. 5 scatter, 1 in the binned map)")
    ap.add_argument("--xlim", type=float, nargs=2, default=None, metavar=("LO", "HI"), help="x range of the phase space scatter (default: fits the points)")
    ap.add_argument("--ylim", type=float, nargs=2, default=None, metavar=("LO", "HI"))
    ap.add_argument("--keep_cold", action="store_true", help="keep the cold start of the chains (link 0): by default it is only the starting mesh")
    ap.add_argument("--only_finished", action="store_true")
    args = ap.parse_args()
    np.seterr(divide="ignore", invalid="ignore")          # the volume centre of an open (flat) mesh is 0/0 in the mesh helpers; nothing here uses it

    batch = os.path.normpath(args.batch)
    args.out = args.out or os.path.join(batch, "one_bead_analysis")
    tags = args.chains or sorted(d for d in os.listdir(batch) if os.path.exists(os.path.join(batch, d, "chain_summary.csv")))
    if not tags:
        raise SystemExit("no <tag>/chain_summary.csv in %s (the batch is the folder that holds the chain folders)" % batch)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chains = {}
    for tag in tags:
        chain = os.path.join(batch, tag)
        if not os.path.exists(os.path.join(chain, "chain_summary.csv")):
            print("%s: no chain_summary.csv, skipped" % tag)
            continue
        print("\n%s" % tag)
        points = collect(tag, chain, args)
        if not points:
            print("  no points")
            continue
        low, n_both, n_back = pick_lowest(points)
        print("  %d points (%d values with both legs: the way back has the lower energy in %d)" % (len(points), n_both, n_back))
        print("   link  leg   %-9s  w~       coverage  Total_E          V-V_target" % points[0]["scan_var"])
        for p in points:
            print("   %4d  %-4s  %-9.5g  %-8.4g  %7.4f  %15.8g  %-10.4g %s" % (p["link"], p["leg"], p["scan_value"], p["wt"], p["cov"], p["Total_E"], p["dV"], "*" if p["lowest"] else ""))
        chains[tag] = points
        if "csv" in args.plots:
            write_csv(os.path.join(args.out, "csv", tag + ".csv"), points)
        if "export" in args.plots:
            export_states(tag, points, args)
        if "coverage" in args.plots:
            plot_coverage(tag, points, args, plt)
        if "volume" in args.plots:
            plot_volume(tag, points, args, plt)
    if not chains:
        raise SystemExit("no points in any chain")

    every = [p for pts in chains.values() for p in pts]
    lowest = [p for p in every if p["lowest"]]
    axes = args.axes if args.axes != "auto" else ("fig2" if all(p["planar"] for p in every) else "fig5")
    print("\nphase space: %d points of %d chains, axes %s" % (len(lowest), len(chains), axes))
    if "phase" in args.plots:
        phase_space(lowest, args, plt, axes)
        phase_map(lowest, args, plt, axes)
        write_csv(os.path.join(args.out, "phase_space", "phase_space_lowest.csv"), lowest)
    if "csv" in args.plots:
        write_csv(os.path.join(args.out, "csv", "all_points.csv"), every)
    overview_plots(chains, args, plt, axes)
    print("\nwrote results below %s/" % args.out)


if __name__ == "__main__":
    main()
