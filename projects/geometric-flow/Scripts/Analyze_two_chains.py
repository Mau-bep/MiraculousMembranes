#!/usr/bin/env python3
"""Per-chain analysis of the two-bead runs (Create_subjob_two_vesicle.py / Run_two_vesicle_chain.py), e.g. the shell_KI*_rep*_AreaConstraint_{correct,noexcess} scans.

usage:  python3 Analyze_two_chains.py                                  # every coverage_r_*_bigmem_correctratio_*_excess*_AreaConstraint folder of --root
        python3 Analyze_two_chains.py --ratio 0.7 --excess 0.001       # only some of them (filters: --r --ratio --excess --KI --rep)
        python3 Analyze_two_chains.py --paper --plots energy coverage  # paper curve on the energy plot, skip the slices
        python3 Analyze_two_chains.py --pattern 'shell_KI*_rep*_AreaConstraint_correct' 'shell_KI*_rep*_AreaConstraint_noexcess'   # the Adhesion KI scans
        python3 Analyze_two_chains.py ../Results/Wrapping_two_vesicle/coverage_r_0.5_bigmem_correctratio_0.7_excess0.001_AreaConstraint   # explicit chain folders

Needs numpy + matplotlib (micromamba env mir_mem: ~/micromamba/envs/mir_mem/bin/python).

A chain is one folder with link*/<n>/ run folders (the cold start link00_cold and the Coverage pre-wrap are skipped, as in Analyze_two_vesicle.py).
Per chain it writes, below --out (default two_beads_analysis/):
    final_states_per_chain/<chain>/<link>.png     slice of every final state (plane through the two bead centres and the volume centre)
    final_states_per_chain/<chain>_sheet.png      all the slices of the chain on one sheet, in chain order
    energy_vs_distance/<chain>.png                bending energy u(s) = E_bend(s) - E_bend(far) in kappa against the distance (--paper adds Fig. 4a of van der Wel et al.)
    coverage_scatter/<chain>.png                  coverage of bead 1 (y) against bead 0 (x), one point per simulation, grey dashed line x = y
    area_difference/<chain>.png                   final area minus the target area of the Area_constraint (A - A_target, and in percent of A_target) against the distance
    coverage_scatter/<chain>_cap.png              the same with the plain solid-angle fraction in the shell of the config (cap_cfg0, cap_cfg1)
    final_states_per_chain/<chain>/<link>.obj     the final mesh of every state, Bead_0_data.txt, Bead_1_data.txt (frame, x y z, force; frame = index
                                                  in states.txt) and states.txt (frame -> link, s, d, step, run folder)
    viewer_runs/<chain>/                          a run folder for main_visualize (Input_file.json, membrane_<frame>.obj as hard links, Bead_*_data.txt):
                                                  main_visualize <out>/viewer_runs   (the slider goes through the final states; the folders of the menu are the chains)
    csv/<chain>.csv                               the numbers of every run
"Coverage" (--cov): auto (default) = adhesion if the beads have an adhesion strength W > 0, else coverage;
                    coverage = the cosine-weighted solid-angle fraction in the +-25 % shell, the quantity the Coverage energy constrains (its target of the
                               config is drawn as a dotted line in the scatter);
                    adhesion = f = -E_bead/(4 pi W a^2), the fraction the adhesion energy of that bead sees (0 = bare, 1 = closed sphere);
                    cap = plain solid-angle fraction of the faces facing the bead inside the shell of the config (geometric, includes the neck).
Display instead of saving: --show, with the plot types as flags (--energy --coverage --cap --area_diff --slices; none = the first four) and the leg with
--legs in|out|both, e.g.  python3 Analyze_two_chains.py --ratio 0.7 --excess 0.001 --show --energy --legs in
Without --show the type flags choose what is saved (instead of --plots); --legs filters the saved plots too (the csv and the export keep everything).
Markers: circle = way in (s decreasing along the chain), square = way back; red cross = beads overlap (d < 2a).
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import Analyze_two_vesicle as A2  # noqa: E402  run folders, mesh helpers, energies, the paper's digitized curve

# DEFAULT_ROOT = os.path.join(HERE, "..", "Results", "Wrapping_two_vesicle")
DEFAULT_ROOT = os.path.join(HERE, "..", "Results", "WrappingTwoVesicle")

LINK_RE = re.compile(r"link(\d+)_s([0-9.]+)(?:_w(\d+))?$")


# ----------------------------------------------------------------------------------------------------------------------
PARAMS = {"r": r"coverage_r_([0-9.]+?)_", "ratio": r"(?:correctratio|_coverage)_([0-9.]+?)_excess", "excess": r"excess([0-9.]+?)_", "KI": r"shell_KI([0-9.]+?)_", "rep": r"_rep(?:lica_)?(\d+)"}


def pick_chains(args):
    if args.chains:
        return [os.path.normpath(c) for c in args.chains]
    out = []
    for pattern in args.pattern:
        for path in sorted(glob.glob(os.path.join(args.root, pattern))):
            name = os.path.basename(path)
            ok = True
            for key, rx in PARAMS.items():
                wanted = getattr(args, key)
                if not wanted:
                    continue
                m = re.search(rx, name)
                if not m or not any(abs(float(m.group(1)) - w) < 1e-9 for w in wanted):
                    ok = False
            if ok and os.path.isdir(path) and path not in out:
                out.append(path)
    return out


def last_row(run):
    """(column index of every 'Bead' column, last data row) of Output_data.txt without parsing the whole file."""
    last = None
    with open(os.path.join(run, "Output_data.txt")) as f:
        header = f.readline().split()
        for line in f:
            if line.strip():
                last = line
    return [i for i, n in enumerate(header) if n == "Bead"], [float(x) for x in last.split()]


def analyze(run, link_name, cov_mode):
    import json
    p = A2.analyze_run(run)
    cfg = json.load(open(os.path.join(run, "Input_file.json")))
    bead_cols, row = last_row(run)
    kappa = p["kappa"]
    mesh, _ = A2.last_mesh(run)
    V, F = A2.read_obj(mesh)
    pos = A2.read_bead_positions(run)
    m = LINK_RE.match(link_name)
    p["link"], p["s_target"], p["warm"] = int(m.group(1)), float(m.group(2)), int(m.group(3) or 0)
    p["link_name"] = link_name
    for i, b in enumerate(cfg["Beads"][:2]):
        a, W = b["radius"], b.get("inter_str", 0.0)
        shell, power = b.get("shell_width", 0.25), b.get("shell_power", 1.0)
        side = -1.0 if b.get("outside", 1) < 0 else 1.0
        _, cap, _ = A2.coverage(V, F, pos[i], a, width=shell, power=power, side=side)
        E_code = row[bead_cols[i]] if i < len(bead_cols) else 0.0
        p["E_adh%d" % i] = E_code / kappa
        p["f_ad%d" % i] = -E_code / (4.0 * math.pi * W * a * a) if W > 0 else float("nan")
        p["cap_cfg%d" % i] = cap
        kind = cov_mode
        if kind == "auto" or (kind == "adhesion" and not W > 0):
            kind = "adhesion" if W > 0 else "coverage"
        p["frac%d" % i] = {"adhesion": p["f_ad%d" % i], "coverage": p["cov%d" % i], "cap": cap}[kind]
        p["cov_kind"] = kind
        p["shell_w"] = shell
    area_e = [e["constants"] for e in cfg["Energies"] if e["Name"] == "Area_constraint"]
    p["A_target"] = area_e[0][1] if area_e else float("nan")  # constants: [KA, target area, ...]
    p["dA"] = p["Area"] - p["A_target"]
    p["dA_rel"] = p["dA"] / p["A_target"]
    p["W"] = cfg["Beads"][0].get("inter_str", float("nan"))
    cov_e = [e["constants"] for e in cfg["Energies"] if e["Name"] == "Coverage"]
    p["target0"], p["target1"] = (cov_e[0][1], cov_e[0][2]) if cov_e and len(cov_e[0]) >= 3 else (float("nan"), float("nan"))
    p["pos0"], p["pos1"] = pos[0], pos[1]
    p["mesh"] = (V, F)
    return p


def chain_xaxis(chain, wanted):
    """arc or chord: the distance the chain was run with (column 'distance' of chain_summary.csv), unless --xaxis says."""
    if wanted != "auto":
        return wanted
    path = os.path.join(chain, "chain_summary.csv")
    if os.path.exists(path):
        for row in csv.DictReader(open(path)):
            return row.get("distance") or "arc"
    return "arc"


def collect(chain, args):
    points = []
    for run in A2.find_runs(chain):
        link_name = os.path.basename(os.path.dirname(os.path.normpath(run)))
        if not LINK_RE.match(link_name):
            continue  # not a link of a chain
        try:
            p = analyze(run, link_name, args.cov)
        except Exception as exc:
            print("  skipped %s: %s" % (run, exc), file=sys.stderr)
            continue
        if args.only_finished and not p["finished"]:
            continue
        points.append(p)
    if not points:
        return points
    points.sort(key=lambda p: (p["link"], p["warm"] or 99))  # w1, w2, then the target itself
    # way in until the smallest target separation, way back after it
    link_min = min(points, key=lambda p: (p["s_target"], -p["link"]))["link"]
    for p in points:
        p["leg"] = "in" if p["link"] <= link_min else "out"
    pts_s = [p for p in points if p["s_over_Dp"] > args.s_far]
    key = "s_over_Dp" if args.xaxis_chain == "arc" else "s_chord"
    ref_pts = [p for p in points if p[key] > args.s_far] or [max(points, key=lambda p: p[key])]
    e0 = float(np.mean([p["Ebend_kappa"] for p in ref_pts]))
    e0_tot = float(np.mean([p["Etot_kappa"] for p in ref_pts]))
    for p in points:
        p["u_kappa"] = p["Ebend_kappa"] - e0
        p["u_tot_kappa"] = p["Etot_kappa"] - e0_tot
        p["x"] = p[key]
    return points


# ----------------------------------------------------------------------------------------------------------------------
def finish(fig, args, plt, folder, filename, dpi=150):
    """Save the figure as <out>/<folder>/<filename>, or, with --show, leave it open: main calls plt.show() once at the end."""
    if args.show:
        return
    d = os.path.join(args.out, folder)
    os.makedirs(d, exist_ok=True)
    fig.savefig(os.path.join(d, filename), dpi=dpi)
    plt.close(fig)


def slice_segments(V, F, pos, args):
    """(segments (n,2,2) in the (u, w) plane, bead centres in that plane, origin) of the plane through both beads and the volume centre."""
    _, _, vc = A2.mesh_geometry(V, F)
    u = pos[0] - pos[1]
    u = u / np.linalg.norm(u)
    w = vc - 0.5 * (pos[0] + pos[1])
    w -= u * (w @ u)
    if np.linalg.norm(w) < 1e-9:
        w = np.cross(u, [0.0, 0.0, 1.0])
    w /= np.linalg.norm(w)
    n = np.cross(u, w)
    P = V - vc
    z = P @ n
    zf = z[F]
    sel = (zf.min(1) < 0) & (zf.max(1) > 0)
    Fs, zs = F[sel], zf[sel]
    Q, C = [], []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        c = (zs[:, i] < 0) != (zs[:, j] < 0)
        t = np.where(c, zs[:, i] / np.where(c, zs[:, i] - zs[:, j], 1.0), 0.0)
        Q.append(P[Fs[:, i]] + t[:, None] * (P[Fs[:, j]] - P[Fs[:, i]]))
        C.append(c)
    Q, C = np.stack(Q, 1), np.stack(C, 1)            # (n, 3 edges, 3 xyz), (n, 3)
    order = np.argsort(~C, axis=1, kind="stable")[:, :2]   # the two edges that are crossed
    pts = np.take_along_axis(Q, order[:, :, None], axis=1)
    seg = np.stack([pts @ u, pts @ w], axis=2)       # (n, 2 points, 2 plane coords)
    centres = np.array([[(q - vc) @ u, (q - vc) @ w] for q in pos])
    return seg, centres


def draw_slice(ax, p, half, shell_arg, show_shell, title_size=8):
    from matplotlib.collections import LineCollection
    import matplotlib.pyplot as plt
    V, F = p["mesh"]
    a = p["a"]
    shell_w = shell_arg if shell_arg is not None else p["shell_w"]
    seg, centres = slice_segments(V, F, [p["pos0"], p["pos1"]], None)
    ax.add_collection(LineCollection(seg, colors="k", linewidths=0.7))
    for c in centres:
        ax.add_patch(plt.Circle(c, a, fill=False, color="tab:red", lw=1.1))
        if show_shell:
            for r in (a * (1 - shell_w), a * (1 + shell_w)):
                ax.add_patch(plt.Circle(c, r, fill=False, color="tab:red", lw=0.4, ls=":"))
    ax.set_title("%s  d=%.2f\ncov %.2f | %.2f" % (p["link_name"], p["d"], p["frac0"], p["frac1"]), fontsize=title_size)
    ax.set_aspect("equal")
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_xticks([])
    ax.set_yticks([])


def make_video(chain_name, points, slice_dir, args):
    """mp4 of the slices in chain order (w1, w2, target of every link), every frame shown --frame_time seconds (needs ffmpeg)."""
    import subprocess
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for p in points:
            f.write("file '%s'\nduration %g\n" % (os.path.abspath(os.path.join(slice_dir, p["link_name"] + ".png")), args.frame_time))
        f.write("file '%s'\n" % os.path.abspath(os.path.join(slice_dir, points[-1]["link_name"] + ".png")))  # concat needs the last file twice
        lst = f.name
    path = os.path.join(args.out, "final_states_per_chain", chain_name + ".mp4")
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", lst, "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white,format=yuv420p",
           "-r", "30", "-c:v", "libx264", path]
    subprocess.run(cmd, check=True)
    os.remove(lst)
    print("  wrote %s (%d frames of %g s)" % (path, len(points), args.frame_time))


def plot_slices(chain_name, points, args, plt):
    d = os.path.join(args.out, "final_states_per_chain", chain_name)
    if not args.show:  # --show displays the sheet only
        os.makedirs(d, exist_ok=True)
        for p in points:
            fig, ax = plt.subplots(figsize=(4.4, 4.6))
            draw_slice(ax, p, args.half, args.shell_width, args.shell, title_size=9)
            fig.tight_layout()
            fig.savefig(os.path.join(d, "%s.png" % p["link_name"]), dpi=110)
            plt.close(fig)
        if args.video:
            make_video(chain_name, points, d, args)
        if args.no_sheet:
            return
    n = len(points)
    ncol = min(args.ncols, n)
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.4 * ncol, 2.6 * nrow), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, p in zip(axes.ravel(), points):
        ax.axis("on")
        draw_slice(ax, p, args.half, args.shell_width, args.shell, title_size=6.5)
    fig.suptitle(chain_name, fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    finish(fig, args, plt, "final_states_per_chain", chain_name + "_sheet.png", dpi=100)


# ----------------------------------------------------------------------------------------------------------------------
def plot_energy(chain_name, points, args, plt):
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    if args.paper:
        paper = A2.load_paper()
        if "experiment_floppy" in paper:
            e = paper["experiment_floppy"]
            ax.plot(e[:, 0], e[:, 1], "o", color="tab:red", ms=3.5, alpha=0.45, label="experiment, floppy (paper Fig. 4a)")
        if "simulation" in paper:
            e = paper["simulation"]
            ax.plot(e[:, 0], e[:, 1], "D", color="k", ms=4, label="paper simulation (Fig. 4a)")
    curves = []
    if args.curves in ("bending", "both"):
        curves.append(("u_kappa", "bending", "tab:blue", "-"))
    if args.curves in ("total", "both"):
        curves.append(("u_tot_kappa", "total", "tab:orange", "--"))
    for key, name, color, ls in curves:
        for leg, marker, alpha in (("in", "o", 1.0), ("out", "s", 0.6)):
            sel = sorted([p for p in points if p["leg"] == leg], key=lambda p: p["x"])
            if sel:
                ax.plot([p["x"] for p in sel], [p[key] for p in sel], marker + ls, color=color, ms=4, lw=1.0, alpha=alpha,
                        label="%s, way %s" % (name, leg))
    bad = [p for p in points if p["overlap"]]
    if bad:
        ax.plot([p["x"] for p in bad], [p["u_kappa"] for p in bad], "x", color="tab:red", ms=7, mew=1.5, label="beads overlap (d < 2a)")
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xlabel(r"$s\ [D_p]$  (arc of a sphere of radius $R_v$)" if args.xaxis_chain == "arc" else r"$d/D_p$  (distance between the bead centres)")
    ax.set_ylabel(r"$u\ [\kappa]$,  $\kappa=K_B/2=%.4g$" % points[0]["kappa"])
    if args.xlim:
        ax.set_xlim(*args.xlim)
    if args.ylim:
        ax.set_ylim(*args.ylim)
    sec = ax.secondary_yaxis("right", functions=(lambda y: y * args.kappa_kbt, lambda y: y / args.kappa_kbt))
    sec.set_ylabel(r"$u\ [k_BT]$ for $\kappa=%.1f\,k_BT$" % args.kappa_kbt)
    ax.set_title(chain_name, fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    finish(fig, args, plt, "energy_vs_distance", chain_name + ".png")


def plot_coverage(chain_name, points, args, plt, key="frac", suffix=""):
    """key: 'frac' (the coverage of --cov) or 'cap_cfg' (the plain solid-angle fraction in the shell of the config); suffix goes to the file name."""
    k0, k1 = key + "0", key + "1"
    fig, ax = plt.subplots(figsize=(6.0, 5.4))
    vals = [p[k0] for p in points] + [p[k1] for p in points]
    lo, hi = (args.cov_lim if args.cov_lim else (min(0.0, min(vals) - 0.02), max(1.0, max(vals) + 0.02)))
    ax.plot([lo, hi], [lo, hi], "--", color="0.5", lw=1.0, label=r"$x=y$ (equal wrapping)", zorder=1)
    sc = None
    for leg, marker in (("in", "o"), ("out", "s")):
        sel = [p for p in points if p["leg"] == leg]
        if sel:
            sc = ax.scatter([p[k0] for p in sel], [p[k1] for p in sel], c=[p["x"] for p in sel], cmap="viridis_r",
                            vmin=min(p["x"] for p in points), vmax=max(p["x"] for p in points), marker=marker, s=28, edgecolor="k", linewidth=0.3,
                            label="way " + leg, zorder=2)
    bad = [p for p in points if p["overlap"]]
    if bad:
        ax.plot([p[k0] for p in bad], [p[k1] for p in bad], "x", color="tab:red", ms=7, mew=1.2, label="beads overlap (d < 2a)", zorder=3)
    if sc is not None:
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label(r"$s\ [D_p]$" if args.xaxis_chain == "arc" else r"$d/D_p$")
    kind = points[0]["cov_kind"] if key == "frac" else "cap"
    what = {"adhesion": "f = -E_bead/(4 pi W a^2)", "coverage": "weighted solid angle in the shell", "cap": "solid angle in the shell of the config"}[kind]
    t0, t1 = points[0]["target0"], points[0]["target1"]
    if not math.isnan(t0) and not math.isnan(t1):  # the target c of the Coverage energy of the config, on both scatters
        ax.axvline(t0, color="tab:green", ls=":", lw=1.0, label="target c = %.3g" % t0 + ("" if abs(t0 - t1) < 1e-12 else " | %.3g" % t1), zorder=1)
        ax.axhline(t1, color="tab:green", ls=":", lw=1.0, zorder=1)
    ax.set_xlabel("coverage of bead 0  (%s)" % what)
    ax.set_ylabel("coverage of bead 1")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_title(chain_name, fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    finish(fig, args, plt, "coverage_scatter", chain_name + suffix + ".png")


def plot_area(chain_name, points, args, plt):
    """Final area minus the target area of the Area_constraint, against the distance."""
    if math.isnan(points[0]["A_target"]):
        print("  no Area_constraint in the config: no area plot")
        return
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    for leg, marker, alpha in (("in", "o", 1.0), ("out", "s", 0.6)):
        sel = sorted([p for p in points if p["leg"] == leg], key=lambda p: p["x"])
        if sel:
            ax.plot([p["x"] for p in sel], [p["dA"] for p in sel], marker + "-", color="tab:blue", ms=4, lw=1.0, alpha=alpha, label="way " + leg)
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xlabel(r"$s\ [D_p]$  (arc of a sphere of radius $R_v$)" if args.xaxis_chain == "arc" else r"$d/D_p$  (distance between the bead centres)")
    ax.set_ylabel(r"$A - A_{target}$")
    A0 = points[0]["A_target"]
    sec = ax.secondary_yaxis("right", functions=(lambda y: 100.0 * y / A0, lambda y: y * A0 / 100.0))
    sec.set_ylabel(r"$(A - A_{target})/A_{target}$ [%%], $A_{target}=%.5g$" % A0)
    if args.xlim:
        ax.set_xlim(*args.xlim)
    ax.set_title(chain_name, fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    finish(fig, args, plt, "area_difference", chain_name + ".png")


def export_states(chain_name, points, args):
    """The final meshes, the bead positions and a run folder for main_visualize (the final states as the frames of the slider)."""
    import json
    import shutil
    d = os.path.join(args.out, "final_states_per_chain", chain_name)
    v = os.path.join(args.out, "viewer_runs", chain_name)
    os.makedirs(d, exist_ok=True)
    os.makedirs(v, exist_ok=True)
    for f in os.listdir(v):
        if f.startswith("membrane_") and f.endswith(".obj"):
            os.remove(os.path.join(v, f))
    nbeads = 2
    bead_rows = [[] for _ in range(nbeads)]
    with open(os.path.join(d, "states.txt"), "w") as st:
        st.write("# frame link s_over_Dp d step leg run_folder\n")
        for k, p in enumerate(points):
            mesh, _ = A2.last_mesh(p["run"])
            obj = os.path.join(d, p["link_name"] + ".obj")
            shutil.copyfile(mesh, obj)
            frame = os.path.join(v, "membrane_%d.obj" % k)
            try:
                os.link(obj, frame)  # a hard link, the viewer only reads regular files and this costs no space
            except OSError:
                shutil.copyfile(obj, frame)
            for i in range(nbeads):
                rows = [l.split() for l in open(os.path.join(p["run"], "Bead_%d_data.txt" % i)) if l.strip() and not l.startswith("#")]
                bead_rows[i].append([str(k)] + rows[-1][1:])  # the last row: the position (and the force) of the final state
            st.write("%d %s %.4f %.4f %d %s %s\n" % (k, p["link_name"], p["s_over_Dp"], p["d"], p["step"], p["leg"], os.path.abspath(p["run"])))
    for i in range(nbeads):
        for folder in (d, v):
            with open(os.path.join(folder, "Bead_%d_data.txt" % i), "w") as f:
                f.write("####### one row per final state: frame (index in states.txt), x y z, force fx fy fz of the last step of the run; radius %g \n" % points[0]["a"])
                for r in bead_rows[i]:
                    f.write(" ".join(r) + " \n")
    # Input_file.json of the viewer run: the config of the first state, pointing at the first frame; the viewer only needs the mesh and the beads
    cfg = json.load(open(os.path.join(points[0]["run"], "Input_file.json")))
    cfg["init_file"] = os.path.abspath(os.path.join(v, "membrane_0.obj"))
    cfg["first_dir"] = os.path.abspath(os.path.dirname(v)) + "/"
    cfg.pop("first_dir_subfolder", None)
    for i, b in enumerate(cfg["Beads"][:nbeads]):
        b["Pos"] = [float(x) for x in bead_rows[i][0][1:4]]
    with open(os.path.join(v, "Input_file.json"), "w") as f:
        json.dump(cfg, f, indent=2)
    print("  exported %d states: %s, %s" % (len(points), d, v))


def write_csv(chain_name, points, args):
    d = os.path.join(args.out, "csv")
    os.makedirs(d, exist_ok=True)
    fields = ["link_name", "leg", "s_target", "s_over_Dp", "s_chord", "d", "Ebend_kappa", "Etot_kappa", "u_kappa", "u_tot_kappa",
              "frac0", "frac1", "f_ad0", "f_ad1", "cap_cfg0", "cap_cfg1", "E_adh0", "E_adh1", "Area", "A_target", "dA", "dA_rel", "Volume", "sigma_eff", "nV", "step", "finished", "overlap", "run"]
    with open(os.path.join(d, chain_name + ".csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for p in points:
            w.writerow([p[k] for k in fields])


# ----------------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("chains", nargs="*", help="chain folders (default: the folders of --root that match --pattern)")
    ap.add_argument("--root", default=DEFAULT_ROOT, help="folder with the chains (default: Results/Wrapping_two_vesicle)")
    ap.add_argument("--pattern", nargs="+", default=["TwoInteraction_r_*_coverage_*_excess*_replica_0_AreaConstraint"],
                    help="glob(s) of the chain folders inside --root (quote them)")
    ap.add_argument("--r", type=float, nargs="+", default=None, help="only these bead radii (the r of coverage_r_<r>_...)")
    ap.add_argument("--ratio", "--c", type=float, nargs="+", default=None, help="only these coverage targets c (the correctratio_<c> or coverage_<c> of the folder name)")
    ap.add_argument("--excess", type=float, nargs="+", default=None, help="only these area excesses")
    ap.add_argument("--KI", type=float, nargs="+", default=None, help="only these adhesion strengths (shell_KI<KI> folders)")
    ap.add_argument("--rep", type=float, nargs="+", default=None, help="only these replicas (shell_..._rep<n> folders)")
    ap.add_argument("--out", default="two_beads_analysis", help="output folder")
    ap.add_argument("--plots", nargs="+", choices=["energy", "coverage", "cap", "area", "slices", "csv", "export"],
                    default=["energy", "coverage", "cap", "area", "slices", "csv", "export"], help="what to save (default all; slices are the slow part)")
    ap.add_argument("--energy", action="store_true", help="plot type: bending energy against the distance")
    ap.add_argument("--coverage", action="store_true", help="plot type: coverage scatter of bead 1 against bead 0")
    ap.add_argument("--cap", action="store_true", help="plot type: the same scatter with cap_cfg (solid angle in the shell of the config)")
    ap.add_argument("--area_diff", action="store_true", help="plot type: A - A_target against the distance")
    ap.add_argument("--slices", action="store_true", help="plot type: the slices of the final states (with --show: the sheet of each chain)")
    ap.add_argument("--show", action="store_true",
                    help="display the plots in windows instead of saving them (nothing is written; csv and export are skipped). Without plot-type flags: "
                         "energy, coverage, cap and area_diff")
    ap.add_argument("--legs", choices=["in", "out", "both"], default="both",
                    help="which leg of the chain the plots show: in (way in), out (way back) or both (default). The zero of u still comes from the whole chain")
    ap.add_argument("--cov", choices=["auto", "coverage", "adhesion", "cap"], default="auto", help="coverage definition, see above (default auto)")
    ap.add_argument("--cov_lim", type=float, nargs=2, default=None, metavar=("LO", "HI"), help="axis range of the coverage scatter (default 0..1)")
    ap.add_argument("--paper", action="store_true", help="overlay the digitized Fig. 4a of van der Wel et al. on the energy plot")
    ap.add_argument("--curves", choices=["bending", "total", "both"], default="bending", help="energy difference(s) on the energy plot (default bending)")
    ap.add_argument("--xaxis", choices=["auto", "arc", "chord"], default="auto",
                    help="auto (default): what the chain was run with (the distance column of chain_summary.csv; chains without it: arc); arc: s = R_v phi; chord: d/D_p")
    ap.add_argument("--s_far", type=float, default=4.0, help="runs with s above this (D_p) give the zero of u (default 4.0)")
    ap.add_argument("--xlim", type=float, nargs=2, default=None, metavar=("LO", "HI"))
    ap.add_argument("--ylim", type=float, nargs=2, default=None, metavar=("LO", "HI"))
    ap.add_argument("--kappa_kbt", type=float, default=20.9)
    ap.add_argument("--only_finished", action="store_true")
    ap.add_argument("--half", type=float, default=8.0, help="half width of a slice panel in bead radii")
    ap.add_argument("--shell", action="store_true", help="draw the adhesion shell (radius a(1 +- shell_width)) around the beads in the slices")
    ap.add_argument("--shell_width", type=float, default=None, help="shell half width drawn by --shell (default: that of the run: shell_width of the beads, else 0.25)")
    ap.add_argument("--video", action="store_true", help="also make final_states_per_chain/<chain>.mp4 from the slices (needs ffmpeg)")
    ap.add_argument("--frame_time", type=float, default=0.5, help="seconds per frame of --video (default 0.5)")
    ap.add_argument("--no_sheet", action="store_true", help="no all-in-one sheet of the slices")
    ap.add_argument("--ncols", type=int, default=6, help="columns of the slice sheet")
    args = ap.parse_args()

    chains = pick_chains(args)
    if not chains:
        raise SystemExit("no chain folders found (root %s, patterns %s)" % (args.root, args.pattern))
    print("chains:\n  " + "\n  ".join(chains))

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # plot types: the flags (--energy --coverage --cap --area_diff --slices) if any is given, else --plots (saving) or the four 2D plots (--show)
    flagged = [t for t, on in (("energy", args.energy), ("coverage", args.coverage), ("cap", args.cap), ("area", args.area_diff), ("slices", args.slices)) if on]
    if flagged:
        types = set(flagged)
    elif args.show:
        types = {"energy", "coverage", "cap", "area"}
    else:
        types = set(args.plots)
    if args.show:
        types &= {"energy", "coverage", "cap", "area", "slices"}

    for chain in chains:
        name = os.path.basename(chain)
        print("\n%s" % name)
        args.xaxis_chain = chain_xaxis(chain, args.xaxis)
        points = collect(chain, args)
        if not points:
            print("  no runs")
            continue
        print("  %d runs; reference E_bend from s > %g D_p" % (len(points), args.s_far))
        print("   link             leg   s/Dp    u/k    cov0   cov1   (cov1-cov0)")
        for p in points:
            print("   %-16s %-4s %6.2f %7.3f  %5.3f  %5.3f   %+6.3f%s" % (p["link_name"], p["leg"], p["s_over_Dp"], p["u_kappa"], p["frac0"], p["frac1"],
                                                                     p["frac1"] - p["frac0"], "  OVERLAP" if p["overlap"] else ""))
        if "csv" in types and not args.show:
            write_csv(name, points, args)
        if "export" in types and not args.show:
            export_states(name, points, args)
        shown = [p for p in points if args.legs in ("both", p["leg"])]  # the plots show the chosen leg only (csv and export keep everything)
        if not shown:
            print("  no runs on leg '%s'" % args.legs)
        else:
            if "energy" in types:
                plot_energy(name, shown, args, plt)
            if "coverage" in types:
                plot_coverage(name, shown, args, plt)
            if "cap" in types:
                plot_coverage(name, shown, args, plt, key="cap_cfg", suffix="_cap")
            if "area" in types:
                plot_area(name, shown, args, plt)
            if "slices" in types:
                plot_slices(name, shown, args, plt)
        for p in points:
            del p["mesh"]  # free the meshes
    if args.show:
        print("\n%d figure(s) open; close the windows to finish" % len(plt.get_fignums()))
        plt.show()
    else:
        print("\nwrote results below %s/" % args.out)


if __name__ == "__main__":
    main()
