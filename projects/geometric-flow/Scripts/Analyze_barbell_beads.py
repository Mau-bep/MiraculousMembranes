#!/usr/bin/env python3
"""Tube of the barbell with beads of radius 1 (Create_subjob_barbell_beads.py): tube radius, area and energies against KA and KB.

usage:  python3 Analyze_barbell_beads.py [--out PREFIX] [--group KB|KA] RESULTS_DIR [RESULTS_DIR ...]

Every folder below a RESULTS_DIR that holds an Input_file.json, an Output_data.txt and a mesh (Final_state.obj, else the last membrane_*.obj) is one run. Per run:
    KA, KB       the Surface_tension and Bending constants of the config; d the distance of the bead centres, a the bead radius
    R_tube       radius of the tube: mean distance from the axis of the vertices of the middle third between the beads (R_min: the smallest ring radius there)
    R_theory     sqrt(KB / (4 KA)), the radius of a free tube (minimum of 2 pi R L (KA + KB/(4 R^2)))
    spread       relative spread of the radius over the middle third (std / mean): the tube is uniform when it is small (a few percent)
    L_tube       length of the tube, the x range where the radius stays below 1.5 R_tube between the beads
    Area, E_bend, E_tension, E_total  from the last row of Output_data.txt; steps = the last step; converged = the run ended before the step limit of the config
Writes <out>.csv and <out>.png (left: R_tube against R_theory, one marker per KB, the line is y = x; right: R_tube / R_theory against KA).
--group KA: one series per KA and the right panel against KB (for a scan of KB at fixed surface tension, the Tension batch of Mem_barbell_beads_up.sh).
Only needs json, csv and numpy (matplotlib for the plot).
"""
import argparse
import csv
import glob
import json
import math
import os
import re

import numpy as np


def read_obj(path):
    V, F = [], []
    with open(path) as f:
        for line in f:
            p = line.split()
            if p and p[0] == "v":
                V.append([float(x) for x in p[1:4]])
            elif p and p[0] == "f":
                F.append([int(x.split("/")[0]) - 1 for x in p[1:4]])
    return np.array(V), np.array(F, dtype=int)


def find_runs(root):
    runs = []
    for dirpath, dirnames, filenames in os.walk(root):
        if "Input_file.json" in filenames and "Output_data.txt" in filenames:
            runs.append(dirpath)
            dirnames[:] = []
    return sorted(runs)


def last_mesh(run):
    final = os.path.join(run, "Final_state.obj")
    if os.path.exists(final):
        return final
    meshes = glob.glob(os.path.join(run, "membrane_*.obj"))
    if not meshes:
        return None
    return max(meshes, key=lambda p: int(re.findall(r"membrane_(\d+)", p)[0]))


def analyze_run(run):
    cfg = json.load(open(os.path.join(run, "Input_file.json")))
    mesh = last_mesh(run)
    if mesh is None:
        return None
    energies = {e["Name"]: e["constants"][0] for e in cfg["Energies"]}
    KA = energies.get("Surface_tension", float("nan"))
    KB = [v for k, v in energies.items() if k in ("Bending", "Bending_tan")][0]
    xs = sorted(b["Pos"][0] for b in cfg["Beads"])
    a = cfg["Beads"][0]["radius"]
    d = xs[-1] - xs[0]
    V, F = read_obj(mesh)
    rho = np.hypot(V[:, 1], V[:, 2])
    middle = np.abs(V[:, 0] - 0.5 * (xs[0] + xs[-1])) < d / 6.0
    R_tube = float(rho[middle].mean())
    # smallest ring radius: bin the middle third in x
    bins = np.linspace(-d / 6.0, d / 6.0, 13)
    xm = V[:, 0] - 0.5 * (xs[0] + xs[-1])
    ring = [rho[(xm >= lo) & (xm < hi)].mean() for lo, hi in zip(bins[:-1], bins[1:]) if ((xm >= lo) & (xm < hi)).any()]
    between = (V[:, 0] > xs[0]) & (V[:, 0] < xs[-1])
    narrow = between & (rho < 1.5 * R_tube)
    L_tube = float(V[narrow, 0].max() - V[narrow, 0].min()) if narrow.any() else float("nan")
    with open(os.path.join(run, "Output_data.txt")) as f:
        header = f.readline().split()
        rows = [[float(x) for x in l.split()] for l in f if l.strip()]
    col = {}
    for i, name in enumerate(header):
        col.setdefault(name, i)
    last = rows[-1]
    steps = last[col["step"]]
    return {"run": run, "KA": KA, "KB": KB, "d": d, "a": a, "R_tube": R_tube, "R_min": float(min(ring)), "R_theory": math.sqrt(KB / (4.0 * KA)),
            "spread": float(rho[middle].std() / R_tube), "L_tube": L_tube, "Area": last[col["Area"]],
            "E_bend": last[col.get("Bending", col.get("Bending_tan"))], "E_tension": last[col["Surface_tension"]], "E_total": last[col["Total_E"]],
            "steps": steps, "converged": steps < cfg["timesteps"] - 1}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("roots", nargs="+", help="results folders (a batch folder with run folders below it)")
    p.add_argument("--out", default="barbell_beads", help="prefix of the output files (default barbell_beads: .csv and .png in the current folder)")
    p.add_argument("--group", choices=["KB", "KA"], default="KB", help="series of the plot: one per KB (scan of KA, default) or one per KA (scan of KB)")
    args = p.parse_args()
    other = "KA" if args.group == "KB" else "KB"

    rows = []
    for root in args.roots:
        for run in find_runs(root):
            r = analyze_run(run)
            if r:
                rows.append(r)
    if not rows:
        raise SystemExit("no run folder found")
    rows.sort(key=lambda r: (r["KB"], r["KA"]))
    keys = list(rows[0].keys())
    with open(args.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print("%-10s %-10s %9s %9s %7s %8s %9s %9s  %s" % ("KA", "KB", "R_tube", "R_theory", "ratio", "spread", "L_tube", "E_total", "run"))
    for r in rows:
        print("%-10g %-10g %9.4f %9.4f %7.3f %8.3f %9.2f %9.3f  %s%s" % (r["KA"], r["KB"], r["R_tube"], r["R_theory"], r["R_tube"] / r["R_theory"], r["spread"],
                                                                          r["L_tube"], r["E_total"], r["run"], "" if r["converged"] else "  (hit the step limit)"))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not found: no plot")
        return
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for g in sorted({r[args.group] for r in rows}):
        sel = [r for r in rows if r[args.group] == g]
        ax[0].plot([r["R_theory"] for r in sel], [r["R_tube"] for r in sel], "o", label="%s = %g" % (args.group, g))
        ax[1].plot([r[other] for r in sel], [r["R_tube"] / r["R_theory"] for r in sel], "o", label="%s = %g" % (args.group, g))
    lim = [0.0, 1.05 * max(max(r["R_theory"], r["R_tube"]) for r in rows)]
    ax[0].plot(lim, lim, "k--", lw=1, label=r"$\sqrt{K_B/4K_A}$")
    ax[0].set_xlabel("R_theory = sqrt(KB / 4 KA)")
    ax[0].set_ylabel("R_tube (measured)")
    ax[0].legend()
    ax[1].axhline(1.0, color="k", lw=1, ls="--")
    ax[1].set_xscale("log")
    ax[1].set_xlabel(other)
    ax[1].set_ylabel("R_tube / R_theory")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(args.out + ".png", dpi=150)
    print("wrote %s.csv and %s.png" % (args.out, args.out))


if __name__ == "__main__":
    main()
