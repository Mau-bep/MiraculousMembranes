#!/usr/bin/env python3
"""Cross-sections of the final states of one-bead runs: the plane through the bead centre and the volume centre of the vesicle (the axis of the pocket).

usage:  python3 Slice_one_vesicle.py [--out slices.png] [label=]RUN_FOLDER [[label=]RUN_FOLDER ...]

RUN_FOLDER is a run folder (Final_state.obj or the last membrane_*.obj, Bead_0_data.txt, Input_file.json), e.g.
../Results/Wrapping_one_vesicle/transition_flat/link12_wt4.4/1 . The black line is the membrane, the red circle is the bead (radius of the config); the title gives
the value of the scanned parameter when the folder belongs to a chain (chain_summary.csv). For a flat disk (--planar) the plane contains the bead and the normal of the disk.
Look at these before trusting a curve: the pocket should be axisymmetric around the line from the vesicle centre to the bead, and a partially wrapped bead sits in a shallow dimple.
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from Analyze_two_vesicle import last_mesh, mesh_geometry, read_bead_positions, read_obj  # noqa: E402
from Slice_two_vesicle import section  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="[label=]run folder")
    ap.add_argument("--out", default="slices.png")
    ap.add_argument("--half", type=float, default=8.0, help="half width of each panel in bead radii")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(args.runs), figsize=(4.6 * len(args.runs), 4.8), squeeze=False)
    for ax, item in zip(axes[0], args.runs):
        run = item.rpartition("=")[2] if "=" in item else item
        name = os.path.basename(os.path.dirname(os.path.normpath(run)))
        label = item.rpartition("=")[0] if "=" in item else name
        title = label
        summary = os.path.join(os.path.dirname(os.path.dirname(os.path.normpath(run))), "chain_summary.csv")
        if os.path.exists(summary):
            for r in csv.DictReader(open(summary)):
                if r["run"] == name:
                    title = "%s\n%s = %g" % (label, r["scan_var"], float(r["scan_value"]))
        mesh, _ = last_mesh(run)
        V, F = read_obj(mesh)
        cfg = json.load(open(os.path.join(run, "Input_file.json")))
        a = cfg["Beads"][0]["radius"]
        pos = read_bead_positions(run, 1)[0]
        if cfg.get("boundary"):
            # planar disk: the plane through the bead normal to the disk (x is the normal), the origin is the foot of the bead
            vc = np.array([0.0, pos[1], pos[2]])
            u = np.array([0.0, 1.0, 0.0])
            w = np.array([1.0, 0.0, 0.0])
        else:
            _, _, vc = mesh_geometry(V, F)
            w = pos - vc
            w /= np.linalg.norm(w)
            ez = np.array([0.0, 0.0, 1.0])
            u = ez - (ez @ w) * w
            if np.linalg.norm(u) < 1e-6:
                u = np.array([1.0, 0.0, 0.0]) - w[0] * w
            u /= np.linalg.norm(u)
        n = np.cross(u, w)
        for (p0, p1) in section(V, F, vc, u, w, n):
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "k-", lw=0.7)
        q = pos - vc
        ax.add_patch(plt.Circle((q @ u, q @ w), a, fill=False, color="tab:red", lw=1.2))
        ax.set_title(title, fontsize=9)
        ax.set_aspect("equal")
        ax.set_xlim(-args.half, args.half)
        ax.set_ylim(-args.half, args.half)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
