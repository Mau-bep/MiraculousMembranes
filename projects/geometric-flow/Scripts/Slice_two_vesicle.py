#!/usr/bin/env python3
"""Cross-sections of the final states of two-bead vesicle runs: the plane through the two bead centres and the volume centre.

usage:  python3 Slice_two_vesicle.py [--out slices.png] [label=]RUN_FOLDER [[label=]RUN_FOLDER ...]

RUN_FOLDER is a run folder (Final_state.obj or the last membrane_*.obj, Bead_*_data.txt, Input_file.json), e.g.
../Results/Wrapping_two_vesicle/floppy_patch143/link05_s2.40/1 . The black line is the membrane, the red circles are the beads
(radius of the config); for a planar disk (--planar) the plane through the two bead centres normal to the disk (up is the normal of the disk).
Look at these before trusting a number: in a tensionless vesicle the beads sit at the end of deep necks, with a tension or the paper's area
budget the dimples are shallow, beads inside the vesicle sit in buds outside the main sphere.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Analyze_two_vesicle import last_mesh, mesh_geometry, read_bead_positions, read_obj  # noqa: E402


def section(V, F, vc, u, w, n):
    P = V - vc
    z = P @ n
    segments = []
    for face in F:
        zz = z[face]
        if zz.min() < 0 < zz.max():
            pts = []
            for i, j in ((0, 1), (1, 2), (2, 0)):
                if (zz[i] < 0) != (zz[j] < 0):
                    t = zz[i] / (zz[i] - zz[j])
                    q = P[face[i]] + t * (P[face[j]] - P[face[i]])
                    pts.append((q @ u, q @ w))
            if len(pts) == 2:
                segments.append(pts)
    return segments


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
        label, _, run = item.rpartition("=") if "=" in item else (os.path.basename(os.path.dirname(os.path.normpath(item))), "", item)
        mesh, _ = last_mesh(run)
        V, F = read_obj(mesh)
        cfg = json.load(open(os.path.join(run, "Input_file.json")))
        a = cfg["Beads"][0]["radius"]
        pos = read_bead_positions(run)
        if cfg.get("boundary"):
            # planar disk (--planar): the plane through the bead centres normal to the disk (x is the normal), the origin is the midpoint of the beads
            vc = 0.5 * (pos[0] + pos[1])
            u = pos[0] - pos[1]
            u[0] = 0.0
            u /= np.linalg.norm(u)
            w = np.array([1.0, 0.0, 0.0])
        else:
            _, _, vc = mesh_geometry(V, F)
            u = pos[0] - pos[1]
            u /= np.linalg.norm(u)
            w = vc - 0.5 * (pos[0] + pos[1])
            w -= u * (w @ u)
            w /= np.linalg.norm(w)
        n = np.cross(u, w)
        for (p0, p1) in section(V, F, vc, u, w, n):
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "k-", lw=0.7)
        for p in pos:
            q = p - vc
            ax.add_patch(plt.Circle((q @ u, q @ w), a, fill=False, color="tab:red", lw=1.2))
        d = np.linalg.norm(pos[0] - pos[1])
        ax.set_title("%s\nbead distance %.2f" % (label, d), fontsize=9)
        ax.set_aspect("equal")
        ax.set_xlim(-args.half, args.half)
        ax.set_ylim(-args.half, args.half)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
