#!/usr/bin/env python3
"""Shaded pictures of closed axisymmetric meshes (the shapes of Run_shape_chain.py), standing on their axis, without Blender: matplotlib only.

usage:  micromamba run -n mir_mem python3 Render_shape_obj.py --out pears.png  a.obj:"c0 2.4, v 0.716"  b.obj:"c0 3, v 0.775"  [--cut] [--color "#d62728"]

Every mesh is rotated so that its principal axis is vertical with the big lobe at the bottom (a pear stands on its big end), lit from the upper left, back faces left out. A title after the
colon is optional. --cut adds under every picture the outline of the meridian cut (what the snapshots on the curves of Analyze_shape_chains.py show).
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Run_shape_chain import meridian_outline, read_obj, shape_frame  # noqa: E402


def shaded_faces(V, F, azim_deg=-60.0, elev_deg=10.0):
    """Triangles (in the frame of the shape: z along the axis), their shades in [0, 1] and the camera direction."""
    c0, e, u, _ = shape_frame(V, F)
    P = V - c0
    X = np.stack([P @ u, P @ np.cross(e, u), P @ e], axis=1)
    tri = X[F]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300)
    az, el = math.radians(azim_deg), math.radians(elev_deg)
    cam = np.array([math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)])
    light = cam + np.array([-0.55, 0.45, 0.75])
    light /= np.linalg.norm(light)
    shade = 0.22 + 0.78 * np.clip(n @ light, 0.0, 1.0)
    front = n @ cam > 0.0
    return tri[front], shade[front], X


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--color", default="#d9534f")
    ap.add_argument("--cut", action="store_true", help="the outline of the meridian cut under every picture")
    ap.add_argument("--size", type=float, default=2.6, help="inches per picture")
    ap.add_argument("meshes", nargs="+", help="mesh.obj or mesh.obj:title")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import to_rgb
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    items = []
    for m in args.meshes:
        path, _, title = m.partition(":")
        V, F = read_obj(path)
        items.append((path, title, V, F))
    ext = max(np.ptp(shaded_faces(V, F)[2][:, 2]) for _, _, V, F in items)         # one scale for all, the longest shape sets it
    rows = 2 if args.cut else 1
    fig = plt.figure(figsize=(args.size * len(items), args.size * (1.25 if args.cut else 1.0) * 1.15 * rows / rows + (1.2 if args.cut else 0.0)))
    base = np.array(to_rgb(args.color))
    for k, (path, title, V, F) in enumerate(items):
        tri, shade, X = shaded_faces(V, F)
        ax = fig.add_subplot(rows, len(items), k + 1, projection="3d")
        pc = Poly3DCollection(tri, facecolors=np.clip(base[None, :] * shade[:, None] + 0.04, 0, 1), edgecolors="none", linewidths=0)
        ax.add_collection3d(pc)
        zc = 0.5 * (X[:, 2].max() + X[:, 2].min())
        half = 0.5 * ext * 1.02
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
        ax.set_zlim(zc - half, zc + half)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=10, azim=-60)
        ax.set_axis_off()
        ax.set_title(title or os.path.basename(os.path.dirname(os.path.dirname(path))), fontsize=9, pad=-6)
        if args.cut:
            ax2 = fig.add_subplot(rows, len(items), len(items) + k + 1)
            segs = meridian_outline(V, F)
            ax2.add_collection(LineCollection(segs, colors=args.color, linewidths=1.1))
            ax2.set_xlim(-0.5 * ext, 0.5 * ext)
            ax2.set_ylim(zc - 0.5 * ext, zc + 0.5 * ext)
            ax2.set_aspect("equal")
            ax2.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.93, wspace=0.0, hspace=0.0)
    fig.savefig(args.out, dpi=180)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
