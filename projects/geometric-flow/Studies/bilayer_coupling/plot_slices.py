"""Cross-sections of vesicle meshes (OBJ) through the axis of symmetry, drawn like the shapes of Seifert et al. (Figs. 5, 7), with
the reduced volume v and the reduced area difference da of every mesh as the label.

    python plot_slices.py out.png mesh1.obj [mesh2.obj ...]

The axis of symmetry is the principal axis of the vertex covariance whose eigenvalue differs most from the other two; the
slice is the intersection of the mesh with the plane through the centroid that holds that axis and the longest other axis.
All the panels share one length scale."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def load(path):
    V, F = [], []
    for line in open(path):
        t = line.split()
        if t and t[0] == "v":
            V.append([float(x) for x in t[1:4]])
        elif t and t[0] == "f":
            F.append([int(x.split("/")[0]) - 1 for x in t[1:4]])
    return np.array(V), np.array(F)


def quantities(V, F):
    P = V[F]
    n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    a = 0.5 * np.linalg.norm(n, axis=1)
    n = n / (2 * a[:, None])
    A = a.sum()
    vol = np.einsum("ij,ij->", P[:, 0], np.cross(P[:, 1], P[:, 2])) / 6
    edges = {}
    for fi, f in enumerate(F):
        for k in range(3):
            i, j = f[k], f[(k + 1) % 3]
            edges.setdefault((min(i, j), max(i, j)), []).append((fi, i, j))
    M = 0.0
    for fs in edges.values():
        if len(fs) != 2:
            continue
        (f1, a1, b1), (f2, _, _) = fs
        e = V[b1] - V[a1]
        l = np.linalg.norm(e)
        e = e / l
        M += 0.5 * l * np.arctan2(e @ np.cross(n[f1], n[f2]), n[f1] @ n[f2])
    R0 = np.sqrt(A / (4 * np.pi))
    return vol / (4 * np.pi / 3 * R0 ** 3), M / (4 * np.pi * R0)


def slice_mesh(V, F):
    c = V.mean(0)
    w, U = np.linalg.eigh(np.cov((V - c).T))
    gaps = [abs(w[i] - w[(i + 1) % 3]) + abs(w[i] - w[(i + 2) % 3]) for i in range(3)]
    axis = U[:, int(np.argmax(gaps))]
    others = [U[:, i] for i in range(3) if i != int(np.argmax(gaps))]
    ext = [np.ptp((V - c) @ o) for o in others]
    u = others[int(np.argmax(ext))]
    nrm = np.cross(axis, u)
    d = (V - c) @ nrm
    segs = []
    for f in F:
        dd = d[f]
        pts = []
        for k in range(3):
            i, j = k, (k + 1) % 3
            if dd[i] * dd[j] < 0:
                s = dd[i] / (dd[i] - dd[j])
                pts.append(V[f[i]] + s * (V[f[j]] - V[f[i]]))
        if len(pts) == 2:
            segs.append([((p - c) @ u, (p - c) @ axis) for p in pts])
    return np.array(segs)


if __name__ == "__main__":
    out, files = sys.argv[1], sys.argv[2:]
    data = []
    for p in files:
        V, F = load(p)
        v, da = quantities(V, F)
        data.append((slice_mesh(V, F), v, da, p))
    R = max(np.abs(s).max() for s, *_ in data) * 1.05
    fig, axes = plt.subplots(1, len(data), figsize=(2.1 * len(data), 2.5), squeeze=False)
    for ax, (s, v, da, p) in zip(axes[0], data):
        ax.add_collection(LineCollection(s, colors="k", linewidths=1.2))
        ax.set_xlim(-R, R)
        ax.set_ylim(-R, R)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title("$\\Delta a$ = %.3f\nv = %.3f" % (da, v), fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
