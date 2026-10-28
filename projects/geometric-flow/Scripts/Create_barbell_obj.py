#!/usr/bin/env python3
"""Initial mesh of a barbell: two caps around two beads of radius a, joined by a tube of radius rt, a surface of revolution about the x axis.

usage:  python3 Create_barbell_obj.py --d 8 --radius 1 --tube 0.3 --out ../Config_files/x_init.obj

The beads sit at x = +-d/2; each cap is a sphere of radius a + gap around the bead centre (the membrane keeps the repulsive core of the bead off the
surface: gap 0.15 = a little more than the WCA distance 2^(1/6) - 1 = 0.12), the tube has radius rt, the caps are joined to it by a smooth neck (softplus of width
"neck"). The mesh is made of rings along the profile, each with 2 pi rho / edge vertices, zipped to the next one: edges close to --edge everywhere, a single
vertex at each tip, outward normals. It is only a starting mesh: the remesher of the run adapts it, and the relaxation finds the actual tube radius (sqrt(KB/(4 sigma)))
and the shape of the necks, so a rough --tube is enough (a tube too thick relaxes in a few thousand steps, a tube too thin is stretched).
"""
import argparse
import math
import sys

import numpy as np


def profile(d, a, gap, rt, neck, n=4000):
    """(x, rho) of the generating curve, x from -(d/2 + rb) to d/2 + rb."""
    rb = a + gap
    h = d / 2.0
    x = np.linspace(-h - rb, h + rb, n)

    def cap(xc):
        return np.sqrt(np.clip(rb ** 2 - (x - xc) ** 2, 0.0, None))

    def softplus(u):
        return neck * np.logaddexp(0.0, u / neck)

    rho = rt + softplus(cap(h) - rt) + softplus(cap(-h) - rt)
    # the tube stops where the caps end: the tips must close at rho = 0
    tip = np.clip((h + rb - np.abs(x)) / (0.5 * rb), 0.0, 1.0)
    rho = rho * np.sqrt(tip) if neck > 0 else rho
    rho[0] = rho[-1] = 0.0
    return x, rho


def barbell_mesh(d, a=1.0, gap=0.15, rt=0.3, edge=0.15, neck=0.05):
    """Vertices (N,3) and triangles (M,3) of the barbell, oriented outward."""
    x, rho = profile(d, a, gap, rt, neck)
    seg = np.hypot(np.diff(x), np.diff(rho))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    K = max(int(round(s[-1] / edge)), 8)
    sk = np.linspace(0.0, s[-1], K + 1)
    xk, rk = np.interp(sk, s, x), np.interp(sk, s, rho)
    V = [[xk[0], 0.0, 0.0]]
    rings = [[0]]
    for k in range(1, K):
        n = max(int(round(2.0 * math.pi * rk[k] / edge)), 6)
        ring = []
        for j in range(n):
            ph = 2.0 * math.pi * j / n
            V.append([xk[k], rk[k] * math.cos(ph), rk[k] * math.sin(ph)])
            ring.append(len(V) - 1)
        rings.append(ring)
    V.append([xk[K], 0.0, 0.0])
    rings.append([len(V) - 1])
    V = np.array(V)

    F = []
    for k in range(K):
        A, B = rings[k], rings[k + 1]
        if len(A) == 1:
            F += [[A[0], B[(j + 1) % len(B)], B[j]] for j in range(len(B))]
            continue
        if len(B) == 1:
            F += [[A[j], A[(j + 1) % len(A)], B[0]] for j in range(len(A))]
            continue
        # zipper: both rings start at angle 0; advance along the ring whose next vertex comes first
        na, nb = len(A), len(B)
        i = j = 0
        while i < na or j < nb:
            if i < na and (j >= nb or (i + 1) * nb <= (j + 1) * na):
                F.append([A[i], A[(i + 1) % na], B[j % nb]])
                i += 1
            else:
                F.append([A[i % na], B[(j + 1) % nb], B[j]])
                j += 1
    F = np.array(F, dtype=int)
    # orientation: outward normals (positive volume)
    a_, b_, c_ = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    if (a_ * np.cross(b_, c_)).sum() < 0:
        F = F[:, ::-1]
    return V, F


def check(V, F):
    e = {}
    for f in F:
        for u, v in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            e[(u, v)] = e.get((u, v), 0) + 1
    bad = [k for k, c in e.items() if c != 1 or (k[1], k[0]) not in e]
    chi = len(V) - len(e) // 2 + len(F)
    a_, b_, c_ = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    vol = (a_ * np.cross(b_, c_)).sum() / 6.0
    area = 0.5 * np.linalg.norm(np.cross(b_ - a_, c_ - a_), axis=1).sum()
    return len(bad), chi, vol, area


def write_obj(path, V, F):
    with open(path, "w") as f:
        f.write("# barbell mesh of revolution (Create_barbell_obj.py)\n")
        for v in V:
            f.write("v %.9g %.9g %.9g\n" % tuple(v))
        for t in F:
            f.write("f %d %d %d\n" % (t[0] + 1, t[1] + 1, t[2] + 1))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--d", type=float, required=True, help="distance between the bead centres")
    p.add_argument("--radius", type=float, default=1.0, help="bead radius a")
    p.add_argument("--gap", type=float, default=0.15, help="distance of the caps from the bead surface")
    p.add_argument("--tube", type=float, default=0.3, help="tube radius")
    p.add_argument("--edge", type=float, default=0.15, help="target edge length")
    p.add_argument("--neck", type=float, default=0.05, help="smoothing length of the junction of cap and tube")
    p.add_argument("--out", required=True)
    args = p.parse_args()
    V, F = barbell_mesh(args.d, args.radius, args.gap, args.tube, args.edge, args.neck)
    bad, chi, vol, area = check(V, F)
    print("%d vertices, %d triangles, non-manifold/boundary half-edges %d, Euler characteristic %d, volume %.3f, area %.3f" % (len(V), len(F), bad, chi, vol, area))
    if bad or chi != 2 or vol <= 0:
        sys.exit("the mesh is not a closed outward-oriented sphere")
    write_obj(args.out, V, F)
