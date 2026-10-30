"""Triangle mesh (OBJ) with the topology of a Jorge-Meeks k-noid: a sphere with k circular holes.

The k boundary circles are the mouths of k flared tubes.  By default the axes lie in the xy-plane, so
the circle centres sit on a ring; --layout sphere spreads them evenly over the sphere and --layout
random places them randomly (minimum separation).  Each circle lies in the plane perpendicular to
its radial direction.  It imitates the topology (genus 0, k boundary circles, k-fold symmetry), it is not
a minimal surface.

Construction: a unit sphere body with k circular holes (Delaunay triangulation of points on the
sphere = convex hull, hole facets removed), and one tube per hole built from rings of vertices
at equal 3D spacing (neighbouring rings with different vertex counts are stitched strip-wise).

Usage:  python Create_knoid_obj.py [--k 7] [--edge 0.08] [--length 1.3] [--radius R] [--layout ring|sphere|random] [--out ../Meshes/Knoid_7.obj]

--length is the tube length outside the unit sphere: 0 puts the holes on the sphere (circles of
radius sin(0.34) = 0.33), larger values push them further out along their radial directions.
"""
import argparse
import os
from collections import Counter

import numpy as np
from scipy.spatial import ConvexHull


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


def tube_profile(rb, theta_b, length, a0, a1, n=4000):
    """Meridian (s, rho) of the tube as a function of arc length; the tangent angle goes smoothly
    from a0 to a1 (angle to the axis), starting on the sphere at polar angle theta_b."""
    t = np.linspace(0.0, 1.0, n)
    alpha = a0 + (a1 - a0) * t * t * (3 - 2 * t)
    dl = length / (n - 1)
    s = rb * np.cos(theta_b) + np.concatenate([[0.0], np.cumsum(0.5 * (np.cos(alpha[1:]) + np.cos(alpha[:-1])) * dl)])
    rho = rb * np.sin(theta_b) + np.concatenate([[0.0], np.cumsum(0.5 * (np.sin(alpha[1:]) + np.sin(alpha[:-1])) * dl)])
    return t * length, s, rho


def stitch(A, angA, B, angB):
    """Triangles between two rings (vertex index lists in increasing azimuth, angles in [0, 2pi))."""
    na, nb = len(A), len(B)
    k0 = int(np.argmin(np.abs(np.angle(np.exp(1j * (angB - angA[0]))))))
    B = np.roll(B, -k0)
    angB = np.roll(angB, -k0)
    angB = np.unwrap(angB - angB[0]) + angB[0]
    angA_u = np.append(angA, angA[0] + 2 * np.pi)
    angB_u = np.append(angB, angB[0] + 2 * np.pi)
    # B[0] may sit slightly before A[0]; compare in a common unwrapped frame
    i = j = 0
    tris = []
    while i < na or j < nb:
        if i == na:
            adv_a = False
        elif j == nb:
            adv_a = True
        else:
            adv_a = angA_u[i + 1] <= angB_u[j + 1]
        if adv_a:
            tris.append((A[i % na], A[(i + 1) % na], B[j % nb]))
            i += 1
        else:
            tris.append((A[i % na], B[(j + 1) % nb], B[j % nb]))
            j += 1
    return tris


def orient(F, nv):
    """Make the face orientations consistent (BFS over shared edges)."""
    F = [list(f) for f in F]
    edge_faces = {}
    for fi, f in enumerate(F):
        for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            edge_faces.setdefault((min(a, b), max(a, b)), []).append(fi)
    done = np.zeros(len(F), dtype=bool)
    for start in range(len(F)):
        if done[start]:
            continue
        done[start] = True
        stack = [start]
        while stack:
            fi = stack.pop()
            f = F[fi]
            for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
                for fj in edge_faces[(min(a, b), max(a, b))]:
                    if fj == fi or done[fj]:
                        continue
                    g = F[fj]
                    # same directed edge a->b in the neighbour means it must be flipped
                    if any((g[q], g[(q + 1) % 3]) == (a, b) for q in range(3)):
                        g[1], g[2] = g[2], g[1]
                    done[fj] = True
                    stack.append(fj)
    return np.array(F)


def taubin(X, F, fixed, iters, lam=0.5, mu=-0.53):
    nbr = [set() for _ in range(len(X))]
    for f in F:
        for q in range(3):
            nbr[f[q]].add(f[(q + 1) % 3])
            nbr[f[q]].add(f[(q + 2) % 3])
    nbr = [np.array(sorted(s)) for s in nbr]
    X = X.copy()
    free = np.nonzero(~fixed)[0]
    for it in range(iters):
        for w in (lam, mu):
            D = np.array([X[nbr[v]].mean(axis=0) - X[v] for v in free])
            X[free] += w * D
    return X


def directions(k, layout, seed=0, min_sep=0.0):
    """k unit vectors (hole/tube axes): on the equator ('ring'), spread evenly over the sphere
    ('sphere', Fibonacci start relaxed by 1/r repulsion) or random with a minimum separation."""
    if layout == "ring":
        psi = 2 * np.pi * np.arange(k) / k
        return np.c_[np.cos(psi), np.sin(psi), np.zeros(k)]
    if layout == "sphere":
        p = fibonacci_sphere(k) if k > 2 else np.array([[0.0, 0.0, 1.0], [0.0, 0.1, -1.0]])[:k]
        p = p / np.linalg.norm(p, axis=1)[:, None]
        for _ in range(3000):
            diff = p[:, None, :] - p[None, :, :]
            r = np.linalg.norm(diff, axis=2) + np.eye(k)
            f = (diff / r[:, :, None] ** 3).sum(axis=1)
            f -= np.sum(f * p, axis=1)[:, None] * p
            p = p + 0.05 * f / (1.0 + np.linalg.norm(f, axis=1).max())
            p /= np.linalg.norm(p, axis=1)[:, None]
        return p
    if layout == "random":
        rng = np.random.default_rng(seed)
        sep = max(min_sep, 0.6 * np.sqrt(4 * np.pi / k))
        for _ in range(200):
            pts = []
            for _ in range(2000):
                v = rng.normal(size=3)
                v /= np.linalg.norm(v)
                if all(np.arccos(np.clip(v @ q, -1, 1)) > sep for q in pts):
                    pts.append(v)
                    if len(pts) == k:
                        return np.array(pts)
        raise ValueError("could not place %d random axes with separation %.2f rad; use fewer holes, a shorter --length or --layout sphere" % (k, sep))
    raise ValueError("unknown layout %r" % layout)


def build(k, edge, layout="ring", seed=0, rb=1.0, theta_b=0.34, length=1.3, a0=-0.30, a1=0.40, smooth=20, radius=None, force=False, auto_edge=True, min_edges=2.0):
    def problem(msg):
        if not force:
            raise ValueError(msg + "  (--force skips this check)")
        print("WARNING (ignored because of --force): " + msg)

    if radius is not None:
        if int(round(length / (0.866 * edge))) == 0:
            # no tube: the hole is the cap cut out of the sphere
            if not 0 < radius < rb:
                raise ValueError("--radius must be in (0, %.2f) when there is no tube" % rb)
            theta_b = np.arcsin(radius / rb)
        else:
            # mouth radius grows monotonically with the end angle a1 of the tube profile: bisect on it
            mouth = lambda a: tube_profile(rb, theta_b, length, a0, a)[2][-1]
            lo, hi = -0.6, 1.3
            if not mouth(lo) <= radius <= mouth(hi):
                problem("--radius %.3f not reachable with --length %.2f (reachable: %.3f to %.3f); change the length"
                        % (radius, length, max(mouth(lo), 0.0), mouth(hi)))
                radius = min(max(radius, mouth(lo)), mouth(hi))   # forced: use the closest reachable radius
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if mouth(mid) < radius else (lo, mid)
            a1 = 0.5 * (lo + hi)
    ell, s_prof, rho_prof = tube_profile(rb, theta_b, length, a0, a1)
    no_tube = length < 0.5 * 0.866 * edge   # holes are the base rings on the sphere
    s_chk = np.interp(np.linspace(0, length, 50), ell, s_prof)
    r_chk = np.interp(np.linspace(0, length, 50), ell, rho_prof)
    # axes closer than this and the holes (or the tubes) overlap
    need_hard = 2 * theta_b if no_tube else max(2 * theta_b, 2 * np.arcsin(min(1.0, np.max(r_chk / s_chk))))
    # random layouts keep a safety margin (independent of the automatic edge length below)
    need_random = max(2 * theta_b + 0.15, 2 * np.arcsin(min(1.0, np.max((2 * r_chk + edge) / (2 * s_chk)))))
    d = directions(k, layout, seed, need_random)
    gmin = np.arccos(np.clip(d @ d.T - 2 * np.eye(k), -1, 1)).min() if k > 1 else np.pi
    if gmin <= need_hard:
        problem("axes only %.2f rad apart but more than %.2f needed: the holes/tubes overlap; use fewer holes, a smaller --radius, "
                "a shorter --length or a different --layout" % (gmin, need_hard))
    elif auto_edge and k > 1:
        # edges between neighbouring boundaries: along the sphere between two holes, and (3D) between two tubes
        gap_surface = rb * (gmin - 2 * theta_b)
        gap_tubes = np.inf if no_tube else np.min(2 * s_chk * np.sin(gmin / 2) - 2 * r_chk)
        new_edge = min(edge, gap_surface / min_edges, gap_tubes)
        if new_edge < edge:
            print("edge length reduced from %.4f to %.4f so that at least %.0f edges fit between neighbouring boundaries "
                  "(%.3f between the holes on the sphere, %.3f between tubes)" % (edge, new_edge, min_edges, gap_surface, gap_tubes))
            edge = new_edge
    area = 4 * np.pi * rb ** 2 + k * 2 * np.pi * np.mean(r_chk) * length
    if area / (0.866 * edge ** 2) > 400000:
        problem("edge length %.4f would need about %.0f vertices" % (edge, area / (0.866 * edge ** 2)))
    if rho_prof.min() < 0.5 * edge * 3:
        problem("tube neck too thin for this edge length (min radius %.3f)" % rho_prof.min())
    # frame (nvec, mvec) perpendicular to each axis (for the ring layout: nvec = tangent, mvec = z)
    ref = np.where(np.abs(d[:, 2:3]) > 0.9, np.array([[1.0, 0.0, 0.0]]), np.array([[0.0, 0.0, 1.0]]))
    nvec = np.cross(ref, d)
    nvec /= np.linalg.norm(nvec, axis=1)[:, None]
    mvec = np.cross(d, nvec)
    # length < ~half a ring spacing: no tube, the holes are the base rings on the sphere
    n_rings = 0 if no_tube else max(1, int(round(length / (0.866 * edge))))
    ell_j = np.linspace(0, length, n_rings + 1)
    s_j = np.interp(ell_j, ell, s_prof)
    rho_j = np.interp(ell_j, ell, rho_prof)
    m_j = np.maximum(6, np.round(2 * np.pi * rho_j / edge).astype(int))

    def ring(kk, s, rho, m):
        a = 2 * np.pi * np.arange(m) / m
        return a, s * d[kk] + rho * (np.cos(a)[:, None] * nvec[kk] + np.sin(a)[:, None] * mvec[kk])

    verts = []

    def add(P):
        start = sum(len(v) for v in verts)
        verts.append(P)
        return np.arange(start, start + len(P))

    # base rings on the sphere (shared by the body and the tubes)
    base_idx, base_ang = [], []
    for kk in range(k):
        a = 2 * np.pi * np.arange(m_j[0]) / m_j[0]
        P = rb * (np.cos(theta_b) * d[kk] + np.sin(theta_b) * (np.cos(a)[:, None] * nvec[kk] + np.sin(a)[:, None] * mvec[kk]))
        base_idx.append(add(P))
        base_ang.append(a)
    nb_base = sum(len(i) for i in base_idx)

    # background points on the sphere, away from the holes
    n_bg = int(round(4 * np.pi * rb ** 2 / (0.866 * edge ** 2)))
    bg = fibonacci_sphere(n_bg)
    ang = np.arccos(np.clip(bg @ d.T, -1, 1))
    bg = bg[(ang > theta_b + 0.8 * edge / rb).all(axis=1)]
    bg_idx = add(rb * bg)

    # body = Delaunay on the sphere, hole facets removed
    V0 = np.concatenate(verts)
    hull = ConvexHull(V0 / rb)
    ringid = -np.ones(len(V0), dtype=int)
    for kk in range(k):
        ringid[base_idx[kk]] = kk
    body = []
    for t in hull.simplices:
        r = ringid[t]
        if r[0] >= 0 and r[0] == r[1] == r[2]:
            continue
        body.append(tuple(t))

    # tubes
    tris = list(body)
    for kk in range(k):
        prev_idx, prev_ang = base_idx[kk], base_ang[kk]
        for j in range(1, n_rings + 1):
            a, P = ring(kk, s_j[j], rho_j[j], m_j[j])
            idx = add(P)
            tris += stitch(prev_idx, prev_ang, idx, a)
            prev_idx, prev_ang = idx, a
    X = np.concatenate(verts)
    F = np.array(tris)
    F = orient(F, len(X))

    # boundary vertices (mouth rings) are fixed in the smoothing
    cnt = Counter()
    for f in F:
        for q in range(3):
            a, b = f[q], f[(q + 1) % 3]
            cnt[(min(a, b), max(a, b))] += 1
    fixed = np.zeros(len(X), dtype=bool)
    for (a, b), c in cnt.items():
        if c == 1:
            fixed[a] = fixed[b] = True
    X = taubin(X, F, fixed, smooth)

    # outward normals
    A = X[F]
    nrm = np.cross(A[:, 1] - A[:, 0], A[:, 2] - A[:, 0])
    if np.sum(np.einsum("ij,ij->i", nrm, A.mean(axis=1))) < 0:
        F = F[:, ::-1]
    return X, F


def report(X, F, k):
    cnt = Counter()
    for f in F:
        for q in range(3):
            a, b = f[q], f[(q + 1) % 3]
            cnt[(min(a, b), max(a, b))] += 1
    nonman = sum(1 for c in cnt.values() if c > 2)
    adj = {}
    for (a, b), c in cnt.items():
        if c == 1:
            adj.setdefault(a, []).append(b)
            adj.setdefault(b, []).append(a)
    seen, loops = set(), []
    for s in adj:
        if s in seen:
            continue
        stack, comp = [s], []
        while stack:
            v = stack.pop()
            if v in seen:
                continue
            seen.add(v)
            comp.append(v)
            stack.extend(adj[v])
        loops.append(comp)
    used = np.unique(F)
    chi = len(used) - len(cnt) + len(F)
    # directed-edge consistency = orientable and consistently oriented
    dirs = Counter((f[q], f[(q + 1) % 3]) for f in F for q in range(3))
    inconsistent = sum(1 for c in dirs.values() if c > 1)
    A = X[F]
    e = np.array([np.linalg.norm(A[:, 1] - A[:, 0], axis=1), np.linalg.norm(A[:, 2] - A[:, 1], axis=1), np.linalg.norm(A[:, 0] - A[:, 2], axis=1)])
    cosang = []
    for q in range(3):
        u = A[:, (q + 1) % 3] - A[:, q]
        v = A[:, (q + 2) % 3] - A[:, q]
        cosang.append(np.sum(u * v, axis=1) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)))
    ang = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
    print("vertices %d  faces %d  edges %d  unused vertices %d" % (len(used), len(F), len(cnt), len(X) - len(used)))
    print("Euler characteristic %d (expected %d)  non-manifold edges %d  inconsistently oriented edges %d" % (chi, 2 - k, nonman, inconsistent))
    print("boundary loops %d (expected %d), sizes %s" % (len(loops), k, sorted(len(c) for c in loops)))
    print("edge length min %.3f mean %.3f max %.3f   angles min %.1f max %.1f" % (e.min(), e.mean(), e.max(), ang.min(), ang.max()))
    for comp in loops:
        P = X[comp]
        c = P.mean(axis=0)
        # circle test: planar and equidistant from the centre
        u, sv, vt = np.linalg.svd(P - c)
        r = np.linalg.norm(P - c, axis=1)
        print("  loop: centre (%.3f, %.3f, %.3f)  distance from origin %.3f  circle radius %.3f (spread %.1e)  planarity %.1e"
              % (c[0], c[1], c[2], np.linalg.norm(c), r.mean(), r.max() - r.min(), sv[2]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--k", type=int, default=7, help="number of boundary circles")
    ap.add_argument("--edge", type=float, default=0.08,
                    help="target edge length (body is a unit sphere); reduced automatically if needed, see --min-edges")
    ap.add_argument("--min-edges", type=float, default=2.0,
                    help="at least this many edges between neighbouring boundaries; the edge length is shrunk to make it so")
    ap.add_argument("--no-auto-edge", action="store_true", help="keep --edge exactly, even if holes are too close for it")
    ap.add_argument("--layout", choices=["ring", "sphere", "random"], default="ring",
                    help="hole directions: on the equator, evenly spread over the sphere, or random (min. separation)")
    ap.add_argument("--seed", type=int, default=0, help="seed for --layout random")
    ap.add_argument("--length", type=float, default=1.3, help="tube length outside the sphere (0: holes lie on the sphere)")
    ap.add_argument("--radius", type=float, default=None,
                    help="radius of the boundary circles (default: ~0.40 for --length 1.3; for --length 0 it sets the hole size)")
    ap.add_argument("--force", action="store_true",
                    help="turn the geometry checks (hole/tube clearance, thin neck, unreachable --radius) into warnings and build anyway")
    ap.add_argument("--smooth", type=int, default=20, help="Taubin smoothing iterations (boundary fixed)")
    ap.add_argument("--out", default="../Meshes/Knoid_7.obj")
    args = ap.parse_args()

    X, F = build(args.k, args.edge, layout=args.layout, seed=args.seed, length=args.length, smooth=args.smooth, radius=args.radius, force=args.force,
                 auto_edge=not args.no_auto_edge, min_edges=args.min_edges)
    report(X, F, args.k)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        for p in X:
            fh.write("v %.9g %.9g %.9g\n" % tuple(p))
        for f in F:
            fh.write("f %d %d %d\n" % (f[0] + 1, f[1] + 1, f[2] + 1))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
