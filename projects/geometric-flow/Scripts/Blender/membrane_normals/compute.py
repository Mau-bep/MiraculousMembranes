"""Per-vertex normals and mean curvature of a triangle mesh (numpy only).

The mean curvature follows the simulation's convention (core/src/geometry.cpp):
the "mean curvature normal" of a vertex is the area gradient

    g_i = 1/2 sum_j (cot a_ij + cot b_ij) (x_i - x_j)

which points the way the area grows (outward on a sphere, 2 H A n for a smooth
surface).  Dividing by twice the barycentric dual area A_i gives H, the mean
curvature (k1 + k2) / 2, so a sphere of radius R reads 1/R everywhere.
"""

import numpy as np


def signed_mean_curvature(V, T, N):
    """Mean curvature H at every vertex, signed along the vertex normals N.

    Positive where the area gradient points along N (outward on a sphere with
    outward normals).  Isolated vertices read 0.  Boundary vertices use the
    area of their half ring, which is consistent for the normal component.
    """
    n = len(V)
    a, b, c = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    two_area = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    ok = two_area > 1e-30
    safe = np.where(ok, two_area, 1.0)

    # cotangent of the angle at each corner (0 for degenerate triangles)
    cot = np.stack([
        np.einsum("ij,ij->i", b - a, c - a),
        np.einsum("ij,ij->i", a - b, c - b),
        np.einsum("ij,ij->i", a - c, b - c),
    ], axis=1) / safe[:, None]
    cot[~ok] = 0.0

    # the edge opposite corner k is (k+1, k+2); both ends get 1/2 cot (x_end - x_other)
    grad = np.zeros_like(V)
    for k in range(3):
        i, j = T[:, (k + 1) % 3], T[:, (k + 2) % 3]
        d = 0.5 * cot[:, k, None] * (V[i] - V[j])
        np.add.at(grad, i, d)
        np.add.at(grad, j, -d)

    area = np.bincount(T.ravel(), weights=np.repeat(two_area / 2.0 / 3.0, 3), minlength=n)
    h = np.zeros(n)
    live = area > 1e-30
    h[live] = np.einsum("ij,ij->i", grad[live], N[live]) / (2.0 * area[live])
    return h


def default_range(h, percentile=98.0):
    """(0, high) colour range from |H|: the high end ignores a few spiky vertices."""
    m = np.abs(h)
    if m.size == 0:
        return 0.0, 1.0
    high = float(np.percentile(m, percentile))
    return 0.0, high if high > 0 else max(float(m.max()), 1.0)


def mean_edge_length(V, T):
    e = np.concatenate([V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 1]], V[T[:, 0]] - V[T[:, 2]]])
    return float(np.linalg.norm(e, axis=1).mean()) if len(e) else 1.0
