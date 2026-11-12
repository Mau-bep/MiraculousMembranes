"""Per-vertex mean, Gaussian and principal curvatures of a triangle mesh (numpy only).

The mean curvature follows the simulation's convention (core/src/geometry.cpp):
the "mean curvature normal" of a vertex is the area gradient

    g_i = 1/2 sum_j (cot a_ij + cot b_ij) (x_i - x_j)

which points the way the area grows (outward on a sphere, 2 H A n for a smooth
surface).  Dividing by twice the barycentric dual area A_i gives H, the mean
curvature (k1 + k2) / 2, so a sphere of radius R reads 1/R everywhere.

The Gaussian curvature is K = angle defect / barycentric dual area, as the simulation
(angleDefect in core/src/geometry.cpp) and geometry-central have it, 1/R^2 on a sphere.
The principal curvatures are H +/- sqrt(H^2 - K) (clamped at 0 under the root), the same
construction as the simulation's principalCurvatures and the library's.
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


def gaussian_curvature(V, T):
    """Gaussian curvature K at every vertex: angle defect over the barycentric dual area.

    The defect is 2 pi minus the sum of the triangle angles at the vertex, pi minus it on the boundary
    (a vertex on a boundary edge), with the area of the half ring, as for the mean curvature.
    Isolated vertices read 0.  Needs no normals, so it does not depend on the orientation.
    """
    n = len(V)
    a, b, c = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    cross = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    # angle at each corner, atan2 of |cross| and dot is exact for thin triangles too
    angles = np.stack([
        np.arctan2(np.linalg.norm(np.cross(b - a, c - a), axis=1), np.einsum("ij,ij->i", b - a, c - a)),
        np.arctan2(np.linalg.norm(np.cross(a - b, c - b), axis=1), np.einsum("ij,ij->i", a - b, c - b)),
        np.arctan2(np.linalg.norm(np.cross(a - c, b - c), axis=1), np.einsum("ij,ij->i", a - c, b - c)),
    ], axis=1)
    angle_sum = np.bincount(T.ravel(), weights=angles.ravel(), minlength=n)
    area = np.bincount(T.ravel(), weights=np.repeat(cross / 2.0 / 3.0, 3), minlength=n)

    # a vertex is on the boundary when one of its edges belongs to a single triangle
    edges = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
    keys, counts = np.unique(edges[:, 0] * n + edges[:, 1], return_counts=True)
    lonely = keys[counts == 1]
    boundary = np.zeros(n, dtype=bool)
    boundary[lonely // n] = True
    boundary[lonely % n] = True

    defect = np.where(boundary, np.pi, 2.0 * np.pi) - angle_sum
    k = np.zeros(n)
    live = area > 1e-30
    k[live] = defect[live] / area[live]
    return k


def principal_curvatures(h, k):
    """(k_min, k_max) from the mean curvature h and the Gaussian curvature k: h -/+ sqrt(h^2 - k).

    The discrete h and k do not always agree (h^2 < k at a vertex that is nearly round): the root is
    then 0 and both principal curvatures read h.
    """
    root = np.sqrt(np.maximum(h * h - k, 0.0))
    return h - root, h + root


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
