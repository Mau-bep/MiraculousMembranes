#!/usr/bin/env python3
"""The continuum bending energy of an axisymmetric vesicle, by a 1D minimisation: the reference that the energies of the 3D meshes (Run_shape_chain.py) are compared with.

usage:  micromamba run -n mir_mem python3 Shape_axisym.py check --c0 2.4 --chain ../Results/Shape_chains/c0_2.4/prolate1 [--out table.csv]
        micromamba run -n mir_mem python3 Shape_axisym.py mesh  --c0 2.4 --mesh Final_state.obj

The functional of Seifert, Berndl and Lipowsky (Phys. Rev. A 44, 1182, eq. 3.2), in the units R0 = sqrt(A / 4 pi) = 1: the contour of the vesicle from the north pole to the south pole is
given by its tilt angle psi(s) against the arclength s (dX/ds = cos psi, dZ/ds = -sin psi, X the distance to the axis),

    F / (8 pi kappa) = (1/8) int_0^L X (psi' + sin(psi)/X - c0)^2 ds,   area  int X ds = 2,   volume  int X^2 sin(psi) ds = 4 nu / 3,   closure  X(L) = 0,

with psi(0) = 0 and psi(L) = pi. A sphere gives (1 - c0/2)^2. psi is given at N + 1 nodes of a uniform arclength grid (the length L is a variable too), the integrals are midpoint sums, and
SLSQP minimises F under the three constraints starting from the profile of a mesh, so it stays on the branch of that mesh: the result is the energy of the continuum shape that the mesh
approximates, at the reduced volume and spontaneous curvature of the mesh. Axisymmetric shapes without overhangs only (prolate, pear, dumbbell; not a stomatocyte).
"""
import argparse
import csv
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Run_shape_chain import read_obj, shape_frame, volume_area  # noqa: E402


def forward(psi, L):
    N = len(psi) - 1
    h = L / N
    pm = 0.5 * (psi[:-1] + psi[1:])
    c, s = np.cos(pm), np.sin(pm)
    X = np.concatenate([[0.0], np.cumsum(h * c)])
    Xm = 0.5 * (X[:-1] + X[1:])
    return h, pm, c, s, X, Xm, (psi[1:] - psi[:-1]) / h


def energy(psi, L, c0):
    h, pm, c, s, X, Xm, dpsi = forward(psi, L)
    return h * np.sum((Xm * dpsi + s - c0 * Xm) ** 2 / Xm) / 8.0


def profile_from_mesh(V, F, N=240, bins=160):
    """Initial guess (psi at N + 1 nodes, length L) in units R0 = 1 from the profile r(s) of a mesh (median radius of the vertices of every slab along the principal axis)."""
    vol, area = volume_area(V, F)
    c0, e, u, _ = shape_frame(V, F)
    r0 = math.sqrt(area / (4.0 * math.pi))
    s = (V - c0) @ e / r0
    r = np.linalg.norm((V - c0) - np.outer(s * r0, e), axis=1) / r0
    edges = np.linspace(s.min(), s.max(), bins + 1)
    idx = np.clip(np.digitize(s, edges) - 1, 0, bins - 1)
    zs, xs = [], []
    for j in range(bins):
        sel = idx == j
        if sel.sum() >= 2:
            zs.append(s[sel].mean())
            xs.append(np.median(r[sel]))
    zs, xs = np.array(zs), np.array(xs)
    # contour from the top (largest s) to the bottom, with the poles on the axis
    order = np.argsort(-zs)
    Z = np.concatenate([[s.max()], zs[order], [s.min()]])
    X = np.concatenate([[0.0], xs[order], [0.0]])
    ds = np.hypot(np.diff(X), np.diff(Z))
    arc = np.concatenate([[0.0], np.cumsum(ds)])
    L = arc[-1]
    t = np.linspace(0.0, L, N + 1)
    Xn, Zn = np.interp(t, arc, X), np.interp(t, arc, Z)
    # smooth a little (the medians are noisy), keep the poles
    k = np.ones(5) / 5.0
    Xs = np.convolve(np.pad(Xn, 2, mode="edge"), k, mode="valid")
    Zs = np.convolve(np.pad(Zn, 2, mode="edge"), k, mode="valid")
    Xs[0] = Xs[-1] = 0.0
    dX, dZ = np.diff(Xs), np.diff(Zs)
    psi_mid = np.arctan2(-dZ, dX)
    psi = np.concatenate([[0.0], 0.5 * (psi_mid[:-1] + psi_mid[1:]), [math.pi]])
    return psi, float(np.hypot(dX, dZ).sum()), vol, area


def _arrays(x, N):
    psi = np.concatenate([[0.0], x[:-1], [math.pi]])
    h = x[-1] / N
    pm = 0.5 * (psi[:-1] + psi[1:])
    c, s = np.cos(pm), np.sin(pm)
    Xm = h * (np.cumsum(c) - 0.5 * c)                       # midpoint of the segments: h (sum_{k<i} c_k + c_i / 2)
    return psi, h, c, s, Xm


def _suffix(g):
    """W_k = sum_{i>k} g_i + g_k / 2, the weight of the angle of segment k in sum_i g_i Xm_i."""
    return g.sum() - np.cumsum(g) + 0.5 * g


def objective(x, c0, N):
    """Energy and its gradient with respect to (psi_1 ... psi_{N-1}, L)."""
    psi, h, c, s, Xm = _arrays(x, N)
    d = (psi[1:] - psi[:-1]) / h
    A = Xm * (d - c0) + s
    G = A * A / Xm
    E = h * G.sum() / 8.0
    gX = 2.0 * A * (d - c0) / Xm - A * A / Xm ** 2
    gd, gs = 2.0 * A, 2.0 * A / Xm
    W = _suffix(gX)
    gpsi = (h / 8.0) * (-(h / 2.0) * (W[1:] * s[1:] + W[:-1] * s[:-1]) + (gd[:-1] - gd[1:]) / h + 0.5 * (gs[:-1] * c[:-1] + gs[1:] * c[1:]))
    dEdh = G.sum() / 8.0 + (h / 8.0) * (gX * Xm / h - gd * d / h).sum()
    return E, np.concatenate([gpsi, [dEdh / N]])


def constraints(x, nu, N):
    """The three constraints (area, volume, closure) and their Jacobian."""
    psi, h, c, s, Xm = _arrays(x, N)
    X_end = h * c.sum()
    g = np.array([h * Xm.sum() - 2.0, h * (Xm ** 2 * s).sum() - 4.0 * nu / 3.0, X_end])
    k = np.arange(N)
    W1 = N - k - 0.5                                          # _suffix of ones
    Wq = _suffix(2.0 * Xm * s)
    J = np.zeros((3, N))
    J[0, :N - 1] = h * (-(h / 2.0) * (W1[1:] * s[1:] + W1[:-1] * s[:-1]))
    J[1, :N - 1] = h * (-(h / 2.0) * (Wq[1:] * s[1:] + Wq[:-1] * s[:-1]) + 0.5 * (Xm[1:] ** 2 * c[1:] + Xm[:-1] ** 2 * c[:-1]))
    J[2, :N - 1] = -(h / 2.0) * (s[:-1] + s[1:])
    J[0, N - 1] = 2.0 * Xm.sum() / N
    J[1, N - 1] = 3.0 * (Xm ** 2 * s).sum() / N
    J[2, N - 1] = c.sum() / N
    return g, J


def solve(psi0, L0, nu, c0, maxiter=300):
    """Minimise the energy at the reduced volume nu and spontaneous curvature c0 starting from (psi0, L0). Returns (F, psi, L, success, residual)."""
    N = len(psi0) - 1
    x0 = np.concatenate([psi0[1:-1], [L0]])
    res = minimize(lambda x: objective(x, c0, N), x0, jac=True, method="SLSQP",
                   constraints=[{"type": "eq", "fun": lambda x: constraints(x, nu, N)[0], "jac": lambda x: constraints(x, nu, N)[1]}], options={"maxiter": maxiter, "ftol": 1e-13})
    psi = np.concatenate([[0.0], res.x[:-1], [math.pi]])
    return float(res.fun), psi, float(res.x[-1]), bool(res.success), float(np.abs(constraints(res.x, nu, N)[0]).max())


def check_gradients(N=60, c0=2.4, nu=0.8, seed=1):
    """Analytic gradients against central differences at a random point (the test of the formulas above)."""
    rng = np.random.default_rng(seed)
    psi = np.linspace(0.0, math.pi, N + 1)
    x = np.concatenate([psi[1:-1] + 0.05 * rng.standard_normal(N - 1), [3.0 + 0.1 * rng.random()]])
    E, g = objective(x, c0, N)
    gc, J = constraints(x, nu, N)
    eps = 1e-6
    gn, Jn = np.zeros_like(g), np.zeros_like(J)
    for i in range(len(x)):
        xp, xm = x.copy(), x.copy()
        xp[i] += eps
        xm[i] -= eps
        gn[i] = (objective(xp, c0, N)[0] - objective(xm, c0, N)[0]) / (2 * eps)
        Jn[:, i] = (constraints(xp, nu, N)[0] - constraints(xm, nu, N)[0]) / (2 * eps)
    return float(np.abs(g - gn).max() / np.abs(gn).max()), float(np.abs(J - Jn).max() / np.abs(Jn).max())


def from_mesh(path, c0, nu=None, N=240):
    """(F, psi, L, success, residual, nu of the mesh). A neck that is nearly closed (a radius of a hundredth of R0) breaks the minimiser from the profile of the mesh (X in the denominator of the
    energy): it is then started again from the same profile smoothed over 21 nodes, which opens the neck a little and lets it relax to the continuum value."""
    V, F = read_obj(path)
    psi0, L0, vol, area = profile_from_mesh(V, F, N)
    nu_mesh = 6.0 * math.sqrt(math.pi) * vol / area ** 1.5
    target = nu_mesh if nu is None else nu
    out = solve(psi0, L0, target, c0)
    if not (out[3] and out[4] < 1e-8 and out[0] > 0.0):
        smooth = np.convolve(np.pad(psi0, 10, mode="edge"), np.ones(21) / 21.0, mode="valid")
        smooth[0], smooth[-1] = 0.0, math.pi
        out = solve(smooth, L0, target, c0)
    return out + (nu_mesh,)


def check(args):
    rows = list(csv.DictReader(open(os.path.join(args.chain, "chain_summary.csv"))))
    out = []
    for r in rows:
        if r["kind"] not in ("seed", "point") or r["F"] == "":
            continue
        F_mesh, nu, c0e = float(r["F"]), float(r["nu"]), float(r["c0_eff"])
        F_cont, psi, L, ok, res = from_mesh(r["mesh"], c0e, nu)[:5]
        out.append((r["tag"], nu, c0e, F_mesh, F_cont, (F_mesh - F_cont) / F_cont * 100.0, ok, res))
        print("%-22s nu %.4f  F mesh %.5f  F continuum %.5f  diff %+.2f %%  %s (constraint residual %.1e)" % (r["tag"], nu, F_mesh, F_cont, out[-1][5], "ok" if ok else "NOT CONVERGED", res), flush=True)
    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["tag", "nu", "c0_eff", "F_mesh", "F_continuum", "diff_percent", "converged", "residual"])
            w.writerows(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("--c0", type=float, required=False)
    c.add_argument("--chain", required=True)
    c.add_argument("--out", default=None)
    sub.add_parser("test", help="analytic gradients against central differences")
    m = sub.add_parser("mesh")
    m.add_argument("--c0", type=float, required=True)
    m.add_argument("--mesh", required=True)
    args = ap.parse_args()
    if args.cmd == "test":
        print("relative error of the gradient of the energy %.2e, of the Jacobian of the constraints %.2e" % check_gradients())
    elif args.cmd == "check":
        check(args)
    else:
        F, psi, L, ok, res, nu = from_mesh(args.mesh, args.c0)
        print("nu %.4f  F continuum %.6f  (L %.4f, converged %s, constraint residual %.1e)" % (nu, F, L, ok, res))


if __name__ == "__main__":
    main()
