#!/usr/bin/env python3
"""The continuum bending energy of an axisymmetric vesicle, by a 1D minimisation: the reference that the energies of the 3D meshes (Run_shape_chain.py) are compared with.

usage:  micromamba run -n mir_mem python3 Shape_axisym.py check --c0 2.4 --chain ../Results/Shape_chains/c0_2.4/prolate1 [--out table.csv]
        micromamba run -n mir_mem python3 Shape_axisym.py mesh  --c0 2.4 --mesh Final_state.obj
        micromamba run -n mir_mem python3 Shape_axisym.py saddle --c0 2.4 --chain ../Results/Shape_chains/c0_2.4/prolate1 --out fig11_saddles.csv

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


def symmetrise(psi):
    """The mean of a contour and its up-down mirror image (psi(L - s) = pi - psi(s)): a shape that equals its image has the up-down symmetry of the prolates."""
    return 0.5 * (psi + (math.pi - psi[::-1]))


def shape_numbers(psi, L):
    """(waist, largest radius, asymmetry): the smallest radius between a quarter and three quarters of the arclength, the largest radius, and max |psi - mirror image of psi|."""
    h, pm, c, s, X, Xm, dpsi = forward(psi, L)
    N = len(psi) - 1
    waist = float(Xm[N // 4:3 * N // 4].min())
    return waist, float(Xm.max()), float(np.abs(psi - (math.pi - psi[::-1])).max())


def stability(psi, L, nu, c0, eps=1e-6):
    """(lowest eigenvalue against asymmetric deformations, lowest eigenvalue against symmetric ones, how much the two sectors are coupled) of the second variation of the energy at a symmetric
    stationary shape, with the multipliers of the area, volume and closure constraints: the Hessian of the Lagrangian by central differences of the analytic gradient, in the two sectors of the up-down
    mirror (psi_k -> pi - psi_{N-k}); the symmetric one is restricted to the deformations that keep the constraints. A negative value means a saddle: the first is negative where the symmetric shape
    is unstable against becoming a pear (the paper's C1^pear is where it changes sign), the second where it is unstable within the symmetric shapes (beyond the minimal volume M^pro)."""
    from scipy.linalg import null_space
    N = len(psi) - 1
    x = np.concatenate([psi[1:-1], [L]])
    gE = objective(x, c0, N)[1]
    g, J = constraints(x, nu, N)
    lam = np.linalg.lstsq(J.T, gE, rcond=None)[0]

    def grad_lagrangian(y):
        return objective(y, c0, N)[1] - constraints(y, nu, N)[1].T @ lam

    H = np.empty((N, N))
    for j in range(N):
        h = eps * max(1.0, abs(x[j]))
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        H[:, j] = (grad_lagrangian(xp) - grad_lagrangian(xm)) / (2.0 * h)
    H = 0.5 * (H + H.T)
    # index i of x is the node k = i + 1 of psi, whose mirror image is the node N - k, index N - 2 - i; the last entry of x is L (it does not change under the mirror)
    pairs = [(i, N - 2 - i) for i in range(N // 2 - 1)]
    middle = N // 2 - 1
    BA = np.zeros((N, len(pairs) + 1))
    BS = np.zeros((N, len(pairs) + 1))
    for c, (i, p) in enumerate(pairs):
        BA[i, c] = BA[p, c] = BS[i, c] = 1.0 / math.sqrt(2.0)
        BS[p, c] = -1.0 / math.sqrt(2.0)
    BA[middle, -1] = 1.0
    BS[N - 1, -1] = 1.0
    lam_asym = float(np.linalg.eigvalsh(BA.T @ H @ BA).min())
    Z = null_space(J @ BS)
    lam_sym = float(np.linalg.eigvalsh(Z.T @ BS.T @ H @ BS @ Z).min())
    coupling = float(np.abs(BA.T @ H @ BS).max() / np.abs(H).max())
    return lam_asym, lam_sym, coupling


def saddle(c0, psi, L, nu, nu_to, step=0.0025, min_step=1.5e-4, waypoints=(), log=None):
    """The up-down symmetric branch of the shape (psi, L) at the reduced volume nu followed DOWN to nu_to, with the symmetry enforced at every step. Below the bifurcation C1^pear the symmetric prolate-1
    shapes are saddles (a minimiser, like the relaxation of the 3D meshes, leaves them for the pear), but the stationary shapes are still there and the constrained minimiser stays on them when the mirror
    symmetry is imposed. The step is halved (down to min_step) where the minimisation fails or the branch jumps (the energy or the waist changes by more than 20 or 25 percent in one step); the point that is
    still not reached is the end of the branch, the minimal volume M^pro where it turns back. Returns the rows (nu, F, length, waist, largest radius, the two lowest eigenvalues) and the reason it stopped."""
    F, psi, L, ok, res = solve(psi, L, nu, c0, maxiter=800)
    psi = symmetrise(psi)
    F, psi, L, ok, res = solve(psi, L, nu, c0, maxiter=800)
    waist, rmax, _ = shape_numbers(psi, L)
    la, ls, _ = stability(psi, L, nu, c0)
    rows = [(nu, F, L, waist, rmax, la, ls)]
    pending = sorted([w for w in waypoints if nu_to < w < nu - 1e-9], reverse=True)
    step_now, reason = step, "reached nu_to = %g" % nu_to
    while nu > nu_to + 1e-9:
        target = max(nu - step_now, nu_to)
        if pending and pending[0] >= target - 1e-12:
            target = pending[0]
        with np.errstate(all="ignore"):                       # a trial step that fails can overflow before it is rejected
            Fn, psin, Ln, ok, res = solve(psi, L, target, c0, maxiter=800)
            psin = symmetrise(psin)
            Fn, psin, Ln, ok, res = solve(psin, Ln, target, c0, maxiter=800)
            wn, rn, _ = shape_numbers(psin, Ln)
        good = np.isfinite(Fn) and 0.0 < Fn < 5.0 and res < 1e-7 and abs(Fn - F) < 0.20 * F and abs(wn - waist) < 0.25 * waist
        if not good:
            step_now = 0.5 * (nu - target)
            if step_now < min_step:
                reason = "the branch cannot be followed below nu = %.4f (minimal volume M^pro, or a jump to another branch)" % nu
                break
            continue
        la, ls, _ = stability(psin, Ln, target, c0)
        F, psi, L, nu, waist = Fn, psin, Ln, target, wn
        rows.append((nu, F, L, waist, rn, la, ls))
        if pending and abs(pending[0] - nu) < 1e-9:
            pending.pop(0)
        step_now = min(step, 1.5 * step_now)
        if log:
            log("nu %.4f  F %.5f  waist %.4f  lowest eigenvalue: asymmetric %+.3e, symmetric %+.3e" % (nu, F, waist, la, ls))
    return rows, reason


def chain_start(chain, c0):
    """(psi, L, nu, lowest nu of the chain, tag) for `saddle`: the lowest link of the chain at which the 1D solver converges from the profile of the mesh (the lowest links are the ones next to the
    bifurcation, where the profile is hardest to fit), and the lowest reduced volume that the chain reached."""
    rows = []
    for r in csv.DictReader(open(os.path.join(chain, "chain_summary.csv"))):
        if r["kind"] in ("seed", "point") and r["F"] != "" and r["status"] != "ended":
            if not os.path.isabs(r["mesh"]):
                r["mesh"] = os.path.join(chain, r["mesh"])
            rows.append(r)
    rows.sort(key=lambda r: float(r["nu"]))
    for r in rows:
        F, psi, L, ok, res, nu_mesh = from_mesh(r["mesh"], c0, float(r["nu"]))
        if ok and res < 1e-8 and F > 0.0:
            return psi, L, float(r["nu"]), float(rows[0]["nu"]), r["tag"]
    raise SystemExit("the 1D solver does not converge from any mesh of %s" % chain)


def run_saddle(args):
    psi, L, nu, nu_last, tag = chain_start(args.chain, args.c0)
    print("start: link %s at nu %.4f; the chain's lowest nu is %.4f" % (tag, nu, nu_last), flush=True)
    rows, reason = saddle(args.c0, psi, L, nu, args.nu_to, step=args.step, waypoints=[nu_last], log=lambda s: print(s, flush=True))
    print("stop:", reason)
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["nu", "F", "c0", "length", "waist", "rmax", "lowest_eigenvalue_asymmetric", "lowest_eigenvalue_symmetric", "note"])
            for k, (nu_k, F_k, L_k, waist_k, rmax_k, la_k, ls_k) in enumerate(rows):
                note = "start (profile of link %s)" % tag if k == 0 else ("the last link of the 3D chain" if abs(nu_k - nu_last) < 1e-9 else ("end: " + reason if k == len(rows) - 1 else ""))
                w.writerow(["%.6f" % nu_k, "%.6g" % F_k, "%g" % args.c0, "%.6g" % L_k, "%.6g" % waist_k, "%.6g" % rmax_k, "%.4g" % la_k, "%.4g" % ls_k, note])
        print("wrote", args.out)


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
    sd = sub.add_parser("saddle", help="the up-down symmetric part of a branch below the last link of its 3D chain, where the symmetric shapes are saddles")
    sd.add_argument("--c0", type=float, required=True)
    sd.add_argument("--chain", required=True, help="folder of the chain (chain_summary.csv) of the branch: prolate1")
    sd.add_argument("--out", default=None, help="csv of the points (nu, F, ...)")
    sd.add_argument("--nu_to", type=float, default=0.68, help="where to stop if the branch does not end before (default 0.68)")
    sd.add_argument("--step", type=float, default=0.0025, help="step in nu (halved where the minimisation fails)")
    args = ap.parse_args()
    if args.cmd == "test":
        print("relative error of the gradient of the energy %.2e, of the Jacobian of the constraints %.2e" % check_gradients())
    elif args.cmd == "check":
        check(args)
    elif args.cmd == "saddle":
        run_saddle(args)
    else:
        F, psi, L, ok, res, nu = from_mesh(args.mesh, args.c0)
        print("nu %.4f  F continuum %.6f  (L %.4f, converged %s, constraint residual %.1e)" % (nu, F, L, ok, res))


if __name__ == "__main__":
    main()
