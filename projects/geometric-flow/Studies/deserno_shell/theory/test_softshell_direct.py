#!/usr/bin/env python
"""Tests of softshell_direct.py.   micromamba run -n mir_mem python test_softshell_direct.py [--slow]

fast checks (~1 min): derivatives, quadrature, analytic energies, KKT <-> fixed-w~ consistency, e'(G) = w~.
--slow adds: node convergence of the lines, R_dom convergence, s -> 0 against Deserno (deserno_theory.py).
Every check prints value, reference, tolerance and PASS/FAIL.
"""
import os
import sys
import time
import math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "Scripts"))
import softshell_direct as sd

SLOW = "--slow" in sys.argv
rows = []


def check(name, value, ref, tol, mode="abs"):
    if mode == "abs":
        ok = abs(value - ref) <= tol
    elif mode == "rel":
        ok = abs(value / ref - 1.0) <= tol
    elif mode == "lt":
        ok = value < ref
    else:
        ok = value > ref
    rows.append(ok)
    print("%-62s value %-14.8g ref %-14.8g tol %-9.3g %s" % (name, value, ref, tol, "PASS" if ok else "FAIL"), flush=True)
    return ok


# ---------------------------------------------------------------------------------------------------
def test_derivatives():
    print("\n[1] gradient and Hessian against finite differences")
    rng = np.random.default_rng(2)
    for p in (1.0, 4.0):
        m = sd.Model(sigma=0.5, s=0.3, p=p, R=4.0)
        n1 = 10
        d = 1.6 / n1
        psi = np.concatenate([(np.arange(n1) + 0.5) * d, np.linspace(1.2, 0.1, 6)])
        N = len(psi)
        ell = np.concatenate([np.full(n1, 2 * math.sin(d / 2)), np.full(6, 0.35)])
        ell /= ell.sum()
        st = sd.State(psi + 0.01 * rng.standard_normal(N), -1.0 + 0.01 * rng.standard_normal(), ell, m.R)
        wt = 5.0
        ev = sd.Eval(m, st)
        gr, Hr = ev.hessian(wt)
        z0 = st.z()
        fd = np.zeros_like(gr)
        Hfd = np.zeros_like(Hr)
        for k in range(len(z0)):
            dd = 1e-6
            zp, zm = z0.copy(), z0.copy()
            zp[k] += dd
            zm[k] -= dd
            fd[k] = (sd.energies(m, st.with_z(zp), wt)[0] - sd.energies(m, st.with_z(zm), wt)[0]) / (2 * dd)
            Hfd[:, k] = (sd.Eval(m, st.with_z(zp), hess=False).grad(wt) - sd.Eval(m, st.with_z(zm), hess=False).grad(wt)) / (2 * dd)
        check("p=%g gradient (max abs error / max |g|)" % p, np.abs(fd - gr).max() / np.abs(gr).max(), 0.0, 1e-7)
        check("p=%g Hessian (max abs error / max |H|)" % p, np.abs(Hfd - Hr).max() / np.abs(Hr).max(), 0.0, 1e-6)


def test_quadrature():
    print("\n[2] adhesion quadrature")
    from scipy.integrate import quad
    for p in (1.0, 4.0):
        m = sd.Model(sigma=0.5, s=0.25, p=p, R=10.0)
        for d in (0.7, 0.9, 1.1):
            st = sd.plate_state(m, d)
            Gn = sd.energies(m, st, 0.0)[3]
            f = lambda r: m.weight(np.array([r]))[0] * d / r ** 2
            Gq = quad(f, max(d, 1 - m.s), 1 + m.s, limit=400)[0]
            check("flat plate d=%.1f p=%g  G (numeric vs 1-d quad)" % (d, p), Gn, Gq, 2e-6)
    # sphere: polyline on the unit sphere, G = 1 - cos(alpha), E_bend = 4 (1 - cos alpha)
    m = sd.Model(sigma=0.5, s=0.25, p=1.0, R=10.0)
    al = 1.3
    for n in (40, 80):
        d = al / n
        psi = (np.arange(n) + 0.5) * d
        L = 2 * math.sin(d / 2)
        # attach a straight tail so that the closure rho_N = R holds: only the cap part is checked
        ell = np.full(n, L)
        stc = sd.State(psi, -1.0, ell / ell.sum(), float((L * np.cos(psi)).sum()))     # R chosen so that S = 1
        e, eb, ea, gm = sd.energies(m, stc, 0.0)
        check("unit-sphere cap n=%d: G vs 1-cos(alpha)" % n, gm, 1 - math.cos(al), 4e-3 * (40 / n) ** 2)
        # open end: the cap energy carries the (missing) half cell at the end: 4 z - rho_end*L/2*4
        ref = 4 * (1 - math.cos(al)) - 4 * math.sin(al) * L / 2
        check("unit-sphere cap n=%d: E_bend vs 4z - end cell" % n, eb, ref, 3e-3 * (40 / n) ** 2)


def test_branch_consistency():
    print("\n[3] KKT branch vs fixed-w~ minimisation, e'(G) = w~, analytic limits")
    m = sd.Model(sigma=0.53, s=0.25, p=1.0, R=10.0)
    kw = sd.GRIDS["mid"]
    st = sd.initial_state(m, 2.0, "mid")
    G0 = sd.energies(m, st, 0)[3]
    br = sd.trace_branch(m, st, 2.0, G0, 1.5, dG0=0.02, dG_max=0.03, grid_kw=kw, regrid_every=0)   # frozen grid: smooth
    a = br.arrays()
    # slope of e(G) from finite differences of the data equals the multiplier w~ (midpoints)
    Gm = 0.5 * (a["G"][1:] + a["G"][:-1])
    slope = np.diff(a["e"]) / np.diff(a["G"])
    wm = 0.5 * (a["wt"][1:] + a["wt"][:-1])
    # (the finite difference is O(dG^2): dG = 0.02-0.03 and w~(G) rises steeply at the start of the branch)
    check("de/dG (finite difference) vs multiplier w~, max relative deviation over branch", (np.abs(slope - wm) / wm).max(), 0.0, 6e-3)
    check("e(G_end) - e(G_0) vs trapezoid integral of w~(G) (relative)", (a["e"][-1] - a["e"][0]) / np.trapezoid(a["wt"], a["G"]), 1.0, 2e-3, mode="rel")
    # the constrained state is a stationary point of E - w~ G: fixed-w~ Newton returns it unchanged
    i = len(a["G"]) // 2
    st_i, w_i = br.states[i], br.wt[i]
    st2, info = sd.minimize_wt(m, st_i, w_i)
    check("fixed-w~ Newton from the KKT state: converged", float(info["converged"]), 1.0, 0.0)
    check("   ... coverage unchanged", sd.energies(m, st2, 0)[3], a["G"][i], 1e-7)
    check("   ... E = e - w~ G", sd.energies(m, st2, w_i)[0], a["e"][i] - w_i * a["G"][i], 1e-8)
    check("branch state is a constrained minimum (bordered-Hessian inertia 1)", sd.constraint_inertia(m, st_i, w_i)[0], 1, 0)
    # unbound reference: a flat plate has e = 0, G_flat(d) and the free minimum over d is the maximum of G
    # (E = -w~ G, no bending, no tension)
    m1 = sd.Model(sigma=0.53, s=0.25, p=1.0, R=10.0)
    stp = sd.plate_state(m1, 0.9)
    stp, info = sd.minimize_wt(m1, stp, 1e-3)
    check("plate at w~ = 1e-3 stays flat: bending energy", sd.energies(m1, stp, 1e-3)[1], 0.0, 1e-6)


def test_envelope_energy():
    print("\n[4] fully enveloped state: E vs -2(w~-4) + 4 sigma~ (exact for s -> 0)")
    for s in (0.05, 0.25):
        m = sd.Model(sigma=0.53, s=s, p=1.0, R=10.0)
        st = sd.initial_state(m, 2.0, "mid")
        G0 = sd.energies(m, st, 0)[3]
        br = sd.trace_branch(m, st, 2.0, G0, 1.995, dG0=0.03, dG_max=0.04, grid_kw=sd.GRIDS["mid"])
        a = br.arrays()
        i = int(np.argmin(np.abs(a["wt"] - 12.0) + 1e3 * (a["G"] < 1.8)))
        st_e, info = sd.minimize_wt(m, br.states[i], 12.0)
        e = sd.energies(m, st_e, 12.0)
        ref = -2 * (12.0 - 4) + 4 * 0.53
        print("   s=%.2f: G=%.4f  E=%.4f  (Ebend %.4f, tension %.4f)  ideal %.4f" % (s, e[3], e[0], e[1], e[2] * 0.53, ref))
        # ideal envelope has G = 2 (<W> = 1); the soft shell loses 2 w~ (1 - <W>) -> compare with that estimate
        check("s=%.2f enveloped E - ideal, in units of w~ (2 - G) (neck + weight loss)" % s, e[0] - ref, 0.0, 0.8)


def test_deserno_limit():
    print("\n[5] s -> 0 against Deserno (deserno_theory.py): shifts must vanish, ~ s^(1/2) for E and S2 (negative line tension)")
    import deserno_theory as dt
    sigma = 0.53
    ref = np.array([dt.spinodal_S1(sigma), dt.w_E(sigma), dt.spinodal_S2(sigma)])
    print("   Deserno S1 %.4f E %.4f S2 %.4f" % tuple(ref))
    res = {}
    for s in (0.05, 0.0125):
        out = sd.soft_lines(sigma, s=s, p=1.0, R=10.0, grid_scan="mid", grid_ref="mid")
        res[s] = np.array([out["S1"], out["E"], out["S2"]]) - ref
        print("   s=%.4f: shifts S1 %+.4f E %+.4f S2 %+.4f (%.0f s)" % (s, res[s][0], res[s][1], res[s][2], out["time"]))
    check("S1 shift at s=0.0125 (small, positive side)", res[0.0125][0], 0.0, 0.05)
    check("E shift at s=0.0125 below s=0.05 value by factor > 1.6", res[0.05][1] / res[0.0125][1], 1.6, 0.0, mode="gt")
    check("S2 shift at s=0.0125 below s=0.05 value by factor > 1.6", res[0.05][2] / res[0.0125][2], 1.6, 0.0, mode="gt")
    check("E shift exponent between s=0.05 and 0.0125 (expect ~0.5-0.7)", math.log(res[0.05][1] / res[0.0125][1]) / math.log(4.0), 0.6, 0.2)
    check("S2 shift exponent between s=0.05 and 0.0125 (expect ~0.5)", math.log(res[0.05][2] / res[0.0125][2]) / math.log(4.0), 0.5, 0.2)


def test_convergence():
    print("\n[6] convergence of the lines at sigma~ = 0.53, s = 0.25")
    base = sd.soft_lines(0.53, 0.25, 1.0, 10.0, grid_scan="mid", grid_ref="mid")
    fine = sd.soft_lines(0.53, 0.25, 1.0, 10.0, grid_scan="fine", grid_ref="fine")
    for k in ("S1", "S2", "E"):
        check("%s: mid vs fine grids" % k, base[k], fine[k], 5e-3)
    ng = sd.soft_lines(0.53, 0.25, 1.0, 10.0, grid_scan="mid", grid_ref="mid", ngl=10)
    for k in ("S1", "S2", "E"):
        check("%s: 6 vs 10 Gauss points" % k, base[k], ng[k], 2e-3)
    big = sd.soft_lines(0.53, 0.25, 1.0, 30.0, grid_scan="mid", grid_ref="mid")
    for k in ("S1", "S2", "E"):
        check("%s: R_dom = 10 vs 30" % k, base[k], big[k], 2e-3)


if __name__ == "__main__":
    t0 = time.time()
    test_derivatives()
    test_quadrature()
    test_branch_consistency()
    test_envelope_energy()
    if SLOW:
        test_deserno_limit()
        test_convergence()
    print("\n%d / %d checks passed (%.0f s)" % (sum(rows), len(rows), time.time() - t0))
    sys.exit(0 if all(rows) else 1)
