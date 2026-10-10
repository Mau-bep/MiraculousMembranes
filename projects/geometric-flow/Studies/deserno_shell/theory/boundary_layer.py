#!/usr/bin/env python
"""
boundary_layer.py -- workstream F1: contact-line boundary layer of a finite-range adhesion shell, and what it does to
Deserno's phase lines.  Everything quoted in line_tension_analysis.md is produced here.

    python boundary_layer.py            # all sections (about 6 min on 2 cores), writes figures/ and data/
    python boundary_layer.py bl         # flat boundary layer: Dc, first integral, T_p, virial identity, scaling
    python boundary_layer.py lines      # first-order line-tension model: lines vs s, p; comparison with the old simulations
    python boundary_layer.py disc       # discretisation error of the boundary layer on a regular 1D mesh
    python boundary_layer.py tilt       # tilt-weighted (Y^q) boundary layer, resolution trade-off
    python boundary_layer.py ref        # planar reference state (plane touching the bead)

UNITS.  a = 1, kappa = KB/2, energies in pi*kappa, w~ = 4 KI a^2/KB, sigma~ = 2 KA a^2/KB, s = shell half width / a.
Boundary layer (flat substrate, small slopes): height d(x) of the membrane above the bead surface,
        e = int [ (kappa/2) d''^2 + U(d, d') ] dx,       U = -w W(d/delta) T(d'),      W(x) = [(1+cos(pi x))/2]^p,
T = 1 (code as it is) or T = (1+d'^2)^(-q/2) (tilt weight Y^q).  Boundary conditions: d -> 0 for x -> -inf (bound, flat),
d -> (Dc/2)(x-xc)^2 + h for x -> +inf (outer free membrane, curvature jump Dc).  The Euler-Lagrange equation is the 4th
order ODE  kappa d'''' - (U_d')' + U_d = 0;  the Dc that makes this problem solvable is the unknown parameter of the BVP
(translation invariance is fixed by pinning d'(L2)).
"""
import os
import sys
import json
import time
import numpy as np

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp, simpson, quad
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from scipy.optimize import minimize
from scipy.signal import argrelextrema

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.normpath(os.path.join(HERE, ".."))
FIG = os.path.join(STUDY, "figures")
DATA = os.path.join(STUDY, "data")
CACHE = os.path.join(HERE, "boundary_layer_cache.json")


# =====================================================================================================
# 0. the shell and a tiny cache for the slow parts
# =====================================================================================================
def W(x, p):
    x = np.asarray(x, float)
    c = np.where(np.abs(x) < 1, np.cos(0.5 * np.pi * x), 0.0)
    return c ** (2 * p)


def Wp(x, p):
    x = np.asarray(x, float)
    ins = np.abs(x) < 1
    c = np.where(ins, np.cos(0.5 * np.pi * x), 0.0)
    s = np.where(ins, np.sin(0.5 * np.pi * x), 0.0)
    return -p * np.pi * c ** (2 * p - 1) * s


def _cache_get(key):
    try:
        return json.load(open(CACHE)).get(key)
    except Exception:
        return None


def _cache_put(key, val):
    try:
        d = json.load(open(CACHE))
    except Exception:
        d = {}
    d[key] = val
    json.dump(d, open(CACHE, "w"))


# =====================================================================================================
# 1. flat boundary layer, scaled so that there is no free parameter
# =====================================================================================================
def _left_projector(k4):
    """Two linear conditions that select the solutions of y' = M y (d'''' = -k4 d) decaying for x -> -inf."""
    M = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [-k4, 0, 0, 0]], float)
    ev, V = np.linalg.eig(M.T)
    rows = []
    for i in range(4):
        if ev[i].real < 0:
            rows += [V[:, i].real, V[:, i].imag]
    return np.linalg.svd(np.array(rows))[2][:2]


def solve_scaled(p, L1=30.0, L2=8.0, Q=6.0, n=3000, P0=1.0):
    """Scaled flat BL: d = delta*D, x = sqrt(delta/Dc)*t, kappa = 1, w = kappa Dc^2/2:  D'''' = W'(D)/2, D''(inf) = P.
    P is returned by the solver (the contact condition says P = 1).  Q pins D'(L2)."""
    Lm = _left_projector(p * np.pi ** 2 / 4)

    def f(t, y, P):
        return np.vstack([y[1], y[2], y[3], 0.5 * Wp(y[0], p)])

    def bc(ya, yb, P):
        return np.array([Lm[0] @ ya, Lm[1] @ ya, yb[2] - P[0], yb[3], yb[1] - Q])

    t = np.linspace(-L1, L2, n)
    sp = np.logaddexp(0, (t - (L2 - Q)) * 2) / 2
    d = 0.5 * sp ** 2
    y0 = np.vstack([d, np.gradient(d, t), np.gradient(np.gradient(d, t), t), np.zeros_like(t)])
    return solve_bvp(f, bc, t, y0, p=[P0], tol=1e-8, max_nodes=400000), Q


def analyse_scaled(p, sol, Q):
    t = np.linspace(sol.x[0], sol.x[-1], 400001)
    D, u, v, w3 = sol.sol(t)
    P = sol.p[0]
    C = u * w3 - 0.5 * v ** 2 - 0.5 * W(D, p)                    # first integral (kappa = 1, w = 1/2)
    tc = sol.x[-1] - Q / P                                       # vertex of the outer parabola
    h = D[-1] - 0.5 * P * (sol.x[-1] - tc) ** 2                  # vertex height ("lift-off offset") in units of delta
    ref = np.where(t < tc, -0.5, 0.5 * P ** 2)
    e = 0.5 * v ** 2 - 0.5 * W(D, p)
    Iexc = simpson(e - ref, x=t)
    Tp = -2.0 * Iexc                                             # tau = -Tp * w * sqrt(delta*ell_c)
    # virial identity (exact, from tau ~ sqrt(delta) and Hellmann-Feynman):  tau/2 = int delta dU/d delta dx = w int x W'(x) dx
    # in scaled units this is  T_p = -2 int D W'(D) dt
    Tvir = -2.0 * simpson(D * Wp(D, p), x=t)
    lw = simpson(W(D, p)[t > tc], x=t[t > tc])                   # int W ds over the lift-off side (units sqrt(delta ell_c))
    return dict(P=P, Cmean=C.mean(), Cspread=C.max() - C.min(), tc=tc, h=h, T=Tp, Tvir=Tvir, lw=lw, Dmin=D.min(),
                t=t, D=D, u=u, v=v)


def trial_bound(p):
    """tau <= trial energy of the sharp profile D = t^2/2 inside the soft shell: T_p >= int_0^sqrt2 W(t^2/2) dt."""
    return quad(lambda t: float(W(t * t / 2, p)), 0, np.sqrt(2), limit=200)[0]


def section_bl(show=True):
    out = {}
    for p in (1, 2, 4, 8):
        sol, Q = solve_scaled(p)
        r = analyse_scaled(p, sol, Q)
        # initial-guess independence of Dc
        s2, _ = solve_scaled(p, P0=0.5, n=2000)
        r["P_alt"] = s2.p[0]
        out[p] = r
        if show:
            print("p=%d status %d  Dc=%.8f (start 0.5: %.8f)  first-integral spread %.1e  T_p=%.4f  (trial bound %.4f)  "
                  "h/delta=%.4f  Dmin/delta=%.4f  lw=%.4f  T/lw=%.4f"
                  % (p, sol.status, r["P"], r["P_alt"], r["Cspread"], r["T"], trial_bound(p), r["h"], r["Dmin"], r["lw"],
                     r["T"] / r["lw"]))
    return out


# =====================================================================================================
# 2. boundary layer with tilt weight (physical lengths in units of ell_c = 1/Dc, kappa = 1, w = 1/2)
# =====================================================================================================
def _T(u, q):
    return (1 + u * u) ** (-q / 2)


def _T1(u, q):
    return -q * u * (1 + u * u) ** (-q / 2 - 1)


def _T2(u, q):
    return -q * (1 + u * u) ** (-q / 2 - 1) + q * (q + 2) * u * u * (1 + u * u) ** (-q / 2 - 2)


def solve_tilt(p, dh, q, L2=8.0, Q=6.0, n=3000):
    """dh = delta/ell_c = s sqrt(w~).  Weight W(d/dh) (1+d'^2)^(-q/2).  Returns the solution; q = 0 reproduces solve_scaled."""
    k4 = p * np.pi ** 2 / (4 * dh ** 2)
    L1 = 14 / (k4 ** 0.25 / np.sqrt(2))
    Lm = _left_projector(k4)

    def f(t, y, P):
        d, u, v, w3 = y
        x = d / dh
        Wv, Wd = W(x, p), Wp(x, p)
        rhs = (Wd / (2 * dh)) * _T(u, q) - (Wd / (2 * dh)) * _T1(u, q) * u - 0.5 * Wv * _T2(u, q) * v
        return np.vstack([u, v, w3, rhs])

    def bc(ya, yb, P):
        return np.array([Lm[0] @ ya, Lm[1] @ ya, yb[2] - P[0], yb[3], yb[1] - Q])

    t = np.linspace(-L1, L2, n)
    sp = np.logaddexp(0, (t - (L2 - Q)) * 3) / 3
    d = 0.5 * sp ** 2
    y0 = np.vstack([d, np.gradient(d, t), np.gradient(np.gradient(d, t), t), np.zeros_like(t)])
    return solve_bvp(f, bc, t, y0, p=[1.0], tol=1e-8, max_nodes=300000), Q


def analyse_tilt(sol, p, dh, q, Q=6.0):
    t = np.linspace(sol.x[0], sol.x[-1], 400001)
    d, u, v, w3 = sol.sol(t)
    P = sol.p[0]
    U = -0.5 * W(d / dh, p) * _T(u, q)
    Up = -0.5 * W(d / dh, p) * _T1(u, q)
    C = u * w3 - u * Up - 0.5 * v ** 2 + U
    tc = sol.x[-1] - Q / P
    ref = np.where(t < tc, -0.5, 0.5 * P ** 2)
    tau = simpson(0.5 * v ** 2 + U - ref, x=t)
    wt = W(d / dh, p) * _T(u, q)
    lw = simpson(wt[t > tc], x=t[t > tc])
    k4 = p * np.pi ** 2 / (4 * dh ** 2)
    return dict(P=P, Cspread=C.max() - C.min(), Theta=-2 * tau, lw=lw, lb=np.sqrt(2) / k4 ** 0.25, dmin=d.min())


def section_tilt(show=True, force=False):
    res = None if force else _cache_get("tilt_scan")
    if res is None:
        res = []
        for p in (1, 2, 4, 8):
            for dh in (0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.1):
                for q in (0, 4, 8, 16):
                    sol, Q = solve_tilt(p, dh, q)
                    if sol.status != 0:
                        sol, Q = solve_tilt(p, dh, q, n=6000)
                    if sol.status != 0:
                        continue
                    r = analyse_tilt(sol, p, dh, q, Q)
                    res.append(dict(p=p, dh=dh, q=q, Theta=r["Theta"], lw=r["lw"], lb=r["lb"], P=r["P"]))
        _cache_put("tilt_scan", res)
    if show:
        for q in (0, 4, 8, 16):
            x = np.array([a["Theta"] / a["lw"] for a in res if a["q"] == q])
            print("q=%2d  Theta/ell_w over %d cases: min %.3f  max %.3f  mean %.3f" % (q, len(x), x.min(), x.max(), x.mean()))
    return res


# =====================================================================================================
# 3. first-order line-tension model on top of Deserno's solution
# =====================================================================================================
TP = {1: 1.1982, 2: 1.0289, 4: 0.8746, 8: 0.7395}          # filled/validated by section_bl (values printed there)
_CACHE_NPZ = os.path.join(DATA, "theory_cache.npz")
_C = None


def _tc():
    global _C
    if _C is None:
        _C = np.load(_CACHE_NPZ)
    return _C


class Deserno:
    """F(z), psi0dot(z) of Deserno's problem at sigma~ = 2k/15 (data/theory_cache.npz), spline-interpolated."""

    def __init__(self, k, Nf=8001):
        C = _tc()
        self.sig = float(C["s%02d_sigma" % k])
        z, F, pd = C["z"], C["s%02d_F" % k], C["s%02d_pd" % k]
        ok = np.isfinite(F) & np.isfinite(pd)
        z, F, pd = z[ok], F[ok], pd[ok]
        wD = (1 - pd) ** 2
        self.Fs = CubicHermiteSpline(z, F, wD - 4 - 2 * self.sig * z)      # exact slope identity
        self.wDs = CubicSpline(z, wD)
        self.zg = np.linspace(z[0], z[-1], Nf)
        self.F, self.wD = self.Fs(self.zg), self.wDs(self.zg)

    def g(self, Tp, s):
        """line-tension term, E~ units:  g = 2 a sin(alpha) tau/kappa,  tau = -Tp * mu * sqrt(delta/Dc),  mu = w_D(z)."""
        z = self.zg
        return -Tp * np.sqrt(s) * np.abs(self.wD) ** 0.75 * np.sqrt(np.clip(z * (2 - z), 0, None))

    def E(self, w, Tp=0.0, s=0.25):
        z = self.zg
        return -(w - 4) * z + self.sig * z ** 2 + self.F + (self.g(Tp, s) if Tp else 0.0)

    def wsoft(self, Tp, s):
        return self.wD + (np.gradient(self.g(Tp, s), self.zg) if Tp else 0.0)


def lines(th, Tp, s, wlo=2.0, whi=14.0, nw=80):
    """S1, S2 (max / following min of w_soft(z)), E line and barrier of E_D + g.  Tp = 0: Deserno itself (z = 2 endpoint)."""
    z = th.zg
    ws = th.wsoft(Tp, s)
    mx = argrelextrema(ws, np.greater, order=20)[0]
    mn = argrelextrema(ws, np.less, order=20)[0]
    out = {}
    if len(mx):
        out["S1"], out["zS1"] = ws[mx[0]], z[mx[0]]
        j = mn[mn > mx[0]]
        if len(j):
            out["S2"], out["zS2"] = ws[j[0]], z[j[0]]

    def mins(w):
        E = th.E(w, Tp, s)
        m = argrelextrema(E, np.less, order=20)[0]
        M = argrelextrema(E, np.greater, order=20)[0]
        if not Tp:
            m = np.append(m, len(z) - 1)
        return E, m, M

    def diff(w):
        E, m, M = mins(w)
        if len(m) < 2:
            return np.nan
        return E[m[0]] - E[m[-1]]

    grid = np.linspace(wlo, whi, nw)
    d = np.array([diff(w) for w in grid])
    for i in range(nw - 1):
        if np.isfinite(d[i]) and np.isfinite(d[i + 1]) and d[i] * d[i + 1] < 0:
            a, b = grid[i], grid[i + 1]
            fa = d[i]
            for _ in range(45):
                c = 0.5 * (a + b)
                dc = diff(c)
                if not np.isfinite(dc):
                    break
                if dc * fa <= 0:
                    b = c
                else:
                    a, fa = c, dc
            out["E"] = 0.5 * (a + b)
            E, m, M = mins(out["E"])
            out["zp"], out["ze"] = z[m[0]], z[m[-1]]
            kk = M[(M > m[0]) & (M < m[-1])]
            out["barrier"] = (E[kk[0]] - E[m[-1]]) if len(kk) else 0.0
            if len(kk):
                out["zb"] = z[kk[0]]
            break
    return out


def sim_lines():
    import csv
    f = os.path.join(DATA, "baseline_lines.csv")
    return {round(float(r["sigma_t"]), 3): r for r in csv.DictReader(open(f))}


def _fmt(d, k):
    return ("%.2f" % d[k]) if k in d else "  - "


def section_lines(show=True):
    sim = sim_lines()
    SIGS = (2, 4, 8)          # sigma~ = 0.267, 0.533, 1.067 (and k=1: 0.133)
    rows = []
    if show:
        print("first-order line tension (p=1, s=0.25) vs old simulation (shell 0.25) vs Deserno")
        print(" sigma~ | Deserno E S1 S2 bar | model E S1 S2 bar | sim E S1 S2")
    for k in range(1, 16):
        th = Deserno(k)
        b = lines(th, 0.0, 0.0)
        r = lines(th, TP[1], 0.25)
        sm = sim.get(round(th.sig, 3))
        if show:
            print(" %.3f | %s %s %s %s | %s %s %s %s | %s %s %s" % (
                th.sig, _fmt(b, "E"), _fmt(b, "S1"), _fmt(b, "S2"), _fmt(b, "barrier"), _fmt(r, "E"), _fmt(r, "S1"),
                _fmt(r, "S2"), _fmt(r, "barrier"), sm["E_sim"] if sm else "", sm["S1_sim"] if sm else "",
                sm["S2_sim"] if sm else ""))
    # scan in s and p at four sigma
    for k in (1, 2, 4, 8):
        th = Deserno(k)
        b = lines(th, 0.0, 0.0)
        for p in (1, 2, 4, 8):
            for s in (0.25, 0.15, 0.1, 0.05, 0.025):
                r = lines(th, TP[p], s)
                row = dict(sigma=th.sig, p=p, s=s, ell=np.sqrt(s) * 5.0 ** -0.25)
                for key in ("E", "S1", "S2", "barrier"):
                    row[key] = r.get(key, np.nan)
                    row[key + "_D"] = b.get(key, np.nan)
                    row["d" + key] = row[key] - row[key + "_D"]
                rows.append(row)
    keys = list(rows[0].keys())
    with open(os.path.join(DATA, "line_tension_shifts.csv"), "w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join("%.5g" % r[k] for k in keys) + "\n")
    return rows


def section_extrapolation(show=True):
    """Recover Deserno's lines from soft-shell results at several s (first-order model as the 'experiment')."""
    out = {}
    for k in (4, 8):
        th = Deserno(k)
        b = lines(th, 0.0, 0.0)
        ss = np.array([0.25, 0.16, 0.09, 0.04, 0.0225])
        data = {key: [] for key in ("E", "S1", "S2", "barrier")}
        for s in ss:
            r = lines(th, TP[1], s)
            for key in data:
                data[key].append(r.get(key, np.nan))
        x = np.sqrt(ss)
        res = {}
        for key in data:
            y = np.array(data[key])
            row = {}
            for deg, sel in ((1, slice(2, 5)), (2, slice(1, 5)), (2, slice(0, 5))):
                if np.all(np.isfinite(y[sel])):
                    c = np.polyfit(x[sel], y[sel], deg)
                    row["deg%d_pts%s" % (deg, sel.start)] = float(np.polyval(c, 0.0))
            row["Deserno"] = b.get(key, np.nan)
            row["values"] = [float(v) for v in y]
            res[key] = row
        out[th.sig] = res
        if show:
            print("sigma~ = %.3f  (s = %s)" % (th.sig, ss))
            for key, row in res.items():
                print("   %-8s Deserno %.3f | extrapolations: %s | raw %s" % (
                    key, row["Deserno"], {a: round(v, 3) for a, v in row.items() if a.startswith("deg")},
                    np.round(row["values"], 3)))
    return out


# =====================================================================================================
# 4. planar reference state: plane touching the bead
# =====================================================================================================
def gamma_ref(s, p):
    """covered (weighted) solid-angle fraction z_ref of a plane tangent to the bead: s int_0^1 W(x)/(1+s x)^2 dx; E_ref = -w~ z_ref."""
    return s * quad(lambda x: float(W(x, p)) / (1 + s * x) ** 2, 0, 1)[0]


def section_ref(show=True):
    out = {}
    for s in (0.25, 0.1, 0.05):
        out[s] = {p: gamma_ref(s, p) for p in (1, 2, 4, 8)}
        if show:
            print("s=%.2f  z_ref (= -E_ref/w~):" % s, {p: round(v, 4) for p, v in out[s].items()})
    return out



# =====================================================================================================
# 4b. checks of the first-order model against the old simulations (shell 0.25, p = 1, read-only data)
# =====================================================================================================
OLD = "/home/mrojasve/Documents/DDG/MiraculousMembranes/projects/geometric-flow/Results/Wrapping_planar_excess/"


def _load_cov(fn):
    L = [l.split() for l in open(OLD + fn) if not l.startswith("#")]
    return np.array([[float(x) for x in l[1:]] for l in L if len(l) == 15])


def section_check(show=True):
    """(i) sigma~ = 0 and 0.133: z(w~), E(w~) of the old flat start vs the first-order model;
       (ii) 263 flat-start points of baseline_partial_branch.csv: E_sim - E_first-order and z_ad - z."""
    out = {}
    try:
        fl = _load_cov("Coverage_data.txt")
    except Exception as e:
        print("old data not available:", e)
        fl = None
    zz = np.linspace(1e-6, 2 - 1e-6, 40001)
    if fl is not None:
        for KA, key in ((0.0, "sigma0"), (0.0667, "sigma0.133")):
            A = fl[np.abs(fl[:, 0] - KA) < 1e-3]
            A = A[np.argsort(A[:, 2])]
            rows = []
            for r in A:
                w = 4 * r[2]
                if KA == 0:
                    wD = 4.0 * np.ones_like(zz)
                    E = -(w - 4) * zz - TP[1] * 0.5 * wD ** 0.75 * np.sqrt(zz * (2 - zz))
                    zf, Ef = zz[np.argmin(E)], E.min()
                else:
                    th = Deserno(1)
                    E = th.E(w, TP[1], 0.25)
                    m = argrelextrema(E, np.less, order=20)[0]
                    zf, Ef = (th.zg[m[0]], E[m[0]]) if len(m) else (np.nan, np.nan)
                rows.append((w, 2 * r[6], -r[12] / (2 * np.pi * r[2]), 2 * r[13] / np.pi, zf, Ef))
            out[key] = np.array(rows)
            if show:
                print(key, " w~  z_sim  z_ad,sim  E_sim | z_model  E_model")
                for q in out[key][::2]:
                    print("   %.1f  %.3f  %.3f  %.3f | %.3f  %.3f" % tuple(q))
    import csv
    rows = [r for r in csv.DictReader(open(os.path.join(DATA, "baseline_partial_branch.csv"))) if r["file"] == "flat"]
    ths, res = {}, []
    for r in rows:
        sg = float(r["sigma"])
        k = int(round(sg * 15 / 2))
        if k < 1 or abs(2 * k / 15 - sg) > 1e-3:
            continue
        th = ths.setdefault(k, Deserno(k))
        w, E, zad = float(r["w"]), float(r["E"]), float(r["zad"])
        Et = th.E(w, TP[1], 0.25)
        m = argrelextrema(Et, np.less, order=20)[0]
        if not len(m):
            continue
        j = m[np.argmin(abs(th.zg[m] - zad))]
        res.append((sg, w, zad, E, th.zg[j], Et[j]))
    res = np.array(res)
    sel = res[(res[:, 2] < 1.6) & (res[:, 4] < 1.7)]
    out["partial"] = sel
    if show:
        d = sel[:, 3] - sel[:, 5]
        print("partial branch, %d flat-start points: E_sim - E_model = %.3f +- %.3f (rms %.3f);  z_ad,sim - z_model = %.3f (rms %.3f)"
              % (len(sel), d.mean(), d.std(), np.sqrt((d ** 2).mean()), (sel[:, 2] - sel[:, 4]).mean(),
                 np.sqrt(((sel[:, 2] - sel[:, 4]) ** 2).mean())))
        for lo, hi in ((0, 0.5), (0.5, 0.9), (0.9, 1.3), (1.3, 1.7), (1.7, 2.1)):
            q = sel[(sel[:, 0] >= lo) & (sel[:, 0] < hi)]
            if len(q):
                print("   sigma~ in [%.1f,%.1f): n=%d  E_sim - E_model = %.3f +- %.3f" % (lo, hi, len(q), (q[:, 3] - q[:, 5]).mean(), (q[:, 3] - q[:, 5]).std()))
    return out


def section_vs_t1(show=True):
    """Compare the first-order model with the full continuum soft-shell solution of T1 (data/softshell_lines.csv),
    for every (sigma~, s, p) where both exist (sigma~ matched to the theory cache 2k/15 within 0.01)."""
    import csv
    f = os.path.join(DATA, "softshell_lines.csv")
    if not os.path.exists(f):
        print("no T1 data")
        return []
    out = []
    if show:
        print(" sigma~    s     p |   T1:  S1     E     S2   |  first order:  S1     E     S2   | first-order minus T1")
    seen = set()
    for r in csv.DictReader(open(f)):
        sg, s_, p = float(r["sigma_t"]), float(r["s"]), int(float(r["p"]))
        k = int(round(sg * 15 / 2))
        if k < 1 or abs(2 * k / 15 - sg) > 0.01 or (k, s_, p) in seen:
            continue
        seen.add((k, s_, p))
        if p not in TP:
            continue
        th = Deserno(k)
        m = lines(th, TP[p], s_)
        t1 = [float(r[x]) for x in ("S1", "E", "S2")]
        fo = [m.get(x, np.nan) for x in ("S1", "E", "S2")]
        out.append(dict(sigma=th.sig, s=s_, p=p, t1=t1, fo=fo))
        if show:
            print(" %.3f  %.3f  %d | %6.3f %6.3f %6.3f | %6.3f %6.3f %6.3f | %+.3f %+.3f %+.3f" % (
                th.sig, s_, p, *t1, *fo, *(np.array(fo) - np.array(t1))))
    return out

# =====================================================================================================
# 5. discretisation of the boundary layer on a regular 1D mesh (units: ell = sqrt(delta/Dc), p = 1)
# =====================================================================================================
def disc_T(p, hh, phi, sample="node", L1=24.0, L2=7.0, tc=2.0):
    t0 = -L1
    n = int(round((L1 + L2) / hh))
    t = t0 + hh * np.arange(n + 1)
    tcv = tc + phi * hh
    nu = n + 1 - 4

    def unpack(x):
        d = np.zeros(n + 1)
        d[2:n - 1] = x[:nu]
        d[n - 1:] = (0.5 * (t - tcv) ** 2 + x[nu])[n - 1:]
        return d

    def E(x):
        d = unpack(x)
        c2 = (d[2:] - 2 * d[1:-1] + d[:-2]) / hh ** 2
        Eb = 0.5 * np.sum(c2 ** 2) * hh
        if sample == "node":
            Ep = -0.5 * np.sum(W(d, p)[1:-1]) * hh
        else:
            Ep = -0.5 * np.sum(W(0.5 * (d[1:] + d[:-1]), p)) * hh
        return Eb + Ep

    x0 = np.concatenate([np.where(t[2:n - 1] > tcv, 0.5 * (t[2:n - 1] - tcv) ** 2, 0.0), [0.0]])
    r = minimize(E, x0, method="L-BFGS-B", options=dict(maxiter=20000, maxfun=10 ** 7, ftol=1e-15, gtol=1e-10))
    a, b = t[1] - 0.5 * hh, t[n - 1] + 0.5 * hh
    ref = (-0.5 * (tcv - a) + 0.5 * (b - tcv)) if sample == "node" else (-0.5 * (tcv - t[0]) + 0.5 * (b - tcv))
    return -2 * (r.fun - ref)


def section_disc(show=True, force=False):
    res = None if force else _cache_get("disc")
    if res is None:
        res = []
        for hh in (1.0, 0.5, 0.25, 0.125):
            row = dict(h_over_ell=hh)
            for smp in ("node", "mid"):
                v = [disc_T(1, hh, ph, smp) for ph in (0.0, 0.5)]
                row[smp] = v
            res.append(row)
        _cache_put("disc", res)
    if show:
        print("discrete BL, p=1, continuum T_1 = %.4f" % TP[1])
        for r in res:
            print("  h/ell=%.3f  nodal sampling T = %.4f / %.4f (lattice phase 0 / 1/2)  centroid-like sampling T = %.4f / %.4f"
                  % (r["h_over_ell"], r["node"][0], r["node"][1], r["mid"][0], r["mid"][1]))
    return res


# =====================================================================================================
# 6. figures
# =====================================================================================================
def figures(bl, lines_rows, tilt, disc, ref, extrap=None):
    os.makedirs(FIG, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25})
    # ---- fig 1: boundary layer
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.6))
    for p, c in zip((1, 2, 4, 8), ("C0", "C1", "C2", "C3")):
        r = bl[p]
        t = r["t"] - r["tc"]
        sel = (t > -9) & (t < 3.5)
        ax[0].plot(t[sel], r["D"][sel], color=c, label="p=%d, T_p=%.3f" % (p, r["T"]))
    tt = np.linspace(0, 3.5, 100)
    ax[0].plot(tt, 0.5 * tt ** 2, "k--", lw=0.8, label="sharp parabola")
    ax[0].axhline(1, color="gray", lw=0.7)
    ax[0].set_ylim(-0.3, 3)
    ax[0].set_xlabel(r"$x-x_c$ in units $\sqrt{\delta/\Delta c}$")
    ax[0].set_ylabel(r"$d/\delta$")
    ax[0].set_title("boundary layer (gray: edge of the shell)")
    ax[0].legend(fontsize=7)
    ax[1].plot(tt, 0 * tt, alpha=0)
    for p, c in zip((1, 2, 4, 8), ("C0", "C1", "C2", "C3")):
        r = bl[p]
        t = r["t"] - r["tc"]
        sel = (t > -10) & (t < 4)
        C = r["u"] * 0 + 0
    ps = np.array([1, 2, 4, 8])
    ax[1].plot(ps, [bl[p]["T"] for p in ps], "o-", label=r"$T_p$ (relaxed boundary layer)")
    ax[1].plot(ps, [trial_bound(p) for p in ps], "s--", label="trial bound (sharp profile)")
    ax[1].plot(ps, [bl[1]["T"] * p ** -0.25 for p in ps], "k:", label=r"$T_1\,p^{-1/4}$")
    ax[1].set_xscale("log", base=2)
    ax[1].set_xlim(0.8, 10)
    ax[1].set_xticks([1, 2, 4, 8])
    ax[1].set_xticklabels(["1", "2", "4", "8"])
    ax[1].set_ylim(0.4, 1.35)
    ax[1].set_xlabel("shell power p")
    ax[1].set_ylabel(r"$T_p$  ($\tau=-T_p\,w\sqrt{\delta\,\ell_c}$)")
    ax[1].legend(fontsize=7)
    ax[1].set_title("excess energy of the contact line")
    h = np.array([r["h_over_ell"] for r in disc])
    ax[2].plot(h, [abs(r["node"][0] - bl[1]["T"]) / bl[1]["T"] for r in disc], "o-", label="nodal sampling, phase 0")
    ax[2].plot(h, [abs(r["node"][1] - bl[1]["T"]) / bl[1]["T"] for r in disc], "o--", label="nodal sampling, phase 1/2")
    ax[2].plot(h, [abs(r["mid"][0] - bl[1]["T"]) / bl[1]["T"] for r in disc], "s-", label="centroid-type sampling")
    ax[2].plot(h, 0.075 * h ** 2, "k:", label=r"$0.075\,(h/\ell)^2$")
    ax[2].set_xscale("log")
    ax[2].set_yscale("log")
    ax[2].set_xlabel(r"$h/\ell$, $\ell=\sqrt{\delta/\Delta c}$")
    ax[2].set_ylabel(r"|relative error of $T_1$| (regular 1D mesh)")
    ax[2].legend(fontsize=7)
    ax[2].set_title("discretisation error")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "boundary_layer.png"), dpi=150)
    plt.close(fig)

    # ---- fig 2: lines in (w~, sigma~): Deserno, first-order model, simulation
    sim = sim_lines()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    sg, DE, D1, D2, ME, M1, M2 = [], [], [], [], [], [], []
    for k in range(1, 16):
        th = Deserno(k)
        b = lines(th, 0.0, 0.0)
        r = lines(th, TP[1], 0.25)
        sg.append(th.sig)
        DE.append(b.get("E", np.nan)); D1.append(b.get("S1", np.nan)); D2.append(b.get("S2", np.nan))
        ME.append(r.get("E", np.nan)); M1.append(r.get("S1", np.nan)); M2.append(r.get("S2", np.nan))
    sg = np.array(sg)
    S = [(float(r["sigma_t"]), float(r["E_sim"]), float(r["S1_sim"]), float(r["S2_sim"])) for r in sim.values()]
    S = np.array(S)
    a = ax[0]
    a.plot(DE, sg, "k-", lw=2, label="Deserno E")
    a.plot(D1, sg, "k--", label="Deserno S1")
    a.plot(D2, sg, "k:", label="Deserno S2")
    a.plot(ME, sg, "-", color="C3", lw=2, label="first-order line tension E")
    a.plot(M1, sg, "--", color="C3", label="S1")
    a.plot(M2, sg, ":", color="C3", label="S2")
    a.plot(S[:, 1], S[:, 0], "o", color="C0", ms=5, label="simulation E")
    a.plot(S[:, 2], S[:, 0], "^", color="C0", ms=5, mfc="none", label="simulation S1")
    a.plot(S[:, 3], S[:, 0], "v", color="C0", ms=5, mfc="none", label="simulation S2")
    a.set_xlabel(r"$\tilde w$")
    a.set_ylabel(r"$\tilde\sigma$")
    a.set_xlim(3, 13)
    a.set_title("s = 0.25, p = 1, no fit parameter")
    a.legend(fontsize=7, ncol=2)
    # panel b: shifts of E vs ell
    b_ = ax[1]
    for p, c in zip((1, 2, 4, 8), ("C0", "C1", "C2", "C3")):
        for sg_, ls in ((0.5333, "-"), (1.0667, "--")):
            rr = [r for r in lines_rows if r["p"] == p and abs(r["sigma"] - sg_) < 1e-3]
            rr.sort(key=lambda r: r["s"])
            b_.plot([np.sqrt(r["s"]) for r in rr], [r["dE"] for r in rr], ls, color=c, marker="o", ms=3,
                    label=("p=%d" % p) if ls == "-" else None)
    b_.set_xlabel(r"$\sqrt{s}$  ($\ell/a=\sqrt{s}\,\tilde w^{-1/4}$)")
    b_.set_ylabel(r"shift of the E line, $\Delta\tilde w_E$")
    b_.set_title(r"first-order shift of E (solid $\tilde\sigma$=0.53, dashed 1.07)")
    b_.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "line_tension_lines.png"), dpi=150)
    plt.close(fig)

    # ---- fig 2b: extrapolation protocol
    if extrap:
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
        for a_, k in zip(ax, sorted(extrap)):
            res = extrap[k]
            x = np.sqrt(np.array([0.25, 0.16, 0.09, 0.04, 0.0225]))
            xx = np.linspace(0, 0.52, 100)
            for key, c in (("E", "C0"), ("S1", "C1"), ("S2", "C2")):
                y = np.array(res[key]["values"])
                a_.plot(x, y, "o", color=c, label=key)
                cf = np.polyfit(x[1:], y[1:], 2)
                a_.plot(xx, np.polyval(cf, xx), "-", color=c, lw=0.8)
                a_.plot([0], [res[key]["Deserno"]], "*", color=c, ms=11, mec="k")
            a_.set_xlabel(r"$\sqrt{s}$")
            a_.set_ylabel(r"$\tilde w$")
            a_.set_title(r"$\tilde\sigma$ = %.2f: quadratic fit in $\sqrt{s}$ (stars: Deserno)" % k)
            a_.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "extrapolation_protocol.png"), dpi=150)
        plt.close(fig)

    # ---- fig 3: regularisation trade-off
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    cols = {0: "k", 4: "C0", 8: "C1", 16: "C3"}
    for q in (0, 4, 8, 16):
        pts = [(a_["lw"], a_["Theta"]) for a_ in tilt if a_["q"] == q]
        ax[0].plot(*zip(*pts), "o", color=cols[q], ms=3, label="q=%d" % q)
    xx = np.linspace(0.15, 1.0, 10)
    ax[0].plot(xx, 1.259 * xx, "k-", lw=0.8)
    ax[0].set_xlabel(r"$\ell_w=\int W\,T\,dx$ over the lift-off side (units of $\ell_c$)")
    ax[0].set_ylabel(r"$\Theta=-\tau/(w\,\ell_c)$")
    ax[0].set_title("tilt weight $Y^q$ at equal lift-off width: no gain")
    ax[0].legend(fontsize=7)
    for p, c in zip((1, 2, 4, 8), ("C0", "C1", "C2", "C3")):
        pts = sorted([(a_["dh"], a_["Theta"]) for a_ in tilt if a_["q"] == 0 and a_["p"] == p])
        ax[1].plot(*zip(*pts), "o-", color=c, ms=3, label="p=%d" % p)
    ax[1].set_xlabel(r"$\delta/\ell_c = s\sqrt{\tilde w}$")
    ax[1].set_ylabel(r"$\Theta$ (q=0)")
    ax[1].legend(fontsize=7)
    ax[1].set_title(r"$\Theta=T_p\sqrt{\delta/\ell_c}$")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "regularisation_tradeoff.png"), dpi=150)
    plt.close(fig)


# =====================================================================================================
def main(argv):
    secs = argv[1:] or ["all"]
    allp = "all" in secs
    bl = tilt = disc = rows = ref = extrap = None
    t0 = time.time()
    if allp or "bl" in secs:
        print("== flat boundary layer ==")
        bl = section_bl()
        for p in (1, 2, 4, 8):
            TP[p] = bl[p]["T"]
    if allp or "ref" in secs:
        print("== planar reference state ==")
        ref = section_ref()
    if allp or "disc" in secs:
        print("== discretisation (regular 1D mesh) ==")
        disc = section_disc()
    if allp or "tilt" in secs:
        print("== tilt weight ==")
        tilt = section_tilt()
    if allp or "lines" in secs:
        print("== first-order line-tension model ==")
        rows = section_lines()
        print("== extrapolation test ==")
        extrap = section_extrapolation()
    if allp or "t1" in secs:
        print("== first-order model vs the full continuum soft-shell solver (T1) ==")
        section_vs_t1()
    if allp or "check" in secs:
        print("== check against old simulations ==")
        section_check()
    if allp and bl and tilt and disc and rows:
        figures(bl, rows, tilt, disc, ref, extrap)
    print("done in %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main(sys.argv)
