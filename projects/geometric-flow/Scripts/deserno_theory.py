#!/usr/bin/env python
"""
deserno_theory.py -- numerical solution of the free-membrane problem of

    M. Deserno, "Elastic deformation of a fluid membrane upon colloid binding",
    Phys. Rev. E 69, 031903 (2004) / cond-mat/0303656,

and the phase lines derived from it (E, W, spinodals S1/S2, energy barrier, Fig. 5 curves).
Pure python, numpy + scipy only.

Notation (a = colloid radius, kappa = bending modulus, sigma = tension, w = adhesion energy / area)
    lambda = sqrt(kappa/sigma)         w~ = 2 w a^2 / kappa        sigma~ = sigma a^2 / kappa = (a/lambda)^2
    E~ = E/(pi kappa)                  z = 1 - cos(alpha) in [0, 2] (degree of wrapping)
    E~(z) = -(w~ - 4) z + sigma~ z^2 + F(z; sigma~)           (paper Eq. 4,  F = E~_free)

The free membrane energy (lengths in units of a)
    F(z; sigma~) = int_0^inf ds  r [ (psi' + sin(psi)/r)^2 + 2 sigma~ (1 - cos psi) ],
    r' = cos(psi), r(0) = sin(alpha), psi(0) = alpha, psi(s) -> 0, psi' -> 0 (s -> inf)
is obtained by SHOOTING on the Hamiltonian shape equations (paper Eqs. 6-13).

METHOD (what this file actually does)
  * Hamilton equations (Eq. 9, p_h = 0, H = 0) are integrated INWARD (stable direction) from a far-field
    Bessel (K1) start at psi = psi_m = 1e-5 with an adaptive 8th-order Runge-Kutta (DOP853, rtol 1e-11).
    The one-parameter family of decaying solutions is labelled by r_m (radius where psi = psi_m), and the
    contact condition r = sin(psi) selects the stopping point u*.  The solution CURVE g(r_m, u) = 0
    (g = r - sin psi) is followed by pseudo-arclength continuation, which passes folds of the
    S-shaped (sigma~ > sigma~_c = 4.72) multivalued regime without any difficulty.  Along the curve one
    reads off z = 1 - cos(psi*), the contact curvature psi_dot_0 = q - 1 (q = p_psi/(2r)) and the free energy
    F (accumulated energy density + analytic energy of the linear tail).
  * Exact relation used for interpolation and for the phase lines (derived from the boundary terms of the
    variation, Eqs. 7, 13, 14):        dF/dz = (1 - psi_dot_0)^2 - 4 - 2 sigma~ z,
                                       dE/dz = (1 - psi_dot_0)^2 - w~ .
    So F(z) is interpolated with cubic Hermite splines using the ODE energies AND these exact slopes, and
    w_eq(z) := (1 - psi_dot_0(z))^2 is the 'equation of state' whose
        maximum  on the partially wrapped branch  = spinodal S1,
        minimum  near the neck (z ~ 1.9)            = spinodal S2.
  * Multivalued regime (sigma~ > sigma~_c): the solution curve is S-shaped in (z, psi_dot_0).  The branch
    attached to z = 0 ('branch 1', psi_dot_0 < 0, partially wrapped) ends in a fold at z_max(sigma~) < 2;
    the lowest energy solution beyond it lives on the branch attached to the neck z = 2 ('branch 3',
    obtained by a reverse trace started at 2 - delta).  free_energy() returns the LOWEST energy solution:
    branch 1 for z < z_x, branch 3 for z > z_x (z_x = energy crossing).  Those branch-3 states have
    psi_dot_0 > 1/a for z in (z_x, ~z_x+) and are geometrically inaccessible (paper, Sec. III.C); they never
    enter the E line, the barrier or S1 (all found on branch 1).  For sigma~ > 1e4 the branch-3 sliver
    2 - z < ~3e-4 is bridged by a cubic Hermite connection to F(2) = 0, F'(2) = -4 sigma~ (flagged).
  * Phase lines: E line = w~ where min_z E(z) (partially wrapped minimum, w_eq(z_p) = w~ on the rising part
    of w_eq) equals E(2) = -2(w~-4) + 4 sigma~ (E~_free(2) = 0); barrier = E(z_b) - E(2) where z_b is the
    maximum of E (w_eq(z_b) = w~ on the falling part).

Independent checks (see test_deserno_theory.py): direct variational minimisation of the discretised functional
(prior/deserno_var.py, not needed here) agrees with F to 1e-8..6e-6 relative; the exact slope identity above
holds along the traced curves; Eq. 20/21 small-gradient formulas, Eq. 26, Eqs. 34-37 are reproduced.

PUBLIC API  (table: deserno_table.npz / deserno_table.csv next to this file; sigma~ = 1e-4 ... 1e6, 20 nodes/decade)
    free_energy(z, sigma_t)              lowest-energy E~_free(z; sigma~)   (solves + caches; instantaneous at table
                                          nodes, otherwise ~30-70 s for the first call at a new sigma~)
    contact_curvature(z, sigma_t)        psi_dot_0 (a=1) of the lowest-energy shape
    total_energy(z, w_t, sigma_t)
    w_E(sigma_t), barrier(sigma_t), spinodal_S1(sigma_t), spinodal_S2(sigma_t), z_partial(...), z_barrier(...)
                                          (cubic splines in log10 sigma~ of the table; outside the table range w_E and
                                          barrier fall back to the paper asymptotics, scaled to match, with a warning)
    fig5_curves(a_over_lambda)           -> (w/sigma on E, w/sigma on W)
    plot_fig2_lines(ax, ylim=(0, 1))     draws E (solid) and S1, S2 (short dashed) on a matplotlib axes (w~ horizontal,
                                          sigma~ vertical);  fig2_line_data(ylim) returns the same lines as arrays
    load_table(), interp_logsigma(), TABLE_COLUMNS, export_csv(), refresh_table(), rebuild_rows()
    solve_curve(sigma_t) -> SolutionCurve (all details), compute_curve / compute_point(sigma_t), generate_table(...)
    asymptotics: F_small_gradient (Eq. 20), F_small_z (Eq. 21), w_E_small_gradient (Eq. 26/27),
                 high_tension_wE / _z / _barrier (Eqs. 34-37)
  Table columns: sigma_t, w_E, w_over_sigma (= w~/(2 sigma~) = w/sigma), a_over_lambda, barrier, z_partial (penetration of
  the partially wrapped state on E), z_barrier, S1, S2, zS1, zS2, wE_minus4_over_sigma, s_shaped (fold resolved in
  z < 2 - 1.5e-4), z_fold, z_cross (energy crossing branch 1 / branch 3), branch3_exact, asymptotic_fallback (all 0),
  scan_err, barrier_scan_err (brute-force verification), npts.   S2 = NaN for sigma~ > 1e4 (branch 3 not traced).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")      # small dense ops only; threads would just add contention
import math
import time
import warnings
import numpy as np
from scipy.integrate import ode, solve_ivp
from scipy.special import kve, k0, k1, lambertw
from scipy.optimize import brentq, minimize_scalar
from scipy.interpolate import CubicSpline, CubicHermiteSpline

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE_PATH = os.path.join(HERE, "deserno_table.npz")

GAMMA_E = 0.5772156649015329
PSI_M = 1e-5                      # far-field start amplitude (domain truncation); 1e-4 / 1e-6 change w_E by <1e-7 rel.
Z_END = 2.0 - 1.5e-4              # forward traces stop here (F(2) = 0, F'(2) = -4 sigma~ known exactly)
PROFILE = dict(dzmax=0.02, dpmax=0.02, dFmax=0.05)    # step control of the continuation (see convergence study)
SIGMA_B3_MAX = 1.0e4              # branch 3 (neck branch) traced for sigma~ <= this; bridged above

# values quoted in the paper (used by the tests only)
PAPER = dict(sigma_c=4.721139, z_c=1.86289, A=5.650, wE_1=6.142, barrier_1=1.05,
             S1_1=7.5, S2_1=2.7, minfig5=(4.4, 1.37))


# =================================================================================================
# 1.  Shooting core
# =================================================================================================
class Shooter:
    """Hamilton equations integrated inward in u = S - s.  State y = [psi, q, r, p_r, E],
    q = p_psi/(2r) = psi_dot + sin(psi)/r, E = accumulated energy (incl. analytic linear tail)."""

    def __init__(self, sigma, rtol=1e-11, psi_m=PSI_M):
        self.sig = float(sigma)
        self.k = math.sqrt(sigma)
        self.rtol = rtol
        self.psi_m = psi_m
        self.ntraj = 0

    def far_start(self, rm):
        sig, k, psi_m = self.sig, self.k, self.psi_m
        x = k * rm
        a0 = kve(0, x)
        a1 = kve(1, x)
        q = -psi_m * k * a0 / a1
        sp = math.sin(psi_m)
        cp = math.cos(psi_m)
        omc = 2.0 * math.sin(0.5 * psi_m) ** 2
        pr = (-rm * q * q + 2.0 * q * sp + 2.0 * sig * rm * omc) / cp       # H = 0
        E = psi_m ** 2 * x * a0 / a1                                       # linear-theory energy of the tail
        return [psi_m, q, rm, pr, E]

    def rhs(self, u, y):
        psi, q, r, pr, E = y
        sig = self.sig
        sp = math.sin(psi)
        cp = math.cos(psi)
        omc = 2.0 * math.sin(0.5 * psi) ** 2
        return [-(q - sp / r), -(sig + pr / (2.0 * r)) * sp, -cp,
                -(q * q - 2.0 * q * sp / r + 2.0 * sig * omc), r * q * q + 2.0 * sig * r * omc]

    def state_at(self, rm, u):
        o = ode(self.rhs).set_integrator('dop853', rtol=self.rtol, atol=1e-20, nsteps=500000)
        o.set_initial_value(self.far_start(rm), 0.0)
        y = o.integrate(u)
        self.ntraj += 1
        if not o.successful():
            return None
        return y

    @staticmethod
    def g(y):
        return y[2] - math.sin(y[0])

    @staticmethod
    def gu(y):
        psi, q, r = y[0], y[1], y[2]
        return -math.cos(psi) * (1.0 - (q - math.sin(psi) / r))

    @staticmethod
    def features(y):
        """(z, psi_dot_0, F) at a contact state."""
        psi, q, r = y[0], y[1], y[2]
        return 2.0 * math.sin(psi / 2.0) ** 2, q - math.sin(psi) / r, y[4]

    def newton_u(self, rm, u, tol=1e-13, maxit=40):
        for _ in range(maxit):
            y = self.state_at(rm, u)
            if y is None:
                return None
            du = -self.g(y) / self.gu(y)
            u += du
            if abs(du) < tol * (1.0 + abs(u)):
                return u
        return None

    def start_point(self, z0):
        """(rm, u) on the curve close to the small-gradient (Bessel) profile at wrapping ~ z0."""
        k = self.k
        alpha = math.acos(1.0 - z0)
        rc = math.sin(alpha)
        target = math.log(self.psi_m) + (math.log(kve(1, k * alpha)) - k * alpha) - math.log(alpha)
        x = brentq(lambda x: math.log(kve(1, x)) - x - target, k * alpha, 5000.0)
        rm = x / k
        return rm, self.newton_u(rm, rm - rc)


def trace(sh, z0=1e-4, z_end=Z_END, h0=0.02, max_pts=6000, dzmax=0.04, dpmax=0.15, dFmax=0.2,
          start=None, verbose=False, tol_acc=1e-5, delta_rel=1e-7, stop_at_fold=False,
          orient_down=False, z_stop=None, extra_after_fold=3):
    """Pseudo-arclength continuation of the solution curve g(r_m, u) = 0.

    Returns (pts, info); pts columns: [rm, u, z, psi_dot_0, F, psi, q, r, p_r, dz/dt].
    stop_at_fold: stop a few points after the first z-extremum (sign change of dz/dt);
    orient_down:  begin in the direction of decreasing z (reverse trace from the neck), stop at z_stop.
    """
    if start is None:
        start = sh.start_point(z0)
    x = np.array(start, float)
    delta = delta_rel * max(1.0, x[0])
    tol_abs = tol_acc * min(1.0, 1.0 / sh.k)          # Newton acceptance: error after step ~ dx^2/lambda

    def grad(x):
        y0 = sh.state_at(x[0], x[1])
        y1 = sh.state_at(x[0] + delta, x[1])
        if y0 is None or y1 is None:
            return None, None, None
        gr_ = np.array([(sh.g(y1) - sh.g(y0)) / delta, sh.gu(y0)])
        return y0, gr_, (y1[0] - y0[0]) / delta

    pts = []
    t_prev = np.array([1.0, 0.0])
    h = h0 * max(1.0, 0.1 * x[0])
    y0, gr, dpr = grad(x)
    if y0 is None:
        return np.zeros((0, 10)), dict(ok=False, reason='start failed')
    nafter = None
    while len(pts) < max_pts:
        t = np.array([-gr[1], gr[0]]) / np.hypot(gr[0], gr[1])
        sp0, c0 = math.sin(y0[0]), y0[1] - math.sin(y0[0]) / y0[2]
        if orient_down and not pts:
            if sp0 * (dpr * t[0] - c0 * t[1]) > 0:
                t = -t
        elif t @ t_prev < 0:
            t = -t
        z, pd, F = sh.features(y0)
        dzdt = sp0 * (dpr * t[0] - c0 * t[1])
        pts.append([x[0], x[1], z, pd, F, y0[0], y0[1], y0[2], y0[3], dzdt])
        if verbose and len(pts) % 20 == 0:
            print(f"  pt {len(pts):4d} rm={x[0]:.5f} u={x[1]:.5f} z={z:.6f} pd0={pd:.5f} F={F:.6g} h={h:.2e}")
        if (not orient_down) and z > z_end:
            break
        if z_stop is not None and z < z_stop:
            break
        if stop_at_fold:
            if len(pts) > 1 and pts[-2][9] * dzdt < 0 and nafter is None:
                nafter = 0
            if nafter is not None:
                nafter += 1
                if nafter > extra_after_fold:
                    break
        while True:
            xp = x + h * t
            xn = xp.copy()
            ok = False
            for it in range(14):
                ya, grn, _ = grad(xn)
                if ya is None:
                    break
                J = np.array([[grn[0], grn[1]], [t[0], t[1]]])
                dx = np.linalg.solve(J, np.array([-sh.g(ya), -t @ (xn - xp)]))
                xn = xn + dx
                if np.hypot(*dx) < tol_abs:
                    ok = True
                    break
            if ok:
                yn, grn, dprn = grad(xn)
                if yn is not None:
                    zn, pdn, Fn = sh.features(yn)
                    d = max(abs(zn - z) / dzmax, abs(np.arcsinh(pdn) - np.arcsinh(pd)) / dpmax,
                            abs(math.log(abs(Fn) + 1e-12) - math.log(abs(F) + 1e-12)) / dFmax)
                    if d > 1.5 and h > 1e-8:
                        h *= 0.5
                        continue
                    if d < 0.5:
                        h *= 1.5
                    break
            h *= 0.5
            if h < 1e-10:
                return np.array(pts), dict(ok=False, reason='step too small', last=(x.copy(), z))
        t_prev = t
        x, y0, gr, dpr = xn, yn, grn, dprn
    return np.array(pts), dict(ok=True)


def neck_start(sh, z_t=None, n=70, verbose=False):
    """(rm, u) of a solution with z = z_t ~ 2 (neck branch): scan r_m at fixed contact angle, handling the
    windows in which the trajectory collapses onto the singular r -> 0 solution.  z_t defaults to
    2 - min(2e-4, 2e-3/sigma~) so that the neck radius is small compared with lambda."""
    if z_t is None:
        z_t = 2.0 - min(2e-4, 2e-3 / sh.sig)
    alpha = math.acos(1.0 - z_t)
    sin_a = math.sin(alpha)
    lam = 1.0 / sh.k
    rm_lo = 0.02 * min(1.0, lam)
    rm_hi = 25.0 * max(lam, 0.05) + 3.0

    def ev_psi(u, y):
        return y[0] - alpha
    ev_psi.terminal = True
    ev_psi.direction = 1

    def ev_r(u, y):
        return y[2] - 1e-9
    ev_r.terminal = True
    ev_r.direction = -1

    def ev_q(u, y):
        return 1e5 - abs(y[1])
    ev_q.terminal = True

    def G(rm):
        sol = solve_ivp(sh.rhs, (0.0, 1.5 * rm + 0.5 * lam + 0.1), sh.far_start(rm), method='DOP853',
                        rtol=sh.rtol, atol=1e-20, events=(ev_psi, ev_r, ev_q))
        sh.ntraj += 1
        if len(sol.t_events[0]):
            return sol.y_events[0][0][2] - sin_a, sol.t_events[0][0]
        return -sin_a, np.nan

    grid = np.geomspace(rm_lo, rm_hi, n)
    res = [G(rm) for rm in grid]
    vals = [v[0] for v in res]
    has = [np.isfinite(v[1]) for v in res]
    roots = []

    def accept(rm):
        gv, tu = G(rm)
        if np.isfinite(tu) and abs(gv) < 1e-8:
            roots.append(rm)

    def edge(a, b, valid_left):
        lo, hi = a, b
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if np.isfinite(G(mid)[1]) == valid_left:
                lo = mid
            else:
                hi = mid
        return lo, hi

    for i in range(len(grid) - 1):
        a, b = grid[i], grid[i + 1]
        if has[i] and has[i + 1] and vals[i] * vals[i + 1] < 0:
            accept(brentq(lambda r: G(r)[0], a, b, xtol=1e-14, rtol=1e-13))
        elif has[i] and not has[i + 1]:
            lo, hi = edge(a, b, True)
            if G(lo)[0] * vals[i] < 0:
                accept(brentq(lambda r: G(r)[0], a, lo, xtol=1e-14, rtol=1e-13))
        elif (not has[i]) and has[i + 1]:
            lo, hi = edge(a, b, False)
            if G(hi)[0] * vals[i + 1] < 0:
                accept(brentq(lambda r: G(r)[0], hi, b, xtol=1e-14, rtol=1e-13))
    if verbose:
        print('neck_start roots', roots)
    if not roots:
        return None
    rm = roots[0]
    u = sh.newton_u(rm, G(rm)[1])
    return rm, u


# =================================================================================================
# 2.  Asymptotic / analytic formulas of the paper
# =================================================================================================
def F_small_gradient(z, sigma_t):
    """Eq. (20): small-gradient free energy, k = sin(alpha) = sqrt(z(2-z)), a/lambda = sqrt(sigma~)."""
    z = np.asarray(z, float)
    kk = np.sqrt(z * (2.0 - z))
    x = kk * np.sqrt(sigma_t)
    return np.sqrt(sigma_t) * kk ** 3 / (1.0 - kk ** 2) * k0(x) / k1(x)


def F_small_z(z, sigma_t):
    """Eq. (21): -2 sigma~ z^2 (2 gamma + ln(sigma~ z / 2)) + O(z^3)."""
    z = np.asarray(z, float)
    return -2.0 * sigma_t * z * z * (2.0 * GAMMA_E + np.log(sigma_t * z / 2.0))


def sigma_E_eq26(w):
    """Eq. (26): sigma~ on the E line as a function of w~ (small-gradient, exact for sigma~ -> 0), w < 4.928."""
    x = -(w - 4.0) / 8.0 * np.exp(2.0 * GAMMA_E)
    W = lambertw(x, -1).real
    return (w - 4.0) / 4.0 * (1.0 + np.sqrt(1.0 + 1.0 / (2.0 * W) + 1.0 / (2.0 * W) ** 2))


def w_E_small_gradient(sigma_t, order=26):
    """Invert Eq. (26) (order=26) or the lowest-order Eq. (27) (order=27) for w~_E(sigma~)."""
    sigma_t = np.atleast_1d(np.asarray(sigma_t, float))
    out = np.full(sigma_t.shape, np.nan)
    wmax = 4.0 + 8.0 * math.exp(-1.0 - 2.0 * GAMMA_E)
    for i, s in enumerate(sigma_t):
        if order == 26:
            f = lambda w: sigma_E_eq26(w) - s
        else:
            f = lambda w: (w - 4.0) / 2.0 * (1.0 + 1.0 / (8.0 * (2.0 * GAMMA_E + np.log((w - 4.0) / 8.0)))) - s   # Eq. 27
        try:
            out[i] = brentq(f, 4.0 + 1e-14, wmax - 1e-9 if order == 26 else 5.0, xtol=1e-14, rtol=1e-14)
        except ValueError:
            pass
    return out


def high_tension_wE(sigma_t, A=PAPER['A']):
    """Eq. (35): (w~ - 4)/sigma~ = 4 - 3 A^(2/3) sigma~^(-1/3)  ->  w~_E."""
    s = np.asarray(sigma_t, float)
    return 4.0 + s * (4.0 - 3.0 * A ** (2.0 / 3.0) * s ** (-1.0 / 3.0))


def high_tension_z(sigma_t, A=PAPER['A']):
    """Eq. (34) penetration on the E line and Eq. (36) location of the barrier."""
    s = np.asarray(sigma_t, float)
    return 2.0 - A ** (2.0 / 3.0) * s ** (-1.0 / 3.0), 2.0 - 0.5 * (2.0 - math.sqrt(3.0)) * A ** (2.0 / 3.0) * s ** (-1.0 / 3.0)


def high_tension_barrier(sigma_t, A=PAPER['A']):
    """Eq. (37)."""
    s = np.asarray(sigma_t, float)
    return 0.75 * (2.0 * math.sqrt(3.0) - 3.0) * A ** (4.0 / 3.0) * s ** (1.0 / 3.0)


# =================================================================================================
# 3.  One monotone branch  +  the full solution curve
# =================================================================================================
class Branch:
    """A solution branch single-valued in z (z strictly increasing).

    F(z): cubic Hermite spline through the ODE energies with the exact slopes F' = (1-pd)^2 - 4 - 2 sigma z
    pd(z): cubic spline of the contact curvature psi_dot_0 (a = 1)
    """

    def __init__(self, sigma, z, pd, F):
        self.sigma = float(sigma)
        z = np.asarray(z, float)
        # truncate at the first non-increase (fold)
        bad = np.where(np.diff(z) <= 0)[0]
        n = bad[0] + 1 if len(bad) else len(z)
        self.z = z[:n]
        self.pd = np.asarray(pd, float)[:n]
        self.F = np.asarray(F, float)[:n]
        self.w = (1.0 - self.pd) ** 2
        self.Fp = self.w - 4.0 - 2.0 * self.sigma * self.z
        self.Fh = CubicHermiteSpline(self.z, self.F, self.Fp)
        self.pds = CubicSpline(self.z, self.pd)
        self.zmin, self.zmax = self.z[0], self.z[-1]

    def F_of(self, z):
        return self.Fh(z)

    def pd_of(self, z):
        return self.pds(z)

    def w_of(self, z):
        return (1.0 - self.pds(z)) ** 2

    def w_max(self):
        i = int(np.argmax(self.w))
        lo = self.z[max(i - 1, 0)]
        hi = self.z[min(i + 1, len(self.z) - 1)]
        r = minimize_scalar(lambda z: -self.w_of(z), bounds=(lo, hi), method='bounded', options=dict(xatol=1e-13))
        return r.x, -r.fun

    def w_min_after(self, zfrom):
        """Interior local minimum of w(z) for z > zfrom (spinodal S2), else (None, None)."""
        sel = self.z > zfrom
        if sel.sum() < 3:
            return None, None
        zz, ww = self.z[sel], self.w[sel]
        i = int(np.argmin(ww))
        if i == 0 or i == len(zz) - 1:
            return None, None
        r = minimize_scalar(lambda z: self.w_of(z), bounds=(zz[i - 1], zz[i + 1]), method='bounded',
                            options=dict(xatol=1e-13))
        return r.x, r.fun


class SolutionCurve:
    """Everything about the free-membrane problem at one sigma~ (branch 1 [+ branch 3] + analysis)."""

    def __init__(self, sigma, pts1, pts3=None, meta=None):
        self.sigma = float(sigma)
        self.meta = dict(meta or {})
        self.pts1 = np.asarray(pts1)
        self.b1 = Branch(sigma, self.pts1[:, 2], self.pts1[:, 3], self.pts1[:, 4])
        self.folded = bool(len(self.b1.z) < len(self.pts1))     # points beyond the first z-fold were traced
        self.b3 = None
        if pts3 is not None and len(pts3) > 3:
            p3 = np.asarray(pts3)
            up = np.where(np.diff(p3[:, 2]) >= 0)[0]               # first point where z stops decreasing = fold
            if len(up):
                p3 = p3[:up[0] + 1]
            p = p3[::-1]                                         # reverse trace: z decreasing -> reverse
            self.b3 = Branch(sigma, p[:, 2], p[:, 3], p[:, 4])
            if len(self.b3.z) < 4:
                self.b3 = None
        self._bridge = None
        self.z_cross = None
        self._setup_envelope()

    # --- envelope (lowest energy branch) ----------------------------------------------------
    def _setup_envelope(self):
        s = self.sigma
        b1, b3 = self.b1, self.b3
        self.zx = None
        if b3 is not None:
            lo = max(b1.zmin, b3.zmin)
            hi = min(b1.zmax, b3.zmax)
            if hi > lo:
                f = lambda z: float(b1.F_of(z) - b3.F_of(z))      # > 0: branch 3 has the lower energy
                flo, fhi = f(lo), f(hi)
                if flo * fhi < 0:
                    self.zx = brentq(f, lo, hi, xtol=1e-14, rtol=1e-14)
                elif flo > 0 and fhi > 0:
                    self.zx = lo              # branch 3 lower on the whole overlap
                else:
                    self.zx = hi              # branch 1 lower on the whole overlap
            else:
                self.zx = b1.zmax
            tail = b3
        else:
            tail = b1
        # bridge from the end of the upper-most branch to F(2) = 0, F'(2) = -4 sigma
        ze = tail.zmax
        Fe = float(tail.F_of(ze))
        Fpe = float(tail.w_of(ze) - 4.0 - 2.0 * s * ze)
        self._bridge = CubicHermiteSpline([ze, 2.0], [Fe, 0.0], [Fpe, -4.0 * s])
        self._zbridge = ze
        self.bridged_from_fold = (b3 is None) and (ze < Z_END - 1e-3)

    def _masks(self, z):
        b1, b3 = self.b1, self.b3
        small = z < b1.zmin
        if b3 is not None:
            m1 = (~small) & (z <= self.zx)
            m3 = (z > self.zx) & (z <= b3.zmax)
            mb = z > b3.zmax
        else:
            m1 = (~small) & (z <= b1.zmax)
            m3 = np.zeros(z.shape, bool)
            mb = z > b1.zmax
        return small, m1, m3, mb

    def F(self, z):
        """Lowest-energy E~_free(z) (vectorised); 0 at z <= 0 and z >= 2."""
        z = np.atleast_1d(np.asarray(z, float))
        out = np.zeros(z.shape)
        inside = (z > 0.0) & (z < 2.0)
        small, m1, m3, mb = [m & inside for m in self._masks(z)]
        if small.any():                  # Eq. 21 form scaled to match the first numerical point
            z0 = self.b1.zmin
            out[small] = self.b1.F[0] * F_small_z(z[small], self.sigma) / F_small_z(z0, self.sigma)
        if m1.any():
            out[m1] = self.b1.F_of(z[m1])
        if m3.any():
            out[m3] = self.b3.F_of(z[m3])
        if mb.any():
            out[mb] = self._bridge(z[mb])
        return out

    def pd(self, z):
        """Contact curvature psi_dot_0 (a=1) of the lowest-energy shape."""
        z = np.atleast_1d(np.asarray(z, float))
        out = np.full(z.shape, -1.0)
        inside = (z > 0.0) & (z < 2.0)
        small, m1, m3, mb = [m & inside for m in self._masks(z)]
        if small.any():
            out[small] = -1.0 + (self.b1.pd[0] + 1.0) * z[small] / self.b1.zmin
        if m1.any():
            out[m1] = self.b1.pd_of(z[m1])
        if m3.any():
            out[m3] = self.b3.pd_of(z[m3])
        if mb.any():
            tail = self.b3 if self.b3 is not None else self.b1
            ze = self._zbridge
            out[mb] = tail.pd[-1] + (-1.0 - tail.pd[-1]) * (z[mb] - ze) / (2.0 - ze)
        return out

    def E(self, z, w):
        z = np.asarray(z, float)
        return -(w - 4.0) * z + self.sigma * z * z + self.F(z).reshape(z.shape)

    @staticmethod
    def Efull(w, sigma):
        return -2.0 * (w - 4.0) + 4.0 * sigma

    # --- phase lines --------------------------------------------------------------------------
    def spinodals(self):
        """S1 = max of w_eq(z) = (1-psi_dot_0)^2 on branch 1 (partially wrapped branch loses stability).

        S2 = smallest w~ for which a stable (E'' = dw_eq/dz > 0) near-neck equilibrium exists:
          * single-valued regime: interior local minimum of w_eq beyond S1 (NaN if none);
          * S-shaped regime: the stable neck branch is branch 3 from z = 2 downwards while psi_dot_0 < 1/a;
            it ends where psi_dot_0 = 1/a (w_eq = 0: S2 = 0, enveloped state stable for every w~ > 0) or at its fold
            (then S2 = w_eq at the fold).  NaN if branch 3 is not available."""
        zS1, wS1 = self.b1.w_max()
        if not self.folded:
            zS2, wS2 = self.b1.w_min_after(zS1)
            wS2 = np.nan if wS2 is None else wS2
            zS2 = np.nan if zS2 is None else zS2
        elif self.b3 is not None:
            p, z = self.b3.pd, self.b3.z
            idx = np.where(p >= 1.0)[0]
            if len(idx):                 # first point (from the neck, i.e. descending z) with psi_dot_0 >= 1
                j = idx[-1]              # largest z with p >= 1
                zS2 = float(z[j] + (z[min(j + 1, len(z) - 1)] - z[j]) * (p[j] - 1.0) / max(p[j] - p[min(j + 1, len(z) - 1)], 1e-300)) \
                    if j + 1 < len(z) else float(z[j])
                wS2 = 0.0
            else:
                zS2, wS2 = float(z[0]), float(self.b3.w[0])
        else:
            zS2, wS2 = np.nan, np.nan
        return dict(S1=wS1, zS1=zS1, S2=wS2, zS2=zS2)

    def wE(self):
        """E line, barrier, penetrations.  Returns dict(wE, zp, zb, barrier, zS1, S1)."""
        b1, s = self.b1, self.sigma
        zS1, wS1 = b1.w_max()

        def zp_of(w):
            return brentq(lambda z: float(b1.w_of(z)) - w, b1.zmin, zS1, xtol=1e-14, rtol=1e-14)

        def D(w):
            zp = zp_of(w)
            return float(-(w - 4.0) * zp + s * zp * zp + b1.F_of(zp)) - self.Efull(w, s)
        lo = max(4.0, float(b1.w[0])) + 1e-13
        hi = wS1 * (1.0 - 1e-10)
        wE = brentq(D, lo, hi, xtol=1e-13, rtol=1e-14)
        zp = zp_of(wE)
        zb = brentq(lambda z: float(b1.w_of(z)) - wE, zS1, b1.zmax, xtol=1e-14, rtol=1e-14)
        # barrier from the integral form E(zb) - E(2) = int_zb^2 (w - w_eq) dz is cancellation-free at large
        # sigma~, but needs the full w_eq(z) up to 2; the direct difference is used here and cross-checked
        Eb = float(-(wE - 4.0) * zb + s * zb * zb + b1.F_of(zb))
        barrier = Eb - self.Efull(wE, s)
        return dict(wE=wE, zp=zp, zb=zb, barrier=barrier, zS1=zS1, S1=wS1)

    def scan_check(self, res=None, n=20001):
        """Brute force verification on a fine z grid using the lowest-energy envelope: at w = w_E the global
        minimum of E(z) over (0,2] must coincide with E(2) (and be attained both at z_p and at 2) and the
        maximum between must reproduce the barrier.  Returns dict of discrepancies."""
        res = res or self.wE()
        w, s = res['wE'], self.sigma
        z = np.linspace(res['zp'] * 0.2, 2.0, n)
        z = np.unique(np.concatenate([z, self.b1.z[self.b1.z > res['zp'] * 0.2]]))
        E = self.E(z, w) - self.Efull(w, s)
        i = int(np.argmin(E))
        sel = z >= res['zp']
        scale = max(abs(w - 4.0) * 2.0, 4.0 * s, 1.0)
        return dict(min_E_minus_Efull=float(E.min()) / scale, z_min=float(z[i]),
                    barrier_scan=float(E[sel].max()), barrier_err=float(E[sel].max() - res['barrier']))

    # --- (de)serialisation ------------------------------------------------------------------------
    def summary(self):
        r = self.wE()
        sp = self.spinodals()
        r.update(sp)
        r['sigma'] = self.sigma
        r['s_shaped'] = bool(self.folded)
        r['z_fold'] = float(self.b1.zmax) if self.folded else np.nan
        r['z_cross'] = float(self.zx) if self.zx is not None else np.nan
        r['branch3'] = self.b3 is not None
        return r


# =================================================================================================
# 4.  Solving one sigma~ (used by the table generator and by free_energy())
# =================================================================================================
def compute_curve(sigma_t, branch3=None, profile=None, verbose=False):
    """Trace the solution curve(s) at one sigma~ and return a SolutionCurve."""
    sigma_t = float(sigma_t)
    prof = dict(PROFILE if profile is None else profile)
    sh = Shooter(sigma_t)
    t0 = time.time()
    attempts = [dict(), dict(tol_acc=1e-6, delta_rel=1e-6),
                dict(tol_acc=1e-6, delta_rel=3e-7, dzmax=0.5 * prof['dzmax'], dpmax=0.7 * prof['dpmax'])]
    pts1, info1 = None, None
    for att in attempts:
        kw = dict(prof)
        kw.update(att)
        try:
            pts1, info1 = trace(sh, stop_at_fold=True, **kw)
        except Exception as e:           # pragma: no cover
            pts1, info1 = None, dict(ok=False, reason=repr(e))
        if info1['ok']:
            break
    if pts1 is None or not info1['ok']:
        raise RuntimeError(f"forward trace failed at sigma~={sigma_t}: {info1}")
    meta = dict(ntraj=sh.ntraj, t_fwd=time.time() - t0)
    folded = bool(np.any(np.diff(pts1[:, 2]) <= 0))
    pts3 = None
    if branch3 is None:
        branch3 = sigma_t <= SIGMA_B3_MAX
    if folded and branch3:
        try:
            b1tmp = Branch(sigma_t, pts1[:, 2], pts1[:, 3], pts1[:, 4])
            sliver = 2.0 - b1tmp.zmax
            z_stop = b1tmp.zmax - 2.0 * sliver - 5e-3
            st = neck_start(sh)
            if st is not None and st[1] is not None:
                for att in attempts:
                    kw = dict(prof)
                    kw.update(att)
                    try:
                        pts3, info3 = trace(sh, start=st, orient_down=True, z_stop=z_stop, stop_at_fold=True, **kw)
                    except Exception as e:      # pragma: no cover
                        pts3, info3 = None, dict(ok=False, reason=repr(e))
                    if pts3 is not None and len(pts3) > 5 and (info3['ok'] or pts3[:, 2].min() < b1tmp.zmax - 1e-3):
                        break
                    pts3 = None
        except Exception as e:                  # pragma: no cover
            if verbose:
                print('branch 3 failed:', repr(e))
            pts3 = None
    meta['t_total'] = time.time() - t0
    meta['ntraj'] = sh.ntraj
    return SolutionCurve(sigma_t, pts1, pts3, meta)


_CACHE = {}


def solve_curve(sigma_t, use_table=True, verbose=False):
    """SolutionCurve at sigma~ (memoised; table nodes are rebuilt from the stored curve points, other values
    are solved from scratch, which takes ~10-60 s)."""
    key = float(sigma_t)
    if key in _CACHE:
        return _CACHE[key]
    cur = None
    if use_table and os.path.exists(TABLE_PATH):
        tab = load_table()
        j = np.where(np.abs(np.log(tab['sigma_t'] / key)) < 1e-9)[0]
        if len(j):
            cur = _curve_from_table(tab, int(j[0]))
    if cur is None:
        cur = compute_curve(key, verbose=verbose)
    _CACHE[key] = cur
    return cur


# =================================================================================================
# 5.  Table
# =================================================================================================
TABLE_COLUMNS = ['sigma_t', 'w_E', 'w_over_sigma', 'a_over_lambda', 'barrier', 'z_partial', 'z_barrier',
                 'S1', 'S2', 'zS1', 'zS2', 'wE_minus4_over_sigma', 's_shaped', 'z_fold', 'z_cross',
                 'branch3_exact', 'asymptotic_fallback', 'scan_err', 'barrier_scan_err', 'npts']
_TABLE = {}


def compute_point(sigma_t, verbose=False):
    """Full computation at one sigma~; returns dict with the table row and the curve arrays."""
    cur = compute_curve(sigma_t, verbose=verbose)
    sm = cur.summary()
    chk = cur.scan_check(sm)
    s = float(sigma_t)
    row = dict(sigma_t=s, w_E=sm['wE'], w_over_sigma=sm['wE'] / (2.0 * s), a_over_lambda=math.sqrt(s),
               barrier=sm['barrier'], z_partial=sm['zp'], z_barrier=sm['zb'], S1=sm['S1'], S2=sm['S2'],
               zS1=sm['zS1'], zS2=sm['zS2'], wE_minus4_over_sigma=(sm['wE'] - 4.0) / s,
               s_shaped=float(sm['s_shaped']), z_fold=sm['z_fold'], z_cross=sm['z_cross'],
               branch3_exact=float(sm['branch3']), asymptotic_fallback=0.0,
               scan_err=chk['min_E_minus_Efull'], barrier_scan_err=chk['barrier_err'] / max(sm['barrier'], 1e-300),
               npts=float(len(cur.b1.z)))
    b3 = cur.b3
    arrays = dict(b1=np.column_stack([cur.b1.z, cur.b1.pd, cur.b1.F]),
                  b3=np.column_stack([b3.z, b3.pd, b3.F]) if b3 is not None else np.zeros((0, 3)))
    return dict(row=row, arrays=arrays, meta=cur.meta)


def _worker(s):
    try:
        return s, compute_point(s)
    except Exception as e:           # pragma: no cover
        return s, dict(error=repr(e))


def generate_table(path=None, sigmas=None, nproc=4, verbose=True, per_decade=20, lo=-4.0, hi=6.0, parts_dir=None):
    """Compute the table (default: sigma~ = 1e-4 ... 1e6, 20 points per decade) with a process pool.
    Each finished sigma~ is stored in parts_dir (resumable); the final npz is assembled at the end."""
    import pickle
    from multiprocessing import Pool
    path = path or TABLE_PATH
    parts_dir = parts_dir or os.path.join(os.path.dirname(path), 'table_parts')
    os.makedirs(parts_dir, exist_ok=True)
    if sigmas is None:
        sigmas = 10.0 ** np.linspace(lo, hi, int(round((hi - lo) * per_decade)) + 1)
    sigmas = [float(s) for s in sigmas]
    fn = lambda s: os.path.join(parts_dir, 'sigma_%.12e.pkl' % s)
    t0 = time.time()
    results = {}
    todo = []
    for s in sigmas:
        if os.path.exists(fn(s)):
            with open(fn(s), 'rb') as f:
                results[s] = pickle.load(f)
        else:
            todo.append(s)
    with Pool(nproc) as pool:
        for k, (s, r) in enumerate(pool.imap_unordered(_worker, todo, chunksize=1)):
            results[s] = r
            if 'error' not in r:
                with open(fn(s), 'wb') as f:
                    pickle.dump(r, f)
            if verbose:
                if 'error' in r:
                    print(f"[{k+1}/{len(todo)}] sigma~={s:g} FAILED {r['error']}", flush=True)
                else:
                    rw = r['row']
                    print(f"[{k+1}/{len(todo)}] sigma~={s:10.4g} wE={rw['w_E']:.8g} bar={rw['barrier']:.6g} "
                          f"S1={rw['S1']:.6g} S2={rw['S2']:.6g} t={r['meta']['t_total']:.0f}s "
                          f"elapsed={time.time()-t0:.0f}s", flush=True)
    ok = [s for s in sigmas if 'error' not in results[s]]
    failed = [s for s in sigmas if 'error' in results[s]]
    save_table(path, [results[s] for s in ok])
    return ok, failed


def save_table(path, results):
    results = sorted(results, key=lambda r: r['row']['sigma_t'])
    out = {c: np.array([r['row'][c] for r in results], float) for c in TABLE_COLUMNS}
    for name in ('b1', 'b3'):
        arrs = [r['arrays'][name] for r in results]
        off = np.concatenate(([0], np.cumsum([len(a) for a in arrs]))).astype(np.int64)
        out[name + '_off'] = off
        out[name + '_z'] = np.concatenate([a[:, 0] for a in arrs]) if len(arrs) else np.zeros(0)
        out[name + '_pd'] = np.concatenate([a[:, 1] for a in arrs])
        out[name + '_F'] = np.concatenate([a[:, 2] for a in arrs])
    out['columns'] = np.array(TABLE_COLUMNS)
    np.savez_compressed(path, **out)


def load_table(path=None, reload=False):
    """Load the precomputed table as a dict of numpy arrays (see TABLE_COLUMNS)."""
    path = path or TABLE_PATH
    if path in _TABLE and not reload:
        return _TABLE[path]
    with np.load(path) as d:
        tab = {k: d[k] for k in d.files}
    _TABLE[path] = tab
    return tab


def export_csv(path=None, table_path=None):
    """Write the scalar table columns as plain CSV (no numpy needed to read it)."""
    tab = load_table(table_path)
    path = path or os.path.join(os.path.dirname(table_path or TABLE_PATH), 'deserno_table.csv')
    with open(path, 'w') as f:
        f.write('# Deserno (2003) phase lines; sigma~ = (a/lambda)^2, w~ = 2 w a^2/kappa, energies in units of pi*kappa.\n')
        f.write('# S2 = NaN where not computed (sigma~ > 1e4); asymptotic_fallback = 1 would flag a table point replaced by the paper asymptotics (none).\n')
        f.write(','.join(TABLE_COLUMNS) + '\n')
        for j in range(len(tab['sigma_t'])):
            f.write(','.join('%.10g' % tab[c][j] for c in TABLE_COLUMNS) + '\n')
    return path


def _curve_from_table(tab, j):
    s = float(tab['sigma_t'][j])
    o1 = tab['b1_off']
    sl = slice(o1[j], o1[j + 1])
    p1 = np.column_stack([np.zeros(sl.stop - sl.start), np.zeros(sl.stop - sl.start), tab['b1_z'][sl],
                          tab['b1_pd'][sl], tab['b1_F'][sl]] + [np.zeros(sl.stop - sl.start)] * 5)
    o3 = tab['b3_off']
    sl3 = slice(o3[j], o3[j + 1])
    p3 = None
    if sl3.stop > sl3.start:
        n3 = sl3.stop - sl3.start
        p3 = np.column_stack([np.zeros(n3), np.zeros(n3), tab['b3_z'][sl3], tab['b3_pd'][sl3], tab['b3_F'][sl3]]
                             + [np.zeros(n3)] * 5)[::-1]     # stored ascending; SolutionCurve expects reverse-trace order
    cur = SolutionCurve(s, p1, p3)
    cur.folded = bool(tab['s_shaped'][j])
    return cur


def rebuild_rows(tab=None):
    """Recompute every derived column of the table from the stored curve arrays (no tracing).  Returns a dict of
    columns; used for self-consistency checks and to refresh the table after analysis changes."""
    tab = tab or load_table()
    cols = {c: np.zeros(len(tab['sigma_t'])) for c in TABLE_COLUMNS}
    for j in range(len(tab['sigma_t'])):
        cur = _curve_from_table(tab, j)
        sm = cur.summary()
        chk = cur.scan_check(sm)
        s = cur.sigma
        row = dict(sigma_t=s, w_E=sm['wE'], w_over_sigma=sm['wE'] / (2.0 * s), a_over_lambda=math.sqrt(s),
                   barrier=sm['barrier'], z_partial=sm['zp'], z_barrier=sm['zb'], S1=sm['S1'], S2=sm['S2'],
                   zS1=sm['zS1'], zS2=sm['zS2'], wE_minus4_over_sigma=(sm['wE'] - 4.0) / s,
                   s_shaped=float(sm['s_shaped']), z_fold=sm['z_fold'], z_cross=sm['z_cross'],
                   branch3_exact=float(sm['branch3']), asymptotic_fallback=float(tab['asymptotic_fallback'][j]),
                   scan_err=chk['min_E_minus_Efull'], barrier_scan_err=chk['barrier_err'] / max(sm['barrier'], 1e-300),
                   npts=float(len(cur.b1.z)))
        for c in TABLE_COLUMNS:
            cols[c][j] = row[c]
    return cols


def refresh_table(path=None):
    """Rewrite the npz with all derived columns recomputed from the stored curves."""
    path = path or TABLE_PATH
    tab = dict(load_table(path, reload=True))
    cols = rebuild_rows(tab)
    for c in TABLE_COLUMNS:
        tab[c] = cols[c]
    np.savez_compressed(path, **tab)
    _TABLE.pop(path, None)
    return cols


# --- interpolation in log sigma~ ------------------------------------------------------------------
def _spline_logsigma(tab, ydata, valid=None):
    x = np.log10(tab['sigma_t'])
    if valid is None:
        valid = np.isfinite(ydata)
    return CubicSpline(x[valid], ydata[valid]), x[valid][0], x[valid][-1]


def interp_logsigma(sigma_t, column, table=None, kind='plain'):
    """Interpolate a table column in log10(sigma~) with a cubic spline.

    kind: 'plain'      -> the column itself
          'log'        -> spline of ln(column)             (positive, power-law like: barrier)
          'w_over_sig' -> column = w~; spline of (w~ - 4)/sigma~       (E line, S1)
          'w4_minus'   -> column = w~; spline of (4 - w~)/sigma~       (S2)
    Returns (values, inside) where inside marks sigma~ within the (valid) table range; outside the values are
    NaN (the callers w_E(), barrier(), ... substitute the asymptotic formulas, flagged)."""
    tab = table or load_table()
    s = np.atleast_1d(np.asarray(sigma_t, float))
    y = np.asarray(tab[column], float)
    sig = tab['sigma_t']
    if kind == 'log':
        yy = np.log(y)
    elif kind == 'w_over_sig':
        yy = (y - 4.0) / sig
    elif kind == 'w4_minus':
        yy = (4.0 - y) / sig
    else:
        yy = y
    cs, xlo, xhi = _spline_logsigma(tab, yy)
    x = np.log10(s)
    inside = (x >= xlo - 1e-12) & (x <= xhi + 1e-12)
    v = np.full(s.shape, np.nan)
    v[inside] = cs(np.clip(x[inside], xlo, xhi))
    if kind == 'log':
        v = np.exp(v)
    elif kind == 'w_over_sig':
        v = 4.0 + v * s
    elif kind == 'w4_minus':
        v = 4.0 - v * s
    return v, inside


def _warn_outside(name, s):
    warnings.warn(f"{name}: sigma~ outside the tabulated range, using asymptotic extrapolation", RuntimeWarning,
                  stacklevel=3)


def w_E(sigma_t, return_flag=False):
    """E line w~_E(sigma~) (partially wrapped -> fully enveloped), from the table (cubic spline in log sigma~)."""
    tab = load_table()
    s = np.atleast_1d(np.asarray(sigma_t, float))
    v, inside = interp_logsigma(s, 'w_E', tab, 'w_over_sig')
    if not inside.all():
        _warn_outside('w_E', s[~inside])
        lo, hi = tab['sigma_t'][0], tab['sigma_t'][-1]
        for i in np.where(~inside)[0]:
            if s[i] < lo:
                g_tab = (tab['w_E'][0] - 4.0) / lo
                g26 = (w_E_small_gradient(lo)[0] - 4.0) / lo
                g = (w_E_small_gradient(s[i])[0] - 4.0) / s[i] * g_tab / g26
                v[i] = 4.0 + s[i] * g if np.isfinite(g) else 4.0 + 2.0 * s[i]
            else:
                g_tab = (tab['w_E'][-1] - 4.0) / hi
                g35 = (high_tension_wE(hi) - 4.0) / hi
                v[i] = 4.0 + s[i] * (high_tension_wE(s[i]) - 4.0) / s[i] * g_tab / g35
    return (v, ~inside) if return_flag else (v if np.ndim(sigma_t) else float(v[0]))


def barrier(sigma_t):
    """Energy barrier E~_barrier = E_barrier/(pi kappa) at the E line."""
    tab = load_table()
    s = np.atleast_1d(np.asarray(sigma_t, float))
    v, inside = interp_logsigma(s, 'barrier', tab, 'log')
    if not inside.all():
        _warn_outside('barrier', s[~inside])
        lo, hi = tab['sigma_t'][0], tab['sigma_t'][-1]
        for i in np.where(~inside)[0]:
            if s[i] > hi:
                v[i] = tab['barrier'][-1] * float(high_tension_barrier(s[i]) / high_tension_barrier(hi))
            else:        # power law continued with the local log-slope
                sl = np.log(tab['barrier'][1] / tab['barrier'][0]) / np.log(tab['sigma_t'][1] / tab['sigma_t'][0])
                v[i] = tab['barrier'][0] * (s[i] / lo) ** sl
    return v if np.ndim(sigma_t) else float(v[0])


def spinodal_S1(sigma_t):
    """w~ at which the partially wrapped branch loses stability (max of (1-psi_dot_0)^2)."""
    tab = load_table()
    s = np.atleast_1d(np.asarray(sigma_t, float))
    v, inside = interp_logsigma(s, 'S1', tab, 'w_over_sig')
    return v if np.ndim(sigma_t) else float(v[0])


def spinodal_S2(sigma_t):
    """w~ at which the enveloped (neck) state loses stability.  Single-valued regime (sigma~ < 4.72): local minimum
    of (1-psi_dot_0)^2 near the neck; S-shaped regime: where the neck branch reaches psi_dot_0 = 1/a (S2 = 0) or
    its fold.  Interpolated in log sigma~ (cubic spline of (4 - S2)/sigma~ over the nodes with sigma~ < 4.72, linear
    interpolation of S2 above)."""
    tab = load_table()
    s = np.atleast_1d(np.asarray(sigma_t, float))
    S2 = np.asarray(tab['S2'], float)
    sg = tab['sigma_t']
    x = np.log10(sg)
    xs = np.log10(s)
    v = np.full(s.shape, np.nan)
    low = np.isfinite(S2) & (tab['s_shaped'] < 0.5)
    cs = CubicSpline(x[low], ((4.0 - S2) / sg)[low])
    m = (xs >= x[low][0]) & (xs <= x[low][-1])
    v[m] = 4.0 - cs(xs[m]) * s[m]
    hi = np.isfinite(S2) & (tab['s_shaped'] > 0.5)
    if hi.sum() > 1:
        m2 = (xs > x[low][-1]) & (xs >= x[hi][0]) & (xs <= x[hi][-1])
        v[m2] = np.interp(xs[m2], x[hi], S2[hi])
    return v if np.ndim(sigma_t) else float(v[0])


def z_partial(sigma_t):
    """Penetration of the partially wrapped state at the E transition."""
    v, _ = interp_logsigma(sigma_t, 'z_partial', None, 'plain')
    return v if np.ndim(sigma_t) else float(v[0])


def z_barrier(sigma_t):
    v, _ = interp_logsigma(sigma_t, 'z_barrier', None, 'plain')
    return v if np.ndim(sigma_t) else float(v[0])


def fig5_curves(a_over_lambda):
    """Paper Fig. 5: (w/sigma on E, w/sigma on W) as functions of a/lambda = sqrt(sigma~).
    w/sigma = w~/(2 sigma~);  W line w~ = 4  ->  w/sigma = 2 (lambda/a)^2."""
    a = np.asarray(a_over_lambda, float)
    s = a ** 2
    wE_ = np.asarray(w_E(np.atleast_1d(s)))
    return (wE_ / (2.0 * np.atleast_1d(s))).reshape(a.shape), 2.0 / s


def fig2_line_data(ylim=(0.0, 1.0), n=400):
    """Numerical lines of the paper's Fig. 2 in the (w~, sigma~) plane for sigma~ in ylim (from the table).
    Returns dict with keys 'E', 'S1', 'S2'; each value is a tuple (w_tilde, sigma_tilde) of arrays (NaN-free).
    The lines start at the triple point T = (4, 0) when ylim[0] <= 0."""
    y0, y1 = float(ylim[0]), float(ylim[1])
    smin = 1e-4
    s = np.geomspace(max(y0, smin), y1, n)
    out = {}
    for key, f in (('E', w_E), ('S1', spinodal_S1), ('S2', spinodal_S2)):
        w = np.asarray(f(s), float)
        ok = np.isfinite(w)
        ss, ww = s[ok], w[ok]
        if y0 <= 0.0:                                   # all three lines emanate from the triple point T = (4, 0)
            ss = np.concatenate(([0.0], ss))
            ww = np.concatenate(([4.0], ww))
        out[key] = (ww, ss)
    return out


def plot_fig2_lines(ax, ylim=(0.0, 1.0), color='k', lw_E=2.2, lw_S=1.2, label_prefix='', n=400, **kw):
    """Draw the numerical phase lines of Deserno's Fig. 2 on a matplotlib axes (x = w~, y = sigma~):
    the discontinuous envelopment line E (solid, bold) and the spinodals S1, S2 (short dashed).
    The W line (w~ = 4) and the zero-energy line w~ = 4 + 2 sigma~ are NOT drawn (trivial; add them yourself).
    Does not change axis limits (set them in the caller, e.g. ax.set_xlim(3, 6); ax.set_ylim(*ylim)).
    Returns {'E': Line2D, 'S1': Line2D, 'S2': Line2D}.  matplotlib is imported lazily."""
    d = fig2_line_data(ylim, n)
    artists = {}
    for key, ls, lw in (('E', '-', lw_E), ('S1', (0, (3.0, 2.0)), lw_S), ('S2', (0, (3.0, 2.0)), lw_S)):
        w, s = d[key]
        lab = (label_prefix + key.replace('S1', r'S$_1$').replace('S2', r'S$_2$')) if label_prefix else None
        artists[key], = ax.plot(w, s, ls=ls, lw=lw, color=color, label=lab, **kw)
    return artists


# --- free energy etc. (solve on demand) -------------------------------------------------------------
def free_energy(z, sigma_t):
    """Lowest-energy free-membrane energy E~_free(z; sigma~) (array in z).  See module docstring for the
    handling of the S-shaped regime.  The first call at a sigma~ that is not a table node solves the shape
    equations (about 10-60 s), later calls are instantaneous."""
    r = solve_curve(sigma_t).F(z)
    return r if np.ndim(z) else float(r[0])


def contact_curvature(z, sigma_t):
    """psi_dot_0 (a = 1) of the lowest-energy shape;  w~ = (1 - psi_dot_0)^2 is the adhesion at which z is an
    equilibrium."""
    r = solve_curve(sigma_t).pd(z)
    return r if np.ndim(z) else float(r[0])


def total_energy(z, w_t, sigma_t):
    """E~(z) = -(w~ - 4) z + sigma~ z^2 + E~_free(z; sigma~)."""
    z_ = np.asarray(z, float)
    r = -(w_t - 4.0) * z_ + sigma_t * z_ ** 2 + np.asarray(free_energy(z_, sigma_t))
    return r if np.ndim(z) else float(r)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'table':
        nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        ok, failed = generate_table(nproc=nproc)
        print('done; failed sigma~:', failed)
    else:
        print(__doc__)
