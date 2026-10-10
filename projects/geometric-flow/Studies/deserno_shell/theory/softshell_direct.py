#!/usr/bin/env python
"""
softshell_direct.py -- continuum limit of the energy the simulation minimises (bead + soft adhesion shell +
bending + excess tension), axisymmetric, solved by DIRECT energy minimisation (workstream T1).
Derivation and all approximations: softshell_direct.md.   Tests: test_softshell_direct.py.

UNITS: a = 1, kappa = KB/2, energies in units of pi*kappa, w~ = 2 w a^2/kappa = 4 KI a^2/KB,
sigma~ = sigma a^2/kappa = 2 KA a^2/KB.   Bead centre at the origin, axis = h, rho = distance from the axis.

SHAPE: a polyline with N segments (nodes 0..N), segment i has constant chord angle psi_i (tangent (cos psi, sin psi)
in the (rho, h) half plane), length L_i = S*ell_i with FIXED fractions ell_i (sum 1).  The pole is on the axis,
rho_0 = 0, h_0 = h_p; the outer end is clamped at rho_N = R (R = disk radius of the simulation) which fixes the total
length S = R / sum(ell_i cos psi_i) EXACTLY, so the unknowns are z = (psi_0 .. psi_{N-1}, h_p).  The normal is the left
normal n = (-sin psi, cos psi) and points to the bead on a cap that hugs the bead from below (psi = polar angle from the
bottom pole, c1 = psi' = +1, c2 = sin psi/rho = +1 on a unit sphere).  The free end is a natural boundary (psi free,
h free <=> bead free to move along the axis).

ENERGY / (pi kappa)        E = Ebend + sigma~ * Aex - w~ * G
    Ebend = sum_{i=1}^{N-1} rho_i lbar_i (c1_i + c2_i)^2,  lbar_i = (L_{i-1} + L_i)/2,  c1_i = (psi_i - psi_{i-1})/lbar_i,
            c2_i = sin(psibar_i)/rho_i,  psibar_i = (L_i psi_{i-1} + L_{i-1} psi_i)/(L_{i-1} + L_i)     (2nd order)
    Aex   = (A - pi R^2)/pi = sum_i (rho_i + rho_{i+1}) L_i (1 - cos psi_i)       (exact frustum area minus projection)
    G     = (1/2pi) int W((r-1)/s) |n.rhat|_{n.rhat<0} dOmega  = sum_i int_{seg i} ds rho W(r) b_i^+ / r^3,
            b_i = rho sin psi - h cos psi (constant on a segment), r = |(rho, h)|,  W(x) = [(1+cos(pi x))/2]^p for |x| < 1.
    (G = z_ad of the simulation: -Bead/(2 pi KI);  full sphere -> G = 2;  s -> 0: G -> z = 1 - cos(alpha).)
Gradient and Hessian are EXACT in the sense of complex-step element derivatives chained analytically through the
cumulative-sum map psi -> positions and the closure S(psi); the Hessian (dense, N x N) is complete.

SOLVERS (class Branch / functions at the bottom)
    minimize_wt(state, wt)             damped Newton with Hessian shift: local minimum at fixed w~
    kkt_solve(state, Gamma)            constrained minimum at fixed coverage G (w~ = Lagrange multiplier = de/dG)
    trace_branch(state, G0, G1, dG)    continuation in G: the curve e(G), w~(G) = e'(G); S1 = max, S2 = min of w~(G)
    regrid(state)                      equidistribute nodes: dense where W changes along the curve and where it bends
    lines_from_branch(...)             S1, S2, E, coverage of the states
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import math
import time
import warnings
import numpy as np
import scipy.linalg as sla
from numpy.polynomial.legendre import leggauss

PI = math.pi
warnings.filterwarnings("ignore", category=sla.LinAlgWarning)


# =====================================================================================================
# 1. local element energies (vectorised, complex-step safe)
# =====================================================================================================
class Model:
    """Parameters of the continuum model.  s = shell half width (fraction of a), p = shell power, R = disk radius."""

    def __init__(self, sigma, s=0.25, p=1.0, R=10.0, ngl=6):
        self.sigma = float(sigma)
        self.s = float(s)
        self.p = float(p)
        self.R = float(R)
        self.ngl = int(ngl)
        self.xg, self.wg = leggauss(self.ngl)

    # ---- element energies (arguments are 1-d arrays, real or complex) -----------------------------------
    def e_tens(self, rho, psi, L):
        return (2.0 * rho + L * np.cos(psi)) * L * (2.0 * np.sin(0.5 * psi) ** 2)

    def e_bend(self, rho, psim, psi, Lm, L):
        lb = 0.5 * (Lm + L)
        c1 = (psi - psim) / lb
        pn = (L * psim + Lm * psi) / (Lm + L)
        c2 = np.sin(pn) / rho
        return rho * lb * (c1 + c2) ** 2

    def weight(self, r):
        """shell weight W((r-1)/s) (complex-step safe, selection by the real part)."""
        u = (r - 1.0) / self.s
        inside = np.abs(u.real) < 1.0
        base = 0.5 * (1.0 + np.cos(PI * u))
        return np.where(inside, base ** self.p, 0.0)

    def e_gam(self, rho, h, psi, L):
        c = np.cos(psi)
        sn = np.sin(psi)
        b = rho * sn - h * c
        t = 0.5 * L[:, None] * (self.xg[None, :] + 1.0)
        x = rho[:, None] + t * c[:, None]
        y = h[:, None] + t * sn[:, None]
        r = np.sqrt(x * x + y * y)
        f = x * self.weight(r) / r ** 3
        gam = 0.5 * L * (f * self.wg[None, :]).sum(axis=1)
        bp = np.where(b.real > 0.0, b, 0.0)
        return gam * bp


FAMILIES = {            # local variable layout v = (rho_i, h_i, psi_{i-1}, psi_i, L_{i-1}, L_i)
    "tens": (0, 3, 5),
    "bend": (0, 2, 3, 4, 5),
    "gam": (0, 1, 3, 5),
}


def _call(model, fam, cols):
    if fam == "tens":
        return model.e_tens(*cols)
    if fam == "bend":
        return model.e_bend(*cols)
    return model.e_gam(*cols)


def local_values(model, fam, V):
    idx = FAMILIES[fam]
    sl = slice(1, None) if fam == "bend" else slice(None)
    out = np.zeros(V.shape[0])
    out[sl] = _call(model, fam, [V[sl, k] for k in idx]).real
    if fam == "bend":
        out[0] += 2.0 * V[0, 3] ** 2          # pole cell [0, L_0/2]: int rho (c1+c2)^2 ds with c1 = c2 = 2 psi_0/L_0
    return out


def local_derivs(model, fam, V, hess=True):
    """gradient (M,6) by complex step, Hessian (M,6,6) by central differences of the complex-step gradient."""
    idx = FAMILIES[fam]
    nv = len(idx)
    sl = slice(1, None) if fam == "bend" else slice(None)
    Vs = V[sl]
    M = Vs.shape[0]
    hc = 1e-30

    def grad_at(Vr):
        gg = np.zeros((M, nv))
        for a, k in enumerate(idx):
            Vc = Vr.astype(complex)
            Vc[:, k] += 1j * hc
            gg[:, a] = _call(model, fam, [Vc[:, m] for m in idx]).imag / hc
        return gg

    g = np.zeros((V.shape[0], 6))
    H = np.zeros((V.shape[0], 6, 6))
    g0 = grad_at(Vs)
    for a, k in enumerate(idx):
        g[sl, k] = g0[:, a]
    if hess:
        Hl = np.zeros((M, nv, nv))
        for b, kb in enumerate(idx):
            d = 1e-6 * (1.0 + np.abs(Vs[:, kb]))
            Vp = Vs.copy()
            Vp[:, kb] += d
            Vm = Vs.copy()
            Vm[:, kb] -= d
            Hl[:, :, b] = (grad_at(Vp) - grad_at(Vm)) / (2.0 * d[:, None])
        Hl = 0.5 * (Hl + Hl.transpose(0, 2, 1))
        for a, ka in enumerate(idx):
            for b, kb in enumerate(idx):
                H[sl, ka, kb] = Hl[:, a, b]
    if fam == "bend":
        g[0, 3] += 4.0 * V[0, 3]
        if hess:
            H[0, 3, 3] += 4.0
    return g, H


# =====================================================================================================
# 2. geometry of a state and assembly of reduced derivatives
# =====================================================================================================
class State:
    """psi (N,), hp, ell (N,) fractions.  Lazily computed geometry."""

    def __init__(self, psi, hp, ell, R):
        self.psi = np.array(psi, float)
        self.hp = float(hp)
        self.ell = np.array(ell, float)
        self.R = float(R)
        self.N = len(self.psi)
        self._geo()

    def _geo(self):
        N = self.N
        self.c = np.cos(self.psi)
        self.sn = np.sin(self.psi)
        self.D = float((self.ell * self.c).sum())
        if self.D <= 1e-6:
            raise FloatingPointError("degenerate closure (D<=0)")
        self.S = self.R / self.D
        self.rho = self.S * np.concatenate([[0.0], np.cumsum(self.ell * self.c)])
        self.h = self.hp + self.S * np.concatenate([[0.0], np.cumsum(self.ell * self.sn)])
        self.L = self.S * self.ell
        V = np.empty((N, 6))
        V[:, 0] = self.rho[:N]
        V[:, 1] = self.h[:N]
        V[:, 2] = np.concatenate([[self.psi[0]], self.psi[:-1]])
        V[:, 3] = self.psi
        V[:, 4] = np.concatenate([[self.L[0]], self.L[:-1]])
        V[:, 5] = self.L
        self.V = V
        self.sv = self.S * self.ell * self.sn / self.D          # dS/dpsi_j

    def copy(self):
        return State(self.psi.copy(), self.hp, self.ell.copy(), self.R)

    def z(self):
        return np.concatenate([self.psi, [self.hp]])

    def with_z(self, z):
        return State(z[:-1], z[-1], self.ell, self.R)

    @property
    def rho_min(self):
        """neck radius: the smallest rho after the maximum of rho on the part of the curve within 4 bead radii
        (a nearly closed bud has a neck; otherwise this is simply rho at r = 4)"""
        near = np.hypot(self.rho, self.h) < 4.0
        idx = np.nonzero(near)[0]
        if len(idx) < 2:
            return float(self.rho[-1])
        k = idx[np.argmax(self.rho[idx])]
        return float(self.rho[k:idx[-1] + 1].min())

    @property
    def rho_min_all(self):
        return float(self.rho[1:].min())


def energies(model, st, wt):
    """(E, Ebend, Aex, G) at the state"""
    eb = local_values(model, "bend", st.V).sum()
    ea = local_values(model, "tens", st.V).sum()
    gm = local_values(model, "gam", st.V).sum()
    return eb + model.sigma * ea - wt * gm, eb, ea, gm


def assemble_grad(st, g):
    """reduced gradient (N+1,) from local gradients g (N,6)"""
    N, S, ell, c, sn = st.N, st.S, st.ell, st.c, st.sn
    g0, g1, g2, g3, g4, g5 = g.T
    cs0 = np.cumsum(g0[::-1])[::-1]
    cs1 = np.cumsum(g1[::-1])[::-1]
    Rr = np.concatenate([cs0[1:], [0.0]])
    Rh = np.concatenate([cs1[1:], [0.0]])
    Gpsi = -S * ell * sn * Rr + S * ell * c * Rh + g3
    Gpsi[:-1] += g2[1:]
    ellm = np.concatenate([[0.0], ell[:-1]])
    GS = ((g0 * st.rho[:N] + g1 * (st.h[:N] - st.hp)) / S).sum() + (g4 * ellm).sum() + (g5 * ell).sum()
    Ghp = g1.sum()
    out = np.empty(N + 1)
    out[:N] = Gpsi + GS * st.sv
    out[N] = Ghp
    return out


def assemble_hess(st, g, H):
    """reduced gradient (N+1,) and dense reduced Hessian (N+1, N+1)"""
    N, S, ell, c, sn, D = st.N, st.S, st.ell, st.c, st.sn, st.D
    iS, ih = N, N + 1
    g0, g1, g2, g3, g4, g5 = g.T
    cs0 = np.cumsum(g0[::-1])[::-1]
    cs1 = np.cumsum(g1[::-1])[::-1]
    Rr = np.concatenate([cs0[1:], [0.0]])
    Rh = np.concatenate([cs1[1:], [0.0]])
    ellm = np.concatenate([[0.0], ell[:-1]])
    # ---- extended gradient
    Gy = np.zeros(N + 2)
    Gy[:N] = -S * ell * sn * Rr + S * ell * c * Rh + g3
    Gy[:N - 1] += g2[1:]
    Gy[iS] = ((g0 * st.rho[:N] + g1 * (st.h[:N] - st.hp)) / S).sum() + (g4 * ellm).sum() + (g5 * ell).sum()
    Gy[ih] = g1.sum()
    # ---- dense rows
    a0 = -S * ell * sn
    a1 = S * ell * c
    Vr = np.zeros((N, N + 2))
    Vh = np.zeros((N, N + 2))
    Vr[:, :N] = np.tril(np.broadcast_to(a0, (N, N)), k=-1)
    Vh[:, :N] = np.tril(np.broadcast_to(a1, (N, N)), k=-1)
    Vr[:, iS] = st.rho[:N] / S
    Vh[:, iS] = (st.h[:N] - st.hp) / S
    Vh[:, ih] = 1.0
    Hl = lambda a, b: H[:, a, b]
    P = Hl(0, 0)[:, None] * Vr + Hl(0, 1)[:, None] * Vh
    Q = Hl(0, 1)[:, None] * Vr + Hl(1, 1)[:, None] * Vh
    Hy = Vr.T @ P + Vh.T @ Q
    # ---- dense x sparse cross terms
    K = np.zeros((N + 2, N + 2))
    K[:, :N] += (Vr * Hl(0, 3)[:, None]).T + (Vh * Hl(1, 3)[:, None]).T
    K[:, :N - 1] += (Vr[1:] * Hl(0, 2)[1:, None]).T + (Vh[1:] * Hl(1, 2)[1:, None]).T
    K[:, iS] += Vr.T @ (Hl(0, 4) * ellm + Hl(0, 5) * ell) + Vh.T @ (Hl(1, 4) * ellm + Hl(1, 5) * ell)
    Hy += K + K.T
    # ---- sparse x sparse
    idxN = np.arange(N)
    Hy[idxN, idxN] += Hl(3, 3)
    Hy[idxN[:-1], idxN[:-1]] += Hl(2, 2)[1:]
    off = Hl(2, 3)[1:]
    Hy[idxN[:-1], idxN[1:]] += off
    Hy[idxN[1:], idxN[:-1]] += off
    colS = np.zeros(N + 2)
    colS[:N - 1] += Hl(2, 4)[1:] * ellm[1:] + Hl(2, 5)[1:] * ell[1:]
    colS[:N] += Hl(3, 4) * ellm + Hl(3, 5) * ell
    Hy[:, iS] += colS
    Hy[iS, :] += colS
    Hy[iS, iS] += (Hl(4, 4) * ellm ** 2 + Hl(5, 5) * ell ** 2 + 2.0 * Hl(4, 5) * ellm * ell).sum()
    # ---- second derivatives of the positions
    Tpsi = -S * ell * (c * Rr + sn * Rh)
    TpS = ell * (-sn * Rr + c * Rh)
    Hy[idxN, idxN] += Tpsi
    Hy[idxN, iS] += TpS
    Hy[iS, idxN] += TpS
    # ---- reduction: S = R/D(psi)
    sv = st.sv
    # Hr = Z^T Hy Z with Z = [[I,0],[sv^T,0],[0,1]] written out (O(N^2))
    Hr = np.empty((N + 1, N + 1))
    cS = Hy[:N, iS]
    Hr[:N, :N] = Hy[:N, :N] + np.outer(cS, sv) + np.outer(sv, cS) + Hy[iS, iS] * np.outer(sv, sv)
    Hr[:N, N] = Hy[:N, ih] + sv * Hy[iS, ih]
    Hr[N, :N] = Hr[:N, N]
    Hr[N, N] = Hy[ih, ih]
    ls = ell * sn
    HS = np.diag(S * ell * c / D) + 2.0 * (S / D ** 2) * np.outer(ls, ls)
    Hr[:N, :N] += Gy[iS] * HS
    Hr = 0.5 * (Hr + Hr.T)
    gr = np.empty(N + 1)
    gr[:N] = Gy[:N] + Gy[iS] * sv
    gr[N] = Gy[ih]
    return gr, Hr


class Eval:
    """Cached local derivatives of the three element families at a state; combine for any (sigma, wt)."""

    def __init__(self, model, st, hess=True):
        self.model, self.st = model, st
        self.d = {fam: local_derivs(model, fam, st.V, hess=hess) for fam in ("bend", "tens", "gam")}
        self.hess = hess

    def _comb(self, wb, wt_, wg):
        g = wb * self.d["bend"][0] + wt_ * self.d["tens"][0] + wg * self.d["gam"][0]
        H = None
        if self.hess:
            H = wb * self.d["bend"][1] + wt_ * self.d["tens"][1] + wg * self.d["gam"][1]
        return g, H

    def grad(self, wt):
        g, _ = self._comb(1.0, self.model.sigma, -wt)
        return assemble_grad(self.st, g)

    def grad_gam(self):
        return assemble_grad(self.st, self.d["gam"][0])

    def grad_el(self):
        g, _ = self._comb(1.0, self.model.sigma, 0.0)
        return assemble_grad(self.st, g)

    def hessian(self, wt):
        g, H = self._comb(1.0, self.model.sigma, -wt)
        return assemble_hess(self.st, g, H)

    def hessian_el_gam(self, lam):
        """gradient and Hessian of E_el - lam*G (lam = w~)"""
        return self.hessian(lam)


# =====================================================================================================
# 3. grids: equidistribution of nodes along the curve, plate / sphere initial states
# =====================================================================================================
def dense_curve(st, M=24000):
    """resample the polyline of a state: arclength s, rho, h, psi (psi linearly interpolated between the segment mid points)"""
    sn_ = np.concatenate([[0.0], np.cumsum(st.L)])
    smid = 0.5 * (sn_[1:] + sn_[:-1])
    sd_ = np.linspace(0.0, sn_[-1], M)
    psi_d = np.interp(sd_, smid, st.psi)            # constant extrapolation at both ends
    ds = np.diff(sd_)
    pm = 0.5 * (psi_d[1:] + psi_d[:-1])
    rho_d = np.concatenate([[0.0], np.cumsum(ds * np.cos(pm))])
    h_d = st.hp + np.concatenate([[0.0], np.cumsum(ds * np.sin(pm))])
    return sd_, rho_d, h_d, psi_d


def monitor_spacing(model, sd_, rho_d, h_d, psi_d, h_base=0.4, df_tol=0.04, dpsi_tol=0.10, h_min=2e-4,
                    growth=0.3, neck_tol=0.25):
    """target node spacing h(s) on the dense curve (smooth: relative growth <= `growth`)"""
    ds = sd_[1] - sd_[0]
    r = np.hypot(rho_d, h_d)
    W = model.weight(r)
    q = np.maximum(rho_d * np.sin(psi_d) - h_d * np.cos(psi_d), 0.0) / np.maximum(r, 1e-12)
    f = rho_d * W * q / np.maximum(r, 1e-12) ** 2          # adhesion density per arclength  (G = int f ds)
    dfds = np.abs(np.gradient(f, ds))
    dpsi = np.abs(np.gradient(psi_d, ds))
    # smooth the density gradients lightly so that a single kink does not dominate
    mon = 1.0 / h_base + dfds / df_tol + dpsi / dpsi_tol
    # neck / small rho (away from the pole): c2 = sin(psi)/rho large
    away = sd_ > 0.15
    mon += np.where(away, np.abs(np.sin(psi_d)) / np.maximum(rho_d, 1e-3) / dpsi_tol, 0.0)
    mon += np.where(away, 1.0 / (neck_tol * np.maximum(rho_d, 1e-3)) * (rho_d < 0.5), 0.0)
    hloc = np.clip(1.0 / mon, h_min, h_base)
    # distance transform in spacing space (limit the relative growth of the spacing along the curve)
    for _ in range(2):
        for i in range(1, len(hloc)):
            hloc[i] = min(hloc[i], hloc[i - 1] + growth * ds)
        for i in range(len(hloc) - 2, -1, -1):
            hloc[i] = min(hloc[i], hloc[i + 1] + growth * ds)
    return hloc


def _dt_pass(h, ds, growth):
    # vectorised-ish distance transform: forward/backward minimum with slope `growth`
    n = len(h)
    idx = np.arange(n) * ds
    a = h - growth * idx
    b = np.minimum.accumulate(a)
    h1 = b + growth * idx
    a2 = h1 + growth * idx
    c = np.minimum.accumulate(a2[::-1])[::-1]
    return c - growth * idx


def monitor_spacing_fast(model, sd_, rho_d, h_d, psi_d, **kw):
    h_base = kw.get("h_base", 0.4)
    df_tol = kw.get("df_tol", 0.04)
    dpsi_tol = kw.get("dpsi_tol", 0.10)
    h_min = kw.get("h_min", 2e-4)
    growth = kw.get("growth", 0.3)
    neck_tol = kw.get("neck_tol", 0.25)
    ds = sd_[1] - sd_[0]
    r = np.hypot(rho_d, h_d)
    W = model.weight(r)
    rs = np.maximum(r, 1e-12)
    q = np.maximum(rho_d * np.sin(psi_d) - h_d * np.cos(psi_d), 0.0) / rs
    f = rho_d * W * q / rs ** 2
    dfds = np.abs(np.gradient(f, ds))
    dpsi = np.abs(np.gradient(psi_d, ds))
    mon = 1.0 / h_base + dfds / df_tol + dpsi / dpsi_tol
    away = sd_ > 0.15
    mon += np.where(away, np.abs(np.sin(psi_d)) / np.maximum(rho_d, 1e-3) / dpsi_tol, 0.0)
    mon += np.where(away & (rho_d < 0.5), 1.0 / (neck_tol * np.maximum(rho_d, 1e-3)), 0.0)
    hloc = np.clip(1.0 / mon, h_min, h_base)
    # make the spacing smooth: growth-limited minimum, then a little averaging in log space is not needed
    hloc = _dt_pass(hloc, ds, growth)
    return np.clip(hloc, h_min, h_base)


def new_grid(model, sd_, rho_d, h_d, psi_d, nmax=1400, **kw):
    """returns ell (fractions) and the arclength positions of the new nodes"""
    hloc = monitor_spacing_fast(model, sd_, rho_d, h_d, psi_d, **kw)
    dens = 1.0 / hloc
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(sd_))])
    ntot = cum[-1]
    N = int(min(nmax, max(12, math.ceil(ntot))))
    targets = np.linspace(0.0, ntot, N + 1)
    snew = np.interp(targets, cum, sd_)
    ell = np.diff(snew) / snew[-1]
    return ell, snew


def regrid(model, st, nmax=1400, **kw):
    sd_, rho_d, h_d, psi_d = dense_curve(st)
    ell, snew = new_grid(model, sd_, rho_d, h_d, psi_d, nmax=nmax, **kw)
    smid = 0.5 * (snew[1:] + snew[:-1])
    psi_new = np.interp(smid, sd_, psi_d)
    return State(psi_new, st.hp, ell, st.R)


def plate_state(model, d=0.9, nmax=1400, **kw):
    """flat membrane at height -d below the bead centre"""
    M = 24000
    sd_ = np.linspace(0.0, model.R, M)
    z = np.zeros(M)
    ell, snew = new_grid(model, sd_, sd_, np.full(M, -d), z, nmax=nmax, **kw)
    return State(np.zeros(len(ell)), -d, ell, model.R)


# =====================================================================================================
# 4. solvers
# =====================================================================================================
class SolveFail(Exception):
    pass


def _admissible(st):
    return st.D > 1e-3 * 1.0 and st.rho_min_all > 5e-4


def minimize_wt(model, st, wt, tol=1e-9, maxit=120, verbose=False, trust=0.3, floor=1e-5):
    """Newton with 'absolute value' Hessian modification (negative curvature is followed downhill, eigenvalues
    floored at `floor`), trust radius `trust` on the step (max-norm in psi, h_p) and Armijo line search.
    Returns (state, info)."""
    E0 = energies(model, st, wt)[0]
    info = dict(it=0, converged=False)
    for it in range(maxit):
        ev = Eval(model, st)
        g, H = ev.hessian(wt)
        gn = np.abs(g).max()
        w, V = sla.eigh(H, check_finite=False)
        mu = np.maximum(np.abs(w), floor)
        step = -V @ ((V.T @ g) / mu)
        sm = np.abs(step).max()
        if sm > trust:
            step *= trust / sm
        slope = float(g @ step)
        if gn < tol and abs(slope) < 1e-16:
            info.update(converged=True, it=it, gmax=gn, lam_min=w[0])
            return st, info
        t = 1.0
        z0 = st.z()
        ok = False
        for _ in range(40):
            try:
                stn = st.with_z(z0 + t * step)
                if _admissible(stn):
                    En = energies(model, stn, wt)[0]
                    if En <= E0 + 1e-4 * t * slope + 1e-14 * abs(E0):
                        ok = True
                        break
            except FloatingPointError:
                pass
            t *= 0.5
        if not ok:
            # nothing decreases the energy any more: converged to round-off if the gradient is tiny
            info.update(it=it, gmax=gn, reason="linesearch", lam_min=w[0], converged=bool(gn < 1e-6))
            return st, info
        st, E0 = stn, En
        if verbose:
            print("  it", it, "E", E0, "|g|", gn, "t", t, "slope", slope, "lam_min", w[0])
    info.update(it=maxit, reason="maxit", gmax=np.abs(Eval(model, st, hess=False).grad(wt)).max())
    return st, info


def hessian_min_eig(model, st, wt):
    ev = Eval(model, st)
    g, H = ev.hessian(wt)
    return sla.eigh(H, eigvals_only=True, subset_by_index=[0, 0])[0]


def kkt_solve(model, st, Gt, lam0, tol=1e-9, maxit=40, verbose=False):
    """constrained stationary state at fixed coverage G = Gt; returns (state, lam, info).  lam = w~ = de/dG."""
    lam = lam0
    N = st.N
    for it in range(maxit):
        ev = Eval(model, st)
        gG = ev.grad_gam()
        gr, Hl = ev.hessian(lam)               # gradient of E_el - lam G and its Hessian
        Gval = local_values(model, "gam", st.V).sum()
        res = Gval - Gt
        rn = max(np.abs(gr).max(), abs(res))
        if verbose:
            print("   kkt it", it, "res", res, "|F1|", np.abs(gr).max(), "lam", lam)
        if rn < tol:
            return st, lam, dict(converged=True, it=it, rn=rn, ev=ev, G=Gval)
        M = np.zeros((N + 2, N + 2))
        M[:N + 1, :N + 1] = Hl
        M[:N + 1, N + 1] = -gG
        M[N + 1, :N + 1] = -gG
        rhs = np.concatenate([-gr, [res]])
        try:
            sol = sla.solve(M, rhs, assume_a="sym", check_finite=False)
        except (sla.LinAlgError, np.linalg.LinAlgError, ValueError):
            sol = np.linalg.lstsq(M, rhs, rcond=None)[0]
        dz, dlam = sol[:N + 1], sol[N + 1]
        phi0 = rn
        t = 1.0
        ok = False
        z0 = st.z()
        for _ in range(25):
            try:
                stn = st.with_z(z0 + t * dz)
                if _admissible(stn):
                    evn = Eval(model, stn, hess=False)
                    lamn = lam + t * dlam
                    gn_ = evn.grad(lamn)
                    resn = local_values(model, "gam", stn.V).sum() - Gt
                    phin = max(np.abs(gn_).max(), abs(resn))
                    if phin < (1.0 - 1e-4 * t) * phi0 or phin < tol:
                        ok = True
                        break
            except FloatingPointError:
                pass
            t *= 0.5
        if not ok:
            return st, lam, dict(converged=False, it=it, rn=rn, reason="linesearch")
        st, lam = stn, lamn
    return st, lam, dict(converged=False, it=maxit, rn=rn, reason="maxit")


def tangent(model, st, lam, ev=None):
    """d(z, lam)/dG along the constrained branch"""
    N = st.N
    ev = ev or Eval(model, st)
    gG = ev.grad_gam()
    gr, Hl = ev.hessian(lam)
    M = np.zeros((N + 2, N + 2))
    M[:N + 1, :N + 1] = Hl
    M[:N + 1, N + 1] = -gG
    M[N + 1, :N + 1] = -gG
    rhs = np.zeros(N + 2)
    rhs[N + 1] = -1.0
    sol = sla.solve(M, rhs, assume_a="sym", check_finite=False)
    return sol[:N + 1], sol[N + 1]


def constraint_inertia(model, st, lam):
    """number of negative eigenvalues of the bordered matrix (1 for a constrained local minimum)"""
    ev = Eval(model, st)
    N = st.N
    gG = ev.grad_gam()
    gr, Hl = ev.hessian(lam)
    M = np.zeros((N + 2, N + 2))
    M[:N + 1, :N + 1] = Hl
    M[:N + 1, N + 1] = -gG
    M[N + 1, :N + 1] = -gG
    w = np.linalg.eigvalsh(M)
    return int((w < 0).sum()), w


# =====================================================================================================
# 5. continuation in the coverage G
# =====================================================================================================
class Branch:
    """e(G), w~(G) = e'(G) of the constrained minimum; list of stored shapes"""

    def __init__(self, model):
        self.model = model
        self.G, self.e, self.wt, self.Eb, self.Aex, self.rmin, self.N, self.states = [], [], [], [], [], [], [], []

    def add(self, st, lam, G):
        e, eb, ea, gm = energies(self.model, st, 0.0)
        self.G.append(G)
        self.e.append(e)
        self.wt.append(lam)
        self.Eb.append(eb)
        self.Aex.append(ea)
        self.rmin.append(st.rho_min)
        self.N.append(st.N)
        self.states.append(st)

    def arrays(self):
        return {k: np.array(getattr(self, k)) for k in ("G", "e", "wt", "Eb", "Aex", "rmin", "N")}


def trace_branch(model, st, lam, G0, G1, dG0=0.02, dG_max=0.05, dG_min=1e-5, regrid_every=2, nmax=1400,
                 verbose=False, tol=1e-9, grid_kw=None, stop_rmin=0.02, br=None):
    """Follow the constrained minimum from coverage G0 (state st with multiplier lam is already stationary there)
    to G1 (either direction).  Returns the Branch."""
    grid_kw = grid_kw or {}
    br = br or Branch(model)
    direction = 1.0 if G1 > G0 else -1.0
    G = G0
    dG = dG0
    st, lam, info = kkt_solve(model, st, G, lam, tol=tol)
    if not info["converged"]:
        raise SolveFail("start point does not converge: %s" % info)
    br.add(st, lam, G)
    nacc = 0
    t0 = time.time()
    while (G1 - G) * direction > 1e-12:
        step = min(dG, abs(G1 - G))
        Gn = G + direction * step
        # predictor
        try:
            tz, tl = tangent(model, st, lam, ev=info.get("ev"))
            z_pred = st.z() + tz * (direction * step)
            lam_pred = lam + tl * (direction * step)
            stp = st.with_z(z_pred)
            if not _admissible(stp):
                raise FloatingPointError
        except (FloatingPointError, sla.LinAlgError, np.linalg.LinAlgError):
            stp, lam_pred = st, lam
        stn, lamn, infon = kkt_solve(model, stp, Gn, lam_pred, tol=tol)
        if not infon["converged"] or infon["it"] > 7:
            if not infon["converged"] and step <= dG_min * 1.01:
                if verbose:
                    print("  continuation stopped at G=%.5f: %s" % (G, infon))
                br.stopped = infon
                return br
            dG = max(step * 0.5, dG_min)
            if verbose and not infon["converged"]:
                print("   step rejected G->%.5f, dG=%.2e (%s)" % (Gn, dG, infon.get("reason")))
            continue
        st, lam, G, info = stn, lamn, Gn, infon
        nacc += 1
        # regrid on the new curve and re-converge at the same G
        if regrid_every and nacc % regrid_every == 0:
            try:
                st2 = regrid(model, st, nmax=nmax, **grid_kw)
                st2, lam2, info2 = kkt_solve(model, st2, G, lam, tol=tol)
                if info2["converged"]:
                    st, lam, info = st2, lam2, info2
            except FloatingPointError:
                pass
        br.add(st, lam, G)
        if verbose and nacc % 5 == 0:
            print("  G=%.4f  w~=%.4f  e=%.5f  N=%d  rmin=%.3f  dG=%.3g  t=%.0fs" % (G, lam, br.e[-1], st.N, st.rho_min, dG, time.time() - t0))
        if info["it"] <= 3:
            dG = min(dG * 1.5, dG_max)
        elif info["it"] > 5:
            dG = max(dG * 0.7, dG_min)
        if st.rho_min < stop_rmin:
            br.stopped = dict(reason="neck closed", rmin=st.rho_min)
            return br
    br.stopped = dict(reason="reached G1")
    return br


# =====================================================================================================
# 6. phase lines from the branch  (S1, S2, E) and their refinement
# =====================================================================================================
GRIDS = {
    "coarse": dict(df_tol=0.03, dpsi_tol=0.08, h_base=0.3),
    "mid": dict(df_tol=0.02, dpsi_tol=0.05, h_base=0.2),
    "fine": dict(df_tol=0.01, dpsi_tol=0.03, h_base=0.15),
    "vfine": dict(df_tol=0.005, dpsi_tol=0.015, h_base=0.1),
}


def initial_state(model, wt0=2.0, grid="mid", d0=0.9, nmax=1400):
    kw = GRIDS[grid] if isinstance(grid, str) else grid
    st = plate_state(model, d0, nmax=nmax, **kw)
    st, info = minimize_wt(model, st, wt0)
    for _ in range(2):
        st = regrid(model, st, nmax=nmax, **kw)
        st, info = minimize_wt(model, st, wt0)
    if not info["converged"]:
        raise SolveFail("initial state not converged: %s" % info)
    return st


def initial_state_bound(model, G_target=0.25, grid="mid", w_start=2.0, dw=0.25, w_max=14.0, d0=0.95, nmax=1400):
    """start state on the lower (partially wrapped) branch by fixed-w~ continuation from the flat plate: raise w~ in steps of dw
    (warm start, regrid) until the coverage exceeds G_target.  Used for very narrow shells, where the plate is stable up to
    w~ ~ 3-4 and the constrained trace from the plate stalls.  Returns (state, w~)."""
    kw = GRIDS[grid] if isinstance(grid, str) else grid
    st = plate_state(model, d0, nmax=nmax, **kw)
    st, info = minimize_wt(model, st, w_start)
    good = (st, w_start)
    w = w_start
    while w <= w_max:
        for _ in range(2):
            st = regrid(model, st, nmax=nmax, **kw)
            st, info = minimize_wt(model, st, w)
        if not info["converged"]:
            break
        G = energies(model, st, w)[3]
        if G > 1.2:                       # jumped past the partial branch: keep the last good state
            break
        good = (st, w)
        if G >= G_target:
            break
        w += dw
    return good


def first_extrema(G, w, prom=0.02):
    """index of the first maximum (S1) and the following minimum (S2) of w(G) that have prominence > prom, or None"""
    n = len(w)
    imax = 0
    for i in range(n):
        if w[i] > w[imax]:
            imax = i
        if w[imax] - w[i] > prom and i > imax:
            break
    else:
        return None, None
    # require imax to be a genuine maximum (not the first point)
    if imax == 0:
        return None, None
    imin = imax
    for j in range(imax, n):
        if w[j] < w[imin]:
            imin = j
        if w[j] - w[imin] > prom and j > imin:
            return imax, imin
    return imax, None


def poly_extremum(G, w, i0, half=4, deg=3, kind="max"):
    """extremum of a cubic least-squares fit to w(G) on points i0-half .. i0+half"""
    lo, hi = max(0, i0 - half), min(len(G), i0 + half + 1)
    x = G[lo:hi]
    xc = x.mean()
    c = np.polyfit(x - xc, w[lo:hi], deg)
    dc = np.polyder(c)
    roots = np.roots(dc)
    roots = roots[np.abs(roots.imag) < 1e-12].real
    roots = roots[(roots > (x.min() - xc)) & (roots < (x.max() - xc))]
    if len(roots) == 0:
        return G[i0], w[i0]
    vals = np.polyval(c, roots)
    k = np.argmax(vals) if kind == "max" else np.argmin(vals)
    return roots[k] + xc, vals[k]


def hermite_e(Gp, ep, wp, Gq):
    """cubic Hermite interpolation of e(G) from nodes Gp with values ep and slopes wp (= w~(G))"""
    from scipy.interpolate import CubicHermiteSpline
    return CubicHermiteSpline(Gp, ep, wp)(Gq)


def refine_extremum(model, br, i0, kind, grid="mid", half=0.06, dG=0.01, nmax=1400, iters=3, tolG=4e-3, tolw=2e-3):
    """Locate the extremum of w~(G) near branch point i0: one adaptive trace (grid re-adapted every step; a frozen grid is
    NOT reliable near the neck, where it changed S2 by 0.2 at sigma~ = 1) through [G_c - half, G_c + half] with step dG,
    polynomial fit of the extremum (averages the ~2e-3 regrid noise), re-centre and repeat until G*, w~* settle.
    Returns (G*, w*, G array, w array, list of states) of the last trace."""
    kw = GRIDS[grid] if isinstance(grid, str) else grid
    cand = [(br.G[k], br.wt[k], br.states[k]) for k in range(len(br.G))]
    Gc = br.G[i0]
    prev = None
    for it in range(iters):
        below = [c for c in cand if c[0] <= Gc - half]
        G0, lam, st = max(below, key=lambda c: c[0]) if below else min(cand, key=lambda c: c[0])
        for _ in range(2):
            st = regrid(model, st, nmax=nmax, **kw)
            st, lam, info = kkt_solve(model, st, G0, lam)
            if not info["converged"]:
                raise SolveFail("refine: no convergence at G=%g" % G0)
        tr = trace_branch(model, st, lam, G0, Gc + half, dG0=dG, dG_max=dG, regrid_every=1, grid_kw=kw, nmax=nmax)
        G = np.array(tr.G)
        w = np.array(tr.wt)
        sel = (G >= Gc - 1.5 * half) & (G <= Gc + half + 1e-9)
        i = int(np.argmax(np.where(sel, w if kind == "max" else -w, -np.inf)))
        nfit = 7 if kind == "max" else 5
        Gs, ws = poly_extremum(G[sel], w[sel], int(np.argmin(np.abs(G[sel] - G[i]))), half=nfit, deg=4 if (kind == "max" and sel.sum() > 11) else 3, kind=kind)
        cand = cand + [(G[k], w[k], tr.states[k]) for k in range(len(G))]
        if prev is not None and abs(Gs - prev[0]) < tolG and abs(ws - prev[1]) < tolw:
            break
        prev = (Gs, ws)
        Gc = Gs
        if kind == "min":                          # the minimum is a sharp V near the neck: finer steps, narrower window
            half, dG = min(half, 0.035), min(dG, 0.005)
    return Gs, ws, G, w, tr.states


def e_line(model, br, iS1, iS2, wS1, wS2, GS1, GS2, extra_lower=(), extra_upper=(), grid="fine", nmax=1400, verbose=False):
    """equal-energy w~: partial minimum (lower branch) vs enveloped minimum (upper branch), fixed-w~ Newton on fine grids.
    Start states: the traced branch plus (wt, G, state) triples from the refinement of S1 / S2; a result is accepted only
    if it lies in the right basin (lower: G < G_S1 + 0.01, upper: G > G_S2 - 0.01)."""
    from scipy.optimize import brentq
    kw = GRIDS[grid] if isinstance(grid, str) else grid
    lower = [(br.wt[k], br.G[k], br.states[k]) for k in range(0, iS1 + 1)] + list(extra_lower)
    upper = [(br.wt[k], br.G[k], br.states[k]) for k in range(iS2, len(br.wt))] + list(extra_upper)

    def solve_at(w, cands, side):
        order = sorted(range(len(cands)), key=lambda j: abs(cands[j][0] - w))
        for j in order[:6]:
            st = cands[j][2]
            try:
                for _ in range(2):
                    st = regrid(model, st, nmax=nmax, **kw)
                    st, info = minimize_wt(model, st, w)
            except FloatingPointError:
                continue
            if not info["converged"]:
                continue
            e, eb, ea, gm = energies(model, st, w)
            if (side == "lower" and gm < GS1 + 0.01) or (side == "upper" and gm > max(GS2 - 0.05, GS1 + 0.1)):
                return e, gm, st
        raise SolveFail("no %s-branch state found at w~=%g" % (side, w))

    def dE(w):
        ep, gp, _ = solve_at(w, lower, "lower")
        ee, ge, _ = solve_at(w, upper, "upper")
        if verbose:
            print("   E-line eval w~=%.5f Ep=%.6f (G=%.4f) Ee=%.6f (G=%.4f) diff=%.6f" % (w, ep, gp, ee, ge, ep - ee))
        return ep - ee

    a = wS2 + 0.03 * (wS1 - wS2)
    b = wS1 - 0.03 * (wS1 - wS2)
    fa, fb = dE(a), dE(b)
    if fa * fb > 0:
        return None, dict(reason="no sign change", fa=fa, fb=fb, a=a, b=b)
    wE = brentq(dE, a, b, xtol=1e-3, rtol=1e-8, maxiter=40)
    ep, gp, stp = solve_at(wE, lower, "lower")
    ee, ge, ste = solve_at(wE, upper, "upper")
    return wE, dict(E=0.5 * (ep + ee), Gp=gp, Ge=ge, state_p=stp, state_e=ste)


def soft_lines(sigma, s=0.25, p=1.0, R=10.0, grid_scan="mid", grid_ref="mid", wt0=2.0, G_end=1.985, dG=0.03,
               prom=0.02, ngl=6, verbose=False, nmax=1400, refine=True, do_E=True):
    """S1, S2, E (+ coverages) of the continuum soft-shell model.  See module docstring."""
    model = Model(sigma, s=s, p=p, R=R, ngl=ngl)
    kw = GRIDS[grid_scan] if isinstance(grid_scan, str) else grid_scan
    t0 = time.time()
    retried = False
    br = None
    try:
        st = initial_state(model, wt0, grid_scan, nmax=nmax)
        G0 = energies(model, st, 0.0)[3]
        br = trace_branch(model, st, wt0, G0, G_end, dG0=0.02, dG_max=dG, grid_kw=kw, nmax=nmax, verbose=verbose)
        ok = br.G[-1] > G_end - 0.05
    except SolveFail:
        ok = False
    if not ok:                                   # narrow shell: the flat plate stays stable to w~ ~ 3-4: start on the bound branch
        retried = True
        st0, wt0 = initial_state_bound(model, grid=grid_scan, nmax=nmax)
        G0 = energies(model, st0, 0.0)[3]
        br = trace_branch(model, st0, wt0, G0, G_end, dG0=0.02, dG_max=dG, grid_kw=kw, nmax=nmax, verbose=verbose)
    a = br.arrays()
    out = dict(sigma=sigma, s=s, p=p, R=R, G0=G0, branch=br, stopped=getattr(br, "stopped", None), wt0=wt0, retried=retried)
    imax, imin = first_extrema(a["G"], a["wt"], prom)
    out["imax"], out["imin"] = imax, imin
    if imax is None or imin is None:
        out.update(S1=np.nan, S2=np.nan, E=np.nan, note="no N-shaped w~(G): no barrier", time=time.time() - t0)
        return out
    Gs1, ws1 = poly_extremum(a["G"], a["wt"], imax, kind="max")
    Gs2, ws2 = poly_extremum(a["G"], a["wt"], imin, kind="min")
    out.update(S1_coarse=ws1, S2_coarse=ws2, G_S1_coarse=Gs1, G_S2_coarse=Gs2)
    if refine:
        Gs1, ws1, g1, w1, sts1 = refine_extremum(model, br, imax, "max", grid=grid_ref, nmax=nmax)
        Gs2, ws2, g2, w2, sts2 = refine_extremum(model, br, imin, "min", grid=grid_ref, nmax=nmax)
        out.update(ref_S1=(g1, w1), ref_S2=(g2, w2))
        extra_lower = [(w1[k], g1[k], sts1[k]) for k in range(len(g1)) if g1[k] <= Gs1]
        extra_upper = [(w2[k], g2[k], sts2[k]) for k in range(len(g2)) if g2[k] >= Gs2]
    out.update(S1=ws1, S2=ws2, G_S1=Gs1, G_S2=Gs2)
    if do_E and ws1 > ws2:
        try:
            wE, info = e_line(model, br, imax, imin, ws1, ws2, Gs1, Gs2, extra_lower=extra_lower if refine else (),
                              extra_upper=extra_upper if refine else (), grid=grid_ref, nmax=nmax, verbose=verbose)
        except SolveFail as exc:
            wE, info = None, dict(reason=str(exc))
        if wE is None:
            out.update(E=np.nan, E_note=info.get("reason"))
        else:
            out.update(E=wE, E_G_partial=info["Gp"], E_G_env=info["Ge"], E_energy=info["E"])
    out["time"] = time.time() - t0
    return out
