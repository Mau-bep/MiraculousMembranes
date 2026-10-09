#!/usr/bin/env python
"""Reproduce the checkpoints against the paper (Deserno 2003) from deserno_theory.py / deserno_table.npz.

    micromamba run -n mir_mem python test_deserno_theory.py            # ~1-2 min
    micromamba run -n mir_mem python test_deserno_theory.py --full     # + solves an off-node sigma~ from scratch

Every check prints: paper value, computed value, tolerance used and PASS/FAIL.  The tolerances are written down
next to each check; where the paper only gives a rounded or read-off number the tolerance reflects that.
"""
import os
import sys
import json
import time
import numpy as np
import deserno_theory as dt

HERE = os.path.dirname(os.path.abspath(__file__))
FULL = '--full' in sys.argv
rows = []


def check(name, paper, computed, tol, mode='abs', note=''):
    """mode 'abs': |c-p| <= tol ; 'rel': |c/p-1| <= tol ; 'lt'/'gt': computed < / > paper (paper = bound)."""
    if mode == 'abs':
        ok = abs(computed - paper) <= tol
    elif mode == 'rel':
        ok = abs(computed / paper - 1.0) <= tol
    elif mode == 'lt':
        ok = computed < paper
    else:
        ok = computed > paper
    rows.append((name, paper, computed, tol, mode, ok, note))
    ps = f"{paper:.6g}" if isinstance(paper, (int, float, np.floating)) else str(paper)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}\n        paper={ps}  computed={computed:.8g}  tol={tol:g} ({mode}) {note}")
    return ok


t_start = time.time()
tab = dt.load_table()
sg = tab['sigma_t']
print(f"table: {len(sg)} sigma~ nodes, {sg[0]:g} .. {sg[-1]:g}; fallback flags set: {int(tab['asymptotic_fallback'].sum())}")
print('S-shaped nodes:', int(tab['s_shaped'].sum()), ' with exact branch 3:', int(tab['branch3_exact'].sum()))
print()


def node(s):
    return float(sg[np.argmin(np.abs(np.log(sg / s)))])


# --------------------------------------------------------------------------------------------------
print('=== 1. E_free(z=2) = 0, E_free = 0 at sigma~ = 0, W line w~ = 4 ===')
for s in (0.01, 1.0, 100.0):
    c = dt.solve_curve(node(s))
    r1 = float(c.F(2.0 - 1e-3)[0]) / (4.0 * c.sigma * 1e-3) - 1.0
    r2 = float(c.F(2.0 - 1.5e-4)[0]) / (4.0 * c.sigma * 1.5e-4) - 1.0
    check(f'F(2-d)/(4 sigma~ d) - 1 shrinks as d -> 0 at sigma~={s:g} (exact F\'(2) = -4 sigma~): d=1.5e-4 value', 5e-3, abs(r2), 0.0, 'lt',
          note=f'(d=1e-3: {r1:.2e} -> d=1.5e-4: {r2:.2e}; first correction ~ d ln d)')
    check(f'   ... and is smaller than at d=1e-3', abs(r1), abs(r2), 0.0, 'lt')
check('free_energy(2, 1) == 0', 0.0, dt.free_energy(2.0, node(1.0)), 1e-14, 'abs')
c = dt.solve_curve(node(1e-4))
zz = np.linspace(0.001, 1.999, 4000)
check('max_z F(z; sigma~=1e-4) (sigma~ -> 0: F -> 0 uniformly, O(sigma~))', 0.0, float(c.F(zz).max()), 1e-3, 'abs',
      note=f'(= {c.F(zz).max()/1e-4:.2f} x sigma~)')
# W line: continuous start of partial wrapping at w~ = 4  <=>  w_eq(z->0) = (1-psi_dot_0)^2 -> 4
worst = 0.0
for j in range(0, len(sg), 20):
    cj = dt.solve_curve(sg[j])
    z0 = cj.b1.zmin
    dev = abs(cj.b1.w[0] - 4.0) / (cj.sigma * z0 * (1.0 + abs(np.log(cj.sigma * z0 / 2.0))))
    worst = max(worst, dev)
check('W line: w_eq(z->0) = 4 up to O(sigma~ z ln z)  [max over nodes of |w_eq(z0)-4|/(sigma~ z0 (1+|ln|)) < 10]', 10.0, worst, 0.0,
      'lt', note='(bounded O(1) coefficient => limit is exactly 4)')
check('E line always right of the zero-energy line w~ = 4 + 2 sigma~: min (w_E-4)/(2 sigma~) > 1', 1.0,
      float(np.min((tab['w_E'] - 4.0) / (2.0 * sg))), 0.0, 'gt')
check('(w_E-4)/sigma~ -> 2 as sigma~ -> 0  (value at 1e-4)', 2.0, tab['wE_minus4_over_sigma'][0], 0.03, 'abs')
check('spinodals bracket the W line: S1 > 4 > S2 at sigma~ = 0.01', 1.0,
      float(tab['S1'][np.argmin(abs(sg - 0.01))] > 4.0 and tab['S2'][np.argmin(abs(sg - 0.01))] < 4.0), 0.0, 'abs')

# --------------------------------------------------------------------------------------------------
print('\n=== 2. small-z expansion (Eq. 21) and small-gradient formula (Eq. 20) ===')
for s in (0.01, 1.0):
    c = dt.solve_curve(node(s))
    for z in (2e-3, 1e-2):
        Fn = float(c.F(z)[0])
        r21 = Fn / float(dt.F_small_z(z, c.sigma))
        r20 = Fn / float(dt.F_small_gradient(z, c.sigma))
        print(f'        sigma~={s:g} z={z:g}: F={Fn:.5e}  F/Eq21={r21:.5f}  F/Eq20={r20:.5f}')
    z = 2e-3
    Fn = float(c.F(z)[0])
    check(f'Eq. 21 at sigma~={s:g}, z=2e-3 (ratio -> 1, corrections O(z) relative)', 1.0, Fn / float(dt.F_small_z(z, c.sigma)),
          0.02, 'abs')
    check(f'Eq. 20 at sigma~={s:g}, z=2e-3 (exact for z -> 0)', 1.0, Fn / float(dt.F_small_gradient(z, c.sigma)), 0.01, 'abs')
c = dt.solve_curve(node(0.01))
print('        Eq. 20 deviation grows with z (sigma~=0.01; Eq. 20 has a pole at z = 1):',
      ', '.join(f'z={z:g}: {float(c.F(z)[0]/dt.F_small_gradient(z, c.sigma)):.3f}' for z in (0.01, 0.1, 0.5, 0.9)))

# --------------------------------------------------------------------------------------------------
print('\n=== 3. S-shaped contact curvature for sigma~ > sigma~_c = 4.721139, first at z_c = 1.86289 ===')
i0 = np.where(sg < 4.721139)[0][-1]
check(f'node sigma~={sg[i0]:.4f} (< sigma_c) is single-valued', 0.0, float(tab['s_shaped'][i0]), 0.0, 'abs')
check(f'node sigma~={sg[i0+1]:.4f} (> sigma_c) is S-shaped', 1.0, float(tab['s_shaped'][i0 + 1]), 0.0, 'abs')
cj = os.path.join(HERE, 'cusp_result.json')
if os.path.exists(cj):
    cr = json.load(open(cj))
    check('sigma~_c from reverse-trace fold bisection + (zmax-zmin)^2 extrapolation', dt.PAPER['sigma_c'], cr['sigma_c'], 5e-3, 'abs',
          note='(%s)' % cr.get('note', ''))
    check('z_c (centre of the nascent S)', dt.PAPER['z_c'], cr['z_c'], 5e-3, 'abs')
else:
    print('        (cusp_result.json not found: run cusp_study.py)')

# --------------------------------------------------------------------------------------------------
print('\n=== 4. sigma~ = 1 : E boundary, barrier, spinodals (Fig. 2, 3, 4, 6) ===')
c1 = dt.solve_curve(node(1.0))
sm = c1.summary()
check('(w_E-4)/sigma~ at sigma~=1 (paper Fig. 6 text: 0.142 above 2)', 2.142, (sm['wE'] - 4.0), 5e-4, 'abs')
check('w_E(1) ~ 6.1 (Fig. 4: enveloped state stable at w~ ~ 6.1)', 6.14, sm['wE'], 0.01, 'abs')
check('barrier(1) ~ 1.05 (Fig. 3; paper: 66 kT with kappa = 20 kT)', 1.05, sm['barrier'], 0.02, 'abs',
      note=f'(= {sm["barrier"]*np.pi*20:.1f} kT for kappa = 20 kT; paper 66 kT)')
check('S1(sigma~=1) ~ 7.5 (barrier to enveloping vanishes); exact (1+sqrt(1+2 sigma~))^2 = 7.4641', 7.5, sm['S1'], 0.05, 'abs')
check('S1(1) vs exact (1+sqrt3)^2 (psi_dot_0^2 = 1+2 sigma~ at z = 1, Eq. 13 / ref. [46])', (1 + np.sqrt(3.0)) ** 2, sm['S1'], 1e-6, 'abs')
check('S2(sigma~=1) ~ 2.7 (unbinding barrier vanishes)', 2.7, sm['S2'], 0.05, 'abs')
check('z at maximum of w_eq (S1) at sigma~=1 is exactly z = 1', 1.0, sm['zS1'], 2e-3, 'abs')
# Fig. 2 read-offs
s05 = 0.55
check('S1 passes (w~, sigma~) = (6.0, ~0.55)  [read off Fig. 2]', 6.0, float(dt.spinodal_S1(s05)), 0.15, 'abs')
s072 = 0.72
check('S2 passes (w~, sigma~) = (3.0, ~0.72)  [read off Fig. 2]', 3.0, float(dt.spinodal_S2(s072)), 0.15, 'abs')
check('E passes (w~, sigma~) = (6.0, ~0.93) [read off Fig. 2]', 6.0, float(dt.w_E(0.93)), 0.1, 'abs')

# --------------------------------------------------------------------------------------------------
print('\n=== 5. low-tension asymptotics (Eq. 26/27), Fig. 3 barrier, 0.86 power law (Fig. 9) ===')
for s in (1e-4, 1e-3, 1e-2, 0.1, 0.4):
    wn = float(dt.w_E(s))
    w26 = float(dt.w_E_small_gradient(s)[0])
    w27 = float(dt.w_E_small_gradient(s, order=27)[0])
    print(f'        sigma~={s:g}: numerics w_E-4={wn-4:.6e}  Eq.26: {w26-4:.6e} (ratio {(wn-4)/(w26-4):.5f})  Eq.27: ratio {(wn-4)/(w27-4):.4f}')
for s in (1e-4, 1e-3):
    wn = float(dt.w_E(s))
    check(f'numerics vs Eq. 26 at sigma~={s:g} (small-gradient asymptotics exact as sigma~ -> 0)', 1.0,
          (wn - 4.0) / (float(dt.w_E_small_gradient(s)[0]) - 4.0), 5e-3, 'abs')
check('Fig. 3: barrier(0.22) ~ 0.35 (SFV: 22 kT)', 0.35, float(dt.barrier(0.22)), 0.01, 'abs',
      note=f'(= {float(dt.barrier(0.22))*np.pi*20:.1f} kT; paper ~22 kT)')
slopes = {}
for lo, hi in ((1e-4, 1e-1), (1e-4, 1.0), (1e-3, 1e-1), (1e-2, 1.0)):
    b = dt.barrier(np.array([lo, hi]))
    slopes[(lo, hi)] = np.log(b[1] / b[0]) / np.log(hi / lo)
print('        log-log slope of barrier vs sigma~ over', {k: round(v, 3) for k, v in slopes.items()})
check('Fig. 9 empirical exponent 0.86 (average log-log slope over 1e-4..1; paper fit range unspecified)', 0.86,
      slopes[(1e-4, 1.0)], 0.05, 'abs', note='(local slope varies 0.8-0.9 over the range)')
check('barrier -> 0 as sigma~ -> 0 (value at 1e-4)', 0.0, float(dt.barrier(1e-4)), 1e-3, 'abs')

# --------------------------------------------------------------------------------------------------
print('\n=== 6. high tension: (w~-4)/sigma~ -> 4 - 3 A^(2/3) sigma~^(-1/3), A ~ 5.650 (Eqs. 34-37, Figs. 6, 8, 9) ===')
A = dt.PAPER['A']
print('        sigma~     (w-4)/s   Eq.35     (2-z_E)s^(1/3)  [A^(2/3)=%.4f]  barrier/(Eq.37 prefactor)  [A^(4/3)=%.4f]' % (A ** (2 / 3), A ** (4 / 3)))
for s in (1e3, 1e4, 1e5, 1e6):
    j = np.argmin(abs(sg - s))
    sj = sg[j]
    print('        %-10.3g %.5f  %.5f  %.5f                         %.4f' % (
        sj, tab['wE_minus4_over_sigma'][j], (dt.high_tension_wE(sj) - 4) / sj, (2 - tab['z_partial'][j]) * sj ** (1 / 3),
        tab['barrier'][j] / (0.75 * (2 * np.sqrt(3) - 3) * sj ** (1 / 3))))
j = np.argmin(abs(sg - 1e6))
check('A from (2-z_E) sigma~^(1/3) = A^(2/3) at sigma~=1e6', A, ((2 - tab['z_partial'][j]) * sg[j] ** (1 / 3)) ** 1.5, 0.02, 'rel')
check('A from barrier prefactor (Eq. 37) at sigma~=1e6', A, (tab['barrier'][j] / (0.75 * (2 * np.sqrt(3) - 3) * sg[j] ** (1 / 3))) ** 0.75, 0.01, 'rel')
check('A from z_barrier (Eq. 36) at sigma~=1e6', A, (((2 - tab['z_barrier'][j]) * sg[j] ** (1 / 3)) / (0.5 * (2 - np.sqrt(3)))) ** 1.5, 0.02, 'rel')
check('(w~-4)/sigma~ -> 4 at sigma~=1e6 (Eq. 35 predicts %.4f)' % ((dt.high_tension_wE(1e6) - 4) / 1e6), 3.905, tab['wE_minus4_over_sigma'][j], 0.02, 'abs')
j2 = np.argmin(abs(sg - 3e2))
check('Fig. 6 crossover: (w-4)/sigma~ at sigma~=1 -> 2.14, at 10^2 -> 2.89', 2.89, float(tab['wE_minus4_over_sigma'][np.argmin(abs(sg - 100))]), 0.05, 'abs')

# --------------------------------------------------------------------------------------------------
print('\n=== 7. Fig. 5 : E curve minimum 1.37 at a/lambda = 4.4, -> 2 at large a/lambda ===')
a = np.geomspace(1.0, 1000.0, 20001)
wE_over, W_over = dt.fig5_curves(a)
i = int(np.argmin(wE_over))
check('min over a/lambda of w/sigma on E (paper 1.37)', 1.37, wE_over[i], 0.02, 'abs')
check('location a/lambda of the minimum (paper 4.4)', 4.4, a[i], 0.15, 'abs', note=f'(sigma~ = {a[i]**2:.2f}; paper ~19.4)')
check('w/sigma on E at a/lambda = 1000 -> 2 (approach from below: 1.96)', 2.0, float(wE_over[-1]), 0.06, 'abs')
check('W line is w/sigma = 2 (lambda/a)^2 exactly', 0.0, float(np.max(np.abs(W_over - 2.0 / a ** 2))), 1e-14, 'abs')
check('particles with w/sigma >= 2 are always enveloped: E curve < 2 for a/lambda > 8', 2.0, float(wE_over[a > 8].max()), 0.0, 'lt')
check('at a/lambda=4.4 the "free" line W is w/sigma=0.103 and E is 1.37 (E above W)', 0.0, float(wE_over[i] - 2.0 / a[i] ** 2), 0.0, 'gt')

# --------------------------------------------------------------------------------------------------
print('\n=== 8. internal consistency / independent checks ===')
# (a) energy vs contact curvature: F_{i+1}-F_i = int F'(z) dz with F' = (1-pd)^2-4-2 sigma z, pd from the spline
xg, wg = np.polynomial.legendre.leggauss(6)
for s in (0.01, 1.0, 100.0, 1e4):
    c = dt.solve_curve(node(s))
    b = c.b1
    z, F = b.z, b.F
    zi = 0.5 * (z[1:] + z[:-1])[:, None] + 0.5 * (z[1:] - z[:-1])[:, None] * xg[None, :]
    Fp = (1.0 - b.pds(zi)) ** 2 - 4.0 - 2.0 * c.sigma * zi
    integ = (0.5 * (z[1:] - z[:-1]) * (Fp * wg[None, :]).sum(1))
    err = np.abs(integ - np.diff(F))
    sel = z[1:] > 0.02
    scale = np.maximum(np.abs(np.diff(F)), 1e-3 * np.abs(F[1:]).max())
    check(f'F(z_(i+1)) - F(z_i) = int F\' dz [energy integral vs contact curvature], max rel err, sigma~={s:g}', 0.0,
          float(np.max(err[sel] / scale[sel])), 3e-3, 'abs', note='(two independent outputs of the ODE solve)')
# (b) Eq. 13 and H = 0 on a short trace
sh = dt.Shooter(1.0)
pts, info = dt.trace(sh, max_pts=60, **dt.PROFILE)
z, pd, pr = pts[:, 2], pts[:, 3], pts[:, 8]
sel = np.abs(1 - z) > 0.05
pr13 = np.sqrt(z * (2 - z)) / (1 - z) * (1 + 2 * 1.0 * z - pd ** 2)
check('Eq. 13 p_r(0) = sqrt(z(2-z))/(1-z) {1 + 2 sigma~ z - psi_dot_0^2}: max rel dev over 60 points', 0.0,
      float(np.max(np.abs(pr[sel] / pr13[sel] - 1.0))), 1e-6, 'abs')
# (c) leave-one-out interpolation error of the table splines (every second node removed)
from scipy.interpolate import CubicSpline
x = np.log10(sg)
keep = np.arange(0, len(sg), 2)
g = (tab['w_E'] - 4.0) / sg
cs = CubicSpline(x[keep], g[keep])
odd = np.arange(1, len(sg) - 1, 2)
err_w = np.max(np.abs(cs(x[odd]) - g[odd]) / g[odd])
cs2 = CubicSpline(x[keep], np.log(tab['barrier'][keep]))
err_b = np.max(np.abs(np.exp(cs2(x[odd])) / tab['barrier'][odd] - 1))
check('interpolation error of (w_E-4)/sigma~ when the grid is thinned to 10 nodes/decade (max rel)', 0.0, float(err_w), 1e-4, 'abs',
      note='(the shipped 20 nodes/decade table is better)')
check('interpolation error of barrier at 10 nodes/decade (max rel)', 0.0, float(err_b), 5e-3, 'abs')
# (d) brute-force scan of the lowest-energy envelope reproduces the analytic E line and barrier
check('scan on fine z grid: min_z E(z; w_E) - E(2) (scaled) over all nodes', 0.0, float(np.max(np.abs(tab['scan_err']))), 1e-9, 'abs')
check('scan on fine z grid: |barrier_scan - barrier|/barrier over all nodes', 0.0, float(np.max(np.abs(tab['barrier_scan_err']))), 5e-4, 'abs')
# (e) variational minimisation (independent method) -- optional, needs ./prior
prior = os.path.join(HERE, 'prior')
if os.path.exists(os.path.join(prior, 'engine_dev.py')):
    sys.path.insert(0, prior)
    try:
        from engine_dev import sweep
        for s, zs in ((1.0, (0.5, 1.0, 1.5)), (10.0, (1.0, 1.5)), (0.01, (0.8,))):
            c = dt.solve_curve(node(s))
            allz = list(np.linspace(0.05, zs[0], 6))
            for z in zs[1:]:
                allz += list(np.linspace(allz[-1], z, 5)[1:])
            recs = sweep(s, allz, rel0=0.4, nlev=3)
            for r in recs:
                if any(abs(r['z'] - z) < 1e-9 for z in zs) and r['ok']:
                    Fs = float(c.F(r['z'])[0])
                    check(f'variational minimisation vs shooting F, sigma~={s:g}, z={r["z"]:.2f} (rel)', 0.0, r['F'] / Fs - 1.0, 2e-5, 'abs')
    except Exception as e:
        print('        (variational comparison skipped:', repr(e), ')')
else:
    print('        (prior/ not present: variational comparison skipped)')

# --------------------------------------------------------------------------------------------------
if FULL:
    print('\n=== 9. --full: off-node sigma~ solved from scratch vs table interpolation ===')
    for s in (0.22, 19.4):
        t0 = time.time()
        cur = dt.compute_curve(s)
        r = cur.summary()
        print(f'        sigma~={s}: solved in {time.time()-t0:.0f}s: w_E={r["wE"]:.7f} barrier={r["barrier"]:.6f} S1={r["S1"]:.5f} S2={r["S2"]}')
        check(f'table interpolation of w_E at off-node sigma~={s} vs direct solve (rel (w-4))', 0.0,
              (float(dt.w_E(s)) - 4.0) / (r['wE'] - 4.0) - 1.0, 2e-5, 'abs')
        check(f'table interpolation of barrier at off-node sigma~={s} vs direct solve (rel)', 0.0, float(dt.barrier(s)) / r['barrier'] - 1.0, 5e-4, 'abs')
        check(f'table interpolation of S1 at off-node sigma~={s} vs direct solve (rel (S1-4))', 0.0,
              (float(dt.spinodal_S1(s)) - 4.0) / (r['S1'] - 4.0) - 1.0, 1e-3, 'abs')

n_ok = sum(1 for r in rows if r[5])
print(f'\n==== {n_ok}/{len(rows)} checks passed in {time.time()-t_start:.0f} s ====')
for r in rows:
    if not r[5]:
        print('FAILED:', r[0], '| paper', r[1], '| computed', r[2], '| tol', r[3], r[4])
sys.exit(0 if n_ok == len(rows) else 1)
