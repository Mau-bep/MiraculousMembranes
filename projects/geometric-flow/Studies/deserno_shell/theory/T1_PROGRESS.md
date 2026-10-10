# T1 progress (continuum soft-shell, direct minimisation)

## Milestone 1 (done): solver works
* `softshell_direct.py` written: polyline in chord angles psi_i with fixed arclength fractions, S = R/sum(ell cos psi) closes the clamp
  rho_N = R exactly; exact complex-step element derivatives chained through the cumulative-sum map; dense reduced Hessian verified
  against finite differences (1e-9). Constrained (fixed coverage G) Newton/KKT continuation gives e(G), w~(G) = e'(G).
* Pole artefact found and fixed: node-sum bending has zero weight at rho_0 = 0, so a conical/catenoid pole (psi_0 ~ 0.6) was a spurious
  local minimum (lowered E by ~0.02). Added the pole cell energy 2 psi_0^2.
* First result sigma~ = 0.53, s = 0.25, p = 1, R = 10 (scan grid 'mid', refinement 'fine'): S1 = 6.650, S2 = 6.374, E = 6.550
  (Deserno 5.953 / 3.208 / 5.118; simulation 7.0 / 6.6 / 6.6 with +-0.2). Coverage of partial / enveloped state at E: 1.395 / 1.901.
  => the continuum soft shell already explains most of the offset. (Still to be converged: grid, R, quadrature.)

## TODO
convergence (N, R, ngl), s -> 0 validation vs Deserno, tests file, scans (s, p, sigma~), sim comparison (partial branch energies), figures, md.

## Milestone 2 (done): convergence at sigma~ = 0.53, s = 0.25, p = 1 (S1 / S2 / E)
* grid scan coarse+ref mid: 6.6505 / 6.3745 / 6.5509; mid+fine: 6.6497 / 6.3736 / 6.5502; mid+vfine (N ~ 1000): 6.6494 / 6.3733 / 6.5499
  (converged to 1e-3).  Gauss points 6 -> 10: change < 1e-4.  R_dom = 6, 10, 20, 40: 6.6504, 6.6497, 6.6497, 6.6498 (S1), S2 6.3746, 6.3736, 6.3735, 6.3735,
  E 6.5510, 6.5502, 6.5502, 6.5502: R_dom-independent at sigma~ = 0.53.
* sigma~ = 0.27, s = 0.05 (more demanding): grids coarse/mid/fine/vfine give S1 5.125/5.122/5.120, S2 4.6546/4.6529/4.6523, E 4.9551/4.9549/4.9548; R = 30: same to 5e-4.
  => the large S2 offset at s = 0.05 (Deserno 3.52) is real model behaviour, not discretisation.
* sigma~ = 0.13, s = 0.1: no N-shaped w~(G) (no barrier): smooth crossover, as in the simulation at sigma~ <= 0.13.
* tools: test_softshell_direct.py, run_softshell_scan.py (resumable; sets main / s0 / R / sim), compare_sim.py written.
* running: s0 series (s -> 0, R = 10) -> ../data/softshell_lines.csv

## Milestone 3 (done): robustness fixes found at sigma~ = 1
* A FROZEN fine grid is not reliable for S2 at large sigma~ (changed S2 by 0.2 at sigma~ = 1, s = 0.2); refinement now re-adapts the grid every step
  (dG = 0.01, then 0.005 / +-0.035 in G for the sharp V of the minimum), polynomial fit.  sigma~ = 1, s = 0.2: S1 8.197, S2 6.101, E 7.640 (adaptive fine trace
  gives min w~ = 6.103 at G = 1.92).  Expect S2 uncertainty ~0.02 at large sigma~ (w~(G) is a steep V there, curvature ~500).
* E-line: start states now also come from the refinement traces; result accepted only in the right basin (previous version failed with "no sign change").
* s -> 0 first look (R = 10, p = 1): sigma~ = 0.53, s = 0.2/0.1/0.05/0.025: S1 6.49/6.19/6.05/5.98 (Deserno 5.953), E 6.33/5.86/5.59/5.42 (5.118),
  S2 6.02/5.12/4.52/4.10 (3.208): S2 shift ~ sqrt(s), S1 shift faster (~s^1.3).  sigma~ <= 0.27 and s >= 0.1-0.2: no barrier (smooth crossover).
* All scans were restarted with the fixed code (chain_0/1.sh: sets s0, main, sim, R -> ../data/softshell_lines.csv; logs scan_*.log).
* TODO after the scans: run test_softshell_direct.py (--slow), compare_sim.py, make_softshell_figures.py, fill in softshell_direct.md section 3.

## Milestone 4 (done): scans
* ../data/softshell_lines.csv: sets s0 (p=1, s = 0.2 ... 0.0125, R = 10, + R = 30 spot checks), main (5 s x 4 p x 4 sigma~, R = 10, 80 rows, no failures),
  R (R = 6, 10, 20, 40 at s = 0.25).  R_dom-independence: sigma~ = 0.53 and 1.0, R = 6 ... 40: all lines change by <= 0.001.
* s -> 0: S2 and E shifts from Deserno scale ~ s^0.5-0.65 (negative-line-tension picture: shell gain beyond the contact line ~ sqrt(s)); S1 shift scales faster (~s^0.8-1.6, noisy);
  at s = 0.00625 and sigma~ >= 0.27 the flat plate stays stable to w~ ~ 3-4 and the trace from w~ = 2 failed -> now auto-retry from a bound state (rerun in progress).
* simulation grid sigma~ = 2k/15 at s = 0.25: continuum S1 / E / S2 e.g. sigma~ = 1.2: 9.10 / 8.47 / 6.63 (sim 9.4 / 8.66 / 7.4, Deserno 8.09 / 6.58 / 2.56);
  sigma~ = 2: 11.99 / 10.81 / 6.43 (sim > 12 / 11.09 / 7.4; Deserno 10.56 / 8.39 / 1.93): the sigma~-independent S2 plateau of the simulation is reproduced (6.4-6.6).
* running: compare_sim.py (partial-branch energies, coverage before/after jumps), s0 rerun for s = 0.00625.

## Milestone 5 (done): final
* test_softshell_direct.py --slow: 37 / 37 PASS (222 s).  compare_sim.py done (partial branch: 236 of 245 rows comparable, E_sim - E_cont mean +0.026, rms 0.050).
* s = 0.00625 rows rerun with the bound-start retry (`initial_state_bound`); s = 0.003125 not possible (trace fails).
* collapse of all p onto the p = 1 curve with s_eff = s (I_p/I_1)^2 found (analyze in md 3.1(d)); figures in ../figures/softshell_*.png (7 files).
* softshell_direct.md completed (derivation, numerics, results, limitations).  Nothing committed (orchestrator commits).  Files outside the assignment touched: none.
