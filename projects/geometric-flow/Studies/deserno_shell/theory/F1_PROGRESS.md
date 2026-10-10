# F1 progress (formal analysis of the finite-range shell)

M1 (done): flat boundary layer solved with solve_bvp (scratchpad bl0.py -> to be merged into theory/boundary_layer.py).
 * Exact scale invariance: with w = kappa*Dc^2/2 the flat BL problem has NO free parameter after d = delta*D, s = sqrt(delta/Dc)*t:
   D'''' = (1/2) W'(D), D''(inf) = 1. So tau = -T_p * w * sqrt(delta*ell_c), ell_c = 1/Dc = sqrt(kappa/(2w)), T_p pure number.
 * Numerics: Dc = 1 to 1e-6 for every start guess (contact condition unchanged by the shell), first integral const to 1e-11,
   T_1 = 1.198, T_2 = 1.029, T_4 = 0.874, T_8 = 0.739 (>= trial bound 0.972/0.838/0.714/0.604), lift-off offset h = 0.019 delta, undershoot d_min = -0.13 delta (p=1).
 * Sim units (a=1): g(z) = 2 a sin(alpha) tau/kappa = -T_p sqrt(s) w~^(3/4) sin(alpha) [E~ units], ell/a = sqrt(s) w~^(-1/4) (0.33 at s=.25,w~=5).
 * Key refinement: at fixed covered weight the BL sees the Lagrange multiplier mu = w_D(z) = (1-psi0dot)^2 (not w~): g(z) = -T_p sqrt(s) w_D(z)^(3/4) sin(alpha).
 * Data: Studies/deserno_shell/data/theory_cache.npz has F(z), psi0dot(z) for sigma~=2k/15 (k=1..15) -> use it (no 60 s solves).
TODO: phase-line shifts (2), reference energy + mini axisymmetric minimiser (4b), tilt-weight BL (4), discretisation numbers (3), md + figures.

M2 (done): first-order line-tension model (scratch lt.py -> boundary_layer.py): E_soft(z) = E_D(z) + g(z), g = -T_p sqrt(s) w_D(z)^(3/4) sin(alpha),
 w_D(z) = (1-psi0dot)^2 from the theory cache. s=0.25, p=1, NO fitted parameter:
 * E line: 6.71 / 7.11 / 8.30 / 9.10 / 10.32 / 11.13 at sigma~ 0.53 / 0.67 / 1.07 / 1.33 / 1.73 / 2.0 vs simulated 6.6 / 7.0 / 8.26 / 9.22 / 10.22 / 11.09 (agreement 0.1).
 * S1: 7.14/7.73/9.44 ... vs sim 7.0/7.4/9.0 (overshoot 0.1-1.2); S2: 5.8/5.65/5.08/3.5 vs sim 6.6/7.0/7.4/7.4 (neck regime breaks down).
 * barrier at E drops 0.67 -> 0.13 at sigma~=0.53 (sim: none); partial-branch energies: E_sim - E_first-order = +0.25 +- 0.07 (offset grows with sigma 0.17->0.29), z agrees to 0.003 mean (rms 0.05) over 263 flat-start points.
 * sigma~ = 0: first-order has NO transition, z(w~) smooth crossing z=1 at w~=4, z=1.75 at ~5.9 (sim threshold crossing 5.4) => "delayed binding at 5.4" is the upper edge of a crossover (baseline.md also says smooth).
TODO: sigma=0 z(w~) vs sim data check, reference energy E_ref, tilt-weight BL, discretisation numbers (edge lengths from old runs), md + figures.

M3 (done, numbers in scratch -> to be in boundary_layer.py): 
 * sigma=0 check vs old sim: first-order z(w~) and E(w~) match old flat data to z 0.02-0.1, E 0.05-0.15 (sim z(4.0)=0.96-1.0 vs 1.00): binding at sigma~=0 is a SMOOTH CROSSOVER centred at w~=4, not a delayed jump (5.4 = z>1.75 threshold).
 * E_ref (plane touching bead, p=1, s=.25): Gamma_ref = z_ref = s int W/(1+sx)^2 = 0.109 (p=2: .084, p=4: .063, p=8: .046; s=.1: .047 for p=1); E_ref = -w~ Gamma_ref  (O(s), subleading to the line tension O(sqrt s)).
 * Old run 50 (KI=.9,KA=.6): mean edge near contact 0.31-0.32 a, ell = 0.36 a => h/ell ~ 0.9 (under-resolved), vertex rho mean 1.0095 vs centroid 0.9931: sag h^2/6a=0.017 confirmed.
 * Discrete 1D BL (regular mesh): T_1 error vs h/ell: nodal sampling +0.18/+0.03/+0.006/+0.001 for h/ell=1/.5/.25/.125 (rel ~0.07 (h/ell)^2, lattice-phase variation 2.6% at h=ell); centroid sampling +0.036/-0.008/-0.002.
 * Tilt-weight (Y^q) BL: contact condition unchanged (P=1, C const). tau = -Theta w ell_c, Theta/ell_w = 1.259 (p, delta-hat independent) for q=0 where ell_w = int W ds over lift-off; with tilt weight Theta/ell_w = 1.3-1.7 (WORSE at equal ell_w): no gain at equal resolution; narrowing s or raising p is equivalent.
 * T1 solver (theory/softshell_direct.py) being run in background for validation of g(z) (scratch runt1.py, sigma=1.0667 s=.25 p=1).
TODO: write boundary_layer.py (all of the above, figures), line_tension_analysis.md.

M4: boundary_layer.py written (sections bl, ref, disc, tilt, lines, check; figures boundary_layer.png, line_tension_lines.png, extrapolation_protocol.png, regularisation_tradeoff.png; data/line_tension_shifts.csv; cache boundary_layer_cache.json).
 Extrapolation test (first-order model as experiment, s = .25,.16,.09,.04,.0225): quadratic fit in sqrt(s) over the 4 smallest s recovers Deserno E within 0.01, S1 0.01-0.04, S2 0.004, barrier 0.05-0.07 (linear fit with the 3 smallest s: 0.17-0.34 off).
 sigma~=0: model vs old sim z and E agree to z 0.02-0.1, E 0.05-0.2.  Partial-branch offset E_sim - E_model = +0.17 (sigma<.5) ... +0.29 (sigma~2), z_ad agrees to 0.002 mean.
TODO: write line_tension_analysis.md; check T1 solver result for validation (optional); final report.

M5 (DONE): line_tension_analysis.md written; validated against T1 data/softshell_lines.csv (section t1 in boundary_layer.py): first-order minus exact soft shell is O(s).
 Deliverables: theory/line_tension_analysis.md, theory/boundary_layer.py (+ boundary_layer_cache.json), figures/{boundary_layer,line_tension_lines,extrapolation_protocol,regularisation_tradeoff}.png, data/line_tension_shifts.csv.
