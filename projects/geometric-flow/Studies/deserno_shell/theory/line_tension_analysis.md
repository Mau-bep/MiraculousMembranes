# Why a finite-range adhesion shell shifts Deserno's phase lines (workstream F1)

Code: `theory/boundary_layer.py` (all numbers below come from it unless marked otherwise; `python boundary_layer.py all`, ~8 min).
Figures: `figures/boundary_layer.png`, `figures/line_tension_lines.png`, `figures/extrapolation_protocol.png`,
`figures/regularisation_tradeoff.png`. Data: `data/line_tension_shifts.csv`.
Notation: a = 1 unless written, kappa = KB/2, w = KI, w~ = 2 w a^2/kappa = 4 KI a^2/KB, sigma~ = 2 KA a^2/KB, E~ = E/(pi kappa),
z = 1 - cos(alpha), s = shell half width / a, delta = s a, p = shell power, W(x) = [(1 + cos(pi x))/2]^p for |x| < 1.
w_D(z) = (1 - psi0dot(z))^2 is Deserno's "equation of state" (dE_D/dz = w_D(z) - w~).

## 0. Summary

1. **The contact condition is unchanged at leading order.** For ANY local adhesion weight U(d, d') with U(0,0) = -w and U = 0 away from the
   shell (symmetric or not, with or without hard core, any p, with or without a tilt weight) the boundary layer has the first integral
   `kappa d' d''' - d' U_{d'} - (kappa/2) d''^2 + U = C`, C = -w on the bound side, C = -(kappa/2) Dc^2 on the free side, hence
   `w = (kappa/2) Dc^2`, i.e. `w~ = (1 - psi0dot)^2` (Dc = curvature jump). Numerically: the BVP solved with Dc as the unknown returns
   Dc = 1.00000000 (units sqrt(2w/kappa)) for p = 1, 2, 4, 8 from any starting guess; C constant to 4e-11 (p = 1). Corrections are O(ell/a).
2. **Exact scaling of the excess line energy (flat substrate, small slope):** the boundary layer problem has NO free parameter once lengths are
   scaled with `ell = sqrt(delta/Dc) = sqrt(delta * ell_c)`, `ell_c = 1/Dc = sqrt(kappa/(2w)) = a/sqrt(w~)`. Therefore
   `tau_eff = - T_p * w * ell`, with pure numbers **T_1 = 1.198, T_2 = 1.029, T_4 = 0.875, T_8 = 0.740** (T_p ~ p^(-1/4)).
   tau_eff < 0 always (variational proof: the sharp profile inside the shell is a trial function, T_p >= 0.972, 0.838, 0.714, 0.604).
   It scales as `sqrt(delta)`, NOT as delta, and is NOT proportional to the lift-off length ell_c.
   In simulation variables: `tau_eff = - T_p KI a sqrt(s) w~^(-1/4)`; in E~ units the contact-line term for a contact circle a sin(alpha) is
   `g(z) = 2 a sin(alpha) tau_eff/kappa = - T_p sqrt(s) w_D(z)^(3/4) sin(alpha)`  (see 2.1 for why w_D(z) and not w~).
3. **Consequence for the phase diagram (first order in ell/a, no fit parameter), s = 0.25, p = 1:** E line moves UP by 1.6 (sigma~ = 0.53) ...
   2.0 (sigma~ = 1.07) in w~, the E line is reproduced to 0.1 over the whole range sigma~ = 0.5 ... 2 of the old simulation; the unwrapping barrier
   drops from 0.67 to 0.13 (sigma~ = 0.53); S1 and S2 approach each other; at sigma~ = 0 there is no transition at all but a smooth crossover centred
   at w~ = 4. The full continuum soft-shell solver of T1 confirms this: the first-order model errs by +0.14 (E), +0.41 (S1), -0.50 (S2) at s = 0.2
   and the errors halve when s halves (they are O(ell^2)).
4. **At s = 0.25 the boundary layer is NOT thin:** ell/a = sqrt(s) w~^(-1/4) = 0.35 (w~ = 4), 0.33 (5), 0.31 (7); tilt of the membrane at the outer edge of the
   shell ~47 degrees. First order is a leading asymptotics at best here (the observed accuracy is nevertheless 0.1 in E).
5. **The simulated "delayed binding at sigma~ = 0 (5.4 vs 4)" is an artefact of the definition**: with the shell the transition is a smooth crossover
   with z = 1 at w~ = 4 exactly (simulated z = 0.96-1.0 there) and the quoted S1 = 5.4 is where z crosses the threshold 1.75 (first-order model: 5.9).
6. **A steeper shell or a tilt-weighted shell gives no gain at equal mesh cost**: p = 8 at s is equivalent to p = 1 at s' = 0.38 s (same tau, same
   resolution requirement); the tilt weight Y^q is worse (Theta/ell_w 1.3-1.7 vs 1.26). Recommended path: finite-resolution convergence at h/ell <= 0.3 and
   extrapolation in sqrt(s) (3.4), which recovers Deserno's lines from soft-shell data to 0.005-0.07.

## 1. The contact-line boundary layer

### 1.1 Set-up and first integral

Locally (distance from the contact line << a) the bead is a plane, the contact line is straight, and the membrane is the graph d(x) of its height above
the bead surface along the coordinate x normal to the line. Energy per unit length of the contact line (small slopes; weights counted per projected
(radial) area, which is what `E_ad = -W a^2 sum_f w(r_f) |omega_f|` does, |omega_f| = dA * Y/rho^2):

    e[d] = int dx [ (kappa/2) d''^2 + U(d) ],     U(d) = -w W(d/delta),  W(0) = 1 = max, W symmetric, W'(0) = 0.

Euler-Lagrange: `kappa d'''' + U'(d) = 0`. The Lagrangian has no explicit x, so the Ostrogradsky (higher-derivative Noether) integral is conserved:

    H = d' (dL/dd' - (dL/dd'')') + d'' dL/dd'' - L = - kappa d' d''' + (kappa/2) d''^2 - U(d)
    =>   kappa d' d''' - (kappa/2) d''^2 + U(d) = C.

* Bound side (x -> -inf): d = 0 is an exact solution (U'(0) = 0); d' = d'' = d''' = 0, U = -w: **C = -w**.
* Free side (x -> +inf, d >> delta): U = 0 and the outer membrane has constant curvature on the scale of the layer, d = (Dc/2)(x - xc)^2 + h,
  d'' = Dc, d''' = 0: **C = -(kappa/2) Dc^2**.
* Equating: **w = (kappa/2) Dc^2.** The shell shape W and its width enter neither C nor the condition: the contact condition (Seifert-Lipowsky;
  Deserno's w~ = (1 - psi0dot)^2) holds for the finite-range shell. Same for any U(d, d') (tilt weights): H gains the term -d' U_{d'}, which vanishes in
  both asymptotic regions.
* Hard core: if the membrane cannot enter d < 0 the contact region is d = 0 identically, d' = 0, so the wall reaction does no work in H: unchanged.
  No hard core, symmetric shell (the simulation): the bound solution is d = 0; the layer undershoots (the membrane dips INTO the bead range) by
  d_min/delta = -0.129 (p = 1), -0.091, -0.064, -0.045 (p = 2, 4, 8) at ~0.68 ell (p = 1; 0.42 ell for p = 8) before the vertex; the linearisation
  d'''' = -k^4 d, k^4 = p pi^2/(4 delta^2) (in units kappa Dc^2 / 2 = w, Dc = 1) gives the damped oscillation exp(-k x/sqrt 2) cos(k x/sqrt 2 + ..).
* Numerical confirmation (`section_bl`): BVP on [-30, 8] (scaled), tol 1e-8, Dc as a free parameter, pinned slope d'(L2): Dc = 1.00000000 for p = 1, 2, 4, 8 and for
  starting guesses 0.5, 1, 2; spread of C: 4e-11.

### 1.2 Sphere curvature and the exact statement

On the sphere the bound state is the sphere itself: for the pure bending energy the sphere is a Willmore-critical surface (E = 8 pi kappa
independent of radius), so it exerts no normal force and the bound zone sits at d = 0. A tension sigma pushes it to d0 = 2 sigma/(a U''(0)), U''(0) = w |W''(0)|/delta^2 ->
d0/a = 8 sigma~ s^2 /(p pi^2 w~) ~ 0.005 at sigma~ = 0.5, w~ = 5, s = 0.25, p = 1 (negligible).
The local analysis above is valid up to corrections of relative order ell/a coming from (i) the variation of the second principal curvature c2 across the layer,
(ii) the curvature of the bead surface, d/a ~ s, in c1 = 1/a - d/a^2 - d_ss, (iii) the slope nonlinearity of the bending energy and (iv) force terms (shear force, tension) that
are continuous across the layer and enter H at order (slope) x (force)/w ~ ell/a. None of them is a rigorous bound; their combined size is measured in
2.3 against the exact continuum solver of T1.

The statement that is exact in the continuum soft-shell model (T2/T1): at fixed covered weight Gamma = (1/2 pi) int W Y_+ dOmega the elastic energy e(Gamma) is
minimised, w~ = de/dGamma is the Lagrange multiplier (the analogue of w_D(z)), S1 = max, S2 = min of w~(Gamma), E = Maxwell construction. The "contact condition" is then
the local statement that the multiplier equals (kappa/2) Dc^2 at the outer edge of the layer, to O(ell/a).

### 1.3 Excess energy of the contact line, tau_eff

Definition (no freedom): the outer solution is the parabola d = (Dc/2)(x - xc)^2 + h; xc is its vertex, h the height of the vertex above the bound plane.
The sharp reference with the same outer field has d = d' = 0 at its contact point; this fixes the contact point = xc (continuity of d and d'; any other
choice gives a reference that is not a sharp solution and differs by 2 w Delta x, the jump between the bound density -w and free density +w). Then

    tau_eff = lim int_{-L1}^{L2} [ (kappa/2) d''^2 + U ] dx - [ -w (xc + L1) + (kappa/2) Dc^2 (L2 - xc) ].

It is independent of L1, L2 and of the pinning (checked to 1e-4). The vertex offset is small: h/delta = 0.0188 (p = 1), 0.018, 0.015, 0.011;
its energy cost (shear force x h) is ~3% of tau (estimated), neglected.

**Exact scaling.** With d = delta D, x = ell t, ell = sqrt(delta/Dc), w = kappa Dc^2/2 the problem becomes `D'''' = W'(D)/2`, `D''(inf) = 1`: no parameter
(the quadratic far field ties length^2 to height). Hence tau_eff = -T_p w ell exactly (flat model), in all regimes delta << ell_c and delta >> ell_c
(the two-parameter dimensional analysis f(delta/ell_c) is a pure power sqrt). Values (`section_bl`):

| p | T_p | trial bound int_0^sqrt2 W(t^2/2) dt | ell_w/ell | T_p/ell_w | d_min/delta | h/delta |
|---|-----|------|------|------|------|------|
| 1 | 1.1982 | 0.9716 | 0.9536 | 1.2565 | -0.129 | 0.0188 |
| 2 | 1.0289 | 0.8379 | 0.8171 | 1.2591 | -0.091 | 0.0184 |
| 4 | 0.8746 | 0.7138 | 0.6939 | 1.2605 | -0.064 | 0.0149 |
| 8 | 0.7395 | 0.6042 | 0.5863 | 1.2611 | -0.045 | 0.0113 |

(ell_w = int over the lift-off side of W(d/delta) dx, the width of the lift-off zone that still adheres.) Observations:
* **Sign: negative, rigorous within the model.** Inserting the sharp profile (d = 0, then Dc x^2/2) into the soft functional lowers the energy by
  w int_0^inf W(Dc x^2/2 delta) dx = w ell * 0.9716 (p = 1); the minimiser is lower still: T_p >= trial bound, T_p/trial = 1.23.
* **Scaling:** tau ~ -w sqrt(delta ell_c) (not w delta, not w ell_c). Steeper shells: T_p ~ p^(-1/4) (Gaussian limit: effective width delta ~ 1/sqrt p).
* **Virial identity** (from tau ~ sqrt(delta) + Hellmann-Feynman): `T_p = -2 int D W'(D) dt` (checked: 1.19821 vs 1.19817).
* **T_p/ell_w = 1.259 +- 0.003 for all p and all delta/ell_c** (28 cases): tau = -1.26 w ell_w, a near-universal relation, used in 4.
* tau depends on kappa, w only through w and Dc = sqrt(2w/kappa): tau = -T_p w^(3/4) (kappa/2)^(1/4) delta^(1/2); in
  simulation variables: **tau_eff = - T_p KI a sqrt(s) w~^(-1/4)**, w~^(-1/4): at fixed KI a more strongly curved (stiffer) membrane lifts off over a shorter distance.

### 1.4 Is a functional with negative line tension valid? Thin layer?

* tau_eff < 0 is not an inconsistency: it is the excess of a regularised contact, and the full layer energy is bounded below (the minimiser exists, T_p is finite).
  A sharp-contact functional with a negative line tension would be unbounded below (create contact line), so it is valid only (a) as an effective description at
  length scales >> ell (contact-line curvature radius >> ell: a sin(alpha) >~ 3 ell, and the line-undulation wavelength lambda_u >~ 3 ell), (b) with the tau of the local state:
  tau depends on the local curvature jump, therefore the correct functional is E_D(z) + 2 pi a sin(alpha) tau(w_D(z)) (see 2.1), not a constant tau.
  Stability against non-axisymmetric line undulations was not analysed (a crude estimate, stiffness kappa Dc vs |tau| k, suggests stability for k ell_c <~ 2/Theta, i.e. wavelengths >~ ell).
* Thinness: ell/a = sqrt(s) w~^(-1/4). s = 0.25: 0.38 (w~ = 3), 0.35 (4), 0.33 (5), 0.31 (7), 0.28 (10). The lift-off length itself is ell_c/a = 1/sqrt(w~) = 0.58 (3), 0.45 (5), 0.32 (10).
  Slope at the outer edge of the shell: d' = sqrt(2 delta/ell_c) = sqrt(2 s sqrt(w~)) = 1.06 (47 degrees) at w~ = 5. Requirement ell/a < 0.1 and edge slope < 0.3:
  s < 0.01 sqrt(w~) = 0.022 (w~ = 5) and s < 0.02. At s = 0.25, ell/a ~ 0.33 and the first-order theory is not a controlled asymptotics (but see 2.3).
* Validity in z: the line-tension picture needs the contact circle radius a sin(alpha) > ell, i.e. z in (0.056, 1.944) at s = 0.25, w~ = 5 (and > 2 ell: z in (0.25, 1.75)).

## 2. Consequences for the phase diagram

### 2.1 First-order functional

At fixed covered weight Gamma the constrained minimiser has the multiplier mu = de/dGamma = w_D(z) acting as the effective adhesion in the boundary layer (the layer is a
stationary solution only if w = (kappa/2) Dc^2, and at fixed Gamma this is what the multiplier provides). Hence in the layer w -> mu = kappa w_D/(2 a^2), Dc = sqrt(w_D)/a, and

    E_soft(Gamma; w~) = E_D(Gamma; w~) + g(Gamma),   g = 2 a sin(alpha) tau(mu)/kappa = - T_p sqrt(s) w_D(Gamma)^(3/4) sqrt(z (2 - z))   [pi kappa]

(Gamma is both the soft covered weight and the sharp z up to O(ell); the weight outside the contact line, -mu * (excess weight), is part of tau.)
Properties: g < 0, g(0) = g(2) = 0, g ~ -sqrt(z) at z -> 0 and ~ -sqrt(2 - z) at z -> 2 (infinite slope: the picture fails where a sin(alpha) ~ ell).
At stationary points of E (equilibrium states, w_D = w~) Hellmann-Feynman gives the first-order shift of the energy as g at the unperturbed state, with the equilibrium
tau; this is what controls the E line and the barrier.

Equation of state: w~ = w_D(z) + g'(z) = w_soft(z); S1 = first max, S2 = following min of w_soft; the state is stable where w_soft increases. E line: lowest-z local minimum of E vs
highest-z local minimum (the enveloped state is displaced to z_e < 2: g' diverges at z = 2, the fully closed state is never a minimum: z_e = 1.86-1.99), barrier = E(z_b) - E(z_e).

### 2.2 Signs and sizes (s = 0.25 first; full table below)

* **E line: up.** The partially wrapped state (z_p ~ 1-1.2) is lowered by |g(z_p)| ~ 1-2 while z = 2 is not: Delta w~_E ~ |g(z_p)|/(2 - z_p).
* **S1: up, not down**, although at fixed z < 1, g' < 0 lowers w_soft: the maximum of w_soft is displaced to z_S1 ~ 1.4-1.6 (Deserno: 0.8-1.1) where g' > 0 (g' > 0 for z > 1).
* **S2: up, strongly**: g' grows without bound towards z = 2.
* **Barrier: strongly reduced** (not exactly removed at first order). sigma~ = 0.53, s = 0.25: 0.67 -> 0.13 (p = 1), 0.21 (p = 8).
* **No E line** (no hysteresis loop) at sigma~ <~ 0.3 for s = 0.25 (S1 and S2 coincide within 0.1 at sigma~ = 0.27); at sigma~ = 0.133 the model has no N-shaped w_soft (smooth crossover) for s >= 0.1 (p = 1), a small one appears at s = 0.05.

Shifts Delta = soft - Deserno in w~ (Deserno: sigma~ = 0.267: E 4.556, S1 5.061, S2 3.527, barrier 0.403; 0.533: 5.125 / 5.964 / 3.204 / 0.671;
1.067: 6.289 / 7.674 / 2.678 / 1.089). "bar" is the barrier at the shifted E line. ell/a at w~ = 5: s = 0.25: 0.33, 0.10: 0.21, 0.05: 0.15, 0.025: 0.11.

| sigma~ | p | s=0.25: dE dS1 dS2 bar | s=0.10: dE dS1 dS2 bar | s=0.05: dE dS1 dS2 bar | s=0.025: dE dS1 dS2 bar |
|---|---|---|---|---|---|
| 0.267 | 1 | - 0.92 2.34 - | 0.70 0.32 1.51 0.04 | 0.43 0.12 1.08 0.08 | 0.27 0.03 0.78 0.14 |
| 0.267 | 4 | 0.87 0.46 1.73 0.02 | 0.45 0.13 1.12 0.08 | 0.28 0.03 0.80 0.13 | 0.18 -0.01 0.57 0.19 |
| 0.267 | 8 | 0.68 0.30 1.48 0.04 | 0.36 0.07 0.95 0.10 | 0.22 0.01 0.68 0.16 | 0.15 -0.02 0.49 0.21 |
| 0.533 | 1 | 1.59 1.18 2.59 0.13 | 0.82 0.44 1.64 0.21 | 0.51 0.19 1.16 0.28 | 0.33 0.07 0.82 0.36 |
| 0.533 | 2 | 1.28 0.86 2.23 0.15 | 0.67 0.31 1.41 0.24 | 0.42 0.13 1.00 0.31 | 0.27 0.05 0.71 0.39 |
| 0.533 | 4 | 1.01 0.61 1.89 0.18 | 0.53 0.20 1.20 0.27 | 0.34 0.08 0.85 0.35 | 0.22 0.02 0.60 0.42 |
| 0.533 | 8 | 0.79 0.41 1.60 0.21 | 0.43 0.13 1.02 0.31 | 0.27 0.05 0.72 0.39 | 0.18 0.01 0.51 0.45 |
| 1.067 | 1 | 2.01 1.76 2.40 0.45 | 1.06 0.68 1.52 0.53 | 0.67 0.34 1.08 0.63 | 0.43 0.17 0.76 0.72 |
| 1.067 | 2 | 1.62 1.29 2.06 0.47 | 0.86 0.50 1.31 0.57 | 0.55 0.25 0.93 0.67 | 0.36 0.12 0.66 0.76 |
| 1.067 | 4 | 1.29 0.92 1.75 0.50 | 0.70 0.36 1.11 0.62 | 0.45 0.18 0.79 0.71 | 0.30 0.09 0.56 0.80 |
| 1.067 | 8 | 1.02 0.65 1.49 0.54 | 0.56 0.26 0.94 0.66 | 0.37 0.13 0.67 0.76 | 0.24 0.07 0.47 0.84 |

(sigma~ = 0.133 and the rest of the grid are in `data/line_tension_shifts.csv`; the cache holds sigma~ = 2k/15, so sigma~ = 1.067 instead of 1.0.) The shifts are close to linear in
sqrt(s) (slightly convex, `figures/line_tension_lines.png` right panel); S2 converges slowest (g' is singular at z = 2).

### 2.3 Tests of the first-order model

(a) **Old simulation (shell 0.25, p = 1, no fitted parameter)**: E line from the model vs simulation, sigma~ = 0.53 / 0.67 / 1.07 / 1.33 / 1.73 / 2.0:
6.71 / 7.11 / 8.30 / 9.10 / 10.32 / 11.13 vs 6.6 / 7.0 / 8.26 / 9.22 / 10.22 / 11.09 (differences <= 0.12; Deserno: 5.12 ... 8.39). S1: 7.14 / 7.73 / 9.44 / 10.55 vs 7.0 / 7.4 / 9.0 / 9.8
(model high by 0.1-0.8). S2: model 5.8 / 5.65 / 5.08 / 3.5 vs simulation 6.6 / 7.0 / 7.4 / 7.4: wrong trend at large sigma~ (neck regime, a sin(alpha) ~ 1.3 ell at z_S2 ~ 1.9).
263 flat-start partial-branch states (baseline_partial_branch.csv): z_ad(sim) - z(model) = 0.002 (rms 0.047); E_sim - E_model = +0.25 +- 0.08 (grows from +0.17 at sigma~ < 0.5 to +0.29 at
sigma~ ~ 2; Deserno without shell is off by 1.5-5.7). The +0.2 residual is NOT explained by this analysis (candidates: unrelaxed remeshing energy kick, discrete bending of the cap,
O(ell^2) terms).
(b) **sigma~ = 0**: model z(w~) and E(w~) vs the old flat start (`section_check`): w~ = 2.0 / 3.6 / 5.2 / 6.8 / 10: z_sim 0.42 / 0.77 / 1.69 / 1.92 / 1.98 (z_ad 0.28 / 0.63 / 1.59 / 1.87 / 1.96), model z 0.24 / 0.77 / 1.58 /
1.86 / 1.96; E_sim - E_model = +0.14 / +0.20 / +0.11 / +0.06 / +0.05.
(c) **Exact continuum soft shell (T1, `data/softshell_lines.csv`, snapshot of 10 Oct)**, p = 1: first-order minus T1 (S1 / E / S2):
s = 0.20, sigma~ = 0.533: +0.41 / +0.14 / -0.50; s = 0.10: +0.21 / +0.09 / -0.28; s = 0.05: +0.105 / +0.049 / -0.149; sigma~ = 0.267: s = 0.1: +0.12 / +0.06 / -0.06, s = 0.05: +0.057 / +0.030 / -0.041;
sigma~ = 0.133, s = 0.025: +0.037 / +0.026 / 0.000. The error is proportional to s (halving s halves it, i.e. O(ell^2)): the first-order model is the correct leading asymptotics of the
sphere problem, and at s = 0.25 it is accurate to ~0.15 in E and ~0.5 in S1, S2. T1's exact values at sigma~ = 0.533, s = 0.25: S1 6.65, S2 6.37, E 6.55 (simulation 7.0, 6.6, 6.6): the barrier
window shrinks from 2.76 (Deserno) to 0.28.

### 2.4 Small z, the neck, the reference state

* **Planar reference.** A plane tangent to the bead sits at rho = sqrt(1 + r^2) >= 1 and gains adhesion inside the shell: the weighted covered fraction is
  `z_ref = s int_0^1 W(x) (1 + s x)^(-2) dx`, `E_ref = - w~ z_ref` (pi kappa). Values: p = 1: 0.109 (s = 0.25), 0.047 (0.10), 0.024 (0.05); p = 2: 0.084, 0.036, 0.018; p = 4: 0.063, 0.026, 0.013; p = 8: 0.046, 0.019, 0.010.
  E_ref(w~ = 5, s = .25, p = 1) = -0.55. It scales ~ s (the line tension ~ sqrt(s)): subleading but not negligible at s = 0.25. The membrane then bends towards the bead (second-order gain, not computed).
* The first-order g(z) -> 0 for z -> 0, so it misses E_ref and overestimates the gain at small z: at (sigma~ 0.133, w~ = 2): model z = 0.23, E = -0.64, simulation z_ad = 0.28, E = -0.47.
  There dE/dz(0+) = -infinity in the model: **there is no barrier to leave the unbound state; the unbound state is already adhesive** (z ~ 0.2-0.4 at w~ = 2-3 in the simulation, model 0.23-0.37).
* Where a sin(alpha) <~ ell the boundary layers on the two sides of the contact circle overlap with the pole (small z) or neck (z -> 2); only the full solver (T1/T2) is valid. The model's S2 failure and the
  simulated S2 plateau (7.4) are neck effects; the fully enveloped state keeps z_e = 1.86-1.99 (not 2).
* Effect on S1 / W: S1 of the soft model is set at z_S1 ~ 1.5 (the partial branch lost beyond the equator), not by the small-z patch, so E_ref hardly matters for S1; it shifts the unbound end of the energy
  balance down by E_ref, the enveloped end by whatever the neck gains (unknown; the agreement of the E line with T1 to 0.15 says the net effect is <= 0.15 in w~).

### 2.5 The "delayed binding at sigma~ = 0" (item 4b)

Deserno at sigma~ = 0: E~ = -(w~ - 4) z for all z: w_D = 4, W line = w~ = 4, degenerate. With the shell, w_soft(z) = 4 + g'(z), g' = -A (1 - z)/sqrt(z (2 - z)), A = T_p sqrt(s) 4^(3/4) = 1.69 (p = 1, s = .25):
monotonic, no max/min: continuous binding with z(w~ = 4) = 1 exactly (g'(1) = 0), z = 1.43 (w~ = 4.8), 1.69 (5.6), 1.75 at w~ = 5.9. Simulation at sigma~ = 0: z = 0.96-1.0 at w~ = 4.0, 1.50 at 4.8, 1.78 at 5.6. The baseline's S1_sim = S2_sim = 5.4
(threshold z > 1.75 on a smooth curve, flagged "soft" in baseline.md) is therefore the upper end of a crossover whose centre is w~ = 4, not a delayed transition. The same reading applies for sigma~ = 0.133.

## 3. Discretisation

Triangulated membrane, edge length h, weight sampled at face centroids r_f, faces whose normal points to the bead.

**(a) Centroid sampling.** The weight W is a smooth function along the layer (variation length ell_w ~ 0.6-0.95 ell): evaluating it at centroids is a midpoint rule: relative error O((h/ell)^2) with coefficient 0.03-0.1 (see (d)).

**(b) Sag.** A face of edge h inscribed in the sphere (vertices at rho = R_v) has its centroid at rho_c = sqrt(R_v^2 - h^2/3) = R_v - h^2/(6 R_v). Measured in the old run 50 (final state, KA = 0.6, KI = 0.9; faces with weight > 0.9, mean edge 0.32):
vertex rho = 1.0095, centroid rho = 0.9931, difference 0.0164 vs h^2/6 = 0.0171. Effects: (i) weights at the centroid W(-sag/delta): loss (pi^2 p/4)(sag/delta)^2 = 1% at sag/delta = 0.066, i.e. a shift of d by sag, equivalent to a sphere of radius a + h^2/6a for the
membrane (the adhesion energy is |omega| which depends only on directions: no change in the bound energy; the bending energy of a sphere is radius-invariant) -> no first-order effect on the bound zone; (ii) the bound
curvature seen by the mesh is 1/R_v, not 1/a: Delta c_bound = -h^2/(6 a^3) -> Delta(Dc) = 1/a - 1/R_v = +0.009 (R_v = 1.0095 measured; at most +0.017 = h^2/6a), Delta w~ = 2 w~^(1/2) Delta(Dc) = +0.04 (at most +0.08) for h = 0.32, w~ = 5; scales as h^2. Resolution criterion h^2/(6a) <~ 0.1 delta.
The discrete curvature of Bending_tan on irregular meshes is not exactly 1/R_v (not measured; test: bending energy of the Planar_full.obj bead cap vs 4 pi kappa z).

**(c) Selection by face normal.** The condition Y > 0 (not differentiated) is non-binding inside the layer (tilt reaches 47 degrees at the outer shell edge) and only removes back-side faces; the normal of an inscribed face of a regular triangle passes
through the centre (Y = 1 for bound faces), irregular triangles deviate by O(h skew/a). No contribution at the order considered.

**(d) Edge length vs ell: ideal 1D mesh** (`section_disc`, p = 1, continuum T_1 = 1.1982): discrete bending (second differences) + sampling of W at nodes or at chord midpoints, vertex of the outer parabola at lattice phase 0 or 1/2:

| h/ell | nodal sampling T_1 (phase 0 / 1/2) | rel. error | midpoint sampling T_1 | rel. error |
|---|---|---|---|---|
| 1.0 | 1.374 / 1.343 | +14.7% / +12.1% | 1.234 / 1.234 | +3.0% |
| 0.5 | 1.227 / 1.219 | +2.4% / +1.7% | 1.190 / 1.200 | -0.7% / +0.1% |
| 0.25 | 1.2037 / 1.2027 | +0.5% | 1.196 / 1.197 | -0.2% |
| 0.125 | 1.1992 / 1.1994 | +0.1% | 1.198 | 0.0% |

Error ~ 0.075 (h/ell)^2 (nodal), ~0.03 (h/ell)^2 (midpoint); lattice-phase (pinning) variation 2.6% at h = ell, < 0.2% at h = 0.25 ell. **Criterion: h/ell <= 0.3 gives <= 1% in tau on a regular mesh; allow a factor 3 for a
triangulated, remeshed surface (not measured): h/ell <= 0.2-0.3.** In the old runs the mesh was NOT resolved: mean edge near the contact 0.31-0.32 a (run 50: bound faces 0.32, shell region 0.27-0.35, size_max 0.5) with ell = 0.36 a (w~ = 3.6) -> h/ell ~ 0.9, i.e. an ideal-mesh error of 6-12% in tau
and a mesh with only ~3 edges across the layer. The resulting additional offset in the line positions is of order the error fraction times the shift (~0.1-0.2 in w~) plus unknown irregular-mesh effects.

**Order of convergence.** (1) h -> 0 at fixed s (continuum soft-shell model, which T1 solves directly; the error is ~ (h/ell)^2 per ell), (2) then s -> 0 (the shifts vanish linearly in ell). With a finite budget keep
h/ell fixed while s decreases: the discretisation error is then a fixed fraction of the shift and extrapolates away together with it.

**Extrapolation protocol (3.4).** Choose s_i, ell_i = sqrt(s_i) w~^(-1/4) a, mesh h_i = r0 ell_i with r0 = 0.25 (second series r0 = 0.15 as a check), e.g. s = 0.25, 0.16, 0.09, 0.04, 0.0225 (ell/a = 0.33, 0.27, 0.20, 0.13, 0.10 at w~ = 5; h = 0.08 ... 0.025).
Measure the lines (S1, S2, E, barrier; or the energy of partial/enveloped states at fixed w~) at each (s_i, h_i) and fit y(sqrt s) = y_0 + y_1 sqrt(s) + y_2 s (quadratic through the 4 smallest s; linear fits are not accurate enough). Test (first-order model as pseudo-data, `section_extrapolation`):
sigma~ = 0.533: E 5.116 (Deserno 5.125), S1 5.951 (5.964), S2 3.207 (3.204), barrier 0.61 (0.67); sigma~ = 1.067: E 6.285 (6.289), S1 7.709 (7.674), S2 2.682 (2.678), barrier 1.045 (1.089) (quadratic fit through s = 0.16, 0.09, 0.04, 0.0225).
Test with exact continuum soft-shell data of T1 (sigma~ = 0.533, s = 0.2, 0.1, 0.05; quadratic through the three points in sqrt(s)): S1 5.943 (Deserno 5.953), E 5.123 (5.118), S2 3.14 (3.21): errors 0.01, 0.005, 0.07. A linear fit through the two smallest s misses by 0.25-0.2 (E 4.92, S1 5.70): do not use it.
The protocol also gives the cheap consistency check that the line shifts at fixed r0 are proportional to ell (y_1) with the sign and the size predicted in 2.2.

## 4. Steeper shells and tilt weights

**4.1 p = 2, 4, 8 at fixed s.** tau = -T_p w ell: T_p ratios 1 : 0.86 : 0.73 : 0.62 (p = 1, 2, 4, 8). The shifts of the lines follow: at sigma~ = 0.533, s = 0.25: dE = 1.59 / 1.28 / 1.01 / 0.79 (p = 1 / 2 / 4 / 8),
dS1 = 1.18 / 0.86 / 0.61 / 0.41. The price: the weight varies over ell_w = (0.95, 0.82, 0.69, 0.59) ell, the bound-side oscillation length is ~1.1 ell p^(-1/4): the resolution requirement h <~ 0.3 ell_w grows as p^(1/4)... At equal ell_w/h
the comparison is exact: tau = 1.26 w ell_w for every p, i.e. at fixed ell_w/h tau/(w h) is the same for all p. **A steeper shell is equivalent to a narrower s at the same mesh cost**: ell(p) = ell_w(p)/0.95 * ell(p = 1, s) ~ ell(p = 1, s') with s' = s (ell_w(p)/ell_w(1))^2 = 0.74 s, 0.53 s, 0.38 s (p = 2, 4, 8); check with the table:
p = 8, s = 0.25 gives dE = 0.79 (sigma~ = 0.53), p = 1, s = 0.10 gives 0.82. No free lunch; the choice is only convenience (p large keeps delta large, which helps the sag criterion h^2/6a <~ 0.1 delta, and the planar reference E_ref is smaller by the same factor
(0.046 vs 0.109 at s = 0.25, p = 8 vs 1)).

**4.2 Tilt weight.** Weight W(d/delta) (1 + d'^2)^(-q/2) = W Y^q (Y = cos(tilt) = -n.rhat, additional to the existing solid-angle factor). `solve_tilt` solves the full nonlinear EL equation. The contact condition is unchanged (Dc = 1.000000, C constant to 1e-11, p = 1, delta/ell_c = 0.56, q = 0, 2, 4, 8, 16).
Result (grid: p = 1, 2, 4, 8; delta/ell_c = 0.1-1.1; q = 0, 4, 8, 16; 110 solutions, `figures/regularisation_tradeoff.png`): tau reduces strongly (delta/ell_c = 0.56, p = 1: Theta = 0.90, 0.78, 0.69, 0.57, 0.44 for q = 0, 2, 4, 8, 16; Theta = -tau/(w ell_c)), BUT the lift-off width reduces with it:
**Theta/ell_w = 1.259 (q = 0) vs 1.39 (q = 4), 1.47 (q = 8), 1.56 (q = 16) on average (up to 1.7)**. At equal ell_w (equal resolution need) the tilt weight gives a LARGER tau than the plain distance cutoff. For delta -> infinity with q fixed the weight is limited by the tilt alone (delta/ell_c = 1.1: ell_w = 0.26-0.54 ell_c, Theta = 0.42-0.82, weakly dependent on p).
It also reduces the planar reference gain (cos^q of the tilt r/rho of the touching plane) but needs the derivative of the face normal in the force (not differentiated today). **What I would implement: nothing new in the weight. Use the plain shell (p = 1 or 2), choose s as small as the
mesh allows (h/ell <= 0.3, i.e. h <= 0.3 sqrt(s) w~^(-1/4) a), sample the weight at the centroid (unchanged), and extrapolate in sqrt(s) (3.4).** If one insists on a smarter weight, the boundary-layer analysis says to look for a weight with smaller tau per ell_w than 1.26 w ell_w - none of the local weights tried (distance, power, tilt)
reaches it; the minimum over local weights is not known (conjecture: it needs U with a different shape of the tail, e.g. a negative lobe, which changes the sign structure).

## 5. Honest assessment

**Established (rigorous in the flat-substrate, small-slope boundary-layer model, numerically confirmed):** (1) the first integral and the invariance of the contact condition w = (kappa/2) Dc^2 for any local U(d, d') that vanishes outside the shell, with or without hard core;
(2) the exact scaling tau = -T_p w sqrt(delta/Dc), the numbers T_p and the near-universal relation tau = -1.26 w ell_w; (3) tau < 0 (variational); (4) the regular-mesh discretisation errors; (5) the algebra of the first-order functional E_D + g and the signs: E up, S2 up, S1 up (maximum displaced to z > 1), barrier down.

**First-order / asymptotic (leading order in ell/a):** g(z) = -T_p sqrt(s) w_D^(3/4) sin(alpha); validated against the exact continuum solver of T1 (errors O(s): 0.15 in E, 0.4-0.5 in S1, S2 at s = 0.25 and halving with s), and against the old simulation (E to 0.1, z to 0.05, constant energy offset +0.17-0.29 unexplained). At s = 0.25 (ell/a = 0.33, edge slope 47 degrees) this is
a controlled-looking but not controlled approximation: the sphere-curvature, slope and force corrections to tau were not computed analytically. The extrapolation protocol (3.4) relies on the leading-order structure and has been verified only on first-order pseudo-data and on three-point T1 data.

**Conjecture / not shown:** (i) that irregular triangulated meshes obey the same (h/ell)^2 law with a coefficient within a factor of 3 of the regular 1D mesh; (ii) that the unexplained +0.2 energy offset is a mesh/relaxation effect; (iii) the neck (z -> 2) physics: the plateau of S2_sim at 7.4 is NOT reproduced by the first-order model (it gives S2 falling from 5.8 to 3.5 with sigma~); the continuum solver of T1 gives S2 = 6.37 at sigma~ = 0.53 (simulation 6.6); (iv) non-axisymmetric stability of a negative-tension contact line; (v) that the
planar reference and the neck compensate in the E line.

**Main caveats:** (1) Dc = sqrt(w~)/a = 2.2/a at w~ = 5 is twice the bead curvature: the free membrane is strongly curved at the contact, ell_c = 0.45 a is not small, so the "thin layer" hypothesis fails already for the sharp problem at the relevant w~; (2) small z and the neck: the contact-line picture requires a sin(alpha) >~ ell; g has an infinite slope at z = 0, 2 and the enveloped state of the model sits at z_e < 2; (3) finite size: the sim disc R = 10 only matters for sigma~ -> 0 where
the tail ~1/r is truncated (T1: R-independent at sigma~ = 0.53); (4) the simulation mesh was unresolved (h/ell ~ 0.9); old data therefore combine the physical shell effect (captured to ~0.15 by the continuum model) with discretisation error of unknown size; (5) tau was derived for the sim's weight per projected area, plane geometry.
