# B0 baseline: old planar simulations in Deserno's variables

Script: `sim/extract_lines.py` (functions + CLI, `python sim/extract_lines.py`), theory pre-solved by
`sim/theory_cache_build.py` (`data/theory_cache.npz`, F(z) and psi_dot_0(z) at sigma~ = 2k/15).
Outputs: `data/baseline_lines.csv` (one row per KA), `data/baseline_partial_branch.csv`, figures `baseline_vs_deserno.png`,
`baseline_offsets.png`, `baseline_energy.png`.
Data: flat start `Coverage_data.txt`, full start `Coverage_data_full.txt` (old campaign, shell 0.25, size_min 0.05 / size_max 0.5).
Conversions: w~ = 4 KI (KB = a = 1), sigma~ = 2 KA, E~ = 2 Total_E/(pi KB), z = 2 CoverageUnion.

## 0. Data handling and quality (facts)

* 416 (KA, KI) points per file, 16 KA (0 ... 1, step 1/15; sigma~ 0 ... 2, NOT 0 ... 1 as CONTEXT.md says) x 26 KI (0.5 ... 3.0, step 0.1,
  w~ = 2 ... 12, step 0.4). No duplicates, no non-finite rows, no missing grid point (unfinished runs were already skipped by PostProcessing).
* Only run folder 50 (flat) and 47 (full) still exist in the old checkout, so per-run convergence could not be checked from the logs.
  Substitute checks from the tables themselves:
  * Total_E = Excess_tension + Bending_tan + Bead to 1.4e-4 in both files.
  * Reference area: A - Excess/KA = 314.138 (314.1375 ... 314.1385, scatter from KA rounded to 4 digits) in BOTH files, so Total_E of the two starts is
    directly comparable (same A0, same Excess_tension definition).
  * Same-state scatter between flat and full start (345 pairs with |dz| < 0.2): robust sigma(dz) = 0.004 (95 % of |dz| < 0.022);
    sigma(dE) = 0.027 code units = 0.017 pi*kappa (95 % < 0.12, max 0.27). This is the run-to-run noise of Total_E (chaotic relaxation).
  * Envelope theorem: at a stationary state dTotal_E/dKI = Bead/KI exactly. Comparing finite differences of Total_E over KI steps (same branch)
    with the trapezoid of Bead/KI: rms residual 0.074 (flat) / 0.047 (full) code units per 0.1 step, median relative 2 %, 95th percentile 26-29 %
    (the latter where |dE| per step is small). Consistent with the 0.03-0.07 energy noise above; no gross unconverged run is visible, but single rows can be
    off by ~0.1-0.3 (e.g. sigma~ = 0.53, KI = 1.5: flat is 0.17 higher than full in the same state).
* Resolution: KI step 0.1 -> bracket width 0.4 in w~. Jumps are bracketed (lo, hi); I quote the bracket midpoint with +-0.2.

## 1. Definitions used

* z threshold 1.75 for "enveloped" (enveloped plateau z = 1.94-1.99 in both starts; partially wrapped states have z < ~1.6). S1_sim = midpoint of the bracket
  (last KI below the threshold, lowest KI from which z >= 1.75 for ALL larger KI) of the flat start; S2_sim the same for the full start (the lowest KI at which the full start stays enveloped).
  Columns `S*_lo/hi` give the bracket, `S*_inc*` the bracket of the largest increment of z (agrees with the threshold bracket to one step wherever the jump is sharp),
  `S*_sim_zad` the same with z_ad = -Bead/(2 pi KI) (adhesion-weighted covered solid angle) instead of z.
  Threshold sensitivity (1.5 / 1.75 / 1.9): none for sigma~ >= 0.67 (identical), up to one step for sigma~ <= 0.53.
* `S*_sharp` = largest increment of z between neighbouring KI >= 0.30 and the two brackets agree. At sigma~ <= 0.13 there is no jump (increments 0.30, smooth crossover; open markers),
  sigma~ = 0 is a smooth curve.
* E_sim: only points where the two starts are in different states (z_full - z_flat > 0.2) carry information. dE = E_flat - E_full changes sign from negative
  (partially wrapped lower) to positive (enveloped lower); the crossing is the linear interpolation, the quoted error combines linear-vs-quadratic interpolation and 0.05 energy noise / slope
  (it is 0.03-0.05 in w~, tiny compared to the 0.4 grid, because dE changes by ~1 per KI step). For sigma~ = 0.53, 0.67 there is only one distinct point (no sign change observed): E is only
  bracketed by the hysteresis window (method `bracket`, open markers, +-0.2). For sigma~ <= 0.4 the two starts coincide within noise everywhere (no hysteresis, E undefined).
* Theory: `spinodal_S1/S2`, `w_E` at the actual sigma~ (sigma~ = 0: all lines at w~ = 4).

## 2. The lines (w~; bracket +-0.2 for S1, S2)

| sigma~ | E_th | E_sim | S1_th | S1_sim | S2_th | S2_sim | window sim | window th | E in sim window |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 4.00 | - | 4.00 | 5.4 (soft) | 4.00 | 5.4 (soft) | 0 | 0 | - |
| 0.133 | 4.28 | - | 4.58 | 5.8 (soft) | 3.72 | 5.8 (soft) | 0 | 0.85 | - |
| 0.267 | 4.56 | - | 5.06 | 6.2 | 3.53 | 6.2 | 0 | 1.53 | - |
| 0.400 | 4.84 | - | 5.52 | 6.6 | 3.36 | 6.2 | 0.4 | 2.16 | - |
| 0.533 | 5.13 | 6.6 (bracket) | 5.96 | 7.0 | 3.20 | 6.6 | 0.4 | 2.76 | yes (bracket) |
| 0.667 | 5.41 | 7.0 (bracket) | 6.40 | 7.4 | 3.06 | 7.0 | 0.4 | 3.34 | yes (bracket) |
| 0.800 | 5.70 | 7.50 | 6.83 | 7.8 | 2.93 | 7.0 | 0.8 | 3.90 | yes |
| 0.933 | 6.00 | 7.87 | 7.25 | 8.2 | 2.80 | 7.0 | 1.2 | 4.45 | yes |
| 1.067 | 6.29 | 8.26 | 7.67 | 9.0 | 2.68 | 7.4 | 1.6 | 5.00 | yes |
| 1.200 | 6.58 | 8.66 | 8.09 | 9.4 | 2.56 | 7.4 | 2.0 | 5.53 | yes |
| 1.333 | 6.88 | 9.22 | 8.51 | 9.8 | 2.45 | 7.4 | 2.4 | 6.06 | yes |
| 1.467 | 7.18 | 9.50 | 8.92 | 10.2 | 2.34 | 7.4 | 2.8 | 6.58 | yes |
| 1.600 | 7.48 | 9.86 | 9.33 | 10.6 | 2.23 | 7.4 | 3.2 | 7.10 | yes |
| 1.733 | 7.78 | 10.22 | 9.74 | 11.4 | 2.13 | 7.4 | 4.0 | 7.61 | yes |
| 1.867 | 8.09 | 10.69 | 10.15 | 11.8 | 2.03 | 7.4 | 4.4 | 8.12 | yes |
| 2.000 | 8.39 | 11.09 | 10.56 | > 12.0 (no jump) | 1.93 | 7.4 | > 4.6 | 8.63 | yes (E < 12) |

(full numbers, brackets and flags in `data/baseline_lines.csv`.)

## 3. What the data show about the offsets

All numbers for sigma~ >= 0.5 and sharp/interpolated points only. "rms" is the rms residual of the description; the pure resolution error of S1, S2 is
0.115 (uniform +-0.2), of E about 0.05.

* **S2 (unwrapping of the enveloped state)** is the largest discrepancy: offset +2.7 (sigma~ = 0.27) growing to +5.5 (sigma~ = 2). S2_sim is nearly
  independent of sigma~: 6.2-7.0 for sigma~ 0.27-0.93 and a plateau 7.4 (bracket 7.2-7.6, i.e. KI = 1.8-1.9) for all sigma~ >= 1.07 up to 2.0, while the theory S2 falls from 3.5 to 1.9.
  Constant-offset model rms 0.63, proportional 1.34, "S2_sim independent of sigma~" 0.26 (7.23), linear in sigma~ (6.66 + 0.45 sigma~) 0.15.
  The enveloped state does unwrap, i.e. there is a barrier-free unwrapping only 4-5 units of w~ above the theoretical spinodal.
* **E (equal energy)**: offset +1.8 (sigma~ = 0.8) -> +2.7 (sigma~ = 2). Constant offset: rms 0.29, i.e. rejected by the data; proportional w_sim = 1.32 w_th: rms 0.052
  (ratio 1.31-1.34 at every one of the ten interpolated points); w_sim = w_th + 1.20 + 0.75 sigma~: rms 0.050 (equally good); "4 + c (w_th - 4)" rms 0.34 (rejected).
  E_sim lies inside the simulated hysteresis window wherever both are defined.
* **S1 (binding jump of the flat start)**: offset +1.0 ... +1.7, nearly constant against sigma~ (+1.05 at 0.27, +0.95 at 0.93, +1.3 at 1.07-1.6, +1.65 at 1.73-1.87; the step at
  sigma~ ~ 1.0 is a 0.4 grid step). Constant offset: rms 0.24, proportional w_sim = 1.155 w_th: rms 0.114, w_th + 0.66 + 0.49 sigma~: rms 0.113 (the last two are at resolution; the
  constant offset is marginally rejected). At sigma~ = 2 the flat start has not jumped at w~ = 12 (offset > +1.4).
* The proportionality constants differ between lines (E 1.32, S1 1.155; S2 not proportional), so the disagreement is NOT a single rescaling of w~ (e.g. w_eff = w/1.3).
* The hysteresis window [S2_sim, S1_sim] is much narrower than the theoretical one at small sigma~ (0 at sigma~ <= 0.27, 0.4-1.2 for sigma~ 0.4-0.93, vs 2.2-4.5)
  and comparable to or wider than it at large sigma~ (4.4 vs 8.1 at 1.87: narrower in absolute terms, because S2_sim does not drop). Window width sim/theory ~ 0.5 over sigma~ 1.5-2.
  Since S1_sim > S1_th and S2_sim >> S2_th, the first-order picture (hysteresis, jump, E inside the window) survives qualitatively; the lines are displaced to larger w~ and the lower edge much more.
* Low sigma~ (0 ... 0.13): no jump, flat and full start coincide: smooth crossover with z = 1.75 reached at w~ = 5.4 (sigma~ = 0, theory: continuous at w~ = 4), 5.8 (0.13). Lines from the
  theory (4.0, 4.3-4.6) are off by +1.2 to +1.4 there as well, but the data do not distinguish S1 from S2 from E at this level.
* At small w~ the shell makes the "free" state bound: z_cov = 0.40-0.42 and E~ = -0.46 at w~ = 2, -1.0 to -1.2 at w~ = 4 where the theory has z = 0, E~ = 0.

## 4. Energies besides the lines (facts)

**Enveloped state (z >= 1.9; 340 rows, both starts).** Identity E~ = -2(w~ - 4) + 4 sigma~ (F(z=2) = 0): Total_E of the code (Excess_tension + Bending_tan + Bead, with A0 = 314.138 the flat disc,
a = 1, kappa = KB/2) reproduces it to within a residual E~_sim - E~_th between -0.42 and +0.07 pi*kappa over w~ 6.8-12 and sigma~ 0-2 (mean -0.14 for w~ 7.6-8.8, i.e. 1-2 % of |E~|).
Linear fit residual = -0.50 + 0.033 w~ + 0.112 sigma~ (rms 0.025). So the energy scale (kappa = KB/2, sigma = KA, w = KI, E~ = E/(pi kappa)) is confirmed by the data, independent of the mesh, and the E-line offset cannot
come from the energy of the enveloped state: a residual of -0.14 would move E by only ~0.07 in w~ (slope 2).
The residual is a small difference of larger terms that partly cancel (per term, E~ units, sigma~ 0 ... 2, w~ 12 ... 7):
* adhesion: mean weight <w> = -Bead/(4 pi KI) = 0.972-0.986 (z_ad = 1.94-1.97 instead of 2): +0.2 ... +0.9 (= 2 w~ (1 - <w>)); an equivalent uniform offset of the membrane from the bead of
  0.014-0.04 a, compatible with face centroids lying slightly inside the sphere (not tested here),
* bending: E~_bend = 7.50 on average (6.7-7.8) instead of 8 (the discrete bending of the closed cap + neck is 3-17 % below the ideal 8 pi kappa): -0.2 ... -1.3,
* tension: 2 Excess/(pi KB) - 4 sigma~ = +0.03 ... +0.13 (excess area (A - A0)/pi = 4.03-4.29 instead of 4).

**Partially wrapped branch (flat start, z < 1.85, theoretical partial minimum existing, w~ < S1_th, sigma~ >= 0.5).** The simulated energy is far BELOW the theoretical partial minimum
(or free state), mean E~_sim - E~_th,partial per w~ bin: w~ 2-3: -0.56; 3-4: -0.86; 4-5: -1.23; 5-6: -1.64; 6-7: -2.04; 7-8: -2.42; 8-9: -2.73; 9-10: -3.04 (bin spread 0.1-0.2), i.e. roughly -0.3 w~ and
nearly independent of sigma~. The simulated coverage is z_cov = 0.41 (w~ 2-3) ... 1.08 (9-10) versus the theoretical z_th = 0 ... 0.76 and z_ad = 0.27 ... 0.98. Hence the
simulated partially wrapped state is more strongly bound than in the theory (the soft shell attracts the membrane also before contact), while the enveloped state is (almost) at its theoretical energy
(an extra -0.3 w~ on the partial branch is the same order as the whole E-line distance E~_env - E~_partial at the theoretical E line). The same sign (partial state favoured) is what pushes E to larger w~.
I do not claim the quantitative link (-0.3 w~ shift => 1.32 factor): a rigid downward shift of the partial branch with the theoretical shape cannot reproduce the observed crossing (tried, it never crosses
before the theoretical S1 where the branch ends), so the branch is modified in shape as well.

## 5. Where is the disagreement? (answer, with what is data and what interpretation)

Data: it is in all three lines, with different character.
1. S2 is wrong by far the most (+2.7 ... +5.5) and is almost independent of sigma~ (a plateau at KI = 1.8-1.9). The sigma~ dependence of the theory (neck instability set by tension) is absent.
2. E is displaced by +1.8 ... +2.7, proportional to w~_th (ratio 1.32 +- 0.01) or equivalently linear in sigma~ (1.2 + 0.75 sigma~); a constant offset is excluded.
3. S1 is displaced the least (+1.0 ... +1.7), about proportional (1.155) or constant + weak growth; consistent with a roughly constant offset at fixed grid resolution only marginally.
4. Not a constant offset in w~ for E or S2; E and S1 have a similar relative offset pattern (proportional), S2 does not.

Interpretation (not shown by these data): (a) the partially wrapped branch is energetically too deep (-0.3 w~) because adhesion acts at a distance through the 0.25 a shell, which moves E and
S1 upwards; (b) the sigma~-independent S2 plateau at KI ~ 1.85 suggests a local, mesh/shell-controlled mechanism at the neck (the neck radius goes to the edge length scale) rather than
the tension-controlled continuum instability. Deciding between these needs the planned runs with different shell widths/steepness and mesh sizes (those would show whether the S2 plateau and the
E ratio change). The only geometry available from the old data is the two surviving run folders; no profile (membrane shape, contact-line curvature) analysis was done here.

## 6. Limitations

* Resolution 0.4 in w~ (grid), noise 0.03-0.07 in Total_E code units; convergence per run not verifiable (only 2 folders left).
* sigma~ is a discrete 16-point grid; the lines are brackets, E needs two distinct states.
* The z threshold (1.75) is arbitrary but results do not depend on it for sigma~ >= 0.67.
* CoverageUnion is a post-processing measure of the soft-shell coverage (non-zero z = 0.4 for a non-wrapped bead): z_sim cannot be mapped onto the geometric wrapping angle of the theory.
  z_ad (energy based) is a second measure; both give the same S1/S2 for sigma~ >= 0.67.
