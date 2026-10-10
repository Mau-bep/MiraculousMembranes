# Planar bead wrapping: simulation vs Deserno, and the role of the adhesion shell

Branch `deserno-shell-study`. Question: why do the simulated phase lines of a bead on a planar membrane
(`Adhesion` + `Bending_tan` + `Excess_tension`) lie at larger w~ than Deserno's (arXiv:cond-mat/0303656, Figs. 2 and 5),
can a steeper / narrower shell get closer, and is there a formal way to defend the disagreement?

Variables: kappa = KB/2, sigma = KA, w = KI, a = BeadRadius; w~ = 4 KI a^2/KB, sigma~ = 2 KA a^2/KB; energies in units pi kappa
(Total_E * 2/(pi KB)). See `CONTEXT.md` for the full mapping and the simulated adhesion model.

## Short answer

1. **The disagreement is the finite range of the adhesion shell, not the mesh and not a calibration error.**
   The simulation implements a soft adhesion potential of range delta = s a (s = 0.25 so far); Deserno's theory is the zero-range,
   hard-contact limit. The continuum limit of the simulated energy (`theory/softshell_direct.py`, axisymmetric, no mesh) removes
   81 / 91 / 85 % of the offsets of S1 / E / S2 of the old runs and reproduces the simulated partially wrapped branch to
   0.05 pi kappa rms (236 points).
2. **Fine-mesh simulations sit on that continuum model** (edges ~0.1 near the bead): flat-start partially wrapped states agree
   to ~0.006 pi kappa for the shells p = 1, 2, 4, 8 (`data/campaign_s025.csv`, `figures/campaign_energies.png`); the simulated E crossing
   (5.86, 6.05, 6.30 for p = 8, 4, 2) and S2 brackets agree with the continuum values (5.84, 6.03, 6.27; 5.07, 5.45, 5.89).
3. **Formal argument** (`theory/line_tension_analysis.md`): a finite-range shell leaves the contact condition
   w~ = (1 - psi0dot)^2 exactly as in Deserno (first integral of the contact boundary layer, verified numerically for p = 1..8), so w~
   needs no recalibration. It adds a NEGATIVE effective contact-line tension tau = -T_p w sqrt(delta/Dc), Dc = sqrt(w~)/a
   (T_p = 1.198, 1.029, 0.875, 0.740 for p = 1, 2, 4, 8): lift-off membrane keeps gaining adhesion over a length
   ell = a sqrt(s) w~^(-1/4). In Deserno's units E_soft(z) = E_Deserno(z) - T_p sqrt(s) w~_D(z)^(3/4) sin(alpha) (no fit): it lowers
   partially wrapped states, vanishes at z = 0 and 2, so E and S1 move up, the unwrapping barrier is reduced (0.67 -> 0.13 at
   sigma~ = 0.53) and the E line lands within 0.12 of the old data. First order only (ell/a ~ 0.33 at s = 0.25); fails near the neck
   (S2 plateau needs the exact solver).
4. **Steeper shell = narrower shell**: all p collapse onto s_eff = s (I_p/I_1)^2, factors 0.74, 0.54, 0.39 for p = 2, 4, 8
   (`figures/softshell_collapse_seff.png`). p = 8 at s = 0.25 halves the E and S1 offsets (E 5.85 vs 6.56, Deserno 5.13 at sigma~ = 0.53),
   but the shifts only decay like sqrt(s): matching Deserno to 0.1 needs s_eff ~ 0.003, out of reach for a direct simulation. The
   practical route is the series s = 0.25 ... 0.05 with the mesh scaled to ell and a quadratic extrapolation in sqrt(s).

## Where things are

| what | file |
|---|---|
| shared briefing (mapping, model, rules) | `CONTEXT.md` |
| continuum soft-shell solver, derivation, tests | `theory/softshell_direct.py`, `.md`, `test_softshell_direct.py` (37/37 with --slow) |
| line-tension boundary layer, first-order model, discretisation, extrapolation protocol | `theory/line_tension_analysis.md`, `boundary_layer.py` |
| literature (what supports / does not support the argument) | `literature/notes.md` |
| baseline extraction of the old simulated lines | `sim/baseline.md`, `sim/extract_lines.py`, `data/baseline_*.csv` |
| shell steepness campaign (42 runs) | `sim/run_condition.py` (driver), `sim/compare_campaign.py`, `sim/plot_campaign.py`, `data/campaign_*` |
| lines of the continuum model vs (s, p, sigma~) | `data/softshell_lines.csv`, `figures/softshell_*.png` |
| C++ option | `"shell_power": p` in the bead (commit d1b2b7ba), test `src/Test_shell_power.cpp` |
| cluster scan for the extrapolation | `Scripts/Mem_planar_shell_series_up.sh` |

## Protocol notes (local runs, see also `sim/S0_PROGRESS.md`)

* `--rescale 1.0` always (template default 3.0 makes a radius-30 disk), explicit `--size_min/--size_max`.
* The coarse old protocol (edges 0.05-0.5, switch at 30000) is chaotic: energies move 0.05-0.15 between replicas. With `refine_angle 0.15`,
  `size_max 1.0`, `BFGS_saved_states 60`, plateau switch / 60000 the results replicate to 2e-3 and follow the continuum model;
  a run takes ~35-70 min (p = 1..8, s = 0.25), longer for narrower shells (~5-10x more faces for s = 0.05).
* Edge length near the bead should be h <~ 0.5 ell with ell = sqrt(s_eff) w~^(-1/4) a.
* The flat start jumps (S1) up to ~0.4 later than the continuum spinodal (slow escape near it); full starts below S2 stop in states
  0.3-0.7 pi kappa above the continuum minimum (not understood: pinning of the contact line on the mesh or incomplete unwrapping).
  Use the energy crossing and the enveloped / partial branch energies, not only the jumps.

## Open

* Narrower p = 8 shells (s = 0.15, 0.10) with scaled meshes, run locally: see `sim/` and the end of this file once they have finished.
* Extrapolation to s -> 0 on the cluster (series script), more sigma~, Fig. 5 (vesicle) needs its own calibration of w.
* Two recent Proc. R. Soc. A papers on contact-line bending energy / Tabor parameter of membrane adhesion could not be opened
  (HTTP 403); not cited.
