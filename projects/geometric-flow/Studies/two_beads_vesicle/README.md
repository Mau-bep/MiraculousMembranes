# Two wrapped beads on a vesicle: the pair attraction of van der Wel et al. (Sci. Rep. 6, 32825), Fig. 3e / 4a

Status 2026-10-11: the template, the runner and the analysis work and are tested; the C++ is an opt-in branch. **The attraction of the paper has not been reproduced.**
The pilot (section 5: beads outside, 30-90 % wrapped, noisy) gave repulsion; the second round (section 6: the paper's own geometry as far as it can be read from the paper, beads inside,
98 % wrapped, 5 % area, plus beads outside at 98 % and the flat limit) gives flat or repulsive curves, replicated to 0.015 kappa. Section 7 lists what is left.

## 1. What the paper measures and simulates

Two colloids (D_p = 0.98 um) wrapped by a floppy GUV attract each other. u(s), s the geodesic distance between the points where the two
particles connect to the membrane (the arc of the vesicle between the necks, in the code `s = R_v phi`), has a well of -3.3 kBT = -0.16 kappa
(kappa = 20.9 kBT) at s ~ 1.5 D_p, goes to 0 at ~ 2.5 D_p and is flat beyond (Fig. 3e experiment, Fig. 4a simulation; digitized in
`paper_fig4a_digitized.csv`). Non-wrapped particles do not interact. In the experiment the wrapped particles sit outside the GUV and
"protrude into the interior" (the membrane invaginates around them).

The simulation (Saric & Cacciuto, Methods of the paper): closed triangulated vesicle (5882 vertices, hard core l, D_v = 50 l = 6.25 D_p),
**volume V_t = V_0 and area A_t = 1.05 A_sphere fixed** (K_A = 1e3 kBT/l^2, K_V = 4e3 kBT/(D_p l^2): both within 0.01 %), bending only
(sqrt3 kappa sum(1 - n_i.n_j), kappa = 20.9 kBT), particles D_p = 8 l with an adhesive patch: vertices within the angle theta_Wr of the normal of the active area
feel -eps (l_m/r)^6, l_m = (l + D_p)/2 the contact distance (a hard core), cut at 1.2 l; the active part is "about 90 % of the particles' total area" in the text, and the
caption of Supplementary Fig. S5 shows an inactive cap at the top of the particle and says "in the given example, theta_Wr = 11 pi/12" (165 degrees: the inactive cone is 30 degrees wide, 1.7 % of
the area, the active part 98.3 %; 11/12 = 91.7 % is the fraction of the polar angle, maybe where "about 90 %" comes from; 144 degrees = 4 pi/5 would be 90.5 % of the area). **Which angle was used for Fig. 4 is not stated**:
the tables below keep both, 90 % (143.13 degrees) and 98 % (165 degrees). Monte Carlo annealing of vertices, particles and edge flips. u(s) is the change of the curvature energy;
the adhesion energy and the wrapped area do not change with s (Supplementary Fig. S5b: u_Curv ~ u_T, u_Ad ~ constant); the zero is the energy with the particles at
opposite poles. The paper says the closed fixed-area vesicle has **no long-range repulsion**, in contrast to the flat membrane (Reynwar et al.: a well of
-3 kBT with a barrier at longer range).

Two things in the paper's simulation that the first version of this study did not copy (found while re-reading it, section 6; they do not matter for the experiment, whose beads are outside):
(i) **the particles look like they are inside the vesicle**: in the cross-sections of Fig. 4b,c and in the snapshots of Supplementary Fig. S5b the vesicle is cut open and the particles sit in cups of its inner wall
(outward bulges); the template had them outside, invaginating, like the experiment ("protruding into the interior of the vesicle");
(ii) **the area budget**: in my units (membrane on the bead surface, R_v/a = 6.25) the floppy vesicle with the beads *inside* needs A = 1.0489 A_sphere for both 90 % and 98 % wrapped beads, with the area constraint 1.0498:
the paper's 5 % fits that geometry. Beads *outside* need 1.0605 A_sphere at 98 % and 1.0667 at 90 %, so the 5 % constraint makes them taut (sigma~ = 1.9 at 98 %). That is a plausibility argument for "inside" in the paper's simulation, not a proof:
in the paper the membrane sits 0.5 l above the particle surface (l_m = 4.5 l, 1.125 a), which changes the budget.

## 2. Why the old runs (Wrapping_two_cov / _fixed, TwoBeadsCov, HighBend) could not show it

(bead radius a = 1, kappa = KB/2, 1 kappa/a^2 = 358 nN/m for kappa = 20.9 kBT and a = 0.49 um)

* **tension.** `Surface_tension [1.0]` with KB = 0.121 is sigma~ = sigma a^2/kappa = 16.5 (5.9 uN/m); HighBend (KB = 0.85) is 2.34
  (0.84 uN/m). The floppy vesicles of the experiment have sigma < 10 nN/m, sigma~ < 0.03; the tense ones (> 1 uN/m, sigma~ > 2.8) are
  the blue points of Fig. 4a and attract much less. At sigma~ = 2-16 the membrane is flat beyond a fraction of a bead radius
  (lambda = sqrt(kappa/sigma) = 0.25-0.65 a): your E_bend(d) contours are flat for d > 2.75, as they should be there.
* **size of the signal.** u_min = 0.16 kappa on a bending energy of 67-70 kappa (vesicle 8 pi kappa + two ~90 % caps), 0.2 %. Re-meshing
  a converged state and relaxing it again at the same separation changes E_bend by 0.007-0.04 kappa.
* **soft contact.** The Coverage energy counts faces in a +-25 % shell with a cosine weight and has no hard core: weighted coverage 0.60
  is 0.66 of the solid angle, 0.90 is 0.93, the membrane in the cap sits at r = 0.988 +- 0.010 (inside the bead), and the lift-off
  region is part of the constraint, so the "wrapped area" depends on the neck shape.
* **area budget.** With the paper's 5 % excess two beads can only be wrapped to ~ 65 %: wrapping used 33 a^2 of area at 88 % and
  20.6 a^2 at 60 % (R_v = 6.25, 5 % is 24.5 a^2). Forcing 90 % in the paper's ensemble makes the vesicle taut (`Area_constraint`
  K_A = 2e5: sigma~ = 1.2-1.4). That is the "wrapping around 50-60 %" you suspected, and the area, not the adhesion, can be what sets it.
* **no stable partial wrapping at low tension.** Deserno's theory (`Scripts/deserno_theory.py`) at sigma~ = 0.03-0.1: only unbound or fully
  enveloped, w~_E = 4.07-4.2 (a jump); partially wrapped states appear for sigma~ >~ 0.3. So the degree of wrapping has to be imposed
  (the paper does it with theta_Wr), "wrap by adhesion only" cannot give 50-60 % in the floppy regime.

## 3. What the template does (`Templates/Wrapping_two_vesicle.txt`)

Closed vesicle of radius R_v = 6.25 a (the paper's D_v/D_p), two beads of radius a, **kappa = 1 energy unit (KB = 2)**: u comes out in kappa,
u [kBT] = 20.9 u, length [um] = 0.49 length. Energies: `Bending_tan`, `Volume_constraint` (V = 4/3 pi R_v^3, K_V = 2e5) and one of

| `--area` | what | use |
|---|---|---|
| `tension` (default) | `Surface_tension` sigma = 0.02 kappa/a^2, area free | the experiment's floppy vesicle (sigma~ = 0.02 = 7 nN/m) |
| `constraint` | `Area_constraint`, A_t = (1 + 0.05) 4 pi R_v^2, K_A = 2e5 | the paper's ensemble |
| `excess` | `Excess_tension`: free area up to A_t, then a tension | soft version of the budget |
| `none` | nothing | |

**What "floppy", "taut" and "tense" mean in the runs** (units: bead radius a = 1, kappa = KB/2 = 1; the exact `Energies` of the generated config):
* **floppy** = `--area tension --sigma 0.02`: `Bending_tan [KB = 2.0, H0 = 0]`, `Volume_constraint [K_V = 2e5, V_bar = 4 pi/3 R_v^3 = 1022.65]`, `Surface_tension [0.02]`. The area is **free**, it only costs sigma = 0.02 kappa/a^2
  per unit area (sigma~ = sigma a^2/kappa = 0.02, i.e. 7 nN/m for kappa = 20.9 kBT and a = 0.49 um; the experiment's floppy GUVs have sigma < 10 nN/m). The vesicle uses the area it needs (1.049 to 1.067 A_sphere, depending on the
  wrapping), so it stays nearly spherical (reduced volume 0.93-0.95). A real GUV has a fixed area, usually with more excess than that, and a tension that is only the Lagrange multiplier of that constraint; the two agree when the
  excess area is not binding. It is not a flaccid, strongly fluctuating vesicle (reduced volume 0.7-0.9).
* **taut** = `--area constraint --excess 0.05`: the tension term is replaced by `Area_constraint [K_A = 2e5, A_bar = 1.05 * 4 pi R_v^2 = 515.42, nu = -1, dA = 1e4]`, the paper's fixed area; the tension is then the multiplier
  (the analysis prints it as sigma~: 1.9 for beads outside at 98 %, 0.03-0.06 for beads inside).
* **tense** = `--sigma 3.0 --W 20`: sigma~ = 3 (1.1 uN/m; the experiment's tense GUVs have > 1 uN/m), W = 20 so that the patch stays wrapped against the tension.
* The bead: `Adhesion` W = 6 kappa/a^2 (w~ = 2 W a^2/kappa = 12; the all-or-nothing threshold of a single bead at low tension is w~_E ~ 4.1-4.2), radial shell `shell_width` 0.1 (faces within +-10 % of the bead radius count), patch of half-angle
  theta_Wr with a 5 degree ramp inside its edge, frozen at the chain's positions, the patch axis fixed in space (pointing to the vesicle centre / down on the disk, refreshed at every separation).

The wrapped area is set by an **adhesive patch with a sharp edge** (`--wrap patch`, the paper's theta_Wr, `--patch_angle` in degrees:
90 = 50 % of the bead area, 101.5 = 60 %, 143.13 = 90 %, 165 = 98 %): `Adhesion` (W = 6 kappa/a^2, w~ = 2 W a^2/kappa = 12 >> 4) restricted to the
faces within theta_Wr of an axis, with a 5 degree cosine ramp. The membrane in the patch sits on the sphere (r = 1.003 +- 0.006, against
0.988 +- 0.010 with Coverage), its contact line is pinned at the edge and the lift-off beyond the edge is not counted, so the energy no longer
depends on how the neck is shaped. `--wrap coverage` and `--wrap adhesion` give the old Coverage constraint and the plain adhesion.
The beads are frozen, the vesicle finds its position; each patch axis points to the centre of the vesicle and is refreshed at every separation
(`--free_axis` lets it turn, gain 0.03 kappa at s = 2.7, but it ran into a closed neck in one test).

**The C++ change** (`bead_patch_adhesion.patch`, branch `bead-patch-adhesion`, worktree `../MM-patch`, tested by `Test_patch_adhesion`: forces vs finite
differences 3e-9, area fraction of the patch vs quadrature, free axis stationary, `planar_adhesion_bfgs` benchmark byte-identical to the old binary) adds the
bead keys `patch_angle`, `patch_width`, `patch_axis`, `patch_free` to the `Adhesion` interaction; without them nothing changes. An old build **ignores the keys
silently** and runs a full-sphere adhesion, so `Create_subjob_two_vesicle.py` refuses to write the job unless `../build/bin/main_cluster` contains
"patch_angle". Apply with `git apply Studies/two_beads_vesicle/bead_patch_adhesion.patch` and rebuild.

**Beads inside the vesicle** (`--inside`, bead key `"outside": -1`): the Adhesion (and the Coverage energy) select the faces whose outward normal points to the bead
(the external-bead convention), so a bead on the inner side of the closed membrane would see nothing. `Adhesion::side` = -1 (set from `"outside": -1` of the bead) selects the
faces whose normal points away from the bead instead; `Test_patch_adhesion` checks that the outward-oriented sphere with side -1 gives bit for bit the energy and forces of the
reversed sphere with side +1 and that the wrong side selects nothing; the default path is unchanged (benchmark byte-identical). The Coverage energy does not support it (SimConfig refuses the
combination), so an inside chain has no pre-wrap: the cold start from the tangent contact takes 30000-60000 steps. The template then puts the beads tangent from the inside
(centres R_v - a from the vesicle centre) with the patch axes pointing outwards; the runner points them away from the volume centre at every link.

**Protocol** (`Run_two_vesicle_chain.py`, all in one sbatch job):
1. pre-wrap (beads outside): the stiff `Coverage` energy (K = 2e4, target 0.5) wraps the beads from the tangent contact in 8000 steps (forcing 0.8 crumples the membrane);
   not available for `--inside` (the Coverage energy knows beads outside only) and not used on the flat disk;
2. cold start: the patch ensemble from that mesh, to convergence (~ 8000-30000 steps; from the tangent contact alone the patch needs > 30000 steps just to wrap).
   **`--inside` and `--planar` use free beads in the cold start** (`--free_cold`): the beads move with the membrane that wraps them to the depth they like (inside: 2.3 a
   outwards, the buds end up beyond the main sphere; planar: they sink 1.5 a) and are frozen there from the first link on. With frozen beads the inside cold start needed
   more than 60000 steps (the whole vesicle has to slide away from the beads and retract the lobe left between them), with free beads 15000-30000;
3. the beads are then rotated about the volume centre of the vesicle by +-(phi_new - phi_old)/2 (along the arc, so they stay at their depth), the previous mesh is the
   initial mesh, the volume and area targets are absolute, the beads stay frozen. A link takes 2-5 minutes instead of an hour and the neighbouring
   meshes are related (re-meshing and re-relaxing a converged state changes E_bend by 0.007-0.04 kappa). Steps larger than `--max_step` (0.5 D_p) go
   through waypoints (a jump of 2.4 D_p unwrapped both beads); a link where the integrator breaks is retried in two half steps; "Sim broke", a run of a few steps or an
   energy > 1e6 count as failures (main_cluster returns 0 and writes a Final_state.obj anyway);
4. the default walks in from the far reference (5.8 D_p) to 1.45 and, unless `--no_return`, back out, so every separation appears twice: the difference is the hysteresis
   (0.03-0.07 kappa floppy 90 %, up to 0.3 at 60 %). The pilot of section 5 went out to 5.9 first and then in, with 20 links.

**s** is the arc of a sphere, s = R_v phi, phi the angle between the two bead centres seen from the volume centre of the vesicle (your choice, no geodesics; the paper's s is
the geodesic between the points where the particles connect to the membrane, the same thing). On the flat disk s = d/D_p, the distance between the bead centres. The 3D distance of
the centres d/D_p is the `s_chord` column (`--xaxis chord`): it is larger than s by R_c/R_v, 1.35 for the inside buds, 0.7 for the deep outside pockets.
Beware that the beads have no excluded volume: with deep pockets (floppy, 90 %) the bead centres sit 3.6 a from the vesicle centre instead of 5.3 and the
beads overlap (d < 2a) below s = 1.83 D_p; `Analyze_two_vesicle.py` marks those points.

**Analysis** (`Analyze_two_vesicle.py`): u(s) = E(s) - E(largest s) for the bending energy and for the **total** energy (the bending energy trades against the adhesion
energy when the patch is tilted: 4 degrees move 0.06 kappa between them and leave the total unchanged to 0.003, so the total is the robust curve), table, csv, plot in kappa and
kBT with the paper's digitized Fig. 4a and a +-0.03 kappa restart-noise band. `Slice_two_vesicle.py` draws cross-sections: look at them.

## 4. How to run

```bash
cd projects/geometric-flow/Scripts
# one curve, here (no sbatch): writes the configs and the chain file, then runs the chain (2-4 h). The paper's ensemble as far as it can be read:
# beads inside, theta_Wr = 165 deg (98 %), fixed volume and area A_t = 1.05 A_sphere
python3 Create_subjob_two_vesicle.py --batch_tag Wrapping_two_vesicle/test --unique_tag test_in98 --inside --area constraint --excess 0.05 --patch_angle 165
python3 Run_two_vesicle_chain.py ../Config_files/test_in98_chain.json
python3 Analyze_two_vesicle.py --ref mean_far --s_far 4.3 --curves total --ylim -0.5 1.0 --out ../Studies/two_beads_vesicle/test in98=../Results/Wrapping_two_vesicle/test
python3 Slice_two_vesicle.py ../Results/Wrapping_two_vesicle/test/link12_s2.70/1 --out ../Studies/two_beads_vesicle/test_slice.png
# the other geometries: beads outside (the experiment: drop --inside), the flat disk (--planar --d_cold 14 --R_disk 24 --s_list 7.0,6.0,5.0,...,1.1), floppy (--sigma 0.02, the default)
# on the cluster: edit and run  tcsh Mem_two_vesicle_up.sh   (inside/outside/planar x floppy/area budget x replicas)
```
`python3 Create_subjob_two_vesicle.py --help` lists the knobs (`--inside`, `--planar`, `--free_cold`, `--s_list`, `--patch_angle`, `--free_axis`, `--area`, `--sigma`, `--excess`,
`--Rv`, `--W`, `--replica`, `--wrap`, mesh settings). An inward chain can be continued or reversed: copy its chain file, set `"start_run"` to its last link folder and `"s_list"` to the
new values (outward leg for the hysteresis, or an extension to smaller s), and give both folders to the analysis as one series (`label=dir1,dir2`). Mesh: edges 0.2 near a bead (`refine_angle` 0.3), 0.6 far (`size_max`), 4000-6000 vertices; finer meshes slow everything down more than they help (the
restart noise is already 0.01-0.04 kappa).

## 5. Results of the pilot (R_v = 6.25 a, kappa = 1, KB = 2, W = 6, mesh as in the template)

Chains of 20 links (out to s ~ 6 D_p and back in to ~ 1.4), data in `data/*.csv`, figure `figures/u_of_s_pilot.png`, cross-sections `figures/slices_s2.7.png`
(`make_figures.py`, `Scripts/Slice_two_vesicle.py`). u = change of the total energy in kappa, zero = mean of the relaxed links beyond 4.3 D_p; the
cold start is not a point (it stops 0.3-0.5 kappa above what the next link finds). Frozen patch axes.

**u_tot(s) [kappa], bands where the repeated points disagree by more than 0.15 are starred** (the paper: -0.15 at 1.5, -0.09 at 1.7, -0.05 at 1.9, 0 beyond 2.3):

| ensemble | far E (std) | 3.9 | 3.4 | 3.0 | 2.7 | 2.4 | 2.2 | 2.0 | 1.8 | 1.65 | 1.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| floppy (sigma~ 0.02), 90 % | -55.50 (0.02) | +0.09 | +0.12 | +0.17 | +0.22 | +0.24 | +0.25 | +0.33* | overlap | overlap | overlap |
| area 5 %, 90 % (taut, sigma~ 1.1) | -62.19 (0.02) | 0.00 | -0.01 | -0.01 | +0.08 | +0.12 | +0.08 | +0.13 | +0.18 | +0.31 | +0.43* |
| floppy 60 % | -24.76 (0.13) | +0.39* | +0.38 | +0.48* | +0.45 | +0.66 | +1.16 | +1.95* | +3.58 | +6.9* | +8.4 |
| area 5 %, 60 % | -34.98 (0.08) | +0.22* | +0.54 | +0.77 | +0.64 | +1.05 | +1.61 | +3.3* | +4.9* | +8.7 | +8.0 |
| floppy 50 % | -15.54 | +0.09* | +0.04 | +0.23 | +0.43 | +0.62* | +0.63* | +0.93 | +1.56 | +3.5 | +5.7* |
| floppy 40 % / 30 % | noisy +-0.5 | | | | | | | +0.85 / +0.68 | +1.44 / +0.70 | +2.0 / +0.78 | +2.9* / +1.1* |

**No ensemble shows the well: in every one the energy rises as the beads approach below s ~ 3 D_p.** The plateau beyond is flat to 0.02 kappa for the taut vesicle
(screening, as in the paper) and for the floppy 90 % one beyond 4.3 D_p; the floppy 90 % chain has a repulsive step of 0.1-0.25 kappa between 4 and 2.4 that the
taut one does not have.

What the states look like (`figures/slices_s2.7.png`): with 90 % wrapped, bead centres sit deep, 3.6 a from the vesicle centre (floppy, necks ~ 2 a long, the beads overlap
below s ~ 1.8 D_p) and 4.5 a (area 5 %); with 30 % the beads sit on the surface. Between two pockets the membrane forms a thin lobe of vesicle interior whose tip tightens as the beads
approach: in the floppy 90 % chain the bending energy of the membrane within 3 a of the beads (outside the caps) rises by 0.46 kappa from s = 6 to 2.7, the envelope gives 0.26 back,
the caps stay constant (cotan decomposition, E_bend there 66.2 against the code's 66.8, the two discretizations agree within 1 %).
For 50 % wrapping, in the last 0.15 a of gap (s = 1.3-1.4) the *bending* energy does turn down (u_bend = -1.2 and -2.3: cooperative wrapping) but the patches then lose 12 % of their
coverage (0.49 to 0.43) and u_tot is +6.

**Noise and checks.** Re-meshing and relaxing a converged state at the same separation: E_bend +0.007 and +0.04. Plateau scatter 0.02 (90 %), 0.1 (60 %), 0.3-0.5 (30-50 %, soft modes).
Going out and coming back at the same s: 0.06 (floppy 90 %), up to 0.2 in the default-pipeline check at s = 3-3.4, which is as large as the signal: replicas (`--replica`) and averaging are needed to resolve 0.15 kappa
in the floppy ensembles. Tilting the frozen patch axes by 4 degrees moves 0.06 kappa between E_bend and E_adhesion and changes the total by 0.003; relaxing from tilts of +-8 degrees ended 0.04-0.13 lower (other basins);
the free axis gained 0.026 at s = 2.7. The default pipeline (pre-wrap, cold start, links with waypoints) ran without failure in 26 minutes for 6 links (`figures`, `data/pipeline_check_floppy_90.csv`).

## 6. Second round: the paper's own geometry, budded beads, the flat limit

Re-reading the paper after the pilot (section 1) changed three things, all implemented (`--inside`, `--patch_angle 165`, `--planar`, `--free_cold`):
beads **inside** the vesicle (my reading of the paper's *simulation*; the experiment has them outside, so the rows 'outside' and 'flat disk' are the ones that correspond to it), theta_Wr = 165 degrees (98 % of the bead, narrow necks; the example of Fig. S5, the text says about 90 %) and a flat disk for the large-vesicle limit of the experiment. Same units, same chain protocol
(free-bead cold start, then warm-started links). Data: `data/*.csv` (`out_*` beads outside, `in_*` inside, `planar_*` flat disk),
figure `figures/u_of_s_round2.png`, cross-sections `figures/slices_inside98.png`, `figures/slices_planar_outside98.png`.

**u_tot(s) [kappa]** (zero = mean of the relaxed links beyond 4.3 D_p, beyond 5.5 on the disk; in brackets the spread of the neighbouring runs; the paper: -0.15 at 1.5, -0.09 at 1.7, 0 beyond 2.3).
s is the arc R_v phi of the vesicle (the distance d/D_p of the bead centres on the disk); on a vesicle with beads inside the centres are 1.35 times farther apart than s, so the paper's
centre-to-centre distance 1.5 D_p is s = 1.1 here (`s_chord` column):

| ensemble | sigma~ | 3.0 | 2.4 | 2.0 | 1.65 | 1.45 | 1.2 | 1.0 |
|---|---|---|---|---|---|---|---|---|
| inside, 98 %, area 5 % (my reading of the paper's simulation) | 0.034 | +0.01 | +0.02 | +0.02 | +0.03 | +0.05 | +0.06 | +0.10 |
| inside, 98 %, floppy | 0.02 | +0.01 | +0.02 | +0.03 | +0.06 | +0.06 | +0.07 | +0.13 |
| inside, 98 %, floppy, replica 1 | 0.02 | +0.01 | +0.02 | +0.02 | +0.04 | +0.05 | +0.07 | +0.11 (0.06) |
| inside, 90 %, floppy | 0.02 | +0.04 | +0.05 | +0.09 | +0.11 (0.14) | +0.27 | +0.41 (0.37) | +0.61 (0.39) |
| inside, 98 %, tense (sigma~ = 3, W = 20) | 3 | +0.02 | -0.01 | +0.01 | +0.04 (0.05) | +0.03 (0.06) | -0.03 (0.13) | +0.11 |
| outside, 98 %, floppy | 0.02 | +0.01 | +0.04 | +0.05 | +0.05 |  |  |  |
| outside, 98 %, area 5 % (taut) | 1.9 | -0.03 | +0.02 | +0.01 | +0.06 | +0.08 |  |  |
| outside, 90 %, floppy (pilot) | 0.02 | +0.16 | +0.24 | +0.31 (0.39) |  |  |  |  |
| outside, 90 %, area 5 % (pilot, taut) | 1.1 | +0.00 | +0.12 | +0.12 (0.07) | +0.30 | +0.55 (0.24) |  |  |
| outside, free contact line (W = 3.5), area 5 % | 0.82 | -0.00 | +0.00 | +0.02 | +0.05 | +0.08 |  |  |
| flat disk, 98 % | 0.02 | -0.01 | +0.03 | +0.02 | +0.04 | +0.06 | +0.09 |  |
| flat disk, 90 %, inward | 0.02 | +0.05 | +0.08 | +0.01 | +0.29 | +0.42 (0.06) | +0.75 (0.42) |  |
| flat disk, 90 %, outward | 0.02 | +0.04 | -0.01 | +0.19 | +0.30 | +0.34 | +0.74 (0.13) |  |
| flat disk, 75 % | 0.02 | +0.13 | +0.28 | +0.50 | +0.94 | +1.44 | +3.59 (2.98) |  |
| flat disk, 60 % | 0.02 | +0.14 | +0.28 | +0.66 | +1.05 | +1.46 | +4.12 (3.26) |  |

* **My reading of the paper's simulation** (inside, 98 % or 90 %, area 5 %; `--inside --area constraint --excess 0.05 --patch_angle 165`): flat to +0.03 down to s = 1.9, then +0.05 at 1.45, +0.06 at 1.2, +0.10 at 1.0, +0.14 at 0.9. **No well.**
  The two replicas of the floppy version agree within 0.015 kappa along the whole curve, so the noise is not the issue (it was in the pilot, where the beads were only partially wrapped).
* floppy (sigma~ = 0.02) and area 5 % give the same curve: with the beads inside and 98 % wrapped the floppy vesicle settles at A = 1.0489 A_sphere (514.9 a^2 for R_v = 6.25 a) and with
  the constraint at 1.0498, and the same for 90 % wrapped beads inside (1.0488): **the paper's 5 % fits two wrapped beads inside the vesicle** (the buds take volume from the main sphere, whose area is 0.55 % smaller); outside the same beads need 1.0605 A_sphere
  (the 5 % constraint then makes the vesicle taut, sigma~ = 1.9) and at 90 % wrapped 1.0667. A plausibility argument for 'inside' in the paper's simulation (section 1), and the reason my pilot (beads outside, 90 %, 5 %) was a taut vesicle.
* beads outside, 98 %: the same, flat to +0.02 for s > 2 and +0.05-0.09 at 1.45-1.65 (floppy and taut); the beads touch below s = 1.55.
* 90 % wrapped (open necks, radius 0.6 a): stronger repulsion, +0.27 at 1.45 and +0.6 at 1.0 for beads inside, +0.1-0.3 at 2.0 for beads outside (pilot, noisy).
* **flat disk, 90 %** (R = 24 a, clamped edge, sigma~ = 0.02, the beads sink 1.52 a): a smooth repulsion from s = 3 (+0.05) through +0.1 (2.2), +0.24 (1.8), +0.4 (1.5) to +1.0 at contact (1.1); the inward and the outward
  leg agree within 0.03 except at isolated s (2.0, 1.4, 2.4: 0.1-0.17, another basin, see below).
* **flat disk, other wrapping fractions**: 98 % (budded) is flat to +0.03 down to s = 2.0, +0.06 at 1.45, +0.09 at 1.2, the same weak repulsion as the budded beads on the vesicle. 75 % and 60 % (wide-open necks) repel
  strongly and from far away: +0.14 at s = 3, +0.28 at 2.4, +0.5-0.7 at 2.0, +1.0 at 1.65, +1.45 at 1.45; below 1.3 the patches start to lose coverage (60 %: 0.586 to 0.543 at 1.1) and u_tot jumps by several kappa.
  **The idea that the attraction lives at 50-60 % wrapping is the opposite of what this model gives on the disk: the less wrapped, the more repulsive** (the pilot's 30-60 % chains on the vesicle said the same, noisily).
* **free contact line** (plain `Adhesion` W = 3.5, area 5 %, beads outside, d = 3 start): one bead went to a closed neck (cov 0.99), the other stayed at 0.79-0.83, and u_tot is flat, -0.02 +- 0.02 for s = 4.6-2.2 and
  +0.07 at 1.4 (the bending energy rises by 0.67 there but the adhesion pays it back).
* **tense closed vesicle** (sigma~ = 3, W = 20, beads inside, 98 %; the blue points of Fig. 4a): the tension screens the interaction (lambda = sqrt(kappa/sigma) = 0.58 a), u_tot = 0 +- 0.04 down to s = 1.3 and +0.11 at 1.0.
  The version with the beads outside broke in its cold start (the integrator) and was not repeated.

**What the states look like.** With the beads inside (`figures/slices_inside98.png`) the vesicle is a sphere with two buds hanging from short necks of radius 0.26 a on its outside (bead centres 8.4 a = 1.35 R_v from the centre); as s goes from 5 to 0.9
the two necks merge into one stalk and the mother between them forms a small arch. With the beads outside and on the disk (`figures/slices_planar_outside98.png`) each bead sits in its own pocket and the membrane between two pockets is a thin finger
of vesicle interior (a tongue below the beads at s = 1.55; a small bubble above the disk at s = 1.1) that tightens as they approach.

**Why it is repulsive here.** A pair of wrapped beads is a pair of strongly curved necks connected through the mother membrane; the energy that changes with s is the bending energy of the membrane between them (the caps are constant, the adhesion too).
That membrane gets thinner and more curved when the beads approach: a finger or a stalk (a tube of radius r and length L costs about pi kappa L / r). The paper's own words are that the linear theory predicts a bending-mediated repulsion and that a nonlinear one is needed for the attraction;
the nonlinear continuum solution of this paper's geometry, which the chain computes, has none down to a bead gap of 0.4 a (inside) or contact (outside, disk).

**Noise and basins.** Replicas of the 98 % chains: 0.015 kappa; plateau scatter on the disk 0.01 (98 %), 0.03 (90 %), 0.05-0.1 (60-75 %). The partially wrapped chains of the pilot on the vesicle: 0.1-0.5. On the disk at 90 %: the inward chain has an isolated point 0.17 kappa below the trend (s = 2.0, coverage 0.889 instead of 0.884, a basin
with more wrapped area and more bending) that the outward chain does not find: basins 0.1-0.2 kappa apart exist at this wrapping, as large as the paper's well, and a single chain gives an upper bound of the minimum energy. The simulated annealing of the paper
finds the lower ones; the pilot's replicas and the in/out legs are the poor man's version (take the minimum of the legs at each s: it is still above zero everywhere here).

**What did not work.** A frozen-bead cold start for beads inside (the vesicle has to slide away from the beads; 60000 steps and not converged, free beads: 5000 steps to the wrapped state); free patch axes on the disk (the chain left its wrapped state at the third waypoint, E +9 kappa) and in the vesicle (closed neck, section 5);
a tense disk (sigma~ = 3, W = 20, R = 12) as a control: the clamped disk domes and the frozen beads feel a vertical force of -5 kappa/a, the "attraction" it shows (-0.5 kappa, smooth over 5 D_p) is most likely that artifact, the run was discarded.

## 7. What is left to try (cluster: `Scripts/Mem_two_vesicle_up.sh`)

The continuum elastic model with a pinned contact line does not give the attraction of the paper in any of the geometries above, replicated where it matters (inside 98 %: 0.015 kappa) and cross-checked (patch vs free contact line, floppy vs 5 % area, in/out legs, vesicle vs disk, inside vs outside).
What differs from the paper's simulation and from the experiment, in order of how likely each seems:

1. **The paper's discretization.** Its curvature energy sqrt3 kappa sum(1 - n_i.n_j) saturates for sharp folds (a hinge costs at most 2 sqrt3 kappa = 3.5 kappa) and its vertices have a hard core (minimum edge l, "two layers of vertices between the particles" at 1.5 D_p); the code's `Bending_tan` grows without bound, so the finger
   or stalk between two wrapped beads is presumably cheaper there than here (not tested). The necks of 98 % wrapped beads are one vertex wide (0.26 a = l) in the paper: their energy is set by the mesh. Testing it needs a saturating hinge energy (not in the code: add `Bending` with (1 - cos theta)) and a vertex hard core, or the authors' code.
2. **Topology.** Two pockets (or buds) that touch can merge into one pocket holding both beads, with one neck instead of two and no finger: a state that a mesh cannot reach by deformation and whose energy I did not compute (it saves the finger, 1 kappa at contact on the disk, and a neck). The experiment's beads at s ~ 1.0-1.4 D_p could be in it.
3. **Sampling.** Basins 0.1-0.2 kappa apart (section 6): several starts per separation, or an annealing/basin-hopping step, then the minimum. Replicas (`--replica`) and the outward leg are in the launcher.
4. **Their unspecified details.** epsilon, K_A, K_V in units of the bead, whether the particles rotate (a free patch axis gains 0.03 kappa at s = 2.7 in the pilot, `--free_axis`, but it is fragile), the annealing schedule. The supplementary material only has the potential and theta_Wr. Worth an email to the authors, with Fig. 4b,c in hand: are the particles inside?
5. A larger closed vesicle (`--Rv 12.5`, 4 times the faces: the experiment's GUVs are 14-40 D_p), the tense vesicle with the beads outside, outward legs of the inside chains for the hysteresis (the planar 90 % legs differ by up to 0.17 at single separations).

Things that did not work and are not worth repeating: Coverage as the only wrapping control (soft shell, roughness), forcing 0.8 coverage from the tangent contact (crumples), moving the beads by more than 0.5 D_p in one link (unwraps them), the free patch axis as default (closed neck in one test, left the state on the disk),
frozen beads in the cold start of an inside run, a clamped disk with a high tension.
