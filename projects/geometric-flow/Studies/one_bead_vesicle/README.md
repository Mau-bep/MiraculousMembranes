# One bead on a vesicle: chains that scan a control parameter (`Wrapping_one_vesicle`)

The sister of the two-bead pipeline of [`../two_beads_vesicle/README.md`](../two_beads_vesicle/README.md). There the two beads move along the vesicle in small steps and every
separation is relaxed from the previous one; here the bead stays where it is and **one number of the config changes** from link to link: the interaction strength, the
tension, the wrapped fraction of an adhesive patch, the excess area, the volume or the Coverage target. Each state is warm-started from the previous membrane, so a
link takes minutes where a cold start takes much longer, and neighbouring states share their mesh history (the noise of a fresh mesh drops out of the differences).

**Defaults: the bead interacts through the `Frenkel_Normal_nopush` potential cut at rc = 1.2 bead radii, and the patch is off** (free contact line). The patch (`--wrap patch`) and the
cosine-shell `Adhesion` (`--wrap adhesion`) are options; the patch exists only for `Adhesion`, there is no patch logic for the Frenkel potential.

| file | what |
|---|---|
| `Templates/Wrapping_one_vesicle.txt` | the config of the cold start (Jinja): bead, vesicle or flat disk, energies, remesher, stopping |
| `Scripts/Create_subjob_one_vesicle.py` | the generator: writes the config, the pre-wrap config, the chain file and the sbatch job (writes, does not run) |
| `Scripts/Run_one_vesicle_chain.py` | what the sbatch job runs: cold start, then one `main_cluster` run per value of the list |
| `Scripts/Analyze_one_vesicle.py` | energies, wrapped fraction and shape against the parameter; the coverage figure |
| `Scripts/Slice_one_vesicle.py` | cross-sections of the final states |
| `Scripts/Mem_one_vesicle_up.sh` | the launcher (tcsh): the transition chain active, the rest as commented examples |

## 1. Units and the potential

Bead radius a = 1, `KB` the bending constant of the code (energy `KB H^2`, `H = (c1+c2)/2`), **kappa = KB/2**. `--KB 42` is the kappa = 21 kBT of the paper: the energies come out in kBT.
The only physics is dimensionless, so the generator takes it that way too:

* **w~ = 4 W a^2/KB = 2 W a^2/kappa**, W the depth of the potential per unit area of membrane in contact (`inter_str`): `--w_tilde`, `--scan wt`. For KB = 42, W = 10.5 w~ kBT/a^2.
* **sigma~ = sigma a^2/kappa = 2 sigma a^2/KB** (`Surface_tension` sigma): `--sigma_tilde`, `--scan st`. sigma~ = 0.02 is 7 nN/m for kappa = 21 kBT and a = 0.49 um (the floppy vesicles of the experiment have sigma~ < 0.03).
* The stiffness of the volume and area constraints (2e5 kappa) scales with kappa by default, so the problem at kappa = 21 is the one at kappa = 1 with every energy 21 times larger
  (checked: same coverage, energies 21 times larger, see section 6; the `Timestep` key of the template is not read by the code).
* The vesicle is 6.25 bead diameters across (`--Rv 6.25`, the paper's D_v/D_p = 50 l / 8 l).
* **The potential** (`--wrap frenkel`, the default; `--rc` is the cutoff in bead radii, default 1.2). The bead's interaction is `Frenkel_Normal_nopush`, the one of the earlier wrapping studies: the faces whose normal
  looks at the bead (the others feel nothing) have the energy `eps V(r)` times their area times the cosine between their normal and the direction to the bead (the mean over their vertices of `V(r_v) Q/r_v`, `Q` the normal dotted with the vector to the bead),
  with `V(r) = alpha (a^2/r^2 - 1)(rc^2/r^2 - 1)^2` of the distance `r` from the bead centre and `alpha = 2 (rc/a)^2 [3/(2((rc/a)^2 - 1))]^3` (114.1 for rc = 1.2). It is zero at the bead surface `r = a`, positive (a repulsive core) inside it,
  has its minimum `-1` at `r_min = a sqrt(3 rc^2/(a^2 + 2 rc^2))` = **1.055 a** and is zero with zero slope at `rc`. For rc = 1.2, V/eps is +100 at 0.8 a, +16 at 0.9, +4.4 at 0.95, 0 at 1.0, -0.65 at 1.02, -1 at 1.055, -0.72 at 1.10, -0.22 at 1.15 and 0 at 1.2:
  half the depth is reached between 1.014 and 1.121 a, the well is 0.2 a wide.
* **W = eps**, the energy per unit area of membrane lying in the minimum, so `w~ = 4 eps a^2/KB` is still the control variable. The membrane wraps a sphere of radius `r_min`, not `a`, so Deserno's sharp contact (zero tension) jumps at `w~ = 4 (a/r_min)^2` = **3.59**
  (3.77 for rc = 1.1, 3.52 for 1.25, 3.26 for 1.5; `Analyze_one_vesicle.py` draws that line). There is no continuum theory of the closed vesicle with this potential.
* **Coverage.** The wrapped fraction of the plots is the geometric one for this potential: the solid angle of the bead covered by membrane faces that look at it with their centre at `r` within `rc - a` of the radius (`cap_cfg`, the well). The analysis also gives
  the fraction the bead energy measures, `f_ad = -E_bead/(4 pi eps r_min^2)`: equal to the geometric one for a membrane lying in the minimum, lower where it does not; and the plain fraction in a shell of +-25 percent (`cap25`, `cov25`, what the Coverage energy constrains).
  For `--wrap adhesion` (the `Adhesion` interaction) `G = -E_adhesion/(2 pi W a^2) = (1/2pi) sum_f w(r_f)|omega_f|` runs from 0 to 2 (a closed sphere), `f = G/2` is the wrapped fraction of the plots (the `z_ad` of the planar comparisons and the `G` of `Studies/deserno_shell/theory/softshell_direct.py`).

## 2. The scanned variables (`--scan`, `--scan_list`; default `wt`)

| `--scan` | what changes | needs | default list | largest step |
|---|---|---|---|---|
| `wt` | w~, the strength of the bead | `--wrap frenkel`, `adhesion` or `patch` | 30 values, 2 ... 12 (section 4) | none |
| `W` | W in the energy unit of the config per a^2 | same | 1,2,3,4,5,6,8 | none |
| `KB` | the bending constant KB (kappa = KB/2); W, sigma, excess and volume stay absolute, so w~ and sigma~ change along the chain (the stiffness of the constraints stays at 2e5 kappa of the first value) | none | 1,2,4,8,16,32 | none |
| `st` | sigma~ | `--area tension` or `excess` | 0.5 ... 0.005 | none |
| `sigma` | sigma in the energy unit of the config per a^2 | same | 3 ... 0.01 | none |
| `f` | the wrapped area fraction of the patch, (1 - cos theta_Wr)/2 | `--wrap patch` | 0.3 ... 0.98 | 0.05 |
| `theta` | the patch angle theta_Wr, degrees | `--wrap patch` | 60 ... 165 | 6 |
| `excess` | the relative excess area of the fixed-area vesicle | `--area constraint` or `excess` | 0.2 ... 0.05 | 0.02 |
| `volume` | the volume, fraction of the sphere of radius R_v | closed vesicle | 1.0 ... 0.8 | 0.05 |
| `cov` | the target of the Coverage energy | `--wrap coverage` | 0.2 ... 0.8 | 0.05 |

A decreasing list starts at its end: the cold start is at `--cold_value` (default the first value) and the chain walks the list from there. A change larger than the largest step is split into waypoints
(`--max_step`, 0 turns them off): the contact line of a patch advances by a dtheta per relaxation, a free contact line with a soft potential follows a change of w~ of 0.1 without help. The generator checks that the
options fit the variable (for instance `--scan st` with `--area constraint` is refused: there the tension is the multiplier of the constraint; `--scan f` needs `--wrap patch`, which is off by default).

## 3. The protocol

1. **Pre-wrap** (a bead outside a closed vesicle; default for `--wrap patch`, and for any mode when `--prewrap_cov` is given): the stiff Coverage energy wraps the bead from the tangent contact (`--prewrap_steps` 8000, target `--prewrap_cov`, default the patch fraction - 0.2
   between 0.1 and 0.5 for a patch, 0.3 otherwise), because a patch run from the tangent contact needs 30000 steps just to wrap. Without it (`--prewrap none`, the default of the other modes) the chain starts from the tangent contact.
2. **Cold start** (`link00_cold`): the full config at the cold value, to convergence. It is not a data point.
3. **Links**: for every value of the list a `main_cluster` run from the previous `Final_state.obj`: the bead stays frozen where it was, the patch axis (if any) points to the centre of the vesicle again, the volume and area
   targets are absolute (also for `Excess_tension`, whose target the loader reads relative to the initial mesh: the runner sets the factor of every link from the area of its own mesh), the value goes into the config.
   `Switch_times` 100, 1500. A failed link (the integrator broke, no 'Final state saved', fewer than 50 steps) is retried once through the middle value, then the chain stops.
4. **Return leg** (default on in the generator, `--no_return` turns it off, the launcher does): the list backwards, ending at the cold value: the hysteresis check.

Output in `../Results/<batch_tag>/`: `link00a_prewrap/`, `link00_cold/`, `link01_wt1.6/1/` ... (each a normal run folder), waypoints `_w1, ...`, `chain_summary.csv` (one line per run: link, leg in/out, variable, value, energies, steps, seconds) and `chain.log`.

**Customisation** (`python3 Create_subjob_one_vesicle.py --help`): `--wrap frenkel|adhesion|patch|coverage` and `--rc` (Frenkel), `--shell_width/--shell_power/--shell_profile` (the cosine shell and the profiles of `Adhesion`, not of the Frenkel potential),
`--patch_angle/--patch_width/--free_axis` (patch only; the free axis needs commit 8ced2182), `--planar` (flat disk, `--R_disk`), `--inside` (patch or adhesion), `--free_cold`,
`--W/--w_tilde`, `--area tension|constraint|excess|none`, `--sigma/--sigma_tilde`, `--excess`, `--volume`,
`--KB`, `--bending`, `--KV/--KA/--Kcov`, the mesh (`--size_min/--size_max/--refine_angle`), how long a run is (`--timesteps`, `--window`, `--patience`, `--bfgs_tol`, `--normal_tol`, `--stall_steps`, and `--link_timesteps/--link_bfgs_tol/--link_normal_tol`
for the links only; `--normal_tol 1e-7` shortens a link by 20-30 percent for 0.0004 kappa), the pre-wrap, `--replica`, `--time`.

## 4. The transition chain (kappa = 21, low tension, soft potential)

The plan: at zero tension Deserno's transition from no wrapping to full wrapping is discontinuous for a sharp contact (a first-order jump at w~ = 4), but with a finite, soft potential the wrapped
fraction rises continuously with the strength (Studies/deserno_shell: no N-shaped w~(G), no barrier), so partial wrapping, 50-60 percent included, is a stable state of a free contact line. The chain scans w~ and measures
the coverage.

* **Setup** (`Scripts/Mem_one_vesicle_up.sh`, section 1): `--scan wt --wrap frenkel --rc 1.2 --KB 42 --sigma_tilde 0.02 --Rv 6.25`, outside, the bead frozen, the vesicle free, no patch. `flat` (from weak adhesion up) is active, `full` (the list down from a fully wrapped cold start) is one line away.
  The command in the launcher carries `--refine_angle 0.30` and `--prewrap_cov 0.3`, set by hand after the first version; the mesh numbers below say why 0.10 is the value for this potential (the generator's default when the flag is absent; the examples in the launcher use 0.10).
* **Where the transition is** is not known for this potential. Sharp contact at the minimum: w~ = 3.59. The `Adhesion` cosine shell has a continuum theory (planar disk, sigma~ = 0, s = 0.25, p = 1, R = 10; `Studies/deserno_shell/data/softshell_sigma0.csv`):
  wrapped fraction 0.23 at w~ = 3.0, 0.38 at 3.5, **0.50 at 3.86, 0.55 at 4.0**, 0.69 at 4.5, 0.79 at 5, 0.90 at 6, 0.94 at 7, 0.97 at 9, and **a local pilot of that shell on this vesicle** (kappa = 21, sigma~ = 0.02, R_v = 6.25, cosine shell 0.25, refine_angle 0.10, forward sweep
  from weak adhesion, 4 points, no return leg) gave f = 0.11 at w~ = 2, 0.14 at 3, 0.20 at 4 and **0.94 at 5**: the closed vesicle moved the transition to the right of the planar theory, into 4-5. For the Frenkel potential there is one cold start
  at w~ = 4.5 (6000 steps, not converged): 56 percent of the bead covered in the well, 51 percent by the energy measure. That is all that is known; it puts the rise in the window below.
* **The list**: 30 values, `2, 2.6, 3, 3.3, 3.5, 3.7` below, **`3.8 ... 5.2` in steps of 0.1 (15 values)**, `5.4, 5.6, 5.8, 6, 6.5, 7, 8, 10, 12` above (a closed bead needs a strong adhesion to close the neck). If the first chain shows the rise
  elsewhere, rerun with a denser `--scan_list` there (`--cold_value` and `start_run` continue from a finished run).
* **Mesh.** The well is 0.2 a wide and the potential is evaluated at the vertices, so the edges near the bead must be about 0.1 once the membrane lies on it. The generator takes `--refine_angle 0.10` for this potential when none is
  given (the template's 0.3 for the other modes). A 5000-step test (not converged, kappa = 21, sigma~ = 0.02, 30 percent pre-wrap, a link at w~ = 2, then a link at w~ = 4.5), only the angle changed:

  | refine_angle | state | vertices with r < 1.2 a | edges there | bead covered in the well | E_bead [kappa] |
  |---|---|---|---|---|---|
  | 0.10 | w~ = 2 | 89 | 0.15 | 0.09 | -0.9 |
  | 0.30 | w~ = 2 | 14 | 0.40 | 0.09 | -1.3 |
  | **0.10** | w~ = 4.5 | **701** | **0.11** | **0.56** | -16.4 |
  | **0.30** | w~ = 4.5 | **51** | **0.24** | **0.16** | -4.1 |

  A cold start at w~ = 4.5 from the tangent contact without pre-wrap (0.10, 6000 steps) gave 793 vertices, 0.10 edges, 0.56 covered: the 0.10 state does not depend on the path, the 0.30 one has no dimple (the mesh cannot follow it). So with 0.30 the coverage curve is pushed to
  the right by the mesh; use 0.10. The pre-wrap (`--prewrap_cov 0.3`) leaves a fine mesh at its end (edges 0.06 with 0.10, 0.12 with 0.30) but the first, weak link unwraps the bead and the mesh coarsens again. The Coverage energy has no repulsive core, so the pre-wrapped bead
  sits inside the Frenkel core (845 vertices at r < a with 0.10) and the first relaxation has to push them out (slivers were polished twice, it ended at the unwrapped state); the pre-wrap costs 6 min with 0.10 (7700 vertices) and the cold start 1.7 min.
  Without a pre-wrap the unwrapped state at w~ = 2 is reached in 30 s. In the unwrapped regime the mesh at the contact stays coarser than 0.1 (the membrane is nearly flat there): the first points of the curve are the baseline, not the transition.
* **Analysis**: `python3 Analyze_one_vesicle.py --out ../Studies/one_bead_vesicle/transition flat=../Results/Wrapping_one_vesicle/transition_flat` writes `transition_coverage.png`: the wrapped fraction against w~ (and W on
  the top axis) with the sharp-contact line of Deserno and a band at 50-60 percent (and, for the cosine shell, the soft-shell theory of the planar disk), plus `transition.png` (energies, elastic energy against the wrapped fraction, tension, depth of the pocket) and `transition.csv`.

## 5. What to look at, and the limits

* Slices (`Slice_one_vesicle.py`) before trusting a number: the pocket should be axisymmetric around the line vesicle centre - bead.
* Hysteresis: the flat and the full chain (or the return leg) must agree where the transition is continuous; where they do not, the states are metastable and the sweep direction matters.
* Say which coverage you plot: for the Frenkel potential the plot shows the geometric one in the well and the energy-based one as open markers; where the membrane leaves the bead gently the two differ.
* The bead is frozen and the vesicle moves to it, as in the two-bead chains: at a strong interaction the first relaxations are violent (the pocket forms), the energies of those runs are not noise but the cost of wrapping from the previous state.
* The well is 0.2 a wide and the potential is evaluated at the vertices: the edges near the bead must be about 0.1 or less where the membrane lies on the bead (section 4).
* The restart noise of a relaxation is 0.01-0.04 kappa for a 98 percent patch pair and larger where the state has several basins (90 percent: 0.1-0.5 kappa); for a coverage curve use replicas (`--replica`) where it matters.

## 6. Checks done (2026-10-16, locally, in a scratch copy of the tree; nothing was written into the repository except the new files)

* **Generator**: the default call (`--wrap frenkel`, no options) writes the Frenkel bead (`Frenkel_Normal_nopush`, `inter_str` = w~ KB/4, `rc` 1.2, no patch keys, no pre-wrap), 30 values, refine_angle 0.10; `--rc 1.25 --planar`, `--wrap adhesion` (cosine shell 0.25), `--wrap patch --scan f` (patch keys,
  pre-wrap file) and the launcher itself (tcsh, sbatch replaced by a script that echoes) write their configs; the refusals are `--shell_width` (and the other Adhesion options) with the Frenkel potential, `--rc` without it, `--inside` with it, `--scan f` and `--free_axis` without `--wrap patch`,
  `--rc 0.9`, `--scan wt` with `--wrap coverage`. Earlier checks (cosine shell): every scan variable, `--planar`, `--inside`, `apply_scan`, the `Excess_tension` factor.
* **Binary**: 400 and 6000 steps of the Frenkel cold start at w~ = 4.5 (finite energies, no vertex in the core, the mesh refines around the bead: 793 vertices in the well, 0.56 covered); a chain of two values (4.5, 5, 400 steps per run, with the return leg) ran through (4 runs,
  `chain_summary.csv`); the two 5000-step chains of section 4 (pre-wrap 0.3, cold start at w~ = 2, link at 2, link at 4.5, refine_angle 0.30 and 0.10) ran through; the analysis wrote the csv, the overview and the coverage figure (with the Deserno line at 3.59 and no theory) for them.
* **kappa = 21** (cosine shell, earlier): `--KB 42` against `--KB 2` (all constraints scaled), w~ = 3.5 and 4: E_total/kappa 32.47 and 31.91 against 32.47 and 31.91, E_adhesion/kappa -3.51 and -5.14 against -3.48 and -5.07: equal within the restart noise (0.01-0.06 kappa).
  The `Timestep` key of the template is not read by `main_cluster`.
* **Not done**: the 30 links (every test run is 400-6000 steps and not converged), a long return leg, the fully wrapped cold start of the `full` chain (w~ = 12: expect it to be slow and the Frenkel core to be stiff, eps = 126 kBT/a^2), `--planar`, `--inside` and the patch scans with a binary,
  the other scan variables (configs and the runner's mapping only). Nothing in the repository was changed except the new files (no C++).
