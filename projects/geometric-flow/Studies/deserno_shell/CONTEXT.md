# Context for the Deserno / adhesion-shell study (read this first)

Worktree: `/home/mrojasve/Documents/DDG/MM-deserno` (git branch `deserno-shell-study`, based on `main` cf153b18).
Project folder: `P = /home/mrojasve/Documents/DDG/MM-deserno/projects/geometric-flow`.
Study folder (tracked, put everything worth keeping here): `P/Studies/deserno_shell/`
(`theory/`, `sim/`, `figures/`, `literature/`, `data/`).
The original checkout `/home/mrojasve/Documents/DDG/MiraculousMembranes` is shared with other live sessions:
**never edit, build in, or run git commands that change anything in it.** It is read-only for you; it only
holds the old simulation data (`projects/geometric-flow/Results/...`, gitignored, so absent from the worktree).

## The question

One bead (radius a) adhering to a planar membrane, simulated with `main_cluster` (triangulated membrane,
`Adhesion` + `Bending_tan` + `Excess_tension`), is compared with the continuum theory of
M. Deserno, *Elastic deformation of a fluid membrane upon colloid binding*, arXiv:cond-mat/0303656
(Fig. 2: phase diagram in (w~, sigma~); Fig. 5: a/lambda vs w/sigma). The simulated phase lines do not match:
the binding jump (flat start) and the line E (partially wrapped <-> fully enveloped, found by comparing the total
energy of a flat start and a full start) lie at larger w~ than the theory, and the enveloped state unwraps at
w~ ~ 7 although the paper's second spinodal S2 is ~2-3 (no unwrapping barrier in the simulation).
Mau's working hypothesis: the soft adhesion shell and the mesh resolution (vertices are not exactly on the
sphere, the contact-line curvature of the discrete mesh is not exactly 1/r) cause it. Mau wants (1) to get
closer to the theory by changing the shell (steeper shells such as cos^2 / cos^4 of the current one, narrower,
finer mesh), and (2) to understand why they disagree and whether a formal argument defends the disagreement.

## Mapping between code and Deserno's variables (audited against the C++)

* kappa = KB/2 (code energy is KB * H^2 with H = (c1+c2)/2, i.e. (kappa/2)(2H)^2), sigma = KA, a = BeadRadius,
  w = KI (adhesion energy per area, `Adhesion` only; E_ad = -KI a^2 sum_f w(r_f) |omega_f|).
* w~ = 2 w a^2 / kappa = 4 KI a^2 / KB, sigma~ = sigma a^2 / kappa = 2 KA a^2 / KB, lambda = sqrt(kappa/sigma),
  a/lambda = a sqrt(2 KA/KB), w/sigma = KI/KA, z = 1 - cos(alpha) in [0,2] = 2 * CoverageUnion (CoverageUnion is
  the post-processed covered fraction of the bead surface).
* Deserno's energy in his units, E~ = E/(pi kappa), lengths in units of a:
  E~(z) = -(w~ - 4) z + sigma~ z^2 + F(z; sigma~), F = int ds r [ (psi' + sin(psi)/r)^2 + 2 sigma~ (1 - cos psi) ] over the
  free (unbound) membrane (profile r(s), tangent angle psi, r' = cos psi, contact at r = sin(alpha), psi(0) = alpha, z = 1 - cos alpha);
  the cap contributes -(w~-4) z + sigma~ z^2 (cap bending 4 pi kappa z, adhesion, tension). See the docstring of
  `Scripts/deserno_theory.py` (read it) for the validated normalisation and API.
  Exact identities (validated numerically): dF/dz = (1 - psi0dot)^2 - 4 - 2 sigma~ z, dE/dz = (1 - psi0dot)^2 - w~,
  with psi0dot the free-membrane curvature at the contact line in units 1/a; contact condition
  w~ = (1 - psi0dot)^2 (Seifert-Lipowsky: w = (kappa/2) (curvature jump)^2).
  Phase lines: W (w~ = 4, small sigma~), E (equal energy partial <-> full), spinodals S1 (binding), S2 (unwrapping).
  At sigma~ = 0.53 the paper has w~_E = 5.1. Planar-run data covers sigma~ in [0,1] (KA in [0,0.5], KB = 1, a = 1).

## The simulated adhesion model (what actually is minimised)

`src/BeadGeometry.cpp` (`coverageShellWeight`, `evaluateFaceCoverage`), `src/Interaction.cpp` (`Adhesion::Face_Energy`,
`Adhesion::Add_Face_Force`), `include/Interaction.h`:

    E_ad = - W a^2 sum_f  w(r_f) |omega_f|        (faces f whose outward normal points towards the bead)
    w(r) = 1/2 (1 + cos(pi x)),  x = (r - a) / (s a),  |x| < 1, else 0       s = shell half-width fraction = 0.25 default
    r_f  = distance bead centre -> face centroid;  omega_f = exact solid angle of the triangle seen from the bead
    selection: dot(radial direction, face normal) < 0 (not differentiated); force on the bead included.

Note w(r) is cos^2(pi x/2): symmetric around r = a, there is NO hard core (the membrane may sit inside r < a),
the weight is 1 only exactly at r = a. `"shell_width"` (Adhesion constants[3], `BeadSpec::shell_width`, parsed in
`SimConfig.cpp`) sets s; default 0.25 is bit-identical to the old hard-coded value.
Other energies: `Excess_tension` constants [KA, factor]: E = KA (A - A0) if A >= A0 else 0, A0 = factor * initial mesh
area (flat start: factor 1; the full start uses Planar_full.obj with factor = area(flat disk)/area(full) = 0.9598874783696941
so A0 = 314.138, the flat disk). `Bending_tan` constants [KB, 0]. Membrane = free-boundary disk, radius ~10, `"boundary": true`.
The first Bending_tan radial cut (rho^2 > 1.6 with boundary) was removed earlier because energy and force disagreed.

Simulation set-up (template `P/Templates/Wrapping_planar.txt`, generator `P/Scripts/Create_subjob_planar.py`):
flat start `input/Big_planar_mem.obj` (disk in the yz plane, bead centre at x = a above it, soft potential so no contact
initially) or full start `input/Planar_full.obj` (membrane already wrapped around the bead at (-1.82316, 0.0744237, 0.00342679),
valid for a = 1). Options: `--inter_str` (KI) `--radius --KA --KB --tension excess --init flat|full --shell_width
--size_min --size_max --batch_tag --unique_tag`. Integrators: BFGS then BFGS-Normal, `adapt_remesh: "quality"`, remesher
`size_max` (default 0.5 here) / `size_min` (0.02 here) bound the edge lengths, `stopping` block (window 50, patience 4,
bfgs_tol_E 5e-5, normal_tol_E 1e-8, ...). A typical coarse run is ~30000 steps / ~3 min on one core; fine meshes are much slower.
Observed facts: runs with finer meshes were not converged after 12000 steps (energy still dropping by up to ~0.9 in the last
2000 steps); a narrow shell (0.1) on 0.25 edges is barely resolved; a pilot at sigma~ = 0.53 with shell 0.10 + fine mesh moved
the E line to w~ 4.8-6.0 (paper 5.1) and restored metastability of the enveloped state, but the flat-start binding jump
did not move towards S1.

## Existing tools and data

* `P/Scripts/deserno_theory.py` + `deserno_table.npz` (numerical solver of Deserno's model, validated 53/53 checks by
  `P/Scripts/test_deserno_theory.py`; API: `free_energy, contact_curvature, total_energy, w_E, barrier, spinodal_S1,
  spinodal_S2, fig5_curves, plot_fig2_lines(ax, ylim), fig2_line_data, interp_logsigma, load_table`). Use it as THE reference
  for the zero-range / hard-contact limit. Read its docstrings. `Scripts/Plots.py`: `planar_phase_deserno` (simulation vs
  theory in Deserno's axes), `deserno_fig5`, `planar_phase_comp`.
* Old simulation data (read-only, original checkout, absolute path prefix
  `/home/mrojasve/Documents/DDG/MiraculousMembranes/projects/geometric-flow/Results/`):
  `Wrapping_planar_excess/Coverage_data.txt` (flat start) and `Coverage_data_full.txt` (full start) = output of `PostProcessing`
  (header `#### DIR KA KB KI BeadRadius rc CoveredArea CoverageUnion MultilayerFrac BeadX Area Excess_tension Bending_tan Bead Total_E`),
  416 (KA,KI) points each, KI 0.5-3.0 step 0.1, KA 0-1 in 16 values, KB = 1, a = 1, default shell 0.25.
  `Wrapping_planar/` is the same with plain surface tension (not for this study). Per-run folders hold `Output_data.txt`
  (time step Volume Area <energies> Total_E ...), `Final_state.obj`, `membrane_*.obj`, `Input_file.json`.
* `PostProcessing`: `build/bin/PostProcessing config.json` where config.json = `{"first_dir": "../Results/<batch>/"}`; it scans the
  numbered sub-folders, **appends** to `<first_dir>Coverage_data.txt` (delete the old file first), skips unfinished runs.

## Environment and rules

* Python: `micromamba run -n mir_mem python ...` (numpy 2.3, scipy 1.17, matplotlib 3.10 with LaTeX, jinja2). System python
  has no matplotlib. Use `matplotlib.use("Agg")` in scripts.
* Build (worktree): `cd P/build && make -j4 main_cluster PostProcessing`. A frozen copy of the unmodified binaries lives in
  `P/build/baseline/` (use these for any run that must not change under you; never overwrite them). `build/bin/*` may be rebuilt by
  the agent that owns the C++ change, so copy a binary to your own location before launching a long campaign with it.
* No `sbatch` on this machine. Run a simulation as: `cd P/Scripts && <binary> ../Config_files/<name>.json` (configs written by
  `Create_subjob_planar.py` go to `P/Config_files/`, results to `P/Results/<batch_tag>/<number>/`; `make_numbered_dir` picks the
  number: give every concurrent run its own batch_tag or start them a few seconds apart).
* **CPU**: 8 cores shared with other sessions. Every simulation must be started through
  `P/Studies/deserno_shell/slot_run.sh <command>` (global semaphore, 6 slots, one thread per run, low priority).
  Never start more than a handful of your own waiting jobs at a time, run them with `nohup ... &` in the background and poll with
  short `sleep` loops (a single tool call is limited to 10 minutes). Pure Python numerics: keep to <= 2 cores
  (`OMP_NUM_THREADS=2`, `OPENBLAS_NUM_THREADS=2`).
* Use the scratchpad dir given in your environment for throw-away files. Results worth keeping: small data files (< 5 MB total
  each), scripts, figures (png/pdf) and a short markdown note go to `P/Studies/deserno_shell/<subfolder>/`. Large run folders
  stay in `P/Results/` (gitignored).
* Git: do **not** commit, checkout, stash, reset or change branches (the orchestrator commits). Do not touch files outside
  your assignment unless you must; if you must edit a shared file, say so in your report. Existing files in `P/Scripts/`
  (deserno_theory.py, Plots.py, ...) may be imported but should not be changed except where your task says so.
* Be honest: report what failed, what is unconverged, and what you could not verify. Numbers in your report must come from runs
  you did; do not extrapolate from memory. Mau is skeptical of claims that a finer mesh or narrower shell "should" fix things:
  state evidence, not expectations.

## Update after the first phase (read this, it corrects the text above)

* A usage limit interrupted the first batch of agents; some files in `theory/` and the tree may be partial. If you are a restarted
  workstream: look at what exists first (`git status`, your files), continue from it, and **write progress notes to your own
  `*_PROGRESS.md` file regularly** so that another interruption does not lose your work.
* Planar data cover sigma~ = 2 KA r^2/KB in [0, 2] (KA in [0, 1]), not [0, 1].
* `shell_power` is implemented, tested and committed on this branch (commit d1b2b7ba): bead option `"shell_power": p`
  (p >= 1, weight [1/2 (1+cos pi x)]^p, p = 1 bit-identical to before), `Create_subjob_planar.py --shell_power p`.
  Frozen binaries with it: `P/build/shellpower/{main_cluster,PostProcessing}` (use these for all new runs; `build/baseline` is the
  pre-change binary). A p = 4 run of 6000 steps takes ~2.5 min on one core with the old mesh settings.
* **Template trap**: commit 762c8cd0 changed `Templates/Wrapping_planar.txt` to `"rescale": 3.0` (scales the MESH about its centre
  of mass by 3, not the bead), default edge bounds `size_max 0.3 / size_min 0.005` and default init `Planar_mem.obj`.
  `Create_subjob_planar.py` always passes Big_planar_mem.obj (disk radius 10, area 314.138) or Planar_full.obj, so with the current
  template the flat start becomes a radius-30 disk (area 2827) and the full start is rescaled around a bead that is not.
  The OLD data used `rescale 1.0, size_max 0.5, size_min 0.05, Switch_times[1] = 30000`. **Always pass `--rescale 1.0`** (new option)
  and explicit `--size_min/--size_max` in this study. (Mau has been told.)
* Baseline analysis (B0, `sim/baseline.md`, `data/baseline_lines.csv`, figures `baseline_*.png`; old data, shell 0.25, mesh 0.05-0.5):
  in w~ the simulated lines lie above the theory: S2_sim ~ 5.4 at sigma~ = 0 rising to a plateau 7.4 for sigma~ >= 1.07
  (theory S2 falls 4 -> 1.9: no unwrapping barrier); E_sim = 1.32 * E_theory (E_sim 6.6 at sigma~ 0.53 where the paper has 5.1,
  11.1 at sigma~ 2 vs 8.4); S1_sim ~ 1.155 * S1_theory (S1_sim = S2_sim = 5.4 even at sigma~ = 0 where the theory says 4,
  i.e. the whole transition is shifted up, not only E). The energy of the fully enveloped state agrees with -2(w~-4) + 4 sigma~ to
  ~0.14 (pi kappa), so the energy scale is right. The simulated partially wrapped state lies ~0.3 w~ (pi kappa units) BELOW the
  theoretical partial branch (-0.56 at w~ 2-3, -3.0 at w~ 9-10), nearly independent of sigma~. For sigma~ <= 0.4 the two starts
  coincide, E cannot be defined. CoverageUnion of an unbound bead is already ~0.4 (soft shell).
* Physical bookkeeping the theory work should keep in mind (not yet verified numerically!): a finite shell changes THREE things:
  (i) the planar membrane within the shell range of the touching bead already gains adhesion (lowers the 'unbound' reference,
  which alone would delay binding), (ii) the membrane lifting off the bead beyond the contact line keeps gaining weight
  (lowers partially wrapped states, vanishing for z -> 0, 2, like a negative contact-line tension), (iii) mesh effects.
  A negative line tension alone cannot explain why binding at sigma~ = 0 moves from 4 to 5.4.
