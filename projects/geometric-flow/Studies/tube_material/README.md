# Material properties from the tube a bead pulls out of the planar membrane (`Mem_tube_material.sh`)

The chain version of `Mem_tube_up.sh`. There every final position XF is a separate job: a bead is pulled by a spring (a constant force of 200) from x = -r through the big planar membrane
to XF, 100 jobs for 2 values of KA. Here **one job per (KA, KB) walks a frozen bead** along x, every position relaxed from the previous one, 0.2 -> 6.0 -> 0.2 -> 0.0, and every relaxation is a
saved run (the way back in its own folder `way_back/`): the force on the bead, the energy and the shape of the membrane against x, and, from the way back, the hysteresis of the tube.

| file | what |
|---|---|
| `Templates/Tube_start.txt` | the start state: the old method (manual bead, spring, `Tube_eq.txt` numerics) up to x = 0.2 |
| `Templates/Tube_material.txt` | the frozen bead: the planar template of the cluster runs (BFGS, BFGS-Normal, stopping block) with the physics of `Tube_eq.txt` |
| `Scripts/Create_subjob_tube_material.py` | the generator: writes the configs, the chain file and the sbatch job (writes, does not run) |
| `Scripts/Run_tube_material_chain.py` | what the sbatch job runs: start state, link 0, the way out, the way back; checks every result; resumes |
| `Scripts/Analyze_tube_material.py` | force, energy, radius and area against x with the infinite-tube theory and the old scan |
| `Scripts/Mem_tube_material.sh` | the launcher (tcsh): KA 32 and 16, KB 4, radius 0.3, 0.2 to 6.0 in steps of 0.1 and back to 0.0 |

## 1. Physics and units (those of `Mem_tube_up.sh`)

`Surface_tension` KA (E = KA area, sigma = KA), the regular `Bending` KB (E = KB H^2, kappa = KB/2), the disk `input/Big_planar_mem.obj` (radius 10, plane x = 0, clamped edge, "boundary"), one bead of radius
0.3 with the shifted LJ of the tube runs (shift = strength = 10, cut at 2^(1/6) radius = 0.337 from the vertices; it acts on the VERTICES of the mesh), the remesher of `Tube_eq.txt`. Energies are in the energy unit of
the config (no division by kappa). The infinite tube has the force **f0 = 2 pi sqrt(KB KA)** (50.27 for KA 16, 71.09 for KA 32) and the radius **R0 = sqrt(KB/KA)/2** (0.25 and 0.177), and its energy grows by f0 per unit of length.
The force on the bead is the last row of `Bead_0_data.txt` (negative = the membrane pulls the bead back; the frozen bead has no spring, the force is the membrane's).

## 2. The start state (the old method) and why

A frozen bead cannot be moved into the flat sheet by the chain: the vertices of the planar mesh (edges 0.2, none within 0.11 of the axis) slide sideways and open a hole that the bead pokes through, since the LJ
only acts on vertices (a step of 0.1 from x = -r, and also steps of 0.025, ended in a bead through a face, a run that polishes for ever, or an energy of 1e5). The bead of the old method is moved a little at a time by the
spring, with remeshing, and carries the membrane with it. So the start state is made that way (`Tube_start.txt`, 250-600 steps with the Tube_eq mesh, ~40 s): the bead (state "manual", Velocity 0.5, a Lineal bond of 200 to a frozen
anchor at x = 6) is pulled from x = -0.3 and freezes by itself at FinalPos = **0.2** (`line_search_bead_step` 0.05 radii so that a fast bead cannot skip the acceptance ball), and relaxes for the rest of `--start_timesteps` (6000, a cap).
Its `Final_state.obj` is the mesh of the chain, the bead is then frozen at exactly x = 0.2 (link 0, **a data point**) and the spring is gone (the links are built from `Tube_material.txt`, never from the start run's config, so no
spring energy or force gets into the data). Below x = 0.2 the force is linear in x (F = 9.5 + 28.8 x for KA 16, 16.5 + 53.3 x for KA 32 in the old scan), so the way OUT does not go there; the way back does (0.2 -> 0.0, `--x_back`), from the
state of the way back. `--start_run <run folder or .obj>` takes a mesh of your own instead of the pull; `--x_start -0.3` is the flat start without a pull (does not work with this mesh).

**The mesh.** `Tube_eq.txt` has `refine_angle` 0.5, which leaves 4-8 vertices under the bead at x = 0.2 (edges 0.2 on the dimple): the vertex-only LJ then lets the bead sit inside a face (0.2 from the membrane after a 0.025 move, vertices
0.325 away, the force going DOWN with x) and the first relaxation of the chain stalled (BFGS-Normal) or was rejected. **`--refine_angle 0.25` is the default here** (the pull and the chain): the same chain ran smoothly to x = 0.6, the force
increasing monotonically (-15.8, -18.7, -21.9, -24.8, -27.7 at x = 0.2 ... 0.6, slope 29 against the 28.8 of the old scan), the bead 1.00-1.02 radii from the membrane, ~18 s per relaxation. 0.15 moves the pulled bead too slowly (x = 0.08 after 6000
steps); `--size_max 0.2` with 0.25 doubles the vertices (17k) and the time of a relaxation without a visible gain.

## 3. The chain

1. start state (above), link 0 at x = 0.2 (frozen, relaxed to convergence with the stopping block of the planar runs);
2. the way out: 0.3, 0.4, ... 6.0 (`--step 0.1`), every position reached through waypoints of at most `--max_step` 0.025 (four relaxations per position; waypoints are saved, kind "waypoint", and left out of the curves);
3. the way back: 5.9, ... 0.2, then **0.1, 0.0** (`--x_back`, default 0.0; 60 positions), each from the state of the way back (kind point, leg "out"), **saved in their own folder `way_back/`** inside the results folder
   (runs, configs and logs; `--back_dir` renames it; `chain_summary.csv` stays one file, its `rundir` column says where each run is). Defaults: 59 + 60 positions, about 473 relaxations, time limit computed from that
   (3 minutes each, at least 24 h, at most 168 h; one relaxation takes ~20 s at small x).

**Every relaxation is checked** (`check_run`): no 'nan line search', no BFGS-Normal that stalled and was ended, no 'small timestep' ending, 'Final state saved' present, the total energy did not RISE during the run (by more than
5 or 0.3 percent; a minimisation that ends high is not a result), finite, the force below 1e4, and the **bead rests on the membrane** (`bead_problem`: the nearest vertex at least 0.85 radius from its centre, its centre not closer to a face
than 0.7 of the distance to the nearest vertex (a bead inside a face of the mesh), the membrane within 1.4 radius of the centre and its tip at least 0.5 radius beyond it; limits from 78 healthy relaxed states of the old runs: vertex distance
0.91-1.1 radius, face distance 0.94-1.08, 0.78 at the lowest on a flat mesh with faces as large as the bead). Unreadable output counts as a failure. A relaxation that fails is retried in two half moves through the middle position; if that fails too the chain stops with a message in `chain.log`. A link that reaches its step cap (`--link_timesteps` 60000) is marked
`capped` (not converged) and the chain goes on. `Switch_times` of the links are 100, 1500; `remesh_rollback max_rise` is 1e9 (a remesh may raise the bead energy: the bead pokes between the vertices of a coarse mesh and only a vertex
inserted inside it feels it; with the default 10 the run polishes for ever).

**Resume.** `chain_summary.csv` (one line per run: link, kind point/waypoint/mid, leg in/out, x, run and folder, step, seconds, vertices, x_tip, d_bead, area, energies, Fx, Fy, Fz, capped) is written atomically after every run.
When the job runs out of time or stops, `sbatch` the same subjob again: the chain continues after the last finished run from the folder named in that line (a pull that was validated is reused, an unvalidated one is redone). The definition of the chain
(positions, way back, max_step, radius, KA, KB, kind of start) is stored in `chain_spec.json` and a resume with another definition is refused. Delete the folder to start over.

## 4. Looking at the results

```
python3 Analyze_tube_material.py --out ../Studies/tube_material/tube a=../Results/Tube_material/KA_32.0_KB_4.0_r_0.3 b=../Results/Tube_material/KA_16.0_KB_4.0_r_0.3
```
prints the table (with the hysteresis of the way back against the way out), writes `tube.csv` and `tube.png`: the force on the bead (dotted -f0), E_total - E_total(first), the radius of the tube (dotted R0), the area, circles the way out and
squares the way back, red crosses for capped links, grey crosses the old scan (`Results/Pulling_and_relaxing_2/Force_data.csv`, a spring pulled the bead to each x) for the same KA. Reference values of that scan: **F(0.2) = 15.6 (KA 16) and 27.6 (KA 32)**,
linear up to x about 1, a peak at the snap of the dimple into a tube (61.7 at x = 2.1 for KA 16, 97.6 at 1.8 for KA 32) and then the plateau 50.1 and 70.4 (99 percent of f0) beyond x = 3.


