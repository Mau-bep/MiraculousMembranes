#!/bin/tcsh
#$ -cwd
#$ -j y

# THE ENERGY OF EACH BRANCH OF VESICLE SHAPES against the reduced volume v at a fixed spontaneous curvature c0: Seifert, Berndl and Lipowsky, Phys. Rev. A 44, 1182 (1991), Fig. 11 (c0 = 2.4) and Fig. 13
# (c0 = 3). ONE CHAIN PER BRANCH: a chain relaxes a seed of the branch at its first v and walks v up (the target area shrinks, the stable direction of the flow) and/or down in steps of 0.01
# (0.005 for the pears' way down), every v relaxed from the previous shape, a step halved where the shape changes class, until the branch stops being a minimum (Run_shape_chain.py).
# Every chain is its own sbatch job on 1 cpu and starts from its own seed, so all of them run at the same time; a chain that stops or hits the time limit resumes from its finished relaxations when
# its subjob is sbatched again (../Subjobs/<Prefix>_c0_<c0>_<branch>_subjob). Remesher: size_min 0.005 (Templates/Shape_chain.txt), refine_angle 0.2, size_max 0.15 (Studies/shape_chains/README.md).
#
#   chain (folder in ../Results/Shape_chains/c0_<c0>/)  what it walks
#   prolate1   a capsule relaxed at v = 0.74 (c0 2.4) or 0.85 (c0 3), up to 0.99 and down to where the symmetric dumbbell becomes a pear
#   prolate2   a long capsule (v 0.55 / 0.40) relaxed at v = 0.60 (c0 2.4) or 0.50 (c0 3) and walked up until it loses its symmetry (c0 2.4: 0.60 to 0.72; c0 3: 0.50 to 0.62)
#   prolate2b  (c0 3 only) the rest of the prolate-2 branch, from a capsule of v 0.55 relaxed at 0.62 and walked up to 0.72, where the neck closes
#   pear       two spheres joined by a smooth throat (c0 2.4: radii 1.3, neck 0.2; c0 3: radii 2, neck 0.3) relaxed at v = 0.72 / 0.78 and walked UP to 0.85 / 0.92 (M^pear)
#   pearb      the same seed, walked DOWN in steps of 0.005 until the pear becomes symmetric (c0 2.4: C^pear near 0.706) or its neck closes (c0 3: near the budding limit v = 0.7458)
#
# WHERE THE DATA IS (all relative to this Scripts folder, like every Mem_*.sh): ../Results/Shape_chains/c0_<c0>/<chain>/
#   chain_summary.csv     one line per relaxation: tag, from_tag, leg, kind (seed, point, jump, rejected, failed), nu_target, nu (reached), F (= F_b/8 pi kappa = Bending/(4 pi KB)), E_total, c0_eff, H0, steps,
#                         status (converged, capped = the step cap, ended), seconds, nV, the shape numbers (length, rmax, r_big, r_small, rneck, neck_ratio, lobe_ratio, skew, nlobes), mesh, note
#   <tag>/1/              the run of one relaxation: Final_state.obj (the shape), Output_data.txt (energies), Convergence_log.txt, Input_file.json; <tag>.json (its config) and <tag>.log
#   seed.obj, chain.json  the seed mesh and the chain file; chain.log the decisions (steps halved, branch ends)
#   ../Results/Shape_chains/c0_<c0>/<chain>/<tag>_output.output  the slurm output (copied there when the job ends); the configs of the jobs are ../Config_files/<Prefix>_c0_<c0>_<chain>_{chain,ConfigFile}.json
# Results/, Config_files/, Subjobs/ and Outputs/ are not in git: copy ../Results/Shape_chains/ back from the cluster (rsync -av cluster:<path>/Results/Shape_chains/ ../Results/Shape_chains/).
#
# AFTER THE JOBS RAN, on your machine (matplotlib, scipy: micromamba run -n mir_mem), one figure per c0 with the shapes on the curves and the numbers behind it:
#   python3 Analyze_shape_chains.py --c0 2.4 --root ../Results/Shape_chains/c0_2.4 --out ../Studies/shape_chains/fig11 --data_dir ../Studies/shape_chains/fig11_data --all_meshes --paper --continuum
#   python3 Analyze_shape_chains.py --c0 3.0 --root ../Results/Shape_chains/c0_3.0 --out ../Studies/shape_chains/fig13 --data_dir ../Studies/shape_chains/fig13_data --all_meshes --paper --continuum
# (<out>.png/.svg the figure, <out>_table.csv and <out>_landmarks.md; <data_dir>/points.csv every point of every curve, snapshots.csv the shapes on the curves with their place, snapshots/ and meshes/
#  the OBJ files, to redraw the snapshots as renders; --continuum the check against the 1D axisymmetric solver, Shape_axisym.py)
# On this machine one chain takes 20 minutes (prolate1) to 2 hours (pear); locally without sbatch: python3 Run_shape_chain.py ../Config_files/<Prefix>_c0_<c0>_<chain>_chain.json
# Your own c0, seeds, ranges and mesh settings: python3 Create_subjob_shape_chain.py --help
#
# Options through the environment (setenv before running): Prefix (names of the jobs and configs, default Shape), Batch_Tag (folder in ../Results/, default Shape_chains), Binary (main_cluster, default
# ../build/bin/main_cluster), Submit (default sbatch; "echo" only writes the jobs and prints what would be submitted), Mesh_Study (1: also the study of the remesher settings, below).

if ( ! $?Prefix ) set Prefix = "Shape"
if ( ! $?Batch_Tag ) set Batch_Tag = "Shape_chains"
if ( ! $?Binary ) set Binary = "../build/bin/main_cluster"
if ( ! $?Submit ) set Submit = "sbatch"
if ( ! $?Mesh_Study ) set Mesh_Study = 0

set Make = ( python3 Create_subjob_shape_chain.py --binary ${Binary} )

# ------------------------------------------------------------------------------------------------ c0 = 2.4 (Fig. 11)
$Make --c0 2.4 --branch prolate1 --batch_tag ${Batch_Tag}/c0_2.4/prolate1 --unique_tag ${Prefix}_c0_2.4_prolate1
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_c0_2.4_prolate1_subjob

$Make --c0 2.4 --branch prolate2 --batch_tag ${Batch_Tag}/c0_2.4/prolate2 --unique_tag ${Prefix}_c0_2.4_prolate2
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_c0_2.4_prolate2_subjob

$Make --c0 2.4 --branch pear --nu_down none --batch_tag ${Batch_Tag}/c0_2.4/pear --unique_tag ${Prefix}_c0_2.4_pear
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_c0_2.4_pear_subjob

$Make --c0 2.4 --branch pear --nu_up none --batch_tag ${Batch_Tag}/c0_2.4/pearb --unique_tag ${Prefix}_c0_2.4_pearb
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_c0_2.4_pearb_subjob

# ------------------------------------------------------------------------------------------------ c0 = 3 (Fig. 13)
$Make --c0 3.0 --branch prolate1 --batch_tag ${Batch_Tag}/c0_3.0/prolate1 --unique_tag ${Prefix}_c0_3.0_prolate1
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_c0_3.0_prolate1_subjob

$Make --c0 3.0 --branch prolate2 --nu_up 0.62 --step 0.01 --batch_tag ${Batch_Tag}/c0_3.0/prolate2 --unique_tag ${Prefix}_c0_3.0_prolate2
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_c0_3.0_prolate2_subjob

$Make --c0 3.0 --branch prolate2 --seed_capsule 0.55 --nu0 0.62 --nu_up 0.72 --step 0.01 --batch_tag ${Batch_Tag}/c0_3.0/prolate2b --unique_tag ${Prefix}_c0_3.0_prolate2b
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_c0_3.0_prolate2b_subjob

$Make --c0 3.0 --branch pear --nu_down none --batch_tag ${Batch_Tag}/c0_3.0/pear --unique_tag ${Prefix}_c0_3.0_pear
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_c0_3.0_pear_subjob

$Make --c0 3.0 --branch pear --nu_up none --batch_tag ${Batch_Tag}/c0_3.0/pearb --unique_tag ${Prefix}_c0_3.0_pearb
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_c0_3.0_pearb_subjob

# ------------------------------------------------------------------------------------------------ the study of the remesher settings (optional: setenv Mesh_Study 1)
# about 50 relaxations of 5 shapes (sphere, prolate-1, prolate-2 and two pears) with other refine_angle, size_max, size_min, aspect_min and integrators, 8 at a time, then the remesh stress test:
# Studies/shape_chains/README.md, section 5. A day on 8 cpus.
if ( $Mesh_Study == 1 ) then
python3 Create_subjob_shape_study.py --unique_tag ${Prefix}_mesh_study --workers 8 --binary ${Binary}
if ( $status != 0 ) exit 1
$Submit ../Subjobs/${Prefix}_mesh_study_subjob
endif
