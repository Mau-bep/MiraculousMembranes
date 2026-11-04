#!/bin/tcsh
#$ -cwd
#$ -j y

# The energy of each BRANCH of vesicle shapes against the reduced volume at a fixed spontaneous curvature (Seifert, Berndl and Lipowsky, Phys. Rev. A 44, 1182 (1991), Fig. 11 for c0 = 2.4 and
# Fig. 13 for c0 = 3): ONE chain per (c0, branch). A chain relaxes a seed of the branch at nu0 and walks the reduced volume nu up (the target area shrinks: the stable direction of the
# flow) and down in steps of 0.01 (0.005 where a branch ends), every nu relaxed from the previous shape (Run_shape_chain.py); a step is halved where the shape changes class, and the chain
# of a branch stops where the branch stops being a minimum (prolate-1 at the neck collapse, the pear when it becomes symmetric or jumps to prolate-1, prolate-2 where it loses its symmetry).
#   prolate1  a capsule (long cylinder with caps) relaxed at nu0, the symmetric dumbbell/prolate branch from the sphere (nu = 1) down to where its neck collapses
#   prolate2  a long capsule (nu 0.55, c0 = 2.4; 0.40, c0 = 3) relaxed at the lowest nu and walked up: two lobes joined by a neck, up to the pear instability
#   pear      two overlapping spheres (radii ratio q, smoothed neck): the asymmetric branch the symmetric seeds never reach
# Results: ../Results/Shape_chains/c0_<c0>/<branch>/ (chain_summary.csv, one folder per relaxation, chain.log). A chain that stopped or hit the time limit resumes from its finished
# relaxations when its subjob is sbatched again. After the jobs ran, on your machine (matplotlib: micromamba run -n mir_mem):
#   python3 Analyze_shape_chains.py --c0 2.4 --root ../Results/Shape_chains/c0_2.4 --out ../Studies/shape_chains/fig11 --paper
#   python3 Analyze_shape_chains.py --c0 3.0 --root ../Results/Shape_chains/c0_3.0 --out ../Studies/shape_chains/fig13 --paper
# Locally (one core each, 10 to 40 minutes): python3 Run_shape_chain.py ../Config_files/Shape_c0_<c0>_<branch>_chain.json
# Seeds, ranges, mesh (--refine_angle, --size_min, ...) and a c0 of your own: python3 Create_subjob_shape_chain.py --help ; ../Studies/shape_chains/README.md

set Batch_Tag = "Shape_chains"

foreach c0 ( 2.4 3.0 )
foreach Branch ( prolate1 prolate2 pear )

set Tag = "c0_${c0}/${Branch}"
python3 Create_subjob_shape_chain.py --c0 ${c0} --branch ${Branch} --batch_tag ${Batch_Tag}/${Tag} --unique_tag Shape_c0_${c0}_${Branch}
if ( $status != 0 ) exit 1
sbatch ../Subjobs/Shape_c0_${c0}_${Branch}_subjob

end
end
