#!/bin/tcsh
#$ -cwd
#$ -j y

# Two beads wrapped from the OUTSIDE on a planar membrane, with the Coverage energy: the large-vesicle limit of the outside batch of Mem_two_vesicle_up.sh (section 4,
# Surface_tension). The membrane is the big disk of input/Big_planar_mem.obj (radius 10 = 20 bead radii, clamped edge, normal along x, no volume), the beads start tangent above it.
# --bead_normal: the beads are frozen through gradient descent and BFGS, so the distance between them is exactly the one of the chain, and are released along the normal of the
# disk inside BFGS-Normal (Constraint Normal + release_in_normal): the depth they sit at is found by the energy instead of being chosen (a free bead is pulled through the
# membrane by the stiff Coverage force within 100 steps, a bead frozen at one radius above the disk has a depth that nothing decided).
# One sbatch job per curve u(s): a cold start at a bead distance of 6 D_p (D_p = 2 a = 1), then the straight distance between the bead centres goes 5.8 ... 1.0 and back
# (Run_two_vesicle_chain.py), warm-starting every separation from the previous one.
# Results: ../Results/WrappingTwoPlanarMembrane/<chain>/link*/ ; afterwards
#   python3 Analyze_two_chains.py --root ../Results/WrappingTwoPlanarMembrane --pattern 'TwoInteraction_*' --tag SurfaceTension
# Needs a build of main_cluster with release_in_normal (the generator checks it). Units: bead radius 0.5, kappa = KB/2 = 21. The disk radius 10 is 20 bead radii: with no tension
# (sigma = 0) the deformation reaches the clamped edge, tension 0.2 and 1.0 decay over sqrt(kappa/sigma) = 10 and 4.6; a bigger disk is --R_disk 20 (the mesh is rescaled by R_disk / 10).

set Batch_Tag = "WrappingTwoPlanarMembrane"

foreach rep ( 0 )
foreach cov ( 0.3 0.7 0.9 0.95 )
foreach sigma ( 0.0 0.2 1.0 )

set Tag = "TwoInteraction_r_0.5_coverage_${cov}_tension${sigma}_replica_${rep}_SurfaceTension"
python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --planar --bead_normal --area tension --sigma ${sigma} --replica ${rep} --switch "100, 25000" --wrap coverage --cov ${cov} --KB 42 --radius 0.5 --Kcov 5000 --R_disk 10 --d_cold 6.0 --s_list 5.8,5.0,4.3,3.8,3.4,3.0,2.7,2.4,2.2,2.05,1.9,1.78,1.66,1.55,1.45,1.35,1.25,1.15,1.07,1.0
sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
end
end
end
