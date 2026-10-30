#!/bin/tcsh
#$ -cwd
#$ -j y

# Two beads wrapped from the INSIDE of a closed vesicle, with the Coverage energy: the same chain as the outside batch of Mem_two_vesicle_up.sh (section 4, Surface_tension),
# but each bead sits inside the membrane ("outside": -1) and the faces that count are the ones whose outward normal points away from it.
# One sbatch job per curve u(s): a cold start at a bead distance of 6 D_p (D_p = 2 a = 1) with the beads frozen, then the straight (chord) distance between the bead centres goes
# 5.8 ... 1.0 and back (Run_two_vesicle_chain.py, the beads frozen), warm-starting every separation from the previous one.
# Results: ../Results/WrappingTwoInsideVesicle/<chain>/link*/ ; afterwards
#   python3 Analyze_two_chains.py --root ../Results/WrappingTwoInsideVesicle --pattern 'TwoInteraction_*' --tag SurfaceTension
# (the analysis folder goes next to the chains: ../Results/WrappingTwoInsideVesicle/two_beads_analysis_SurfaceTension)
# Needs a build of main_cluster with the Coverage energy for beads inside (Bead::side); an older build refuses the config when it loads.
# Units: bead radius 0.5, kappa = KB/2 = 21; the vesicle has R_v = 6.25.

set Batch_Tag = "WrappingTwoInsideVesicle"

foreach rep ( 0 )
foreach cov ( 0.3 0.7 0.9 0.95 )
foreach sigma ( 0.0 0.2 1.0 )

set Tag = "TwoInteraction_r_0.5_coverage_${cov}_tension${sigma}_replica_${rep}_SurfaceTension"
python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --inside --area tension --sigma ${sigma} --replica ${rep} --switch "100, 25000" --wrap coverage --cov ${cov} --KB 42 --KV 4200000.0 --radius 0.5 --Kcov 5000 --Rv 6.25 --distance chord --d_cold 6.0
sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
end
end
end
