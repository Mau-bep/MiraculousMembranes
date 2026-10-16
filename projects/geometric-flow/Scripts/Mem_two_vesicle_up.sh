#!/bin/tcsh
#$ -cwd
#$ -j y

# Two wrapped beads: the pair interaction u(s) of van der Wel et al., Sci. Rep. 6, 32825 (2016), Fig. 3e / 4a.
# One sbatch job per curve u(s): it wraps the two beads, relaxes them with the adhesive patch and then walks the beads towards and
# away from each other, warm-starting every separation from the previous one (Run_two_vesicle_chain.py, 2-5 hours).
# Results: ../Results/Wrapping_two_vesicle/<curve>/link*/ ; afterwards, for one or several curves (replicas of one curve: same label prefix),
#   python3 Analyze_two_vesicle.py --ref mean_far --s_far 4.3 --curves total --ylim -0.5 1.0 --out ../Studies/two_beads_vesicle/u_of_s \
#           inside98=../Results/Wrapping_two_vesicle/inside_area5_98_rep0 floppy98=../Results/Wrapping_two_vesicle/inside_floppy_98_rep0
#   python3 Slice_two_vesicle.py ../Results/Wrapping_two_vesicle/inside_area5_98_rep0/link12_s2.70/1 ...    (look at the states)
# Needs a build of main from the merge of bead-patch-adhesion (68af27c6) or later (the adhesive patch, the inside option); the generator checks it.
# Units: bead radius 1, kappa = KB/2 = 1 (KB = 2), so u in kappa; u [kBT] = 20.9 u, length [um] = 0.49 length.
# Physics, choices, results (no attraction found so far) and what to scan next: ../Studies/two_beads_vesicle/README.md

set Batch_Tag = "Wrapping_two_vesicle"

# --- 1. the paper's simulation as far as it can be read from the paper: the beads INSIDE the closed vesicle (Fig. 4b,c), theta_Wr = 165 deg
#        = 98 percent of the bead (Fig. S5), volume and area fixed with A_t = 1.05 A_sphere (the floppy chain of this geometry has A = 1.0498 A_sphere
#        without any constraint: the budget is the area of two budded beads). Replicas: a negligible change of the start (the chaotic mesh history),
#        the noise is 0.03-0.1 kappa, the signal 0.15.
# foreach rep ( 0 1 2 )
# set Tag = "inside_area5_98_rep${rep}"
# python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --inside --area constraint --excess 0.05 --patch_angle 165 --replica ${rep}
# sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
# end

# --- 2. the same with a free area (floppy, sigma~ = sigma a^2/kappa = 0.02, the experiment's sigma < 10 nN/m), beads inside and outside
# foreach rep ( 0 1 )
# # set Tag = "inside_floppy_98_rep${rep}"
# # python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --inside --sigma 0.02 --patch_angle 165 --replica ${rep}
# # sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob

# set Tag = "outside_floppy_98_rep${rep}"
# python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --sigma 0.02 --patch_angle 165 --replica ${rep}
# sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
# end

# --- 3. the flat limit of the experiment's large vesicles (D_v/D_p >= 14): a clamped disk of radius 24 a, floppy, wrapping 60, 75, 90, 98 percent
#        (theta_Wr = 101.5, 120, 143.13, 165 degrees); no closure, no volume, the beads reach contact (s = 1)
# foreach angle ( 101.5 120 143.13 165 )
# set Tag = "planar_floppy_${angle}"
# python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --planar --d_cold 14 --R_disk 24 --s_list 7.0,6.0,5.0,4.2,3.6,3.2,2.9,2.6,2.4,2.2,2.0,1.8,1.65,1.5,1.4,1.3,1.2,1.1 --sigma 0.02 --patch_angle ${angle}
# sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
# end

# --- More to try (uncomment): the pilot's ensembles (beads outside, 50 and 90 percent, floppy and area 5 percent), a tense closed vesicle (--sigma 3.0 --W 20: the blue points of Fig. 4a;
#     not on the planar disk: with a tension the clamped disk domes and the frozen bead heights are 5 kappa/a off their equilibrium), a larger vesicle (R_v = 12.5 a: 4 times the faces, give it more time),
#     free patch axes (the beads can tilt: --free_axis; it once ran into a closed neck), a plain adhesion without the patch (--wrap adhesion --W 3.5: the contact line is free)
# foreach angle ( 90 143.13 )
# set Tag = "outside_area5_${angle}"
# python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --area constraint --excess 0.05 --patch_angle ${angle}
# sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
# end
# python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/inside_floppy_98_Rv12 --unique_tag ${Batch_Tag}_inside_floppy_98_Rv12 --inside --sigma 0.02 --patch_angle 165 --Rv 12.5 --time 120:00:00
# sbatch ../Subjobs/${Batch_Tag}_inside_floppy_98_Rv12_subjob


#== 4 . Mau experiment

foreach rep ( 0 1 )
# set Tag = "inside_floppy_98_rep${rep}"
# python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --inside --sigma 0.02 --patch_angle 165 --replica ${rep}
# sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob

set Tag = "shell_nopatch_rep${rep}"
python3 Create_subjob_two_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --area tension --sigma 0.42 --replica ${rep} --prewrap_cov 0.6 --prewrap_steps 40000 --switch "100, 25000" --wrap adhesion --W 48.51 --KB 42 --KV 4200000.0
sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
end