#!/bin/tcsh
#$ -cwd
#$ -j y

# One bead wrapped by a closed vesicle: chains that scan ONE control parameter, every state warm-started from the previous one (Run_one_vesicle_chain.py).
# One sbatch job per chain. Results: ../Results/Wrapping_one_vesicle/<tag>/linkNN_<variable><value>/1/ and chain_summary.csv; afterwards, on your machine,
#   python3 Analyze_one_vesicle.py --out ../Studies/one_bead_vesicle/transition flat=../Results/Wrapping_one_vesicle/transition_flat
#   (several chains: label=folder for each, label=folder1,folder2 for the legs of one series; the coverage figure transition_coverage.png needs matplotlib: micromamba run -n mir_mem)
#   python3 Slice_one_vesicle.py ../Results/Wrapping_one_vesicle/transition_flat/link12_wt4.4/1 ...    (look at the pockets)
# The potential of the bead is the Frenkel_Normal_nopush of the earlier wrapping studies (every build has it), cut at rc = 1.2 bead radii, and the PATCH IS OFF: the contact line is free
# (--wrap frenkel is the default of the generator; --wrap patch turns the patch on, it needs a build of main from the merge of bead-patch-adhesion, 68af27c6 or later, and is not available for
# the Frenkel potential: see the examples at the end). The generator checks the build for what it needs.
# Everything about the pipeline, the units, the scanned variables and the options: ../Studies/one_bead_vesicle/README.md ; python3 Create_subjob_one_vesicle.py --help

# set Batch_Tag =   "Wrapping_one_flat"

# --- 1. the wrapping transition at low tension with a soft potential: the coverage as a function of the interaction strength ----------------------------------------------
#        kappa = 21 kBT (--KB 42, the energies come out in kBT), the vesicle 6.25 bead diameters across (the paper's D_v/D_p, --Rv 6.25 bead radii), sigma~ = sigma a^2/kappa = 0.02 (7 nN/m,
#        the floppy vesicles of the experiment), no patch. The potential of the bead is Frenkel_Normal_nopush with rc = 1.2 (--rc, in bead radii): V(r) = eps alpha (a^2/r^2 - 1)(rc^2/r^2 - 1)^2 of the
#        distance r from the bead centre to the membrane, zero at the bead surface r = a, repulsive inside it, minimum -eps at r = 1.055 a, zero value and slope at 1.2 a, felt by the faces that look
#        at the bead. eps = W is the energy per area of membrane in the minimum, and the scanned variable is the dimensionless w~ = 4 W a^2/KB = 2 W a^2/kappa of the Deserno mapping
#        (W = 10.5 w~ kBT/a^2 for KB = 42). The membrane lies in the minimum, 1.055 a from the bead centre, so Deserno's sharp contact (zero tension) jumps at w~ = 4 (a/r_min)^2 = 3.59 here; the soft potential
#        and the closed vesicle smooth and move it, which is what the chain measures (no theory for this potential: the figure has the sharp line only).
#        30 values: a few below, steps of 0.1 from 3.8 to 5.2, a few up to 12. Where the rise is, is a guess: a pilot of the cosine shell (--wrap adhesion, below) went from 20 to 94 percent wrapped between
#        w~ = 4 and 5, and one cold start of this potential at 4.5 was 56 percent wrapped after 6000 steps and still going. If the first chain puts it elsewhere, give the chain a denser --scan_list there.
#        MESH: the well is only 0.2 a wide and the potential is evaluated at the vertices, so the edges near the bead must be about 0.1 once the membrane lies on it. The generator takes --refine_angle 0.10 for
#        this potential when none is given. In a 5000-step test at w~ = 4.5, 0.30 left 51 vertices in the well (edges 0.24) and 16 percent of the bead covered, 0.10 gave 793 vertices (edges 0.10) and 56 percent:
#        with 0.30 the mesh pushes the coverage curve to the right. --prewrap_cov 0.3 wraps the bead to 30 percent with the Coverage energy first (the mesh near the bead is fine then, 0.12); the first weak link
#        unwraps it again and the mesh coarsens (edges 0.3-0.4 at w~ = 2 with refine_angle 0.30), so it does not replace the fine refine_angle.
#        flat: from weak adhesion upwards (the cold start is the weakest); full: the same list downwards from a fully wrapped cold start (slower, untested: the cold start has to wrap the bead), for the
#        hysteresis. The way back along the list is switched off here (--no_return): leave it out to walk the list back too (doubles the runs).
#        Time: about 35 s per 5000 steps at 3000-6000 vertices, a link takes minutes, the whole chain a few hours; the time limit is generous.
# set ListUp = "2,2.6,3,3.3,3.5,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5.1,5.2,5.4,5.6,5.8,6,6.5,7,8,10,12"
# set ListUp = "2,2.5,3,3.5,4,4.2,4.4,4.3,4.4,4.6,4.62,4.64,4.66,4.68,4.7,4.72,4.74,4.76,4.78,4.8"
# set ListDown = `python3 -c "print(','.join(reversed('${ListUp}'.split(','))))"`
# foreach s_tilde (0.02 0.04)
# foreach Init ( flat )
# # foreach Init ( flat full )
# if ( ${Init} == "flat" ) set List = ${ListUp}
# if ( ${Init} == "full" ) set List = ${ListDown}
# set Tag = "detailed_transition_02_shell_${Init}"
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --shell_width 0.2 --unique_tag ${Batch_Tag}_${Tag} --scan wt --scan_list ${List} --wrap adhesion --KB 42 --sigma_tilde ${s_tilde} --Rv 6.25 --refine_angle 0.20 --time 24:00:00 --prewrap_cov 0.3
# sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
# end
# end


set Batch_Tag = "Phase_one_flat"

set ListUp = "0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5"
set List = ${ListUp}

foreach KA (0 0.06667 0.1333 0.2 0.2667 0.3333 0.4 0.4667 0.5333 0.6 0.6667 0.7333 0.8 0.8667 0.9333 1)


set Tag = "KA_${KA}"
python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --wrap adhesion --shell_width 0.2 --sigma ${KA} --unique_tag ${Batch_Tag}_${Tag} --scan W --scan_list ${List} --wrap adhesion --KB 1 --refine_angle 0.20 --planar --area excess --excess 0.0 --size_min 0.02 --size_max 0.8 --time 24:00:00
sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
end


set ListUp = "2.915,2.144,1.577,1.16,0.8536,0.628,0.462,0.3398,0.25,0.1839,0.1353,0.09953,0.07322,0.05386,0.03962,0.02915,0.02144,0.01577,0.0116,0.008536,0.00628,0.00462,0.003398,0.0025"
set List = ${ListUp}

foreach KI (0 0.3333 0.6667 1 1.333 1.667 2 2.333 2.667 3 3.333 3.667 4 4.333 4.667 5)


set Batch_Tag = "Phase_one_vesicle"
set Tag = "KI_${KI}"
python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --wrap adhesion --shell_width 0.2 --W  ${KI} --sigma 1.0 --unique_tag ${Batch_Tag}_${Tag} --scan KB --scan_list ${List} --wrap adhesion --refine_angle 0.20 --area excess --excess 0.0 --size_min 0.01 --size_max 0.6 --time 24:00:00
sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
end


# --- More to try (uncomment) ---------------------------------------------------------------------------------------------------------------------------------------------------
# the same at other tensions (sigma~ = 0 is a free area, 0.1 and 0.3 are stiffer), and replicas for error bars
# foreach St ( 0 0.1 0.3 )
# foreach rep ( 0 1 )
# set Tag = "transition_flat_st${St}_rep${rep}"
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --scan wt --scan_list ${ListUp} --wrap frenkel --rc 1.2 --KB 42 --sigma_tilde ${St} --Rv 6.25 --refine_angle 0.10 --no_return --replica ${rep} --time 24:00:00
# sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
# end
# end

# other cutoffs of the Frenkel potential (1.25 is the rc of the earlier vesicle studies; a shorter range is closer to the sharp contact, a longer one softer; the well is rc - 1 bead radii wide, keep the edges near the bead below a fifth of it)
# foreach Rc ( 1.1 1.25 1.5 )
# set Tag = "transition_flat_rc${Rc}"
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --scan wt --scan_list ${ListUp} --wrap frenkel --rc ${Rc} --KB 42 --sigma_tilde 0.02 --Rv 6.25 --refine_angle 0.10 --no_return --time 24:00:00
# sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
# end

# the Adhesion interaction with its cosine shell instead (no patch, free contact line; shell_width 0.25 is the soft shell of the Deserno study and of the pilot above, --shell_power 2 or 4 steepens it;
# --shell_profile frenkel, lj or morse are the repulsive-core profiles of Adhesion, a different thing from the Frenkel_Normal_nopush of --wrap frenkel, they need commit 7c4576c2)
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/transition_cosine --unique_tag ${Batch_Tag}_transition_cosine --scan wt --scan_list ${ListUp} --wrap adhesion --shell_width 0.25 --KB 42 --sigma_tilde 0.02 --Rv 6.25 --refine_angle 0.10 --no_return --time 24:00:00
# sbatch ../Subjobs/${Batch_Tag}_transition_cosine_subjob
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/transition_p4 --unique_tag ${Batch_Tag}_transition_p4 --scan wt --scan_list ${ListUp} --wrap adhesion --shell_width 0.25 --shell_power 4 --KB 42 --sigma_tilde 0.02 --Rv 6.25 --refine_angle 0.10 --no_return

# --- The patch (OFF above): Adhesion limited to an angular patch, the pinned contact line of the paper (--wrap patch); the Frenkel potential has no patch ---------------------------------------
# the wrapped fraction is then set by the patch and the chain scans it, 30 to 98 percent, floppy and with the paper's 5 percent excess area (cosine shell of half width 0.1, W = 6 kappa/a^2: w~ = 12)
# foreach Area ( tension constraint )
# set Tag = "fraction_${Area}"
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --wrap patch --scan f --area ${Area} --KB 42 --sigma_tilde 0.02 --excess 0.05 --time 48:00:00
# sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob
# end
# the transition with a 90 percent patch: the coverage rises with the strength and stops at the patch (--patch_angle 143.13; 165 is 98 percent); --prewrap none starts from the tangent contact
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/transition_patch90 --unique_tag ${Batch_Tag}_transition_patch90 --wrap patch --patch_angle 143.13 --scan wt --scan_list ${ListUp} --prewrap none --KB 42 --sigma_tilde 0.02 --Rv 6.25 --refine_angle 0.10 --no_return --time 24:00:00
# sbatch ../Subjobs/${Batch_Tag}_transition_patch90_subjob

# --- Other scans: the tension of a free contact line (Deserno's wrapping against tension), the excess area of a fixed-area vesicle, its volume --------------------------------------
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/tension --unique_tag ${Batch_Tag}_tension --scan st --wrap frenkel --w_tilde 5 --KB 42 --refine_angle 0.10
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/excess --unique_tag ${Batch_Tag}_excess --scan excess --area constraint --wrap patch --patch_angle 143.13 --KB 42
# python3 Create_subjob_one_vesicle.py --batch_tag ${Batch_Tag}/volume --unique_tag ${Batch_Tag}_volume --scan volume --area constraint --wrap patch --patch_angle 143.13 --KB 42
