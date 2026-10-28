#!/bin/tcsh
#$ -cwd
#$ -j y

# Barbell with beads of radius 1: two frozen beads inside a membrane with bending and surface tension and NO volume constraint, joined by a tube of radius
# R_t = sqrt(KB / (4 KA)). One sbatch job per (KA, KB), a few minutes each from the default start mesh. The old pipeline (Mem_barbell_up.sh) had beads of radius 0.2.
# Two batches (comment one out to run the other):
#   Bending: KB fixed per series (1, 2, 4), KA scanned (2..50)       -> ../Results/Barbell_beads/Bending_KB_<KB>_KA_<KA>/1/
#   Tension: KA = surface tension fixed per series (1, 2, 5), KB scanned (0.25..4) -> ../Results/Barbell_beads/Tension_KB_<KB>_KA_<KA>/1/
# (one folder per job: jobs that start together would otherwise race for the same numbered run folder); afterwards, per batch
#   python3 Analyze_barbell_beads.py --out ../Results/Barbell_beads/bending ../Results/Barbell_beads/Bending_*
#   python3 Analyze_barbell_beads.py --group KA --out ../Results/Barbell_beads/tension ../Results/Barbell_beads/Tension_*
# (tube radius, measured against sqrt(KB/(4 KA)), area and energies for every run: barbell.csv and barbell.png).
# Start mesh: --init mesh (default) = a barbell of revolution, relaxes in minutes; --init shrink = a sphere around both beads that shrinks onto them (slow, add --time 48:00:00).
# Beads at x = +-d/2, d = 10 by default (--d). The tube is thinner than a bead (R_t < 1) for KA > KB/4; mesh resolution: size_min is set from R_t by the generator.
# Smoke test of the mechanics on one run first: python3 Create_subjob_barbell_beads.py --KA 5 --KB 2 --batch_tag Barbell_beads_test, then run the config by hand.

set Batch_Tag = "Barbell_beads"
set Nsim = 1

# --- Bending batch: KB fixed, KA scanned. KA: 12 values from 2 to 50 equally spaced in the tube radius, i.e. in KA^(-1/2) (R_t = sqrt(KB/(4 KA)) changes fastest at small KA): dense there, sparse at large KA
foreach KB ( 1 2 4 )
foreach KA ( 2 2.3 2.7 3.3 4 5 6.3 8.3 11.5 17 27 50 )

set Tag = "Bending_KB_${KB}_KA_${KA}"
python3 Create_subjob_barbell_beads.py --KA ${KA} --KB ${KB} --Nsim ${Nsim} --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Tag}_Nsim_${Nsim} --d 10
sbatch ../Subjobs/${Tag}_Nsim_${Nsim}_subjob

end
end

# --- Tension batch: KA (the surface tension) fixed, KB scanned. KB: 12 values from 0.25 to 4 equally spaced in sqrt(KB), i.e. in the tube radius at fixed KA: dense at small KB
foreach KA ( 1 2 5 )
foreach KB ( 0.25 0.4 0.6 0.83 1.1 1.4 1.7 2.1 2.5 3 3.5 4 )

set Tag = "Tension_KB_${KB}_KA_${KA}"
python3 Create_subjob_barbell_beads.py --KA ${KA} --KB ${KB} --Nsim ${Nsim} --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Tag}_Nsim_${Nsim} --d 10
sbatch ../Subjobs/${Tag}_Nsim_${Nsim}_subjob

end
end

# Mau's start: beads at -5 and 5 inside a sphere of radius 6.5 that shrinks (no volume constraint)
# foreach KA ( 2 8 )
# python3 Create_subjob_barbell_beads.py --KA ${KA} --KB 2 --Nsim ${Nsim} --batch_tag ${Batch_Tag}_shrink/Bending_KB_2_KA_${KA} --unique_tag Shrink_KB_2_KA_${KA}_Nsim_${Nsim} --init shrink --d 10
# sbatch ../Subjobs/Shrink_KB_2_KA_${KA}_Nsim_${Nsim}_subjob
# end
