#!/bin/tcsh
#$ -cwd
#$ -j y


# Pilot: does a narrower adhesion shell (and a finer mesh) move the simulated E boundary
# (partially wrapped <-> fully enveloped) towards the Deserno one (arXiv:cond-mat/0303656)?
# One surface tension, sigma~ = 2 KA r^2/KB = 0.53 (the paper has w~_E = 5.1 there, the first
# scans had the jump around w~ = 7), a scan of KI = w~ KB/(4 r^2), both starts.
#   A: shell 0.25 (code default), template mesh     B: shell 0.25, fine mesh
#   C: shell 0.10, fine mesh                         D: shell 0.05, finer mesh
# B separates the effect of the mesh from the one of the shell. Compare the coverage and the
# Total_E of flat and full starts at each KI with PostProcessing / planar_phase_deserno.
set Nsim=1
set KB = 1
set radius = 1.0
set KA = 0.2667

foreach Cond ( A B C D )
if ( ${Cond} == "A" ) set Extra = ( )
if ( ${Cond} == "B" ) set Extra = ( --size_min 0.02 --size_max 0.25 )
if ( ${Cond} == "C" ) set Extra = ( --shell_width 0.1 --size_min 0.02 --size_max 0.25 )
if ( ${Cond} == "D" ) set Extra = ( --shell_width 0.05 --size_min 0.01 --size_max 0.15 )

foreach Init ( flat full )
set Batch_Tag = "Shell_pilot_${Cond}_${Init}"

# w~ = 4.8, 5.6, ... 8.0
foreach Strg ( 1.2 1.4 1.6 1.8 2.0 )

set Unique_Tag = "${Batch_Tag}_Strg_${Strg}_r_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}"

python3 Create_subjob_planar.py --inter_str ${Strg} --radius ${radius} --KA ${KA} --KB ${KB} --tension excess --init ${Init} ${Extra} --batch_tag ${Batch_Tag} --unique_tag ${Unique_Tag}
sbatch ../Subjobs/${Unique_Tag}_subjob

end
end
end
