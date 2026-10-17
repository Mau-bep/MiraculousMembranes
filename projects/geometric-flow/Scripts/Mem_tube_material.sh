#!/bin/tcsh
#$ -cwd
#$ -j y

# Material properties of the membrane from the tube a bead pulls out of the big planar membrane: ONE chain per (KA, KB). The START STATE is made the old way (Mem_tube_up.sh): a bead
# of radius 0.3 is pulled by a spring (the Lineal bond of Tube_eq.txt, force 200) through the flat disk of input/Big_planar_mem.obj and freezes at x = 0.2 (start_x0.2/ in the results
# folder, about 250-600 steps; the force is linear in x below 0.2 and needs no chain). From that mesh the bead is FROZEN at x = 0.2, 0.3, ... 6.0 (every position relaxed from the
# previous one, in waypoints of 0.025), then walked back 5.9, 5.8, ... 0.2, every run saved (../Results/Tube_material/KA_<KA>_KB_<KB>_r_<radius>/linkNNN_x<x>/1/, chain_summary.csv with
# the force on the bead and the distance of the bead to the mesh, chain.log). The way back shows the hysteresis (the tube that holds on, snaps, retracts) and Mem_tube_up's force-extension
# (Results/Pulling_and_relaxing_2/Force_data.csv: 15.6 at x = 0.2 for KA 16, 27.6 for KA 32, a plateau 50.1 and 70.4 = 99 percent of 2 pi sqrt(KB KA) beyond x = 3) is the check.
# Physics: Surface_tension KA, the regular Bending KB (kappa = KB/2), shifted LJ of Mem_tube_up (shift = strength = 10), the remesher of Tube_eq.txt (size_max 0.3, refine_angle 0.5).
# Every relaxation is checked (no nan line search, the energy did not rise, the bead is not through a face of the mesh); one that fails is retried through the middle position, and the
# chain stops with a message otherwise. A chain that stopped, hit the time limit or failed resumes after its last finished run when its subjob is sbatched again.
# After the jobs ran, on your machine (matplotlib: micromamba run -n mir_mem):
#   python3 Analyze_tube_material.py --out ../Studies/tube_material/tube a=../Results/Tube_material/KA_32.0_KB_4.0_r_0.3 b=../Results/Tube_material/KA_16.0_KB_4.0_r_0.3
# Options (range, step, waypoints, radius, mesh, own start state with --start_run, ...): python3 Create_subjob_tube_material.py --help ; ../Studies/tube_material/README.md

set Batch_Tag = "Tube_material"
set radius = 0.3
set XS = 0.2
set XF = 6.0
set Step = 0.05

foreach KA ( 32.0 24 16.0 )
foreach KB ( 4.0 )

set Tag = "KA_${KA}_KB_${KB}_r_${radius}"
python3 Create_subjob_tube_material.py --KA ${KA} --KB ${KB} --radius ${radius} --x_start ${XS} --x_end ${XF} --step ${Step} --batch_tag ${Batch_Tag}/${Tag} --unique_tag ${Batch_Tag}_${Tag} --x_back 0.0
if ( $status != 0 ) exit 1
sbatch ../Subjobs/${Batch_Tag}_${Tag}_subjob

end
end
