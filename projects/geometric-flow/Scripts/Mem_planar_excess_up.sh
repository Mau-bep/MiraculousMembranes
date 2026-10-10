#!/bin/tcsh
#$ -cwd
#$ -j y


# Excess_tension instead of Surface_tension: KA (A - A0) only while the area A is above
# the initial area A0 (the flat disk's, for both starts).
# Dense scan of the window of Fig. 2 of Deserno (arXiv:cond-mat/0303656), with KB = r = 1:
#   w~ = 4 KI r^2/KB in [3, 6]      -> KI in [0.75, 1.5]  (steps of 0.05, w~ steps of 0.2)
#   sigma~ = 2 KA r^2/KB in [0, 1]  -> KA in [0, 0.5]
# and further to the right, up to w~ = 10 (KI = 2.5, steps of 0.1): the soft adhesion shell of the
# simulation delays the wrapping transitions (binding jump at w~ about 5 for sigma~ = 0.13, about 9 for
# sigma~ = 1 in the first scans, Deserno's S1 is at 4.6 and 7.7), so the comparison needs that range.
# Both starts run on the same grid (the lower energy of the two is the equilibrium one):
#   flat: Big_planar_mem.obj with the bead above the plane
#   full: input/Planar_full.obj, the membrane already wrapped around the bead
# The results go to ../Results/<Base_Tag> (flat) and ../Results/<Base_Tag>_full
set Nsim=1
set Base_Tag = "Wrapping_planar_excess_thinShell_p4"
set KB = 1
set radius = 1.0

# Optional overrides of the template defaults, for example a narrower adhesion shell on a finer mesh:
#   set Extra = ( --shell_width 0.1 --size_max 0.25 )
# set Extra = ( --shell_width 0.05 )
set Extra = ( --shell_width 0.05 --shell_power 4 --rescale 4.0 --size_min 0.005 --size_max 0.25--switch_normal 30000 )

# foreach Init ( flat full )
foreach Init ( flat )

set Batch_Tag = "${Base_Tag}"
if ( ${Init} == "full" ) set Batch_Tag = "${Batch_Tag}_full"

# 26 values: 0.75 ... 1.5 in steps of 0.05 (the window of the paper), 1.6 ... 2.5 in steps of 0.1
foreach Strg ( 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 )
# 16 values, linspace(0, 0.5, 16)
foreach KA ( 0 0.03333 0.06667 0.1 0.1333 0.1667 0.2 0.2333 0.2667 0.3 0.3333 0.3667 0.4 0.4333 0.4667 0.5 )

set Unique_Tag = "${Batch_Tag}_Strg_${Strg}_r_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}"

python3 Create_subjob_planar.py --inter_str ${Strg} --radius ${radius} --KA ${KA} --KB ${KB} --tension excess --init ${Init} ${Extra} --batch_tag ${Batch_Tag} --unique_tag ${Unique_Tag}
sbatch ../Subjobs/${Unique_Tag}_subjob

end
end
end
