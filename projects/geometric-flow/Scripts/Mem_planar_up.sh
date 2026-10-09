#!/bin/tcsh
#$ -cwd
#$ -j y


set Nsim=1
# "full": start from input/Planar_full.obj (membrane already wrapped around the bead),
# "flat": start from Big_planar_mem.obj with the bead above the plane
set Init = "full"
set Batch_Tag = "Wrapping_planar"
if ( ${Init} == "full" ) set Batch_Tag = "${Batch_Tag}_full"
set KB = 1
set radius = 1.0

# 16 values, linspace(1.5, 3, 16) not anymore sorry 
foreach Strg ( 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 )
# 16 values, linspace(0, 1, 16)
foreach KA ( 0 0.06667 0.1333 0.2 0.2667 0.3333 0.4 0.4667 0.5333 0.6 0.6667 0.7333 0.8 0.8667 0.9333 1 )

set Unique_Tag = "${Batch_Tag}_Strg_${Strg}_r_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}"

python3 Create_subjob_planar.py --inter_str ${Strg} --radius ${radius} --KA ${KA} --KB ${KB} --init ${Init} --batch_tag ${Batch_Tag} --unique_tag ${Unique_Tag}
sbatch ../Subjobs/${Unique_Tag}_subjob

end
end
