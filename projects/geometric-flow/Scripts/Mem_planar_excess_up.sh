#!/bin/tcsh
#$ -cwd
#$ -j y


# Same scan as Mem_planar_up.sh, with Excess_tension instead of Surface_tension:
# KA (A - A0) only while the area A is above the initial area A0
set Nsim=1
set Batch_Tag = "Wrapping_planar_excess"
set KB = 1
set radius = 1.0

# 16 values, linspace(1.5, 3, 16)
foreach Strg ( 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 )
# 16 values, linspace(0, 1, 16)
foreach KA ( 0 0.06667 0.1333 0.2 0.2667 0.3333 0.4 0.4667 0.5333 0.6 0.6667 0.7333 0.8 0.8667 0.9333 1 )

set Unique_Tag = "${Batch_Tag}_Strg_${Strg}_r_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}"

python3 Create_subjob_planar.py --inter_str ${Strg} --radius ${radius} --KA ${KA} --KB ${KB} --tension excess --batch_tag ${Batch_Tag} --unique_tag ${Unique_Tag}
sbatch ../Subjobs/${Unique_Tag}_subjob

end
end
