#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
# foreach Strg ( `seq 20 5 140`)
# foreach Strg ( 400.0 500.0 600.0 700.0 )
foreach theta ( 0.3 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_two_beads.py ${theta} ${Strg} ${radius} ${KA} ${KB} ${Nsim}
sbatch ../Subjobs/subjob_serial_two_beads_theta_${theta}_Strg_${Strg}_radius_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}

end
