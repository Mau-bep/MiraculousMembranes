#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
# set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
# set radius = 0.3
# set KA = 4
# set KB = 20
# foreach Strg ( `seq 20 5 140`)
# foreach Strg ( 225.0  )
foreach theta ( 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.7 1.9 2.0 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

# python3 Create_subjob_two_beads.py ${theta} ${Strg} ${radius} ${KA} ${KB} ${Nsim}
# sbatch ../Subjobs/subjob_serial_two_beads_theta_${theta}_Strg_${Strg}_radius_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}

# python3 Create_subjob_two_beads.py ${theta} -1 -1 
# sbatch ../Subjobs/subjob_two_beads_without_spring_theta_${theta}_inside_inside_finer

# python3 Create_subjob_two_beads.py ${theta} -1 1 
# sbatch ../Subjobs/subjob_two_beads_without_spring_theta_${theta}_inside_outside_finer

python3 Create_subjob_two_beads.py ${theta} 1 1 
sbatch ../Subjobs/subjob_two_bead_theta_${theta}_outside_outside_BFGS_M2


end
