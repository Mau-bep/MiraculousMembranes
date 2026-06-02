#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
# set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
set radius = 0.25
# set KA = 5
# set KB = 20
# foreach Strg ( `seq 20 5 140`)
# foreach Strg ( 225.0  )
# foreach KA ( `seq 0 5 15` )
foreach KA ( 0 10 50 )
foreach theta ( 0.6 0.625 0.65 0.675 0.7 0.725 0.75 0.775 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 )
# foreach theta ( 2.6  )

#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

# python3 Create_subjob_two_beads.py ${theta} ${Strg} ${radius} ${KA} ${KB} ${Nsim}
# sbatch ../Subjobs/subjob_serial_two_beads_theta_${theta}_Strg_${Strg}_radius_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}

# python3 Create_subjob_two_beads.py ${theta} -1 -1 ${radius}
# sbatch ../Subjobs/subjob_two_bead_r_${radius}_theta_${theta}_inside_inside_BFGS_Fixed_3

# python3 Create_subjob_two_beads.py ${theta} -1 1 ${radius}
# sbatch ../Subjobs/subjob_two_bead_r_${radius}_theta_${theta}_inside_outside_BFGS_Fixed_3

python3 Create_subjob_two_beads.py ${theta} 1 1 ${radius} ${KA}
sbatch ../Subjobs/subjob_two_bead_r_${radius}_theta_${theta}_outside_outside_BFGS_ST_${KA}_May


end

end