#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!

set Nsim=1
set Init_cond=2


#  4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.5 9.0 10.0
# foreach Init_cond (1 3)
foreach radius (0.2)
foreach Strg ( 0.1 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5 5.75 6.0 )
foreach KB ( 1.0 )
foreach KA ( 100000 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads_pull.py ${radius} ${Strg} ${Init_cond} ${Nsim} ${KB} ${KA}
sbatch ../Subjobs/subjob_serial_bead_pull_radius_${radius}_KA_${KA}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KB_${KB}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
end 
end


set Nsim2=100
foreach radius (0.2)
foreach Strg ( 10.0 )
foreach KB ( 1.0 2.0 5.0 10.0 15.0 20.0 25.0  )
foreach KA ( 100000)
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads_pull.py ${radius} ${Strg} ${Init_cond} ${Nsim2} ${KB} ${KA}
sbatch ../Subjobs/subjob_serial_bead_pull_radius_${radius}_KA_${KA}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim2}_KB_${KB}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
end 
end