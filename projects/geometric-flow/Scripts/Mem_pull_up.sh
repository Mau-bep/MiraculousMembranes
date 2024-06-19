#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!

set Nsim=1
set Init_cond=1


#  4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.5 9.0 10.0
# foreach Init_cond (1 3)
foreach rc (2.0)
foreach Strg ( 10.0 )
foreach KB ( 0.1 0.5 1.0 2.0 4.0 6.0 8.0 10.0 14.0 20.0 )
foreach KA ( 100000 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads_pull.py ${rc} ${Strg} ${Init_cond} ${Nsim} ${KB} ${KA}
sbatch ../Subjobs/subjob_serial_bead_pull_rc_${rc}_KA_${KA}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KB_${KB}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
end 
end


foreach rc (2.0)
foreach Strg ( 10.0 )
foreach KB ( 1.0  )
foreach KA ( 100 500 1000 5000 10000 50000)
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads_pull.py ${rc} ${Strg} ${Init_cond} ${Nsim} ${KB} ${KA}
sbatch ../Subjobs/subjob_serial_bead_pull_rc_${rc}_KA_${KA}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KB_${KB}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
end 
end