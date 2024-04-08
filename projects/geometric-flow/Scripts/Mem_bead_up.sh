#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
set Nsim=1
# set Init_cond=1
foreach v (0.0 100.0 1000.0)
foreach Init_cond ( 1 2)
foreach Strg ( 1.0 5.0 10.0 15.0 )
foreach KA ( 100000 )
foreach Minrel (0.001 0.0001)
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_beads.py ${v} ${Strg} ${Init_cond} ${Nsim} ${KA} ${Minrel} 
sbatch ../Subjobs/subjob_serial_bead_v_${v}_minrel_${Minrel}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KA_${KA}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
end
end 
end