#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
set Nsim=1
# set Init_cond=1
foreach v (0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0)
foreach Init_cond ( 1 3)
foreach Strg ( 4.0)
foreach KA ( 500)
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads.py ${v} ${Strg} ${Init_cond} ${Nsim} ${KA}
sbatch ../Subjobs/subjob_serial_bead_v_${v}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KA_${KA}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
end 
end