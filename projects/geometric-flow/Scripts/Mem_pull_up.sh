#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!

set Nsim=100
set Init_cond=3
#  4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.5 9.0 10.0
foreach PullF (1000.0 5000.0 10000.0)
foreach Strg (0.1 0.25 0.5 0.75 1.0 1.25 1.5 )
foreach KA ( 500 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads_pull.py ${PullF} ${Strg} ${Init_cond} ${Nsim} ${KA}
sbatch ../Subjobs/subjob_serial_bead_pull_PullF_${PullF}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KA_${KA}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
end 