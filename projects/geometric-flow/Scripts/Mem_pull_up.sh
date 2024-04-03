#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!

set Nsim=100
set Init_cond=3
# foreach Strg ( 5.0 1.0 0.1 20.0 0.05 0.0005 0.0001 0.00001 )
foreach PullF (1000.0 5000.0 10000.0)
foreach Strg ( 0.001 0.01 0.1 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.5 3.0 3.5 4.0 5.0 6.0)
foreach KA ( 500 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads_pull.py ${PullF} ${Strg} ${Init_cond} ${Nsim} ${KA}
sbatch ../Subjobs/subjob_serial_bead_pull_PullF_${PullF}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KA_${KA}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
end 