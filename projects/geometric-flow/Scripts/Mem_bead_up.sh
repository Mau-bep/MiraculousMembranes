#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
set v=1.0
set Nsim=15
set Init_cond=100
# foreach Strg ( 5.0 1.0 0.1 20.0 0.05 0.0005 0.0001 0.00001 )
foreach Strg ( 0.001 0.01 0.1 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.5 3.0 3.5 4.0)
foreach KA (10 100 200 500 1000)
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads.py ${v} ${Strg} ${Init_cond} ${Nsim} ${KA}
sbatch ../Subjobs/subjob_serial_bead_v_${v}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KA_${KA}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end