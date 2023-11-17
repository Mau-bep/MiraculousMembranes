#!/bin/tcsh
#$ -cwd
#$ -j y



set v=1.0
set Nsim=1
set Init_cond=1
# foreach Strg ( 5.0 1.0 0.1 20.0 0.05 0.0005 0.0001 0.00001 )
foreach Strg ( 0.1 1.0 0.01 0.001 10.0 )

#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads.py ${v} ${Strg} ${Init_cond} ${Nsim}
sbatch ../Subjobs/subjob_serial_bead_v_${v}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
