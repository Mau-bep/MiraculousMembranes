#!/bin/tcsh
#$ -cwd
#$ -j y



set KB=0.005
set Nsim=2
foreach v ( 0.2 0.3 0.4 0.5 0.6 0.63 0.7 0.8 1.0 )
foreach Init_cond ( 1 2 )

#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_serial.py ${v} ${Init_cond} ${Nsim}
sbatch ../Subjobs/subjob_serial_correct_v_${v}_init_cond_${Init_cond}_Nsim_${Nsim}

#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
