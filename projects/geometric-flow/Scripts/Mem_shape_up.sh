#!/bin/tcsh
#$ -cwd
#$ -j y


# set KA=10.0
# set c0=0.0
set KB=0.005
# foreach KB ( 0.005)
foreach v ( 0.4 )
foreach Init_cond ( 1 2 )

#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_serial.py ${v} ${Init_cond} ${Nsim}
sbatch ../Subjobs/subjob_serial_correct_v_${v}_init_cond_${Init_cond}_Nsim_${Nsim}

#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
