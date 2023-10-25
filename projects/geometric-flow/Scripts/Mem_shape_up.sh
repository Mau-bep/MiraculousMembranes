#!/bin/tcsh
#$ -cwd
#$ -j y


set KA=10.0
set c0=0.0
foreach KB ( 0.005)
#foreach v ( 0.4 )
foreach v ( 0.25 )

#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_serial.py ${v} ${c0} ${KA} ${KB}
sbatch subjob_serial_correct_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
