#!/bin/tcsh
#$ -cwd
#$ -j y


set KA=10.0
set c0=0.0

set KB=0.005

set v=1.0
foreach Strg ( 5.0 1.0 0.1 20.0 0.05 0.0005 0.0001 0.00001 )
# foreach Strg ( 0.1 20.0  )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads.py ${v} ${c0} ${KA} ${KB} ${Strg}
sbatch ../Subjobs/subjob_serial_bead_frenkel_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}_Strg_${Strg}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
