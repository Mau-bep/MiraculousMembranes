#!/bin/tcsh
#$ -cwd
#$ -j y


set KA=10.0
set c0=0.0
foreach c0 (-12.0 12.0)
foreach KB ( 0.0 )
# foreach v ( 0.6 0.65 0.7 0.8 1.0 )
foreach v ( 1.0 )

python3 Create_subjob_test.py ${v} ${c0} ${KA} ${KB}

sbatch ../Subjobs/subjob_test_th_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}

end
end
end