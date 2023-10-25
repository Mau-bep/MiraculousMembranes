#!/bin/tcsh
#$ -cwd
#$ -j y


module load python/3.9.7
EXPORT PATH=/nfs/scistore16/wojtgrp/mrojasve/.local/bin:$PATH

#foreach v ( 0.1 0.3 0.58 0.6 0.65 0.75 0.8 1.0 )
foreach v ( 0.4 0.62 0.8 1.0 )
foreach KB ( 0.01 0.005)

python3 Analyze_test_cluster.py ${v} ${KB}
# python3 Seing_meshes.py

end
end
