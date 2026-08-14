#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
# foreach Strg ( 400.0  )
foreach KA ( 32.0  16.0 )

foreach radius ( 0.1 0.3 )
foreach KB ( 1.0  4.0 )

foreach XF ( 0.0 0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 )

python3 Create_subjob_tube.py ${KA} ${KB} ${radius} ${XF}
sbatch ../Subjobs/subjob_tube_KA_${KA}_KB_${KB}_r_{radius}_XF_{XF}

end
end
end
end 
