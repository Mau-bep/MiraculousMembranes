#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
# foreach Strg ( 400.0  )
foreach KA ( 32.0  16.0 )

foreach radius ( 0.1 0.3 )
foreach KB ( 4.0 )

foreach XF ( 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 1.05 1.1 1.15 1.2 1.3 1.4 1.5 1.6 1.8 2.0 )

python3 Create_subjob_tube.py ${KA} ${KB} ${radius} ${XF}
sbatch ../Subjobs/subjob_tube_KA_${KA}_KB_${KB}_r_${radius}_XF_${XF}

end
end
end
end 
