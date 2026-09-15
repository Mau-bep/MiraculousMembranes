#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
# foreach Strg ( 400.0  )
foreach KA ( 32.0  16.0 )

foreach radius ( 0.3 )
foreach KB ( 4.0 )

# foreach XF ( 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1.0 1.05 1.1 1.15 1.2 1.3 1.4 1.5 1.6 1.8 2.0 )
foreach XF ( 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 )

python3 Create_subjob_tube.py ${KA} ${KB} ${radius} ${XF}
sbatch ../Subjobs/subjob_tube_KA_${KA}_KB_${KB}_r_${radius}_XF_${XF}

end
end
end
end 
