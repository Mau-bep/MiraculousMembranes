#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
# foreach Strg ( 400.0  )
set direction = 1
foreach KA ( 32.0  16.0 )

foreach radius ( 0.3 )
foreach KB ( 4.0 )

foreach XF ( 1.9 1.95 2.0 2.05 2.1 2.15 2.2 2.25)

python3 Create_subjob_tube.py ${KA} ${KB} ${radius} ${XF} ${direction}
sbatch ../Subjobs/subjob_tube_KA_${KA}_KB_${KB}_r_${radius}_XF_${XF}_direction_${direction}

end
end
end
end

set otherdirection = -1
foreach KA ( 32.0  16.0 )

foreach radius ( 0.3 )
foreach KB ( 4.0 )

foreach XF ( 1.7 1.75 1.8 1.85 1.9 1.95 2.0 )

python3 Create_subjob_tube.py ${KA} ${KB} ${radius} ${XF} ${otherdirection}
sbatch ../Subjobs/subjob_tube_KA_${KA}_KB_${KB}_r_${radius}_XF_${XF}_direction_${otherdirection}

end
end
end
end
