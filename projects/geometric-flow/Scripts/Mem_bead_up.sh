#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
foreach Strg ( `seq 100 50 600`)
# foreach Strg ( 400.0  )
foreach KA ( 0 5 10 15 20 25 30  )
# foreach KA ( 0..045 50 )

foreach radius ( 0.3 )
foreach KB ( 20.0 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_beads.py ${Strg} ${radius} ${KA} ${KB} ${Nsim}
sbatch ../Subjobs/subjob_BFGS_wrapping_Strg_${Strg}_radius_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}

end
end
end
end 
