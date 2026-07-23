#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
foreach Strg ( 1 1.5 2 2.5 3 3.5 4 4.5 5  )
# foreach Strg ( 400.0  )
foreach KA ( 1   )
# foreach KA ( 0..045 50 )

foreach radius ( 1.0 )
foreach KB ( 0.005 0.01 0.05 0.1 0.5 1.0   )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_beads.py ${Strg} ${radius} ${KA} ${KB} ${Nsim}
sbatch ../Subjobs/subjob_LBFGS_wrapping_Strg_${Strg}_r_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}

end
end
end
end 
