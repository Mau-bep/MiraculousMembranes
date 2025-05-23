#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
set radius = 0.2
set KA = 10
set KB = 20
foreach Strg ( `seq 400 25 700`)
# foreach Strg (   )
# foreach theta ( 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_fields.py ${Strg} ${KA} ${KB} ${Nsim}
sbatch ../Subjobs/subjob_serial_fields_Strg_${Strg}_KA_${KA}_KB_${KB}_Nsim_${Nsim}


end
end
