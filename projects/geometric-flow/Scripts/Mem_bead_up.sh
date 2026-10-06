#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
# 16 values, linspace(0, 5, 16)
foreach Strg ( 0 0.3333 0.6667 1 1.333 1.667 2 2.333 2.667 3 3.333 3.667 4 4.333 4.667 5 )
# foreach Strg ( 400.0  )
foreach KA ( 1   )
# foreach KA ( 0..045 50 )

foreach radius ( 1.0 )
# 16 values, logspace from 1/4 to 1/400
foreach KB ( 2.915 2.144 1.577 1.16 0.8536 0.628 0.462 0.3398 )
# foreach KB ( 0.25 0.1839 0.1353 0.09953 0.07322 0.05386 0.03962 0.02915 0.02144 0.01577 0.0116 0.008536 0.00628 0.00462 0.003398 0.0025 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_beads.py ${Strg} ${radius} ${KA} ${KB} ${Nsim}
sbatch ../Subjobs/subjob_WrapPSLog_Strg_${Strg}_r_${radius}_KA_${KA}_KB_${KB}_Nsim_${Nsim}

end
end
end
end 
