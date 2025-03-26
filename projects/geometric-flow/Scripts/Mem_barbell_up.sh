#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
# foreach Strg ( 10 50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000)
# foreach Strg ( 400.0 500.0 600.0 700.0 )
foreach KA (  `seq 40 40 800` )
# foreach radius ( 0.2 0.3 0.4 )

# foreach KB ( `seq 2 2 50` )
foreach KB ( 5 10  )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_barbell.py ${KA} ${KB} ${Nsim}
sbatch ../Subjobs/subjob_serial_barbell_KA_${KA}_KB_${KB}_Nsim_${Nsim}



end
end 
