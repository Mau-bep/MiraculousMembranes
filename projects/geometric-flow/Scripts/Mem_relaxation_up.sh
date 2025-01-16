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
foreach KA ( 0.05)
foreach KB (10.0)
# foreach radius ( 0.2 0.3 0.4 )

foreach relaxation_step ( `seq 0 2 20` )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_relaxation.py ${KA} ${KB} ${relaxation_step}
sbatch ../Subjobs/subjob_tube_relaxation_KA_${KA}_KB_${KB}_Nsim_${relaxation_step}

end

foreach relaxation_step ( 30 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_relaxation.py ${KA} ${KB} ${relaxation_step}
sbatch ../Subjobs/subjob_tube_relaxation_KA_${KA}_KB_${KB}_Nsim_${relaxation_step}

end


foreach relaxation_step ( `seq 180 2 220` )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_relaxation.py ${KA} ${KB} ${relaxation_step}
sbatch ../Subjobs/subjob_tube_relaxation_KA_${KA}_KB_${KB}_Nsim_${relaxation_step}

end


foreach relaxation_step ( `seq 220 5 280` )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_relaxation.py ${KA} ${KB} ${relaxation_step}
sbatch ../Subjobs/subjob_tube_relaxation_KA_${KA}_KB_${KB}_Nsim_${relaxation_step}

end



end
end 
