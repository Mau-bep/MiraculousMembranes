#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
# set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
set radius = 1.00
set Batch_Tag = "Wrapping_two_coverage"
# set KA = 5
# set KB = 20
# foreach Strg ( `seq 20 5 140`)
# foreach Strg ( 225.0  )
# foreach KA ( `seq 0 5 15` )
# foreach KA ( 0 10 50 )
# foreach theta ( 0.6 0.625 0.65 0.675 0.7 0.725 0.75 0.775 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 )

# foreach theta ( 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0)
# foreach theta ( 2.5  )
foreach theta ( 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0)
foreach Cov ( 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ) 
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

set Unique_Tag = "${Batch_Tag}_r_${radius}_theta_${theta}_Cov_${Cov}"

python3 Create_subjob_two_beads.py --angle ${theta} --outside1 1 --outside2 1 --radius ${radius} --batch_tag ${Batch_Tag} --unique_tag ${Unique_Tag} --Cov ${Cov}
sbatch ../Subjobs/${Unique_Tag}_subjob


end

end