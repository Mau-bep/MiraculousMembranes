#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
# set Nsim=1
# set Init_cond=1
# foreach v ( 1.0 )
# foreach Init_cond ( 2 )
foreach theta ( 0.6 0.625 0.65 0.675 0.7 0.725 0.75 0.775 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 )

python3 Create_subjob_two_beads_simple.py ${theta}
sbatch ../Subjobs/subjob_two_${theta}_rescaleArea

end