#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!
# set v=1.0
set Nsim=2
# set Init_cond=1
foreach v (0.0 )
foreach Init_cond ( 2)
foreach Strg ( 10.0 50.0 )
foreach KA ( 100000 )
foreach radius ( 1.0 )
foreach KB ( 1.0 )
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}

python3 Create_subjob_beads.py ${v} ${Strg} ${Init_cond} ${Nsim} ${KA} ${radius} ${KB}
sbatch ../Subjobs/subjob_serial_bead_curvadap_${v}_radius_${radius}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KA_${KA}_KB_${KB}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
# end
end
end
end
end 
end
end 