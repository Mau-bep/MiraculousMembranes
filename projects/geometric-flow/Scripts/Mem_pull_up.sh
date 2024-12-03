#!/bin/tcsh
#$ -cwd
#$ -j y


# #   SARIConGPT143!

set Nsim=1
set Init_cond=2



# #  4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.5 9.0 10.0
# # foreach Init_cond (1 3)
# foreach Nsim (0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 45 50 55 60 65 70)
# foreach radius (0.2)
# foreach Strg ( 0.0)
# foreach KB ( 1.0 )
# foreach KA ( 100000 )
# #python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
# python3 Create_subjob_beads_pull.py ${radius} ${Strg} ${Init_cond} ${Nsim} ${KB} ${KA}
# sbatch ../Subjobs/subjob_serial_bead_pull_radius_${radius}_KA_${KA}_Strg_${Strg}_init_cond_${Init_cond}_Nsim_${Nsim}_KB_${KB}
# #sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
# end
# end
# end 
# end

set Init_cond2=6

set Nsim2=1
foreach radius (0.2)
foreach Strg ( 100.0 )
foreach KB ( 2.0 4.0 6.0 8.0 10.0 12.0 14.0 16.0 18.0 20.0 22.0 24.0 26.0 28.0 30.0 32.0 34.0 36.0 38.0 40.0 42.0 44.0 46.0 48.0 50.0 )
foreach KA ( 0.05)
#python3 Create_subjob.py ${v} ${c0} ${KA} ${KB}
python3 Create_subjob_beads_pull.py ${radius} ${Strg} ${Init_cond2} ${Nsim2} ${KB} ${KA}
sbatch ../Subjobs/subjob_serial_bead_pull_radius_${radius}_KA_${KA}_Strg_${Strg}_init_cond_${Init_cond2}_Nsim_${Nsim2}_KB_${KB}
#sbatch subjob_parallel_memshape_v_${v}_c0_${c0}_KA_${KA}_KB_${KB}
end
end
end 
end