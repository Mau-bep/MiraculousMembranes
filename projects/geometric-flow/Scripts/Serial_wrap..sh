#!/bin/tcsh
#$ -cwd
#$ -j y




# foreach angle ( 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.4 1.5 1.6 )
foreach angle ( 0.2 ) 
    python3 Create_subjob_two_beads.py ${angle} 1 1
    ../build/main_cluster ../Config_files/Wrapping_two_${angle}_outside_outside.json 1
    # sbatch ../Subjobs/subjob_serial_two_beads_theta_${angle}_inside_outside
    # sbatch ../Subjobs/subjob_serial_two_beads_theta_${angle}_outside_outside
end



