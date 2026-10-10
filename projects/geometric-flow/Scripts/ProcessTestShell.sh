#!/bin/tcsh
#$ -cwd
#$ -j y


set Nsim=1
set KB = 1
set radius = 1.0
set KA = 0

set Rescale = 1.0
grep -q "^FLAT_MESH = '../../../input/Planar_mem.obj'" Create_subjob_planar.py
if ( $status == 0 ) set Rescale = 5.0
echo "rescale ${Rescale}"

foreach Pw ( 1 2 4 )
foreach Sw ( 0.25 0.1 0.05 )

# refine_angle (edge length near the bead), BFGS -> BFGS-Normal step and slurm time limit for this shell
if ( ${Sw} == "0.25" ) then
set Refine = 0.15
set Switch = 60000
set Time = 12:00:00
endif
if ( ${Sw} == "0.1" ) then
if ( ${Pw} == "1" ) set Refine = 0.15
if ( ${Pw} == "2" ) set Refine = 0.14
if ( ${Pw} == "4" ) set Refine = 0.12
set Switch = 100000
set Time = 24:00:00
endif
if ( ${Sw} == "0.05" ) then
if ( ${Pw} == "1" ) set Refine = 0.11
if ( ${Pw} == "2" ) set Refine = 0.10
if ( ${Pw} == "4" ) set Refine = 0.08
set Switch = 120000
set Time = 48:00:00
endif

set Batch_Tag = "Shell_test_p${Pw}_s${Sw}"

bin/PostProcessing ../Config_files/${Batch_Tag}_Strg_0.5_r_1.0_KA_0_KB_1_Nsim_1_ConfigFile.json


end
end