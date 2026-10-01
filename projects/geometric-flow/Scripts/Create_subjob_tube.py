import sys 
import os
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 

# Lets do a rework here
import json
from jinja2 import Environment, FileSystemLoader
import numpy as np
# We need jinja and json here





# Nsim=int(sys.argv[1])
# ini_config=int(sys.argv[2])
# fin_config=int(sys.argv[3])
# Target_val=float(sys.argv[4])

KA=sys.argv[1]
KB=sys.argv[2]
radius = sys.argv[3]
finalX = sys.argv[4]
Nsim = 1

# How do i do this, cause there may be more/less inputs 

direction = sys.argv[5]
# Now direction can be +1 or -1 or 0. 
# If direction is 0, then we just do the first call



def Create_json_tube_eq():
    # theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))


    template = env.get_template('Tube_eq.txt')
    
    
    location = [1,"outside","inside"]

    dir = '"../Results/Pulling_and_relaxing_Sept_final/"'.format(KA,KB,radius,finalX)

    ka = float(KA)
    kb = float(KB)
    # I need the force 
    
    output_from_parsed_template = template.render(Dir = dir,KA=KA,KB=KB,r=radius,Xini = -1*float(radius),Xfinal=finalX)

    # print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)


    # print("something\n")
    Config_path = '../Config_files/Tube_to_eq_{0}_{1}_{2}_{3}.json'.format(KA,KB,radius,finalX) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path



def Create_json_tube_continuation(direction):
    # theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))


    template = env.get_template('Tube_eq_continuation.txt')

    if(direction == 1):
        obj = "BeforeSnap.obj"
    else:
        obj = "AfterSnap.obj"
    
    location = [1,"outside","inside"]

    dir = '"../Results/Pulling_and_relaxing_Sept_continuation/"'.format(KA,KB,radius,finalX)
    
    ka = float(KA)
    kb = float(KB)
    # I need the force 
    x_ini = 1.9 
    if(direction == -1):
        x_ini = 2.0


    vx = 0.5*direction
    output_from_parsed_template = template.render(Dir = dir,OBJ = obj,KA=KA,KB=KB,r=radius,Xini = x_ini,Xfinal=finalX,vx = vx)

    # print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)

    # print("something\n")
    Config_path = '../Config_files/Tube_to_eq_{0}_{1}_{2}_{3}_direction_{4}.json'.format(KA,KB,radius,finalX,direction) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path






os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)



# Config_path, sim_path = Create_json_wrapping_two(KA,KB,radius,Strength,angle)
# # Hopefully this works
# Config_path, sim_path = Create_json_wrapping_two_outside(angle,outside1,outside2)

if(int(direction) == 0):
    Config_path, sim_path = Create_json_tube_eq()
    Output_name = 'output_tube_KA_{0}_KB_{1}_r_{2}_XF_{3}.output'.format(KA,KB,radius,finalX)
    f=open('../Subjobs/subjob_tube_KA_{0}_KB_{1}_r_{2}_XF_{3}'.format(KA,KB,radius,finalX),mode='w+')
else:
    Config_path, sim_path = Create_json_tube_continuation(int(direction))
    Output_name = 'output_tube_KA_{0}_KB_{1}_r_{2}_XF_{3}_direction_{4}.output'.format(KA,KB,radius,finalX,direction)
    f=open('../Subjobs/subjob_tube_KA_{0}_KB_{1}_r_{2}_XF_{3}_direction_{4}'.format(KA,KB,radius,finalX,direction),mode='w+')


Output_path = '../Outputs/'+Output_name

f.write('#!/bin/bash \n')
f.write('# \n')

f.write('#SBATCH --job-name=Tube\n')
f.write('#SBATCH --output={}'.format(Output_path))
f.write('\n#\n')

# f.write('module load boost\n')

f.write('#number of CPUs to be used\n')
f.write('#SBATCH --ntasks=1\n')
f.write('#Define the number of hours the job should run. \n')
f.write('#Maximum runtime is limited to 10 days, ie. 240 hours\n')
f.write('#SBATCH --time=6:01:20\n')

f.write('#\n')
f.write('#Define the amount of system RAM used by your job in GigaBytes\n')
f.write('#SBATCH --mem=3G\n')
f.write('#\n')

#f.write('#Send emails when a job starts, it is finished or it exits\n')
#f.write('#SBATCH --mail-user=mrojasve@ist.ac.at\n')
#f.write('#SBATCH --mail-type=ALL\n')
#f.write('#\n')


f.write('#SBATCH --no-requeue\n')
f.write('#\n')


f.write('\n')
f.write('#Do not export the local environment to the compute nodes\n')
f.write('#SBATCH --export=NONE\n')
f.write('\n')
# f.write('#SBATCH --error=%x_%j.err \n')
f.write('unset SLURM_EXPORT_ENV\n')
f.write('#for single-CPU jobs make sure that they use a single thread\n')
f.write('export OMP_NUM_THREADS=1\n')
f.write('#SBATCH --nodes=1\n')
f.write('#SBATCH --cpus-per-task=1\n')

f.write('\n')


# f.write('source /nfs/scistore16/wojtgrp/mrojasve/.bashrc\n')
f.write('export PATH="/nfs/scistore16/wojtgrp/mrojasve/.local/bin:$PATH"\n')
f.write('echo $PATH\n')


f.write('module load conda\n')
f.write('conda activate mir_membranes\n')

f.write('pwd\n')

f.write('date\n')
f.write('srun time -v ../build/bin/main_cluster {}\n'.format(Config_path))
f.write('date\n')
#  Here we can tell the script to move the output file

f.write('cp {} {}/{} \n'.format(Output_path,sim_path,Output_name) )
# I need to acces the data in the config file.

f.write('\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | head -n 1\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | tail -n 1\n')
f.write('\n')

f.close()
