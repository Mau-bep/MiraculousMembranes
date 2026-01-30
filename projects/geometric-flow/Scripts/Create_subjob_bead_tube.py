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

# v=float(sys.argv[1])
# c0=float(sys.argv[2])
# KA=float(sys.argv[3])
# KB=float(sys.argv[4])
angle = sys.argv[1]
Nsim = 1
strg = 300

def Create_json_wrapping_with_tube(angle, strg):
    theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))


    template = env.get_template('Wrapping_with_tube.txt')
    
   
    
    x1 = 2.0*np.cos(theta) 
    z1 = 2.0*np.sin(theta)
    # We should do 
    
    # Leq = np.sqrt( (x1-x2)**2 + y2**2 )


    
    output_from_parsed_template = template.render( theta = theta,x1 = x1,z1 = z1, strg = strg)

    # print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)


    # print("something\n")
    Config_path = '../Config_files/Wrapping_with_tube_{}_{}.json'.format(angle,strg) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path



os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)



# Config_path, sim_path = Create_json_wrapping_two(KA,KB,radius,Strength,angle)
# # Hopefully this works
Config_path, sim_path = Create_json_wrapping_with_tube(angle,strg)


# def main():
Output_name = 'output_bead_tube_theta_{}_{}.output'.format(angle,strg)
Output_path = '../Outputs/'+Output_name

f=open('../Subjobs/subjob_bead_tube_theta_{}_{}'.format(angle,strg),'w')

f.write('#!/bin/bash \n')
f.write('# \n')

f.write('#SBATCH --job-name=Wrap\n')
f.write('#SBATCH --output={}'.format(Output_path))
f.write('\n#\n')

# f.write('module load boost\n')

f.write('#number of CPUs to be used\n')
f.write('#SBATCH --ntasks=1\n')
f.write('#Define the number of hours the job should run. \n')
f.write('#Maximum runtime is limited to 10 days, ie. 240 hours\n')
f.write('#SBATCH --time=10:01:20\n')

f.write('#\n')
f.write('#Define the amount of system RAM used by your job in GigaBytes\n')
f.write('#SBATCH --mem=5G\n')
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
f.write('#SBATCH --error=%x_%j.err \n')
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
f.write('srun time -v ../build/bin/main_cluster {} {}\n'.format(Config_path,Nsim))
f.write('date\n')
#  Here we can tell the script to move the output file

f.write('cp {}.output {}/{}.txt \n'.format(Output_path,sim_path,Output_name) )
# I need to acces the data in the config file.

f.write('\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | head -n 1\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | tail -n 1\n')
f.write('\n')

f.close()
