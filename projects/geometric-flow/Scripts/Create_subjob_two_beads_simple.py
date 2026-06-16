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
finaldist = sys.argv[1]
# outside1 = int(sys.argv[2])
# outside2 = int(sys.argv[3])
# radius = float(sys.argv[4])
# ka = sys.argv[5]
Nsim = 1




# Strength=sys.argv[2]


# KA = sys.argv[4]
# KB = sys.argv[5]
# KE = 1
# # Init_cond=sys.argv[3]
# Nsim=sys.argv[6]

location = ["unavailable", "outside", "inside"]




def Create_json_wrapping_two(dist):
    # theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))


    template = env.get_template('ContinueTwoBeads.txt')
    V1x = -10.0
    L0 = 0.4
    if(float(dist) > 1.4):
        V1x = 10.0
        # I also need to change the L0
        L0 = 2.4
    FinalX1 = float(dist)/2.0
    # We should do  
    output_from_parsed_template = template.render(Finaldist = dist, Finalx1 = FinalX1, Finalx2 = -FinalX1, V1x = V1x, V2x = -V1x,L0 = L0)

    print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)


    # print("something\n")
    Config_path = '../Config_files/Wrapping_two_{0}_simplestart.json'.format(dist) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path



os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)



# Config_path, sim_path = Create_json_wrapping_two(KA,KB,radius,Strength,angle)
# # Hopefully this works
# Config_path, sim_path = Create_json_wrapping_two_outside(angle,outside1,outside2)

Config_path, sim_path = Create_json_wrapping_two(finaldist)





# # def main():
Output_name = 'output_two_{0}_simplestart.output'.format(finaldist)
Output_path = '../Outputs/'+Output_name

f=open('../Subjobs/subjob_two_{0}_simplestart'.format(finaldist),'w')

f.write('#!/bin/bash \n')
f.write('# \n')

f.write('#SBATCH --job-name=TwoSimple\n')
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
f.write('#SBATCH --mem=2G\n')
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
f.write('srun time -v ../build/bin/main_cluster {} {}\n'.format(Config_path,Nsim))
f.write('date\n')
#  Here we can tell the script to move the output file

f.write('cp {} {}/{} \n'.format(Output_path,sim_path,Output_name) )
# I need to acces the data in the config file.

f.write('\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | head -n 1\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | tail -n 1\n')
f.write('\n')

f.close()
