import sys 
import os
import argparse
from cli_args import number_str
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 

# Lets do a rework here
import json
from jinja2 import Environment, FileSystemLoader
import numpy as np
# We need jinja and json here






parser = argparse.ArgumentParser(description="Writes the config (if any) and the sbatch subjob for one run.")
parser.add_argument("--inter_str", type=number_str, required=True, help="bead-membrane interaction strength")
parser.add_argument("--KA", type=number_str, required=True, help="surface tension constant KA")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB")
parser.add_argument("--Nsim", type=number_str, required=True, help="simulation number (part of the file names)")
args = parser.parse_args()

Strength = args.inter_str

KA = args.KA
KB = args.KB

Nsim = args.Nsim





def Create_json_fields(ka,kb,inter_str):
    
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Linear_field.txt')
    
    # Radius of the position of the beads is R_v-2*r_b


    output_from_parsed_template = template.render(KA = ka, KB = kb,interaction=inter_str)
    data = json.loads(output_from_parsed_template)
    Config_path = '../Config_files/Linear_field_strg_{}_KA_{}_KB_{}.json'.format(inter_str,ka,kb) 
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path



os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)

Config_path = Create_json_fields(KA,KB,Strength)


f=open('../Subjobs/subjob_serial_fields_Strg_{}_KA_{}_KB_{}_Nsim_{}'.format(Strength,KA,KB,Nsim),'w')

f.write('#!/bin/bash \n')
f.write('# \n')

f.write('#SBATCH --job-name=Mem3DGpa\n')
f.write('#SBATCH --output=../Outputs/output_serial_fields_Strg_{}_KA_{}_KB_{}_Nsim_{}'.format(Strength,KA,KB,Nsim))
f.write('#\n')
f.write('#number of CPUs to be used\n')
f.write('#SBATCH --ntasks=1\n')
f.write('#Define the number of hours the job should run. \n')
f.write('#Maximum runtime is limited to 10 days, ie. 240 hours\n')
f.write('#SBATCH --time=40:00:00\n')

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

f.write('unset SLURM_EXPORT_ENV\n')
f.write('#for single-CPU jobs make sure that they use a single thread\n')
f.write('export OMP_NUM_THREADS=1\n')
f.write('#SBATCH --nodes=1\n')
f.write('#SBATCH --cpus-per-task=1\n')

f.write('\n')


# f.write('source /nfs/scistore16/wojtgrp/mrojasve/.bashrc\n')
f.write('export PATH="/nfs/scistore16/wojtgrp/mrojasve/.local/bin:$PATH"\n')
f.write('echo $PATH\n')


f.write('pwd\n')

f.write('srun time -v ../build/bin/main_cluster {}\n'.format(Config_path))



f.write('\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | head -n 1\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | tail -n 1\n')
f.write('\n')

f.close()
