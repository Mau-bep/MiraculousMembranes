import sys 
import os
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 

# Lets do a rework here
import json
from jinja2 import Environment, FileSystemLoader

# We need jinja and json here





# Nsim=int(sys.argv[1])
# ini_config=int(sys.argv[2])
# fin_config=int(sys.argv[3])
# Target_val=float(sys.argv[4])

# v=float(sys.argv[1])
# c0=float(sys.argv[2])
# KA=float(sys.argv[3])
# KB=float(sys.argv[4])
Strength=sys.argv[1]
radius = float(sys.argv[2])

KA = sys.argv[3]
KB = sys.argv[4]

# Init_cond=sys.argv[3]
Nsim=sys.argv[5]






def Create_json_wrapping(ka,kb,r,inter_str):

    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Wrapping.txt')
    output_from_parsed_template = template.render(KA = ka, KB = kb,radius = r,xpos = r*1.95 ,interaction=inter_str)

    data = json.loads(output_from_parsed_template)

    Config_path = '../Config_files/Wrapping_strg_{}_radius_{}_KA_{}_KB_{}.json'.format(inter_str,r,ka,kb) 
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path




os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)

Config_path = Create_json_wrapping(KA,KB,radius,Strength)


f=open('../Subjobs/subjob_serial_wrapping_Strg_{}_radius_{}_KA_{}_KB_{}_Nsim_{}'.format(Strength,radius,KA,KB,Nsim),'w')

f.write('#!/bin/bash \n')
f.write('# \n')

f.write('#SBATCH --job-name=Mem3DGpa\n')
f.write('#SBATCH --output=../Outputs/output_serial_wrapping_Strg_{}_radius_{}_KA_{}_KB_{}_Nsim_{}'.format(Strength,radius,KA,KB,Nsim))
f.write('#\n')
f.write('#number of CPUs to be used\n')
f.write('#SBATCH --ntasks=1\n')
f.write('#Define the number of hours the job should run. \n')
f.write('#Maximum runtime is limited to 10 days, ie. 240 hours\n')
f.write('#SBATCH --time=40:00:00\n')

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

f.write('srun time -v ../build/bin/main_cluster {} {}\n'.format(Config_path,Nsim))



f.write('\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | head -n 1\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | tail -n 1\n')
f.write('\n')

f.close()
