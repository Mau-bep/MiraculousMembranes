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
Strength=sys.argv[2]
radius = float(sys.argv[3])

KA = sys.argv[4]
KB = sys.argv[5]

# Init_cond=sys.argv[3]
Nsim=sys.argv[6]





def Create_json_wrapping_two(ka,kb,r,inter_str,angle):
    theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Two_beads.txt')
    
    # Radius of the position of the beads is R_v-2*r_b
    Rpos_beads = 1.5-0.2*2
    xpos  = Rpos_beads*np.cos(theta)
    ypos1 = Rpos_beads*np.sin(theta)
    ypos2 = -Rpos_beads*np.sin(theta)


    disp = -1*Rpos_beads*np.cos(theta)+ 0.5*(np.sqrt(1.6*1.6-Rpos_beads*Rpos_beads*np.sin(theta)**2 ) + np.sqrt(1.5*1.5-Rpos_beads*Rpos_beads*np.sin(theta)**2)) 

    output_from_parsed_template = template.render(KA = ka, KB = kb,radius = r,xdisp = disp,xpos1 = xpos,xpos2 =xpos, ypos1= ypos1, ypos2 = ypos2 ,interaction=inter_str, theta = theta)
    data = json.loads(output_from_parsed_template)
    Config_path = '../Config_files/Wrapping_two_{}_strg_{}_radius_{}_KA_{}_KB_{}.json'.format(angle,inter_str,r,ka,kb) 
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path



os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)

Config_path = Create_json_wrapping_two(KA,KB,radius,Strength,angle)


f=open('../Subjobs/subjob_serial_two_beads_theta_{}_Strg_{}_radius_{}_KA_{}_KB_{}_Nsim_{}'.format(angle,Strength,radius,KA,KB,Nsim),'w')

f.write('#!/bin/bash \n')
f.write('# \n')

f.write('#SBATCH --job-name=Mem3DGpa\n')
f.write('#SBATCH --output=../Outputs/output_serial_two_beads_theta_{}_Strg_{}_radius_{}_KA_{}_KB_{}_Nsim_{}'.format(angle,Strength,radius,KA,KB,Nsim))
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

f.write('srun time -v ../build/bin/main_cluster {} {}\n'.format(Config_path,Nsim))



f.write('\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | head -n 1\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | tail -n 1\n')
f.write('\n')

f.close()
