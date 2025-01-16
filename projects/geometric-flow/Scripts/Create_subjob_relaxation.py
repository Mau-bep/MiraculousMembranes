import sys 
import os
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 

# Lets do a rework here
import json
from jinja2 import Environment, FileSystemLoader

# We need jinja and json here




KA = sys.argv[1]
KB = sys.argv[2]
relaxation_step = sys.argv[3]
# Nsim=sys.argv[4]




def Create_json_relaxation(KA,KB,relaxation_step):
    os.makedirs("../Config_files/",exist_ok = True)

    env = Environment(loader=FileSystemLoader('../Templates/'))
    template = env.get_template('Tube_relaxation.txt')
    # I need to get the xpos ypos zpos

    folderpath = "../Results/Tube_pulling_on_plane/Surface_tension_0.0500_Bending_20.0000_Bead_radius_0.2000_str_10.0000_Bead_radius_0.4000_str_0.0000_Bonds_Lineal_1000.0000_Lineal_1000.0000_Nsim_4/"

    [x,y,z] = get_bead_pos(folderpath,relaxation_step)
    filename = folderpath +"membrane_{}.obj".format(relaxation_step*100)

    output_from_parsed_template = template.render(KA = KA, KB = KB, xpos = x, ypos = y, zpos = z,init_file = filename )
    data = json.loads(output_from_parsed_template)

    Config_path = '../Config_files/Tube_relaxation_step_{}_KA_{}_KB_{}.json'.format(relaxation_step,KA,KB) 
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path




os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)

Config_path = Create_json_relaxation(KA,KB,relaxation_step)



f=open('../Subjobs/subjob_tube_relaxation_{}_KA_{}_KB_{}_Nsim_{}'.format(relaxation_step,KA,KB,relaxation_step),'w')

f.write('#!/bin/bash \n')
f.write('# \n')

f.write('#SBATCH --job-name=Mem3DGpa\n')
f.write('#SBATCH --output=../Outputs/output_tube_relaxation_{}_KA_{}_KB_{}_Nsim_{}'.format(relaxation_step,KA,KB,relaxation_step))
f.write('#\n')
f.write('#number of CPUs to be used\n')
f.write('#SBATCH --ntasks=1\n')
f.write('#Define the number of hours the job should run. \n')
f.write('#Maximum runtime is limited to 10 days, ie. 240 hours\n')
f.write('#SBATCH --time=40:00:00\n')

f.write('#\n')
f.write('#Define the amount of system RAM used by your job in GigaBytes\n')
f.write('#SBATCH --mem=4G\n')
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

f.write('srun time -v ../build/bin/main_cluster {} {}\n'.format(Config_path,relaxation_step))



f.write('\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | head -n 1\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | tail -n 1\n')
f.write('\n')

f.close()
