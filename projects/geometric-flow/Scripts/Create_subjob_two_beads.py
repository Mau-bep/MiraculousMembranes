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
outside1 = int(sys.argv[2])
outside2 = int(sys.argv[3])
Nsim = 1

# Strength=sys.argv[2]
# radius = float(sys.argv[3])

# KA = sys.argv[4]
# KB = sys.argv[5]
# KE = 1
# # Init_cond=sys.argv[3]
# Nsim=sys.argv[6]

location = ["unavailable", "outside", "inside"]




def Create_json_wrapping_two(ka,kb,r,inter_str,angle):
    theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Two_beads.txt')
    
    # Radius of the position of the beads is R_v-2*r_b
    R_vesicle = 3.0
    r_bead = 0.3
    Rpos_beads = R_vesicle-r_bead*2
    xpos  = Rpos_beads*np.cos(theta)
    ypos1 = Rpos_beads*np.sin(theta)
    ypos2 = -Rpos_beads*np.sin(theta)


    disp = -1*Rpos_beads*np.cos(theta)+ 0.5*(np.sqrt((R_vesicle+r_bead)*(R_vesicle+r_bead)-Rpos_beads*Rpos_beads*np.sin(theta)*np.sin(theta) ) + np.sqrt(R_vesicle*R_vesicle-Rpos_beads*Rpos_beads*np.sin(theta)*np.sin(theta))) 
    d1 = -1*Rpos_beads*np.cos(theta)+ np.sqrt((R_vesicle+r_bead)*(R_vesicle+r_bead)-Rpos_beads*Rpos_beads*np.sin(theta)*np.sin(theta) )
    d2 = -1*Rpos_beads*np.cos(theta)+ np.sqrt(R_vesicle*R_vesicle-Rpos_beads*Rpos_beads*np.sin(theta)*np.sin(theta))
    disp = -1*0.5*(d1+d2)
    print("d1 is {}, d2 is {}".format(d1,d2))
    
    print("We should have then that this is  0 ? {}".format( Rpos_beads*Rpos_beads+d1*d1+2*d1*Rpos_beads*np.cos(theta) - (R_vesicle+r_bead)*(R_vesicle+r_bead)  ))
    print("We should have then that this is  0 ? {}".format( Rpos_beads*Rpos_beads+d2*d2+2*d2*Rpos_beads*np.cos(theta) - R_vesicle*R_vesicle  ))

    print("disp is {}".format(disp))
    
    print("D1 gets me a distance of {}".format( np.sqrt( (d1+Rpos_beads*np.cos(theta) )*(d1+Rpos_beads*np.cos(theta) )+ Rpos_beads*np.sin(theta)*Rpos_beads*np.sin(theta)  ) ))
    print("D2 gets me a distance of {}".format( np.sqrt( (d2+Rpos_beads*np.cos(theta) )*(d2+Rpos_beads*np.cos(theta) )+ Rpos_beads*np.sin(theta)*Rpos_beads*np.sin(theta)  ) ))
    print("disp gets me a distance of {}".format( np.sqrt( (disp+Rpos_beads*np.cos(theta) )*(disp+Rpos_beads*np.cos(theta) )+ Rpos_beads*np.sin(theta)*Rpos_beads*np.sin(theta)  ) ))


    disp= 0.0
    xpos = (R_vesicle+r_bead)*np.cos(theta)
    ypos1 = (R_vesicle+r_bead)*np.sin(theta)
    ypos2 = -(R_vesicle+r_bead)*np.sin(theta)
    output_from_parsed_template = template.render(KA = ka, KB = kb,radius = r,xdisp = disp,xpos1 = xpos,xpos2 =xpos, ypos1= ypos1, ypos2 = ypos2, vx1= -5*xpos, vy1 = -5*ypos1, vx2 = -5*xpos, vy2 = -5*ypos2 ,interaction=inter_str, theta = theta,KE  = KE)
    data = json.loads(output_from_parsed_template)
    Config_path = '../Config_files/Wrapping_two_{}_strg_{}_radius_{}_KA_{}_KB_{}.json'.format(angle,inter_str,r,ka,kb) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path


def Create_json_wrapping_two_outside(angle, outside1, outside2):
    theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))


    template = env.get_template('Wrapping_two_spring.txt')
    
    # Radius of the position of the beads is R_v-2*r_b
    R_vesicle = 2.0
    r_bead = 0.2
    
    location = [1,"outside","inside"]

    dir = '"../Results/Two_beads_without_spring_{}_{}/"'.format(location[outside1],location[outside2])

    v1x = 10.0*(outside1*-1)
    x1 = 2.0 + 0.3*outside1 

    r2 = 2.0 + 0.3*outside2
    # We should do 
    x2 = r2*np.cos(theta)
    y2 = r2*np.sin(theta)
    Leq = np.sqrt( (x1-x2)**2 + y2**2 )

    v2x = 10*np.cos(theta)*outside2*-1
    v2y = 10*np.sin(theta)*outside2*-1

    
    output_from_parsed_template = template.render(Dir = dir,theta =theta, outside1 = outside1, v1x = v1x, x1 = x1,L0 = Leq, outside2 = outside2, v2x = v2x, v2y = v2y,x2 = x2, y2 = y2 )

    # print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)


    # print("something\n")
    Config_path = '../Config_files/Wrapping_two_without_{}_{}_{}.json'.format(angle,location[outside1],location[outside2]) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path




# os.makedirs('../Subjobs/',exist_ok=True)
# os.makedirs('../Outputs/',exist_ok=True)



# Config_path, sim_path = Create_json_wrapping_two(KA,KB,radius,Strength,angle)
# # Hopefully this works
Config_path, sim_path = Create_json_wrapping_two_outside(angle,outside1,outside2)






# # def main():
Output_name = 'output_two_beads_without_spring_theta_{}_{}_{}.output'.format(angle,location[outside1],location[outside2])
Output_path = '../Outputs/'+Output_name

f=open('../Subjobs/subjob_serial_two_beads_without_spring_theta_{}_{}_{}'.format(angle,location[outside1],location[outside2]),'w')

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
