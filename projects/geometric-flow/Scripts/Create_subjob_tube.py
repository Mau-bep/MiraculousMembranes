import sys 
import os
import argparse
from cli_args import number_str
from subjob import write_subjob
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 

# Lets do a rework here
import json
from jinja2 import Environment, FileSystemLoader
import numpy as np
# We need jinja and json here






parser = argparse.ArgumentParser(description="Writes the config (if any) and the sbatch subjob for one run.")
parser.add_argument("--KA", type=number_str, required=True, help="surface tension constant KA")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB")
parser.add_argument("--radius", type=number_str, required=True, help="bead radius")
parser.add_argument("--finalX", type=number_str, required=True, help="final x of the tube end")
parser.add_argument("--direction", type=int, required=True, help="pulling direction: 0 (first call), 1 or -1")
args = parser.parse_args()

KA = args.KA
KB = args.KB
radius = args.radius
finalX = args.finalX
Nsim = 1

# How do i do this, cause there may be more/less inputs 

direction = args.direction
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


if(int(direction) == 0):
    Config_path, sim_path = Create_json_tube_eq()
    Subjob_name = 'subjob_tube_KA_{0}_KB_{1}_r_{2}_XF_{3}'.format(KA,KB,radius,finalX)
    Output_name = 'output_tube_KA_{0}_KB_{1}_r_{2}_XF_{3}.output'.format(KA,KB,radius,finalX)
else:
    Config_path, sim_path = Create_json_tube_continuation(int(direction))
    Subjob_name = 'subjob_tube_KA_{0}_KB_{1}_r_{2}_XF_{3}_direction_{4}'.format(KA,KB,radius,finalX,direction)
    Output_name = 'output_tube_KA_{0}_KB_{1}_r_{2}_XF_{3}_direction_{4}.output'.format(KA,KB,radius,finalX,direction)

write_subjob(Subjob_name, Output_name,
             '../build/bin/main_cluster {}'.format(Config_path),
             job_name='Tube', time='6:01:20', copy_output_to=sim_path)
