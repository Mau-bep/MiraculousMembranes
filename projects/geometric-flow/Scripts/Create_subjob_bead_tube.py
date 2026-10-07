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
parser.add_argument("--angle", type=number_str, required=True, help="angle (as in the file names)")
args = parser.parse_args()

angle = args.angle
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


Output_name = 'output_bead_tube_theta_{}_{}.output'.format(angle,strg)

write_subjob('subjob_bead_tube_theta_{}_{}'.format(angle,strg),
             Output_name,
             '../build/bin/main_cluster {}'.format(Config_path),
             job_name='Wrap', time='10:01:20', mem='5G', error='%x_%j.err', copy_output_to=sim_path)
