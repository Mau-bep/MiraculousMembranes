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

# We need jinja and json here

import numpy as np





parser = argparse.ArgumentParser(description="Writes the config (if any) and the sbatch subjob for one run.")
parser.add_argument("--inter_str", type=number_str, required=True, help="bead-membrane interaction strength")
parser.add_argument("--radius", type=number_str, required=True, help="bead radius")
parser.add_argument("--KA", type=number_str, required=True, help="surface tension constant KA")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB")
parser.add_argument("--batch_tag", type=str, required=True, help="tag of the batch (also the results folder name)")
parser.add_argument("--unique_tag", type=str, required=True, help="unique tag of this run (prefix of the config, subjob and output files)")
args = parser.parse_args()

Strength = args.inter_str
radius = args.radius

KA = args.KA
KB = args.KB

Batch_tag = args.batch_tag
Unique_tag = args.unique_tag

KE = 1.0

# Each batch gets its own results folder
Batch_dir = '../Results/{}/'.format(Batch_tag)






def Create_json_wrapping(ka,kb,r,inter_str):

    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Wrapping.txt')
    output_from_parsed_template = template.render(KA = ka, KB = kb,radius = r,xpos = 2+r ,interaction=inter_str)

    data = json.loads(output_from_parsed_template)

    Config_path = '../Config_files/Wrapping_strg_{}_radius_{}_KA_{}_KB_{}.json'.format(inter_str,r,ka,kb) 
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path

def Create_json_wrapping_vesicle(ka,kb,r,inter_str):

    rad = float(r)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Wrapping_vesicle.txt')
    output_from_parsed_template = template.render(KA = ka, KB = kb,radius = r,rc=rad*1.25,xpos = 7.0+rad*1.1 ,interaction=inter_str)

    data = json.loads(output_from_parsed_template)

    data['first_dir'] = Batch_dir

    Config_path = '../Config_files/{0}_ConfigFile.json'.format(Unique_tag)
    
    sim_path = data['first_dir']

    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path, sim_path




os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)

Config_path, sim_path = Create_json_wrapping_vesicle(KA,KB,radius,Strength)


Output_name = '{0}_output.output'.format(Unique_tag)

write_subjob('{0}_subjob'.format(Unique_tag),
             Output_name,
             '../build/bin/main_cluster {}'.format(Config_path),
             job_name=Batch_tag, mem='16G', copy_output_to=sim_path)
