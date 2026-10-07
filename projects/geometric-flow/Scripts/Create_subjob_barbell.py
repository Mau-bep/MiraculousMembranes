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




parser = argparse.ArgumentParser(description="Writes the config (if any) and the sbatch subjob for one run.")
parser.add_argument("--KA", type=number_str, required=True, help="surface tension constant KA")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB")
parser.add_argument("--Nsim", type=number_str, required=True, help="simulation number (part of the file names)")
args = parser.parse_args()

KA = args.KA
KB = args.KB
Nsim = args.Nsim



def Create_json_barbell(ka,kb):

    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Barbell.txt')
    output_from_parsed_template = template.render(KA = ka, KB = kb)

    data = json.loads(output_from_parsed_template)

    Config_path = '../Config_files/Barbell_KA_{}_KB_{}.json'.format(ka,kb) 
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path




os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)

Config_path = Create_json_barbell(KA,KB)


write_subjob('subjob_serial_barbell_KA_{}_KB_{}_Nsim_{}'.format(KA,KB,Nsim),
             'output_serial_barbell_KA_{}_KB_{}_Nsim_{}'.format(KA,KB,Nsim),
             '../build/bin/main_cluster {}'.format(Config_path),
             time='40:00:00', mem='4G')
