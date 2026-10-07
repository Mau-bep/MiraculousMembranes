import sys 
import os
import argparse
from cli_args import number_str
from subjob import write_subjob
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 
import numpy as np
# Lets do a rework here
import json

from jinja2 import Environment, FileSystemLoader

# We need jinja and json here




parser = argparse.ArgumentParser(description="Writes the config (if any) and the sbatch subjob for one run.")
parser.add_argument("--KA", type=number_str, required=True, help="surface tension constant KA")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB")
parser.add_argument("--relaxation_step", type=int, required=True, help="saved step to relax from")
args = parser.parse_args()

KA = args.KA
KB = args.KB
relaxation_step = args.relaxation_step


def get_bead_pos(folderpath,relaxation_step):
    filename = folderpath + "Bead_0_data.txt"
    data = np.loadtxt(filename,skiprows=1)
    # Ok now i need the position of the bead 
    x = data[relaxation_step,0]
    y = data[relaxation_step,1]
    z = data[relaxation_step,2]
    return [x,y,z]

def Create_json_relaxation(KA,KB,relaxation_step):
    os.makedirs("../Config_files/",exist_ok = True)

    env = Environment(loader=FileSystemLoader('../Templates/'))
    template = env.get_template('Tube_relaxation.txt')
    # I need to get the xpos ypos zpos

    # folderpath = "../Results/Tube_pulling_on_plane/Surface_tension_0.0500_Bending_20.0000_Bead_radius_0.2000_str_10.0000_Bead_radius_0.4000_str_0.0000_Bonds_Lineal_1000.0000_Lineal_1000.0000_Nsim_4/"
    folderpath = "../Results/Tube_pulling_on_plane/Surface_tension_0.0050_Bending_30.0000_Bead_radius_0.1000_str_10.0000_Bead_radius_0.4000_str_0.0000_Bonds_Lineal_1500.0000_Lineal_1500.0000_Nsim_1002/"
    
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


write_subjob('subjob_tube_relaxation_KA_{}_KB_{}_Nsim_{}'.format(KA,KB,relaxation_step),
             'output_tube_relaxation_{}_KA_{}_KB_{}_Nsim_{}'.format(relaxation_step,KA,KB,relaxation_step),
             '../build/bin/main_cluster {}'.format(Config_path),
             time='1-4:00:00', mem='4G')
