import numpy
import matplotlib.pyplot as plt 
import json
from jinja2 import Environment, FileSystemLoader

import os


def Create_json_wrapping(ka,kb,r,inter_str):

    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Wrapping.txt')
    output_from_parsed_template = template.render(KA = ka, KB = kb,radius = r,xpos = r*1.9 ,interaction=inter_str)

    data = json.loads(output_from_parsed_template)

    Config_path = '../Config_files/Wrapping_strg_{}_radius_{}_KA_{}_KB_{}.json'.format(inter_str,r,ka,kb) 
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path





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



def get_bead_pos(folderpath,relaxation_step):
    filename = folderpath + "Bead_0_data.txt"
    data = numpy.loadtxt(filename,skiprows=1)
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

    folderpath = "../Results/Tube_pulling_on_plane/Surface_tension_0.0500_Bending_20.0000_Bead_radius_0.2000_str_10.0000_Bead_radius_0.4000_str_0.0000_Bonds_Lineal_1000.0000_Lineal_1000.0000_Nsim_4/"

    [x,y,z] = get_bead_pos(folderpath,relaxation_step)
    filename = folderpath +"membrane_{}.obj".format(relaxation_step*100)

    output_from_parsed_template = template.render(KA = KA, KB = KB, xpos = x, ypos = y, zpos = z,init_file = filename )
    data = json.loads(output_from_parsed_template)

    Config_path = '../Config_files/Tube_relaxation_step_{}_KA_{}_KB_{}.json'.format(relaxation_step,KA,KB) 
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path


# Create_json_barbell(0.05,14.0)
# Create_json_wrapping(0.005,10.0,0.4,200)

# Create_json_relaxation(0.05,10.0,61)


import numpy as np


def Func(A):
    # val = A.shape
    B = np.zeros((3*A.shape[0],3*A.shape[1]))
    for i in range(2):
        for j in range(2):
            B[3*i,3*j] = A[i,j]
            B[3*i+1,3*j+1] = A[i,j]
            B[3*i+2,3*j+2] = A[i,j]
            

    print(B)


A = np.array([[1,2],[3,4]])


Func(A)