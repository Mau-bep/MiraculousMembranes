import numpy
import matplotlib.pyplot as plt 
import json
from jinja2 import Environment, FileSystemLoader
import numpy as np
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



def Create_json_wrapping_two(ka,kb,r,inter_str,angle):
    theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Two_beads.txt')
    
    # Radius of the position of the beads is R_v-2*r_b
    Rpos_beads = 2.0-0.2*2
    xpos  = Rpos_beads*np.cos(theta)
    ypos1 = Rpos_beads*np.sin(theta)
    ypos2 = -Rpos_beads*np.sin(theta)


    disp = -1*Rpos_beads*np.cos(theta)+ 0.5*(np.sqrt(2.1*2.1-Rpos_beads*Rpos_beads*np.sin(theta)**2 ) + np.sqrt(2.0*2.0-Rpos_beads*Rpos_beads*np.sin(theta)**2)) 

    output_from_parsed_template = template.render(KA = ka, KB = kb,radius = r,xdisp = disp,xpos1 = xpos,xpos2 =xpos, ypos1= ypos1, ypos2 = ypos2 ,interaction=inter_str, theta = theta)
    data = json.loads(output_from_parsed_template)
    Config_path = '../Config_files/Wrapping_two_{}_strg_{}_radius_{}_KA_{}_KB_{}.json'.format(angle,inter_str,r,ka,kb) 
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path




Create_json_wrapping_two(10,20,0.2,600,0.5)


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


# Create_json_barbell(800,10.0)
# Create_json_wrapping(0.005,10.0,0.4,400)

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