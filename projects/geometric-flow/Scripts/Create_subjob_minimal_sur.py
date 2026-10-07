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
parser.add_argument("--finaldist", type=number_str, required=True, help="final distance between the beads")
args = parser.parse_args()

finaldist = args.finaldist
Nsim = 1






# KE = 1

location = ["unavailable", "outside", "inside"]



def RescaleMesh(Targetdist):
    


    scaleFactor = float(Targetdist)
    # Now i need to load the mesh
    fRead = open("../../../input/Opencil.obj",'r')
    InputFileDir =FirstDir+"Init_{}.obj".format(Targetdist) 
    print("File at {}".format(InputFileDir))
    fWrite = open(InputFileDir,'w+')
    for line in fRead:
        text = line.strip()
        splittedText = text.split(' ')
        if(splittedText[0]=='v'):
            # Then i need to rescale
            x = float(splittedText[1])
            if( x < 0.15 and x> -0.15):
                xN = x*scaleFactor
            else:
                if(x<0):
                    xN = x + displacement
                if(x>0):
                    xN = x - displacement
            
            fWrite.write("v {} {} {}\n".format(xN,splittedText[2],splittedText[3]))
        else:
            fWrite.write(line)
    
    return InputFileDir






def Create_json_wrapping_two(dist):
    # theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))


    template = env.get_template('ContinueTwoBeads.txt')
    V1x = -10.0
    L0 = 0.4
    if(float(dist) > 1.4):
        V1x = 10.0
        # I also need to change the L0
        L0 = 2.4
    FinalX1 = float(dist)/2.0
    # We should do  
    output_from_parsed_template = template.render(Finaldist = dist, Finalx1 = FinalX1, Finalx2 = -FinalX1, V1x = V1x, V2x = -V1x,L0 = L0)

    # print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)


    # print("something\n")
    Config_path = '../Config_files/Wrapping_two_{0}_Frenkel_2.json'.format(dist) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path



def Create_json_wrapping_two(dist,InputFileDir):
    # theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))


    template = env.get_template('ContinueTwoBeadsRescaleArea.txt')
    
    # W
    x1 = float(dist)/2.0
    x2 = -x1


    output_from_parsed_template = template.render(InputFile = InputFileDir,FirstDir=FirstDir ,Finaldist = dist, x1 = x1, x2 = x2)

    # print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)


    # print("something\n")
    Config_path = '../Config_files/Wrapping_two_{0}_rescaleArea.json'.format(dist) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path




os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)



# Config_path, sim_path = Create_json_wrapping_two(KA,KB,radius,Strength,angle)
# # Hopefully this works
# Config_path, sim_path = Create_json_wrapping_two_outside(angle,outside1,outside2)
FirstDir = "../Results/TwoBeadsManualRescaleArea/"
os.makedirs(FirstDir,exist_ok=True)
InputFileDir = RescaleMesh(finaldist)
Config_path, sim_path = Create_json_wrapping_two(finaldist,InputFileDir)


Output_name = 'output_two_{0}_rescaleArea.output'.format(finaldist)

write_subjob('subjob_two_{0}_rescaleArea'.format(finaldist),
             Output_name,
             '../build/bin/main_cluster {}'.format(Config_path),
             job_name='TwoSimple', time='6:01:20', mem='2G', copy_output_to=sim_path)
