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
parser.add_argument("--outside1", type=int, required=True, help="bead 1 location: 1 outside, 2 or -1 inside")
parser.add_argument("--outside2", type=int, required=True, help="bead 2 location: 1 outside, 2 or -1 inside")
parser.add_argument("--radius", type=float, required=True, help="bead radius")
parser.add_argument("--batch_tag", type=str, required=True, help="tag of the batch (prefix of the file names)")
parser.add_argument("--unique_tag", type=str, required=True, help="unique tag for the simulations ")
parser.add_argument("--Cov", type=float, required=False, help="Coverage for the two beads")
args = parser.parse_args()

angle = args.angle
outside1 = args.outside1
outside2 = args.outside2
radius = args.radius
Batch_tag = args.batch_tag
Unique_tag = args.unique_tag
Cov = args.Cov
ka = 1.0
Nsim = 1






# KE = 1

location = ["unavailable", "outside", "inside"]


Batch_dir = '../Results/{}/'.format(Batch_tag)
# Batch_tag = 'Wrap2'

def Create_json_wrapping_two(ka,kb,r,inter_str,angle):
    theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Two_beads.txt')
    
    # Radius of the position of the beads is R_v-2*r_b
    R_vesicle = 7.0
    r_bead = 1.0
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
    R_vesicle = 7.0
    r_bead = 1.0
    
    location = [1,"outside","inside"]

    dir = '"../Results/Two_beads_{}_{}_BFGS_MAR3/"'.format(location[outside1],location[outside2])

    v1x = 1.0*(outside1*-1)
    x1 = R_vesicle + 1.1*r_bead*outside1 

    r2 = R_vesicle + 1.1*r_bead*outside2
    # We should do 
    x2 = r2*np.cos(theta)
    y2 = r2*np.sin(theta)
    Leq = (R_vesicle)*theta

    v2x = 1*np.cos(theta)*outside2*-1
    v2y = 1*np.sin(theta)*outside2*-1

    
    output_from_parsed_template = template.render(Dir = dir,theta =theta, outside1 = outside1, v1x = v1x, x1 = x1,L0 = Leq, outside2 = outside2, v2x = v2x, v2y = v2y,x2 = x2, y2 = y2 )

    # print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)


    # print("something\n")
    Config_path = '../Config_files/Wrapping_two_{}_{}_{}_BFGS_M3.json'.format(angle,location[outside1],location[outside2]) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path





def Create_json_wrapping_two_fixed(dist, outside1, outside2):
    # theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))


    template = env.get_template('Wrapping_two_rigid.txt')
    
    # Radius of the position of the beads is R_v-2*r_b
    R_vesicle = 7.0
    r_bead = radius*0.9
    
    location = [1,"outside","inside"]

    dir = '"../Results/{0}_r_{1:.2f}_{2}_{3}/"'.format(Batch_tag,radius,location[outside1],location[outside2])

    x1 = float(dist)/2.0
    x2 = -float(dist)/2.0 


    disp = -1*np.sqrt(  (R_vesicle+r_bead)**2 -(float(dist)/2.0)**2)
    disp2 = 0.0

    print("Disp is {}".format(disp))


    if(outside1<0 and outside2 <0 ):
        disp =-1*np.sqrt(  (R_vesicle-r_bead)**2 -(float(dist)/2.0)**2)

    if(outside1*outside2<0):
        # In this case we need to  do math
        disp2 = ((R_vesicle+r_bead)**2-(R_vesicle-r_bead)**2)/(2*float(dist))
        disp = np.sqrt(  (R_vesicle+r_bead)**2 - (disp2+float(dist)/2)**2      ) 

    
    # return 
    # We should do  
    output_from_parsed_template = template.render(Dir = dir,r=radius,rc = radius*1.25,dist = dist, outside1 = outside1,disp = disp,disp2 = disp2, x1 = x1, outside2 = outside2, x2 = x2 )



    # print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)
    data['first_dir'] = Batch_dir

    # print("something\n")
    Config_path = '../Config_files/{0}_{1:.2f}_{2}_{3}_{4}_ST_{5}.json'.format(Batch_tag,radius,angle,location[outside1],location[outside2],ka) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path



def Create_json_wrapping_two_cov(dist, outside1, outside2):
    # theta = float(angle)
    os.makedirs("../Config_files/",exist_ok = True)
    env = Environment(loader=FileSystemLoader('../Templates/'))


    template = env.get_template('Wrapping_two_cov.txt')
    
    # Radius of the position of the beads is R_v-2*r_b
    R_vesicle = 7.0
    r_bead = radius*0.9
    
    location = [1,"outside","inside"]

    dir =Batch_dir
    x1 = float(dist)/2.0
    x2 = -float(dist)/2.0 


    disp = -1*np.sqrt(  (R_vesicle+r_bead)**2 -(float(dist)/2.0)**2)
    disp2 = 0.0

    print("Disp is {}".format(disp))


    if(outside1<0 and outside2 <0 ):
        disp =-1*np.sqrt(  (R_vesicle-r_bead)**2 -(float(dist)/2.0)**2)

    if(outside1*outside2<0):
        # In this case we need to  do math
        disp2 = ((R_vesicle+r_bead)**2-(R_vesicle-r_bead)**2)/(2*float(dist))
        disp = np.sqrt(  (R_vesicle+r_bead)**2 - (disp2+float(dist)/2)**2      ) 

    
    # return 
    # We should do  
    output_from_parsed_template = template.render(Dir = dir,dist = dist,disp = disp,disp2 = disp2, x1 = x1, x2 = x2, cov1 = Cov, cov2 = Cov )



    # print(output_from_parsed_template)
    data = json.loads(output_from_parsed_template)
    data['first_dir'] = Batch_dir

    # print("something\n")
    Config_path = '../Config_files/{0}_ConfigFile.json'.format(Unique_tag) 
    
    sim_path = data['first_dir']
    
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path , sim_path



os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)

Config_path, sim_path = Create_json_wrapping_two_cov(angle,outside1,outside2)

Output_name = '{0}_output.output'.format(Unique_tag)

write_subjob('{0}_subjob'.format(Unique_tag),
             Output_name,
             '../build/bin/main_cluster {}'.format(Config_path),
             job_name='Wrap', time='12:01:20', copy_output_to=sim_path)
