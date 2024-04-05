import numpy as np
import os 
import matplotlib.pyplot as plt 

import sys

import matplotlib.cm as cm



def read_output(filename):
    N_vert_list=[]
    Time_taken=[]
    Sim_step=[]
    file = open(filename)
    line=file.readline()
    while(line):
        splitted_line = line.split(' ')
        if(len(splitted_line)>4):
            if(splitted_line[3]=="vertices"):
                N_vert=int(splitted_line[5])
            if(splitted_line[1]=="system"):
                step = int(splitted_line[4])
                # if(step>=500000):
                #     break
            if(splitted_line[1]=="thousand"):
                taken = float(splitted_line[4])/1000
                N_vert_list.append(N_vert)
                Sim_step.append(step)
                Time_taken.append(taken)
            # print("Line is" + line + "size {}".format(len(splitted_line)))
        line=file.readline()
    # I need to set the color 
    norm = plt.Normalize()
    
    colors = plt.cm.viridis(norm(Sim_step))
    print(max(Sim_step))
    # plt.plot(N_vert_list,Time_taken,color='purple',alpha=0.5)
    plt.scatter(N_vert_list,Time_taken,c=colors) 



file= "../Outputs/output_sample.txt"

# read_output(file)


file_ex=""
# for nu in [0.525, 0.55,0.575,0.5,0.625,0.65,0.675,0.6,0.725,0.75,0.775,0.7,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0]:
#     file_ex="../Outputs/output_Mem3DG_v_{}_KB_0.01_evol_Init_cond_1_Nsim_11".format(nu)
#     # print(file_ex)
#     read_output(file_ex)

for nu in [0.525, 0.55,0.575,0.5,0.625,0.65,0.675,0.6,0.725,0.7]:
    file_ex="../Outputs/output_Mem3DG_v_{}_KB_0.01_evol_Init_cond_1_Nsim_11".format(nu)
    # print(file_ex)
    read_output(file_ex)

for nu in [0.45, 0.475, 0.525, 0.55,0.575,0.5,0.625,0.65,0.675,0.6,0.725,0.75,0.7]:
    file_ex="../Outputs/output_Mem3DG_v_{}_KB_0.01_evol_Init_cond_2_Nsim_11".format(nu)
    # print(file_ex)
    read_output(file_ex)

for nu in [0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45, 0.475, 0.525, 0.55,0.575,0.5,0.6]:
    file_ex="../Outputs/output_Mem3DG_v_{}_KB_0.01_evol_Init_cond_3_Nsim_11".format(nu)
    # print(file_ex)
    read_output(file_ex)


plt.ylabel('time [s]')
plt.xlabel('N Vertices ')

    
plt.show()
