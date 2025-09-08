import numpy as np
import matplotlib.pyplot as plt






# We do this function and we can go home




def Create_data_file(folder_dir,KA,KB,KI, theta):
    dir = folder_dir + "Bending_{0:.4f}_Surface_tension_{1:.4f}_Edge_reg_1.0000_Bead_radius_0.3000_Frenkel_Normal_nopush_str_{2:.4f}_theta_const_{3:.4f}_Bead_radius_0.3000_Frenkel_Normal_nopush_str_{2:.4f}_theta_const_{3:.4f}_Switch_Newton_Switch_Free_beads_Nsim_1/".format(KB,KA,KI,theta)


    # Gorgeous i have the precise ffolder

    file_E = open(dir+  "Energy_output.txt","w+")
    output_data = open(dir+ "Output_data.txt","r")
    line = output_data.readline()
    # The first line doesnt contain anything
    line = output_data.readline()

    while(line):
        splitted_line = line.split(' ')
        if(len(splitted_line)>4):
            time = float(splitted_line[0])
            Ebend = float(splitted_line[3])
            Esurf = float(splitted_line[4])
            
            Eedge = float(splitted_line[5])
            Ebead1 = float(splitted_line[6])
            Ebead2 = float(splitted_line[7])
            Etot = float(splitted_line[8])
                # Now we have everything we need
            file_E.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f}\n".format(time,Ebend,Esurf,Eedge,Ebead1,Ebead2, Etot))
        line = output_data.readline()
    
    # That is all girl






folder = "../Results/Two_beads_with_switch_sept_dynamic/"

KA = 15
KB = 10

KI = 400

thetas = [0.3 +0.05*i for i in range(20)]

for i in thetas:
    print("Doing theta = {}".format(i))
    Create_data_file(folder,KA,KB,KI,i)