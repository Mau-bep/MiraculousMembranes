import numpy as np 
import matplotlib.pyplot as plt 
import os


Strengths = [0.1, 0.25,0.5,0.75,1.0,1.25,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,5.75,5.0,5.25,5.5,5.75,6.0]
# Learn to do it with linspace
KB = 1.0
KA = 100000
Nsim = 1
Init_cond = 2




def main():
    dLs = []
    Force = []
    for strength in Strengths:
        folder = "../Results/Mem3DG_Bead_pulling_up_oct_arcsim/nu_1.000_radius_0.200_KA_100000.000_KB_{:.6f}_strength_{:.6f}_Init_cond_2_Nsim_1/".format(KB,strength)
        # I need to open 3 files bead fixed, bead moving and Last frame.

        # 
        filename_bead_moving = folder+ "Bead_0_data.txt"
        f = open(filename_bead_moving)
        Lines = f.readlines()
        line_first_bead = Lines[-1]
        line_first_bead = line_first_bead.split(" ")
        F_moving = np.array([line_first_bead[3], line_first_bead[4], line_first_bead[5] ],dtype=float)

        f.close()

        # Ok so this one has the position and the 
        filename_bead_fixed = folder + "Bead_1_data.txt"
        f = open(filename_bead_fixed)
        Lines = f.readlines()
        line_second_bead = Lines[-1]
        line_second_bead = line_second_bead.split(" ")
        F_fixed = np.array([line_second_bead[3], line_second_bead[4], line_second_bead[5] ],dtype=float)

        f.close()

        # Now i need to find the last  

        files = os.listdir(folder)
        files = [f for f in files if os.path.isfile(folder+'/'+f)]
        higher_index = 0

        for file in files:
            # I want to find the highest number
            line = file.split("_")
            if(line[0][0]=="m"):
                line = line[1].split(".")
                index = int(line[0])
                if(index>higher_index):
                    higher_index = index

        print("The highest index is {}".format(higher_index))
        high, low = find_higher_lowest( folder+"membrane_{}.obj".format(higher_index))
        dL = high-low
        F_tot = F_fixed+ F_moving
        F_tot = np.sqrt(np.sum(F_tot*F_tot))

        dLs.append(dL)
        Force.append(F_tot)
    

    plt.plot(Force,dLs)
    plt.xlabel(r"$\Delta L$")
    plt.ylabel(r"Force")
    plt.savefig(folder+"Pulling_force_plot.jpg",bbox_inches = 'tight')





def find_higher_lowest(filename):
    f = open(filename)
    line = f.readline()

    line = f.readline()
    lowest = 0.0
    highest = 0.0
    while(line):

        x = line.split(" ")[1]
        x = float(x)

        if(x > highest):
            highest = x
        if(x < lowest):
            lowest = x

        line = f.readline()

    return [lowest,highest]

        




main()