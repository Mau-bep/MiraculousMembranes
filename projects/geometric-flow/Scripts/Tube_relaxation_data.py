import numpy as np
import matplotlib.pyplot as plt 


def main(folderpath,Nsims):
    
    dx = []
    fxs = []

    for Nsim in Nsims:

        Sim_path = folderpath + "Surface_tension_0.0500_Bending_10.0000_Bead_radius_0.2000_str_10.0000_Nsim_{}/".format(Nsim)

        Bead_data = Sim_path+"Bead_0_data.txt"

        data = np.loadtxt(Bead_data,skiprows=1)
        # Ok so i have the bead data (I kinda care only about the force in the x direction at the last step)

        xpos = data[-1,0]
        fx = data[-1,3]
        dx.append(xpos)
        fxs.append(fx)
    
    # I should have all the data that i need at this point

    f = open(folderpath + "Tube_relaxation_data.txt",'w')

    for i in range(len(Nsims)):
        f.write("{} {} \n".format(dx[i],fxs[i]))
    
    f.close()


    plt.scatter(dx,fxs,color='black')
    plt.xlabel("dx")
    plt.ylabel("fx")
    plt.savefig(folderpath + "Tube_relaxation_data.png",bbox_inches='tight')
    return 


Nsims = [i for i in range(0,340,20)]

main(folderpath="../Results/Tube_for_relaxation/",Nsims=Nsims)

