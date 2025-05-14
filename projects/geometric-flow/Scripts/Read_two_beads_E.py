import numpy as np 
import matplotlib.pyplot as plt 






thetas = [i/100 for i in range(5,130,5)]
path = "../Results/Two_beads/"



def main(theta):

    folderpath = path +"Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Nsim_20/".format(theta)

    Data = np.loadtxt(folderpath + "Output_data.txt", skiprows=1)


    plt.plot(Data[:,7],color='black',label='Total E')
    plt.xlabel("Time")
    plt.ylabel("Total Energy")
    # plt.yscale("log")

    plt.show()

    return 



for i in range(len(thetas)):
    main(float(thetas[i]))
    # plt.show()
    # plt.clf()
    # plt.close()
