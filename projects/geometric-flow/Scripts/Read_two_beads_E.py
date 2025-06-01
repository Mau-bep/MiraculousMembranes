import numpy as np 
import matplotlib.pyplot as plt 






thetas = [i/100 for i in range(5,130,5)]
path = "../Results/Two_beads/"

path2 = "../Results/Two_beads_with_switch/"

def main(theta):

    folderpath = path +"Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Nsim_20/".format(theta)
    folderpath2 = path2 + "Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Switch_Free_beads_Switch_t_200000_Nsim_2/".format(theta)
    folderpath3 = path2 + "Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Switch_No_remesh_Switch_t_200000_Nsim_1/".format(theta)

    Data = np.loadtxt(folderpath + "Output_data.txt", skiprows=1)
    Data2 = np.loadtxt(folderpath2 + "Output_data.txt", skiprows=1)
    Data3 = np.loadtxt(folderpath3 + "Output_data.txt", skiprows=1)

    plt.plot(Data[:,4]/(20.0*4*np.pi),color='black',label='Total E')
    plt.plot(Data2[:,4]/(20.0*4*np.pi),color='blue',label='Total E with switch')
    plt.plot(Data3[:,4]/(20.0*4*np.pi),color='red',label='Total E with switch no remesh')
    plt.legend()
    plt.title("Theta = {0:.4f}".format(theta))
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
