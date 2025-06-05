import numpy as np 
import matplotlib.pyplot as plt 






thetas = [i/100 for i in range(5,130,5)]
path = "../Results/Two_beads/"

path2 = "../Results/Two_beads_with_switch/"

def main(theta):

    folderpath = path +"Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Nsim_20/".format(theta)
    folderpath2 = path2 + "Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Switch_Free_beads_Switch_t_200000_Nsim_2/".format(theta)
    folderpath3 = path2 + "Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Switch_No_remesh_Switch_t_200000_Nsim_1/".format(theta)

    folderpath4 = path2 + "Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_Shifted_LJ_Normal_nopush_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_Shifted_LJ_Normal_nopush_str_600.0000_theta_const_{0:.4f}_Switch_Free_beads_Switch_t_200000_Nsim_3/".format(theta)

    Data = np.loadtxt(folderpath + "Output_data.txt", skiprows=1)
    Data2 = np.loadtxt(folderpath2 + "Output_data.txt", skiprows=1)
    Data3 = np.loadtxt(folderpath3 + "Output_data.txt", skiprows=1)
    Data4 = np.loadtxt(folderpath4 + "Output_data.txt", skiprows=1)

    # plt.plot(Data[:,7],color='black',label='Total E')
    plt.plot(Data2[:,7],color='blue',label='Total E with switch')
    plt.plot(Data3[:,7],color='red',label='Total E with switch no remesh')
    # plt.plot(Data4[:,7],color='green',label='Total E longersim')
    plt.axvline(x=200,color='black',ls='dashed')
    plt.legend()
    plt.title("Theta = {0:.4f}".format(theta),fontsize=15)
    plt.xlabel("Time",fontsize=15)
    plt.ylabel("Total Energy",fontsize = 15)
    # plt.yscale("log")

    plt.show()

    return 



for i in range(len(thetas)):
    main(float(thetas[i]))
    # plt.show()
    # plt.clf()
    # plt.close()
