import numpy as np

import matplotlib.pyplot as plt 





angles = [i/100.0 for i in range(5,130,5)]

print(angles)

def main(Nsim):

    # So the idea is to read the data
    base_path = "../Results/Two_beads_with_switch/"
    
    Bending_energies = []
    Interaction_1 = []
    Interaction_2 = []
    Surface_tension = []

    aangles =[]
    for i in range(len(angles)):
        # Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_1.2500_Bead_radius_0.2000_str_600.0000_theta_const_1.2500_Switch_No_remesh_Switch_t_200000_Nsim_1
        folder = "Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Switch_No_remesh_Switch_t_200000_Nsim_{1}".format(angles[i],Nsim)
        # Great i have the folder, next thing is the actual bending energy

        Output_data = np.loadtxt(base_path + folder + "/Output_data.txt", skiprows=1)
        print(angles[i])
        print(Output_data.shape)
        if(Output_data.shape[0] == 14):
            print("\t \t No data")
            print("This is for angle = {}".format(angles[i])) 
            continue

        aangles.append(angles[i])
                # I just want the bending energy
        E_B = Output_data[-1,4]/(20.0*4*np.pi)
        # E_I1 = Output_data[-1,5]/(20.0*4*np.pi)

        # E_B = Output_data[-1,5]
        print(E_B)
        Bending_energies.append(E_B)

    
    aangles = np.array(aangles)
    plt.scatter(aangles,Bending_energies,color='purple',label='Switch remesh'.format(Nsim))
    # plt.show("No remesh")
    plt.xlabel("Angle [rad]")
    plt.ylabel("Bending Energy")

    return


def main_2(Nsim):

    # So the idea is to read the data
    base_path = "../Results/Two_beads_with_switch/"
    
    Bending_energies = []
    Interaction_1 = []
    Interaction_2 = []
    Surface_tension = []

    aangles =[]
    for i in range(len(angles)):
        # Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_1.2500_Bead_radius_0.2000_str_600.0000_theta_const_1.2500_Switch_No_remesh_Switch_t_200000_Nsim_1
        folder = "Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Switch_Free_beads_Switch_t_200000_Nsim_{1}".format(angles[i],Nsim)
        # Great i have the folder, next thing is the actual bending energy

        Output_data = np.loadtxt(base_path + folder + "/Output_data.txt", skiprows=1)
        print(angles[i])
        print(Output_data.shape)
        if(Output_data.shape[0] == 14):
            print("\t \t No data")
            print("This is for angle = {}".format(angles[i])) 
            continue

        aangles.append(angles[i])
                # I just want the bending energy
        E_B = Output_data[-1,4]/(20.0*4*np.pi)
        # E_I1 = Output_data[-1,5]/(20.0*4*np.pi)

        # E_B = Output_data[-1,5]
        print(E_B)
        Bending_energies.append(E_B)

    
    aangles = np.array(aangles)
    plt.scatter(aangles,Bending_energies,color='black',label='Switch free bead'.format(Nsim))
    
    plt.xlabel("Angle [rad]")
    plt.ylabel("Bending Energy")

    return


main(1)
plt.legend()
plt.show()

main_2(2)
plt.legend()


plt.show()

main(1)
main_2(2)
plt.legend()
plt.show()