import numpy as np

import matplotlib.pyplot as plt 





angles = [i/100.0 for i in range(5,150,5)]

print(angles)

def main(Nsim):

    # So the idea is to read the data
    base_path = "../Results/Two_beads/"
    
    Bending_energies = []

    for i in range(len(angles)):
        folder = "Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Bead_radius_0.2000_str_600.0000_theta_const_{0:.4f}_Nsim_{1}".format(angles[i],Nsim)
        # Great i have the folder, next thing is the actual bending energy

        Output_data = np.loadtxt(base_path + folder + "/Output_data.txt", skiprows=1)

        # I just want the bending energy
        E_B = Output_data[-1,5]/(20.0*4*np.pi)
        E_B = Output_data[-1,5]
        Bending_energies.append(E_B)

    
    aangles = np.array(angles)
    plt.plot(aangles,Bending_energies,color='black',label='Nsmim = {}'.format(Nsim))
    
    plt.xlabel("Angle [rad]")
    plt.ylabel("Bending Energy")

    return


# main(1)
main(10)
main(20)


plt.show()