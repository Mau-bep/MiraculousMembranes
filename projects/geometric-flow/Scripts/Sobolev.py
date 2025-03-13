import numpy as np 
import matplotlib.pyplot as plt 






def main():

    # We want to see the energy evolution.


    Output_data = np.loadtxt("../Results/Sobolev/nu_0.650_c0_0.000_KA_1.000_KB_1.000000_Nsim3/Output_data.txt",skiprows=1)


    plt.plot(Output_data[:,2]/(4*np.pi*1.0))
    plt.xlabel("Step")
    plt.ylabel("Energy")
    # plt.yscale('log')

    plt.show()

main()