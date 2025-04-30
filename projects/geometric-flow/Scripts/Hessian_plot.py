import numpy as np
import matplotlib.pyplot as plt 






def main():
    filename = "../Results/BFGS_tests/Bending_20.0000_Surface_tension_10.0000_Bead_radius_0.2000_str_600.0000_Nsim_10/Hessian_2.txt"

    data = np.loadtxt(filename, skiprows=0)
    print("The min value is " + str(np.min(data)) + " and the max is " + str(np.max(data)))
    plt.axis("equal")
    plt.pcolormesh(data)
    plt.colorbar()
    plt.show()

main()