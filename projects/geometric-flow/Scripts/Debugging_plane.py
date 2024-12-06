import numpy as np 
import matplotlib.pyplot as plt 







def main():

    Data = np.loadtxt("../Results/Tests/Force_dist.txt",skiprows= 1)


    r = Data[:,3]
    F_norm = Data[:,7]
    plt.xlabel('r')
    plt.ylabel("|F|")
    plt.scatter(r,F_norm,color= "purple")
    plt.yscale('log')
    plt.axvline(1.6,ls='dashed',color = "black")
    plt.show()
    print(len(r))





main()