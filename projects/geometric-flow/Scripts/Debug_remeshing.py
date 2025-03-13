import numpy as  np
import matplotlib.pyplot as plt 




def main():
    # We need the directory     

    parent_dir = "../Results/Debug_remesh/"
    dir = parent_dir+"H2_Bending_1.0000_1.0000_1.0000_0.0000_H2_Bead_radius_0.3000_str_400.0000_Nsim_7/"
    

    dihedralsdir = dir+"dihedrals_evol.txt"
    edgelengthsdir = dir+"edgelengths.txt"

    
    dih = open(dihedralsdir, "r")

    line = dih.readline()
    line = dih.readline()
    counter = 0


    while(line):
        
        if(counter<-400):
            splitted = line.split(" ")[:-1]
        
            distribution = np.array(splitted).astype(np.float)
            distribution = np.abs(distribution)    
            plt.xlabel("Dihedral angle (radians)")
            plt.ylabel("Frequency")
            plt.title("Dihedral angle distribution at step "+str(counter))
            plt.hist(distribution, bins=100)
            plt.yscale("log")
            plt.show()


        plt.clf()
        counter+=1
        line = dih.readline()

    dih.close()

    
    edges = open(edgelengthsdir, "r")
    line = edges.readline()

    counter = 0

    while(line):
        if(counter>=400):
            splitted = line.split(" ")[:-1]
        
            distribution = np.array(splitted).astype(np.float)
            distribution = np.abs(distribution)    
            plt.xlabel("Edgelength")
            plt.ylabel("Frequency")
            plt.title("Edgelength distribution at step "+str(counter))
            plt.hist(distribution, bins=100)
            plt.yscale("log")
            plt.show()


        plt.clf()
        counter+=1
        line = edges.readline()



    return


main()