import matplotlib.pyplot as plt 
import numpy as np 
import matplotlib as mpl



main_dir = "../Results/Mem3DG_Bead_Reciprocal_arcsim_up_Phase/"

def main():
    # i WANT TO READ THE FILE 


    filepath = main_dir + "Coverage_final.txt"

    file = open(filepath)

    line = file.readline()
    line = file.readline()

    X = []
    Y = []
    Z = []

    # I will read the file twice
    Kbs_1 = []
    Strengths_1 = []

    while(line):
        splitted_line = line.split(' ')

        Kbs_1.append(splitted_line[1])
        Strengths_1.append(splitted_line[2])



        line = file.readline()

    # print(np.unique(Kbs))
    # print(np.unique(Strengths))
    file.close()

    Strengths = np.sort(np.array(np.unique(Strengths_1),dtype = float))
    Kbs = np.sort(np.array(np.unique(Kbs_1) ,dtype = float))


    print(Kbs)
    print(Strengths)
    file = open(filepath)

    line = file.readline()
    line = file.readline()
    print(len(Kbs))
    print(len(Strengths))
    Z = np.zeros((len(Strengths),len(Kbs)))
    print(Z.shape)

    while(line):

        splitted_line = line.split(' ')
        KB = float(splitted_line[1])
        STR = float(splitted_line[2])

        coverage = float(splitted_line[3])
        id2 = np.where(abs(Kbs-KB)<1e-5)[0][0]
        id1 = np.where(abs(Strengths-STR)<1e-5)[0][0]
    

        Z[id1][id2]=coverage
        




        line = file.readline()

    print(Z)


    # plt.contour(Kbs,Strengths,Z,10)
    # plt.show()

    x = Kbs  # len = 11
    y = Strengths  # len = 7

    cmap = plt.colormaps['spring']
    norm = mpl.colors.Normalize(vmin = 0, vmax =1.0)
    fig, ax = plt.subplots()
    im = ax.pcolormesh(x, y, Z,cmap=cmap,vmin=0.0,vmax = 1.0) 
    # im = ax.contourf(x,y,Z, norm =norm, cmap = cmap,levels=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],vmin=0.0,vmax = 1.0,extend = 'max')
    ax.set_xlabel("Kb")
    ax.set_ylabel("Interaction strength")
    
    fig.colorbar(im, ax=ax, label = "Coverage")
    plt.show()



main()