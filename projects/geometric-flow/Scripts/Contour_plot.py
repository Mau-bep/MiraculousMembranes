import matplotlib.pyplot as plt 
import numpy as np 
import matplotlib as mpl



main_dir = "../Results/Mem3DG_Bead_Reciprocal_arcsim_up_Phase/"
main_dir = "../Results/Particle_wrapping_on_plane_phase_space_march/"

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

    KB1_plot = []
    EKB1_plot = []
    KB2_plot = []
    EKB2_plot = []
    KB5_plot = []
    EKB5_plot = []

    Coverages = []
    Energies = []

    while(line):
        splitted_line = line.split(' ')

        Kbs_1.append(splitted_line[1])
        Strengths_1.append(splitted_line[2])
        
        Coverages.append( float(splitted_line[3]))
        Energies.append( -1*float(splitted_line[2]))
        # print(Kbs_1[-1])
        if(Kbs_1[-1] == "1"):
            KB1_plot.append(Coverages[-1])
            EKB1_plot.append(Energies[-1])
            # print("This is first true")
        
        if(Kbs_1[-1] == "2"):
            KB2_plot.append(Coverages[-1])
            EKB2_plot.append(Energies[-1])
            # print("This is second true")
        
        if(Kbs_1[-1] == "5"):
            KB5_plot.append(Coverages[-1])
            EKB5_plot.append(Energies[-1])
            # print("This is third true")
        
        
        line = file.readline()

    # print(np.unique(Kbs))
    # print(np.unique(Strengths))
    file.close()
    # plt.clf()
    # print(Energies)
    # plt.scatter(Energies,Coverages)
    # plt.show()
    

    Strengths = np.sort(np.array(np.unique(Strengths_1),dtype = float))
    Kbs = np.sort(np.array(np.unique(Kbs_1) ,dtype = float))


    cmap = plt.colormaps['viridis']
    norm = mpl.colors.Normalize(vmin = 1.0,vmax = np.max(Kbs[4]))

    plt.scatter(EKB1_plot,KB1_plot, label = 'Kb = 1.0', color = cmap(norm(1.0)))
    # plt.plot(EKB1_plot,KB1_plot, ls = 'dashed', color = cmap(norm(1.0)))
    
    plt.scatter(EKB2_plot,KB2_plot, label = 'Kb = 2.0', color = cmap(norm(2.0)))
    # plt.plot(EKB2_plot,KB2_plot, ls = 'dashed', color = cmap(norm(2.0)))
    
    plt.scatter(EKB5_plot,KB5_plot, label = 'Kb = 5.0', color = cmap(norm(5.0)))
    # plt.plot(EKB5_plot,KB5_plot, ls = 'dashed', color = cmap(norm(5.0)))
    plt.axhline(y=1.0,ls='dashed',color='black')
    plt.legend()
    plt.xlabel("Interaction Energy")
    plt.ylabel("Coverage")
    plt.savefig(main_dir+"Coverage_plot.png",bbox_inches = 'tight')
    plt.show()


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

    # print(Z)
    
    cmap = plt.colormaps['viridis']
    norm = mpl.colors.Normalize(vmin = 1.0,vmax = np.max(Kbs[4]))
    print(Kbs)
    print(Strengths)
    plt.scatter(Strengths,Z[:,1],c =cmap(norm(Kbs[1])),label='Kb = {}'.format(Kbs[1]))
    plt.scatter(Strengths,Z[:,2],c =cmap(norm(Kbs[2])),label='Kb = {}'.format(Kbs[2]))
    plt.xlabel("Interaction strength constant")
    plt.ylabel("Coverage")
    plt.scatter(Strengths,Z[:,3],c =cmap(norm(Kbs[3])),label='Kb = {}'.format(Kbs[3]))
    plt.scatter(Strengths,Z[:,4],c =cmap(norm(Kbs[4])),label='Kb = {}'.format(Kbs[4]))
    plt.axhline(1.0,ls ='dashed',c='black')
    plt.legend()
    # plt.show()
    plt.clf()



    # plt.contour(Kbs,Strengths,Z,10)
    # plt.show()

    x = Kbs  # len = 11
    y = Strengths  # len = 7



    # cmap = plt.colormaps['viridis']
    norm = mpl.colors.Normalize(vmin = 0, vmax =1.0)
    
    # x = np.linspace(0,100,3)
    # plt.plot(x,x,c = cmap(0.5))
    # plt.show()
    fig, ax = plt.subplots()
    im = ax.pcolormesh(x, y, Z,cmap=cmap,vmin=0.0,vmax = 1.0) 
    # im = ax.contour(x,y,Z, norm =norm, cmap = cmap,levels=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],vmin=0.0,vmax = 1.0,extend = 'max')
    ax.set_xlabel("Kb")
    ax.set_ylabel("Interaction strength")
    
    fig.colorbar(im, ax=ax, label = "Coverage")
    # plt.show()
    plt.clf()




def main2():
    



main()