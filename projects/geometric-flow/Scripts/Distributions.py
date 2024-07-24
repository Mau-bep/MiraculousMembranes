import numpy as np
import matplotlib.pyplot as plt



# i WANT THE DIRECTORY


def checking_radius_dist():

    nu=1.0
    radius=1.0
    curvadap=0.0
    minrel=0.1
    KA=100000.0
    KB=1.0
    strength = 70.0
    init_cond=1 
    Nsim=12

    basic_name="../Results/Mem3DG_Bead_Reciprocal_arcsim/nu_{:.3f}_radius_{:.3f}_KA_{:.3f}_KB_{:.6f}_strength_{:.6f}_Init_cond_{}_Nsim_{}/".format(nu,radius,KA,KB,strength,init_cond,Nsim)




    radius_vec= np.loadtxt(basic_name+"Radius_distribution_step_"+str(6)+".txt",delimiter=' ')

    print(min(radius_vec))
    print(max(radius_vec))
    plt.hist(radius_vec,bins=200)
    plt.show()
    plt.clf()



    Touching_data=np.loadtxt(basic_name+"Touching_step_"+str(6)+".txt")

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(Touching_data[:,0], Touching_data[:,1], Touching_data[:,2], marker='o')
    plt.axis('equal')
    plt.show()

    # for angle in range(0, 360*4 + 1):
    #     # Normalize the angle to the range [-180, 180] for display
    #     angle_norm = (angle + 180) % 360 - 180

    #     # Cycle through a full rotation of elevation, then azimuth, roll, and all
    #     elev = azim = roll = 0
    #     if angle <= 360:
    #         elev = angle_norm
    #     elif angle <= 360*2:
    #         azim = angle_norm
    #     elif angle <= 360*3:
    #         roll = angle_norm
    #     else:
    #         elev = azim = roll = angle_norm

    #     # Update the axis view and title
    #     ax.view_init(elev, azim, roll)
    #     plt.title('Elevation: %d°, Azimuth: %d°, Roll: %d°' % (elev, azim, roll))

    #     plt.draw()
    #     plt.pause(.001)


def Plotting_coverage():

    # So i will have 3 radius
    # For every radius there will be a curve with a different KA
    # All the KA curves are plotted as a function of Strenght

    file = open("../Results/Mem3DG_Bead_Reciprocal_finemesh/Coverage_final.txt")
    line = file.readline()
    print(line)
    line = file.readline()

    plot_rad_1_x = [[],[],[],[]]
    plot_rad_1_y = [[],[],[],[]]

    plot_rad_2_x = [[],[],[],[]]
    plot_rad_2_y = [[],[],[],[]]

    plot_rad_3_x = [[],[],[],[]]
    plot_rad_3_y = [[],[],[],[]]


    while(line):
        splitted= line.split(' ')

        if(splitted[0])=="1":
            index=int(np.log10(float(splitted[1]))-2)
            plot_rad_1_x[index].append(float(splitted[2]))
            plot_rad_1_y[index].append(float(splitted[3]))

        if(splitted[0])=="2":
            index=int(np.log10(float(splitted[1]))-2)
            plot_rad_2_x[index].append(float(splitted[2]))
            plot_rad_2_y[index].append(float(splitted[3]))
        # I need to decide the plots

        if(splitted[0])=="3":
            index=int(np.log10(float(splitted[1]))-2)
            plot_rad_3_x[index].append(float(splitted[2]))
            plot_rad_3_y[index].append(float(splitted[3]))


        line =file.readline()
        # I need to decide the plots



    # Ok i have the plots ready 

    KAs = [1e2,1e3,1e4,1e5]
    for i in range(len(plot_rad_1_x)):
        plt.plot(plot_rad_1_x[i],plot_rad_1_y[i],label='KA = {}'.format(KAs[i]))
    plt.xlabel("Interacion strenght")
    plt.ylabel(" % Coverage")
    plt.legend()
    plt.ylim(0.0,1.1)
    plt.show()


    for i in range(len(plot_rad_1_x)):
        plt.plot(plot_rad_2_x[i],plot_rad_2_y[i],label='KA = {}'.format(KAs[i]))
    plt.xlabel("Interacion strenght")
    plt.ylabel(" % Coverage")
    plt.ylim(0.0,1.1)
    plt.legend()
    plt.show()

    for i in range(len(plot_rad_1_x)):
        plt.plot(plot_rad_3_x[i],plot_rad_3_y[i],label='KA = {}'.format(KAs[i]))
    plt.xlabel("Interacion strenght")
    plt.ylabel(" % Coverage")
    plt.ylim(0.0,1.1)
    plt.legend()
    plt.show()





def Plotting_coverage_varKB():

    # So i will have 3 radius
    # For every radius there will be a curve with a different KA
    # All the KA curves are plotted as a function of Strenght

    file = open("../Results/Mem3DG_Bead_Reciprocal_finemesh/Coverage_fina_varKB.txt")
    line = file.readline()
    print(line)
    line = file.readline()

    plot_rad_1_x = [[],[],[],[]]
    plot_rad_1_y = [[],[],[],[]]

    plot_rad_1_x_2 = [[],[],[],[]]
    plot_rad_1_y_2 = [[],[],[],[]]

    plot_rad_2_x = [[],[],[],[]]
    plot_rad_2_y = [[],[],[],[]]

    plot_rad_2_x_2 = [[],[],[],[]]
    plot_rad_2_y_2 = [[],[],[],[]]


    while(line):
        splitted= line.split(' ')

        if(splitted[0])=="1":
            index=int(np.log10(float(splitted[1]))-4)
            if(np.log10(float(splitted[4]))==0):
                plot_rad_1_x[index].append(float(splitted[2]))
                plot_rad_1_y[index].append(float(splitted[3]))
            else:
                plot_rad_1_x_2[index].append(float(splitted[2]))
                plot_rad_1_y_2[index].append(float(splitted[3]))
                

        if(splitted[0])=="2":
            index=int(np.log10(float(splitted[1]))-4)
            if(np.log10(float(splitted[4]))==0):
                plot_rad_2_x[index].append(float(splitted[2]))
                plot_rad_2_y[index].append(float(splitted[3]))
            else:
                plot_rad_2_x_2[index].append(float(splitted[2]))
                plot_rad_2_y_2[index].append(float(splitted[3]))
                
        # I need to decide the plots




        line =file.readline()
        # I need to decide the plots


    file.close()

    file = open("../Results/Mem3DG_Bead_Reciprocal_finemesh/Coverage_final.txt")
    line = file.readline()
    print(line)
    line = file.readline()

    plot_rad_1_x_3 = [[],[],[],[]]
    plot_rad_1_y_3 = [[],[],[],[]]

    plot_rad_2_x_3 = [[],[],[],[]]
    plot_rad_2_y_3 = [[],[],[],[]]



    while(line):
        splitted= line.split(' ')

        if(splitted[0])=="1":
            index=int(np.log10(float(splitted[1]))-4)
            if(index<0):
                line=file.readline()
                continue
            plot_rad_1_x_3[index].append(float(splitted[2]))
            plot_rad_1_y_3[index].append(float(splitted[3]))

        if(splitted[0])=="2":
            index=int(np.log10(float(splitted[1]))-4)
            plot_rad_2_x_3[index].append(float(splitted[2]))
            plot_rad_2_y_3[index].append(float(splitted[3]))
        # I need to decide the plots


        line =file.readline()

    linest=['solid','dashed']
    # Ok i have the plots ready 
    print(plot_rad_1_x)
    KAs = [1e4,1e5]
    for i in range(len(KAs)):
        plt.plot(plot_rad_1_x_3[i],plot_rad_1_y_3[i],ls=linest[i],label='KA = {}, KB = {}'.format(KAs[i],0.1),c='purple')

    for i in range(len(KAs)):
        plt.plot(plot_rad_1_x[i],plot_rad_1_y[i],ls=linest[i],label='KA = {}, KB = {}'.format(KAs[i],1.0),c='magenta')
    for i in range(len(KAs)):   
        plt.plot(plot_rad_1_x_2[i],plot_rad_1_y_2[i],ls=linest[i],label='KA = {}, KB = {}'.format(KAs[i],10.0),c='black')
        
    plt.xlabel("Interacion strenght")
    plt.ylabel(" % Coverage")
    plt.legend()
    plt.ylim(0.0,1.1)
    plt.show()

    plt.clf()
    # print(plot_rad_2_x_2)
    # print(plot_rad_2_y_2)
    for i in range(len(KAs)):
        plt.plot(plot_rad_2_x_3[i],plot_rad_2_y_3[i],ls=linest[i],label='KA = {}, KB = {}'.format(KAs[i],0.1),c='purple')
    for i in range(len(KAs)):
        plt.plot(plot_rad_2_x[i],plot_rad_2_y[i],ls=linest[i],label='KA = {}, KB = {}'.format(KAs[i],1.0),c='magenta')
    for i in range(len(KAs)):
        plt.plot(plot_rad_2_x_2[i],plot_rad_2_y_2[i],ls=linest[i],label='KA = {}, KB = {}'.format(KAs[i],10.0),c='black')
        
    plt.xlabel("Interacion strenght")
    plt.ylabel(" % Coverage")
    plt.ylim(0.0,1.1)
    plt.legend()
    plt.show()


# Plotting_coverage()
Plotting_coverage_varKB()



# I want to show the distribution of edgelengths before and after remeshing


def Edge_dist():


    file = open("../Results/Tests_remesher/Edge_dist_before.txt")

    line = file.readline()
    print(line)
    line = file.readline()
    Edges_before=[]
    Edges_after=[]
    while(line):
        Edges_before.append(float(line))
        line = file.readline()
    file.close()

    file = open("../Results/Tests_remesher/Edge_dist_after.txt")
    line = file.readline()
    print(line)
    line = file.readline()
    while(line):
        Edges_after.append(float(line))
        line = file.readline()
    file.close()

    avg_len = np.mean(Edges_before)

    plt.axvline(avg_len*0.5)
    plt.axvline(avg_len*3.0)
    # plt.hist(Edges_after,range=(min(Edges_before),max(Edges_before)),bins=100,label='After')
    # plt.show()
    plt.hist([Edges_before,Edges_after],range=(min(Edges_after),max(Edges_after)),bins=100,label=['Before, {} edges'.format(len(Edges_before)),"After, {} edges".format(len(Edges_after))])
    # plt.hist(Edges_after,range=(min(Edges_before),max(Edges_before)),bins=100,label='After')
    

    plt.legend()

    plt.show()




def Dih_dist():

    file = open("../Results/Tests_remesher/Dihedral_distribution_bad.txt")
    line = file.readline()

    print(line)

    line = file.readline()
    Dihs = []
    while(line):
        Dihs.append(float(line))
        line = file.readline()


    file.close()
    file = open("../Results/Tests_remesher/Dihedral_distribution_before.txt")
    line = file.readline()

    print(line)

    line = file.readline()
    Dihs2 = []
    while(line):
        Dihs2.append(float(line))
        line = file.readline()

    file.close()
    file = open("../Results/Tests_remesher/Dihedral_distribution_after.txt")
    line = file.readline()

    print(line)

    line = file.readline()
    Dihs3 = []
    while(line):
        Dihs3.append(float(line))
        line = file.readline()





    plt.hist([Dihs2,Dihs3],100,label=["Before","After"])
    plt.legend()
    plt.show()

# Edge_dist()
# Dih_dist()
