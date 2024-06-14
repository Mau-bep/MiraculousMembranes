import numpy as np
import matplotlib.pyplot as plt



# i WANT THE DIRECTORY


def checking_radius_dist():

    nu=1.0
    radius=1.0
    curvadap=0.0
    minrel=0.1
    KA=100.0
    KB=0.1
    strength = 10.0
    init_cond=1 
    Nsim=12

    basic_name="../Results/Mem3DG_Bead_Reciprocal_finemesh/nu_{:.3f}_radius_{:.3f}_curvadap_{:.2f}_minrel_{:.4f}_KA_{:.3f}_KB_{:.6f}_strength_{:.6f}_Init_cond_{}_Nsim_{}/".format(nu,radius,curvadap,minrel,KA,KB,strength,init_cond,Nsim)




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


Plotting_coverage()
