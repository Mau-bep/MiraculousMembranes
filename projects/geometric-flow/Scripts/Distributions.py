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
    strength = 100
    init_cond=1 
    Nsim=12

    basic_name="../Results/Mem3DG_Bead_Reciprocal_arcsim_up/".format(nu,radius,KA,KB,strength,init_cond,Nsim)




    radius_vec= np.loadtxt(basic_name+"Radius_distribution_strength_"+"{:.6f}".format(strength)+".txt",delimiter=' ')

    print(min(radius_vec))
    print(max(radius_vec))
    plt.hist(radius_vec,bins=200)
    plt.show()
    plt.clf()

    

    Touching_data=np.loadtxt(basic_name+"Touching_strength_"+"{:.6f}".format(strength)+".txt")

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(Touching_data[:,0], Touching_data[:,1], Touching_data[:,2], marker='o')
    plt.axis('equal')
    plt.show()
    # plt.clf()

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


# checking_radius_dist()

def Ploting_coverage_2():
    file = open("../Results/Mem3DG_Bead_Reciprocal_arcsim_up_oct/Coverage_final.txt")
  
    line = file.readline()
    print(line)
    line = file.readline()

    strengths = []
    coverage = []

    while(line):
        splitted= line.split(' ')
        if(len(splitted)>4):
            break
        strengths.append(float(splitted[2]))
        coverage.append(float(splitted[3]))

        line = file.readline()
    plt.scatter(strengths,coverage,color='magenta',)
    plt.xlabel("Interaction strength")
    plt.ylabel("Coverage ")
    plt.show()


# Ploting_coverage_2()
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
# Plotting_coverage_varKB()







def Plotting_coverage_varKB_2():

    # So i will have 3 radius
    # For every radius there will be a curve with a different KA
    # All the KA curves are plotted as a function of Strenght

    file = open("../Results/Particle_wrapping_on_plane_var_jan/Coverage_final.txt")
    line = file.readline()
    print(line)
    line = file.readline()
    # I have 3 radius
    # For every radius i have 3 values of KB 
    # And there is the interaction strength constant
    # I think i should plot every radius differently


    # 
    # The thin is that for every radius there is one KB
    # [rad][KB]
    # I can make this 

    Coverages_rads= [[[],[],[]], [[],[],[]], [[],[],[]]]
    rads = [[[],[],[]], [[],[],[]], [[],[],[]]]

    Aproximate_coverage_rads = [[[],[],[]], [[],[],[]], [[],[],[]]]
    
    Interaction_strenghts_rads = [[[],[],[]], [[],[],[]], [[],[],[]]]

    # The radius is 0.2 0.3 and 0.4 
    # The Kbs are 10 20 30 
    while(line):
        splitted= line.split(' ')
        
        i = int(float(splitted[0])*10)-2
        j = int(float(splitted[1])/10)-1
        print(" {} and {}".format(i,j))

        Coverages_rads[i][j].append(float(splitted[3]))
        Aproximate_coverage_rads[i][j].append( abs(float(splitted[2])/(float(splitted[4])*(4.0)*np.pi*(float(splitted[0])*1.15 )**2  )))
        Interaction_strenghts_rads[i][j].append( abs(float(splitted[4])))
        


        line =file.readline()
        # I need to decide the plots


    file.close()

    print(Coverages_rads)

    for i in range(3):
        # For every radius
        for j in range(3):
            # For every KB 

            plt.scatter(Interaction_strenghts_rads[i][j],Coverages_rads[i][j],label = "KB = {}".format( (j+1)*10))

            # plt.xlabel(r"$\frac{E_I}{A_{B}K_{I}}$", usetex=True, fontsize = 30)
            plt.xlabel(r"${K_I}$",usetex = True, fontsize = 25)
            plt.ylabel(r" Coverage",usetex = True, fontsize = 25)
            plt.legend()
            plt.ylim(0.0,1.1)

        plt.title("Radius = {}".format((i+2)/10))
        plt.show()
    
        plt.clf()
    
    # plt.show()


    id = np.linspace(0,1,2)
    plt.plot(id,id,ls='dashed',color='black')
    
    for i in range(3):
        # For every radius
        for j in range(3):
            # For every KB 

            plt.scatter(Aproximate_coverage_rads[i][j],Coverages_rads[i][j],label = "KB = {}".format( (j+1)*10))

            plt.xlabel(r"$\frac{E_I}{A_{B}K_{I}}$", usetex=True, fontsize = 30)
            # plt.xlabel(r"$\cfrac{E_I}{AK_{I}}$",usetex = True)
            plt.ylabel(r" Coverage",usetex=True, fontsize = 30)
            plt.legend()
            plt.ylim(0.0,1.1)

        plt.title("Radius = {}".format((i+2)/10))
        # plt.show()
    
        # plt.clf()
    
    # plt.show()



    plt.show()



# I want to show the distribution of edgelengths before and after remeshing

# Plotting_coverage_varKB_2()





def Plotting_phase_space():

    # So i will have 3 radius
    # For every radius there will be a curve with a different KA
    # All the KA curves are plotted as a function of Strenght

    file = open("../Results/Particle_wrapping_on_plane_var_jan/Coverage_final.txt")
    line = file.readline()
    print(line)
    line = file.readline()
    # I have 3 radius
    # For every radius i have 3 values of KB 
    # And there is the interaction strength constant
    # I think i should plot every radius differently


    # 
    # The thin is that for every radius there is one KB
    # [rad][KB]
    # I can make this 
    
    # The variables i have are the radius, the KB and the interaction strength
    # The i can do one plot that is on the x axis the interaction strength 
    # And in the y axis the particle size 

    # I COULD ALSO try to recreate the plot that is being shown 


    # 


    # What will i do then. well i will put the points in their respective area and color them respective to their coverage ofc

    # So there will

    
    x_plot = []
    y_plot = []
    c_plot = []


    Coverages_rads= [[[],[],[]], [[],[],[]], [[],[],[]]]
    rads = [[[],[],[]], [[],[],[]], [[],[],[]]]

    Aproximate_coverage_rads = [[[],[],[]], [[],[],[]], [[],[],[]]]
    
    Interaction_strenghts_rads = [[[],[],[]], [[],[],[]], [[],[],[]]]

    # The radius is 0.2 0.3 and 0.4 
    # The Kbs are 10 20 30 
    while(line):
        splitted= line.split(' ')
        
        i = int(float(splitted[0])*10)-2
        j = int(float(splitted[1])/10)-1
        print(" {} and {}".format(i,j))

        w = float(splitted[4])/10
        a = float(splitted[0])*1.15
        KB = float(splitted[1])/4.0
        
        cov = np.clip(float(splitted[3]),0.0,1.0) 
        wc = KB/(a*a)


        print("w is {} and wc is {}".format(w,wc))
        x_plot.append(w/wc)
        y_plot.append( 0.05*100*a*a/(KB))
        c_plot.append(cov)


        line =file.readline()
        # I need to decide the plots


    file.close()

    ws = np.linspace(0,10.0,10)
    # plt.plot(ws,2*ws-2,color= "magenta")
    # plt.plot(ws,3*ws/5-3/5)
    # plt.xlim(0,4)
    # plt.ylim(0,0.4)
    plt.axvline(1.0,ls = 'dashed',color='black')
    plt.scatter(x_plot,y_plot,c=c_plot)
    # plt.xlabel(r"$\frac{\omega}{\omega_c}$", usetex=True, fontsize = 30)

    # plt.ylabel(r"$(\omega - \omega_c)\frac{a^2}{\kappa}",usetex=True, fontsize = 30)
    plt.legend()
    
    


    plt.show()


Plotting_phase_space()


def checking_radius_dist_var(i):

    nu=1.0
    radius=1.0
    curvadap=0.0
    minrel=0.1
    KA=100000.0
    KB=1.0
    strength = 1800
    init_cond=1 
    Nsim=12

    basic_name="../Results/Particle_wrapping_on_plane_var/"



    radius_vec= np.loadtxt(basic_name+"Radius_{}_distribution_strength_".format(i)+"{:.6f}".format(strength)+".txt",delimiter=' ')

    print(min(radius_vec))
    print(max(radius_vec))
    plt.hist(radius_vec,bins=200)
    plt.show()
    plt.clf()

    


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



    return


# checking_radius_dist_var(0)

# checking_radius_dist_var(1)
# checking_radius_dist_var(2)









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
