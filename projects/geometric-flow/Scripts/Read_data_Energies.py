import numpy as np
import matplotlib.pyplot as plt






# We do this function and we can go home




def Create_data_file(folder_dir,KA,KB,KI, theta):
    dir = folder_dir + "Bending_{0:.4f}_Surface_tension_{1:.4f}_Edge_reg_1.0000_Bead_radius_0.3000_Frenkel_Normal_nopush_str_{2:.4f}_theta_const_{3:.4f}_Bead_radius_0.3000_Frenkel_Normal_nopush_str_{2:.4f}_theta_const_{3:.4f}_Switch_Newton_Switch_Free_beads_Nsim_1/".format(KB,KA,KI,theta)


    # Gorgeous i have the precise ffolder

    file_E = open(folder_dir+  "Energy_output_KA_{:.4f}_KB_{:.4f}_KI_{:.4f}_theta_{:.4f}.txt".format(KA,KB,KI,theta),"w+")
    output_data = open(dir+ "Output_data.txt","r")
    line = output_data.readline()
    # The first line doesnt contain anything
    line = output_data.readline()

    while(line):
        splitted_line = line.split(' ')
        if(len(splitted_line)>4):
            time = float(splitted_line[0])
            Ebend = float(splitted_line[3])
            Esurf = float(splitted_line[4])
            
            Eedge = float(splitted_line[5])
            Ebead1 = float(splitted_line[6])
            Ebead2 = float(splitted_line[7])
            Etot = float(splitted_line[8])
                # Now we have everything we need
            file_E.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f}\n".format(time,Ebend,Esurf,Eedge,Ebead1,Ebead2, Etot))
        line = output_data.readline()
    file_E.close()
    # That is all girl


def Read_data_file(folder_dir,KA,KB,KI,theta):
    file_E = open(folder_dir+  "Energy_output_KA_{:.4f}_KB_{:.4f}_KI_{:.4f}_theta_{:.4f}.txt".format(KA,KB,KI,theta),"r")
    
    data = np.loadtxt(folder_dir+  "Energy_output_KA_{:.4f}_KB_{:.4f}_KI_{:.4f}_theta_{:.4f}.txt".format(KA,KB,KI,theta))

    print("The shape of the data is {}\n".format(data.shape))
    if data.shape[0]<50:
        print("Not enough data points")
        return
    time = data[:,0]

    Ebend = data[:,1]
    Esurf = data[:,2]
    Eedge = data[:,3]
    Ebead1 = data[:,4]
    Ebead2 = data[:,5]
    Etot = data[:,6]
    plt.figure(figsize=(10,8))
    plt.axvline(100000/2000,color='black', linestyle='--')
    plt.axvline(300000/2000,color='black', linestyle='--')
    plt.plot(Etot,label="Total energy", color='black')
    plt.plot(Ebend,label="Bending energy", color='blue')
    plt.plot(Esurf,label="Surface energy", color='red')
    plt.plot(Eedge,label="Edge energy", color='green')
    plt.plot(Ebead1,label="Bead 1 energy", color='orange')
    plt.plot(Ebead2,label="Bead 2 energy", color='purple')
    # plt.yscale('log')
    plt.xlabel("Time")
    plt.ylabel(r"\text{Energy}", fontfamily='serif')
    plt.title("Energies vs time for KA={:.2f}, KB={:.2f}, KI={:.2f}, theta={:.2f}".format(KA,KB,KI,theta))
    plt.legend()
    plt.grid()
    plt.savefig(folder_dir+  "Energy_plot_KA_{:.4f}_KB_{:.4f}_KI_{:.4f}_theta_{:.4f}.png".format(KA,KB,KI,theta))
    
    plt.show()
    plt.close()


folder = "../Results/Two_beads_with_switch_sept_dynamic/"

KA = 10
KB = 10

KI = 350

thetas = [0.3 +0.05*i for i in range(20)]

for i in thetas:
    print("Doing theta = {}".format(i))
    # Create_data_file(folder,KA,KB,KI,i)




for i in thetas:
    print("Doing theta = {}".format(i))
    Read_data_file(folder,KA,KB,KI,i)



    