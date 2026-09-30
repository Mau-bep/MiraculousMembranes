import numpy as np
import matplotlib.pyplot as plt 
import os 
import json
def main(folderpath,Nsims):
    
    dx = []
    fxs = []

    for Nsim in Nsims:

        Sim_path = folderpath + "Surface_tension_0.0500_Bending_10.0000_Bead_radius_0.2000_str_10.0000_Nsim_{}/".format(Nsim)

        Bead_data = Sim_path+"Bead_0_data.txt"

        data = np.loadtxt(Bead_data,skiprows=1)
        # Ok so i have the bead data (I kinda care only about the force in the x direction at the last step)

        xpos = data[-1,0]+0.2
        fx = np.abs(data[-1,3])
        dx.append(xpos)
        fxs.append(fx)
    
    
    # I should have all the data that i need at this point

    f = open(folderpath + "Tube_relaxation_data.txt",'w')

    for i in range(len(Nsims)):
        f.write("{} {} \n".format(dx[i],fxs[i]))
    
    f.close()


    plt.scatter(dx,fxs,color='black')
    plt.xlabel("dx")
    plt.ylabel("fx")
    plt.savefig(folderpath + "Tube_relaxation_data.png",bbox_inches='tight')
    return 


Nsims = [i for i in range(0,340,20)] + [i for i in range(0,20,2)] + [30] + [i for i in range(180,220,2)]
Nsims.sort()
Nsims = np.unique(Nsims)
# print(Nsims)



# main(folderpath="../Results/Tube_for_relaxation/",Nsims=Nsims)

def Read_data(folderpath):
    dirs = os.listdir(folderpath)
    f = open(folderpath+"Force_data.csv","w+")
    for dir in dirs:
        # print(dir)
        split_dir = dir.split(".")
        if( len(split_dir) == 1):
            # Tihs is where the directory is a simulation
            # So now we need to read the file
            Bead_file = folderpath+dir+"/Bead_0_data.txt"
            # print(Bead_file)
            Bead_f = np.loadtxt(Bead_file,skiprows = 1)
            if(len(Bead_f)==7):
                print(dir)
                continue
            # print(Bead_f[-1]
            lastrow = Bead_f[-1]
            # I need the ka and the r 
            with open(folderpath+dir+"/Input_file.json", 'r') as fjson:
                d = json.load(fjson)

            # print(d)
            # print(d["Energies"])
            # print(type(d))
            # res = d.replace("'",'"')
            # print(res)
        
            for Energie in d["Energies"]:
                if(Energie["Name"]=="Surface_tension"):
                    ka = Energie["constants"][0]
            for Bead in d["Beads"]:
                if(Bead["state"]=="manual"):
                    r = Bead["radius"]
            # Entonces d tiene el diccionario completo
            


            f.write("{0} {1} {2} {3} {4} {5} {6}\n".format(dir,ka, r, lastrow[1],lastrow[4],lastrow[5],lastrow[6]))

    f.close()

# Read_data("../Results/Pulling_and_relaxing/")
def Read_data_continuation(folderpath):
    dirs = os.listdir(folderpath)
    f = open(folderpath+"Force_data.csv","w+")
    for dir in dirs:
        # print(dir)
        split_dir = dir.split(".")
        if( len(split_dir) == 1):
            # Tihs is where the directory is a simulation
            # So now we need to read the file
            Bead_file = folderpath+dir+"/Bead_0_data.txt"
            # print(Bead_file)
            Bead_f = np.loadtxt(Bead_file,skiprows = 1)
            if(len(Bead_f)==7):
                print(dir)
                continue
            # print(Bead_f[-1]
            lastrow = Bead_f[-1]
            # I need the ka and the r 
            with open(folderpath+dir+"/Input_file.json", 'r') as fjson:
                d = json.load(fjson)

            # print(d)
            # print(d["Energies"])
            # print(type(d))
            # res = d.replace("'",'"')
            # print(res)
        
            for Energie in d["Energies"]:
                if(Energie["Name"]=="Surface_tension"):
                    ka = Energie["constants"][0]
            for Bead in d["Beads"]:
                if(Bead["state"]=="manual"):
                    r = Bead["radius"]
            
            # I want the sign right so 
            sign = d["Beads"][0]["velocity"][0]
            # Entonces d tiene el diccionario completo
            
            print(sign)

            f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(dir,ka, r, lastrow[1],lastrow[4],lastrow[5],lastrow[6],sign))

    f.close()

def Plot_data(directory):
    # Ok so he datafile would be
    Data = np.loadtxt(directory+"Force_data.csv")
    Distance_x= []
    Force_x = []
    Distance_x_2 = []
    Force_x_2 = []
    for i in range(len(Data)):
        # Ok so thi is gonna go over whats in data
        # print(Data[i,0])
        if(Data[i,2]==0.3 ):
            # print(Data[i,3])
            if(Data[i,1] == 16):
                Distance_x.append(Data[i,3])
                Force_x.append(-1*Data[i,4])
            else:
                Distance_x_2.append(Data[i,3])
                Force_x_2.append(-1*Data[i,4])
                

    plt.axhline(2*np.pi*np.sqrt(16*4),ls='dashed',color='black')
    plt.axhline(2*np.pi*np.sqrt(32*4),ls='dashed',color='black')
    plt.scatter(Distance_x,Force_x)
    plt.xlabel("Bead Final Pos")
    plt.ylabel("Force due to membrane response")
    plt.axvline(2.1)
    plt.axvline(2.3)
    plt.scatter(Distance_x_2,Force_x_2)
    plt.ylim(0,150)
    # print(Distance_x)
    # print(Force_x)
    plt.savefig(directory + "Scatter_plot.png",bbox_inches='tight')
    # plt.show()



def Plot_data_continuation(directory):
    # Ok so he datafile would be
    Data = np.loadtxt(directory+"Force_data.csv")
    Distance_x= []
    Force_x = []

    Distance_x_2= []
    Force_x_2 = []
    
    for i in range(len(Data)):
        # Ok so thi is gonna go over whats in data
        # print(Data[i,0])
        if(Data[i,7]==0.5 ):
            # print(Data[i,7])
            Distance_x.append(Data[i,3])
            Force_x.append(-1*Data[i,4])
        if(Data[i,7]==-0.5 ):
                    # print(Data[i,7])
                    Distance_x_2.append(Data[i,3])
                    Force_x_2.append(-1*Data[i,4])
                
    plt.scatter(Distance_x,Force_x)
    plt.scatter(Distance_x_2,Force_x_2)
    
    plt.axhline(2*np.pi*np.sqrt(16*4),ls='dashed',color='black')
    plt.axhline(2*np.pi*np.sqrt(32*4),ls='dashed',color='black')
    plt.scatter(Distance_x,Force_x)
    plt.xlabel("Bead Final Pos")
    plt.ylabel("Force due to membrane response")
    # plt.axvline(2.1)
    # plt.axvline(2.3)
    plt.scatter(Distance_x_2,Force_x_2)
    # plt.ylim(0,150)
    # print(Distance_x)
    # print(Force_x)
    plt.savefig(directory + "Scatter_plot.png",bbox_inches='tight')
    # plt.show()




Read_data("../Results/Pulling_and_relaxing_Sept_continuation/")
Plot_data("../Results/Pulling_and_relaxing_Sept_continuation/")