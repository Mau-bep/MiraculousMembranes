import numpy as np 
import matplotlib.pyplot as plt 
import os


Strengths = [0.1, 0.25,0.5,0.75,1.0,1.25,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,5.75,5.0,5.25,5.5,5.75,6.0]
# Learn to do it with linspace
KB = 1.0
KA = 100000
Nsim = 1
Init_cond = 2



base ="../Results/Mem3DG_Bead_pulling_up_oct_arcsim/"
def main():
    dLs = []
    Force = []

    for strength in Strengths:
        folder = "../Results/Mem3DG_Bead_pulling_up_oct_arcsim/nu_1.000_radius_0.200_KA_100000.000_KB_{:.6f}_strength_{:.6f}_Init_cond_2_Nsim_1/".format(KB,strength)
        # I need to open 3 files bead fixed, bead moving and Last frame.

        # 
        filename_bead_moving = folder+ "Bead_0_data.txt"
        f = open(filename_bead_moving)
        Lines = f.readlines()
        line_first_bead = Lines[-1]
        line_first_bead = line_first_bead.split(" ")
        F_moving = np.array([line_first_bead[3], line_first_bead[4], line_first_bead[5] ],dtype=float)

        f.close()

        # Ok so this one has the position and the 
        filename_bead_fixed = folder + "Bead_1_data.txt"
        f = open(filename_bead_fixed)
        Lines = f.readlines()
        line_second_bead = Lines[-1]
        line_second_bead = line_second_bead.split(" ")
        F_fixed = np.array([line_second_bead[3], line_second_bead[4], line_second_bead[5] ],dtype=float)

        f.close()

        # Now i need to find the last  

        files = os.listdir(folder)
        files = [f for f in files if os.path.isfile(folder+'/'+f)]
        higher_index = 0

        for file in files:
            # I want to find the highest number
            line = file.split("_")
            if(line[0][0]=="m"):
                line = line[1].split(".")
                index = int(line[0])
                if(index>higher_index):
                    higher_index = index

        # print("The highest index is {}".format(higher_index))
        [low,high] = find_higher_lowest( folder+"membrane_{}.obj".format(higher_index))
        dL = high-low
        F_tot = F_fixed+ F_moving
        F_tot = np.sqrt(np.sum(F_tot*F_tot))
        # print("THen this value turns into {}".format(F_tot))
        
        dLs.append(dL)
        
        Force.append(F_tot)
        # print("DL is {} and Force is {}".format(dL,F_tot))

    plt.scatter(dLs,Force,c="black")
    




def find_higher_lowest(filename):
    f = open(filename)
    line = f.readline()

    line = f.readline()
    lowest = 0.0
    highest = 0.0
    while(line):

        if(line[0]=="v"):
            x = line.split(" ")[1]
        
            x = float(x)

            if(x > highest):
                highest = x
            if(x < lowest):
                lowest = x

        if(line[0]=="e"):
            break
        line = f.readline()

    return [lowest,highest]

        


def intermediate_points(strength):

    folder = "../Results/Mem3DG_Bead_pulling_up_oct_arcsim/nu_1.000_radius_0.200_KA_100000.000_KB_{:.6f}_strength_{:.6f}_Init_cond_2_Nsim_1/".format(KB,strength)
    # I have the folder i now need to iterate

    dLs = []
    Forces = []

    files = os.listdir(folder)
    files = [f for f in files if os.path.isfile(folder+'/'+f)]
    higher_index = 0

    for file in files:
        # I want to find the highest number
        line = file.split("_")
        if(line[0][0]=="m"):
            line = line[1].split(".")
            index = int(line[0])
            if(index>higher_index):
                higher_index = index

    # I have the higher index

    Bead_moving_file = open(folder+ "Bead_0_data.txt")
    Bead_fixed_file = open(folder+"Bead_1_data.txt")

    line1 = Bead_moving_file.readline()
    line2 = Bead_fixed_file.readline()

    for i in range(0,higher_index,200):
        line1 = Bead_moving_file.readline()
        line2 = Bead_fixed_file.readline()
        [low,high] = find_higher_lowest(folder+"membrane_{}.obj".format(i))
        # if(high-low <2.4):
        #     # print("Point {} is not streched enough\n".format(i))
        #     # In this case
        #     continue
        
        
        # So its an extended tube.
        #  
        line_moving_bead = line1.split(" ")
        F_moving = np.array([line_moving_bead[3], line_moving_bead[4], line_moving_bead
        [5] ],dtype=float)
        line_fixed_bead = line2.split(" ")
        F_fixed = np.array([line_fixed_bead[3], line_fixed_bead[4], line_fixed_bead
        [5] ],dtype=float)


        dL = high-low
        F_tot = F_fixed+ F_moving
        F_tot = np.sqrt(np.sum(F_tot*F_tot))

        dLs.append(dL)
        Forces.append(F_tot)
    
    plt.plot(dLs,Forces,color=cmap(strength/6.0),label=str(strength))




    # 

def Get_higher_index(folder):
    files = os.listdir(folder)
    files = [f for f in files if os.path.isfile(folder+'/'+f)]
    higher_index = 0

    for file in files:
        # I want to find the highest number
        line = file.split("_")
        if(line[0][0]=="m"):
            line = line[1].split(".")
            index = int(line[0])
            if(index>higher_index):
                higher_index = index

    return higher_index


def Get_files(strength):

    folder = "../Results/Mem3DG_Bead_pulling_up_oct_arcsim/nu_1.000_radius_0.200_KA_100000.000_KB_{:.6f}_strength_{:.6f}_Init_cond_2_Nsim_1/".format(KB,strength)
    # I have the folder i now need to iterate

    dLs = []
    Forces = []
    
    higher_index = Get_higher_index(folder)

    # I have the higher index

    Bead_moving_file = open(folder+ "Bead_0_data.txt")
    Bead_fixed_file = open(folder+"Bead_1_data.txt")

    line1 = Bead_moving_file.readline()
    line2 = Bead_fixed_file.readline()
    Bead_positions = []
    for i in range(0,higher_index,200):
        line1 = Bead_moving_file.readline()
        line2 = Bead_fixed_file.readline()
        [low,high] = find_higher_lowest(folder+"membrane_{}.obj".format(i))
        # if(high-low <2.4):
        #     # print("Point {} is not streched enough\n".format(i))
        #     # In this case
        #     continue
        
        splitted_line = Bead_fixed_file.split(" ")
        Bead_positions.append( [float(splitted_line[0]),float(splitted_line[1]), float(splitted_line[2])])
        
        
        # So its an extended tube.
        #  
        line_moving_bead = line1.split(" ")
        F_moving = np.array([line_moving_bead[3], line_moving_bead[4], line_moving_bead
        [5] ],dtype=float)
        line_fixed_bead = line2.split(" ")
        F_fixed = np.array([line_fixed_bead[3], line_fixed_bead[4], line_fixed_bead
        [5] ],dtype=float)


        dL = high-low
        F_tot = F_fixed+ F_moving
        F_tot = np.sqrt(np.sum(F_tot*F_tot))

        dLs.append(dL)
        Forces.append(F_tot)
    
    # plt.plot(dLs,Forces,color=cmap(strength/6.0),label=str(strength))










def Tube_growth_data(folder,Kb):
    # So i 
   
    
    # So i have the file
    # The file is folder+"membrane_{}.obj".format(higher_index)
    # We need to read the file and write one that contains only the vertex position data.
    directory = folder+"nu_1.000_radius_0.200_KA_100000.000_KB_{:.6f}_strength_0.000000_Init_cond_5_Nsim_1/".format(Kb)
    higher_index = Get_higher_index(directory)
    print("The higher index is {}".format(higher_index))
    membrane = open(folder+"nu_1.000_radius_0.200_KA_100000.000_KB_{:.6f}_strength_0.000000_Init_cond_5_Nsim_1/".format(Kb)+"membrane_{}.obj".format(higher_index),"r")
    line = membrane.readline()
    # Ok i need to write in another file right away right?

    Output = open(folder+"Vertex_tube_KB_{:.6f}.txt".format(Kb),"w")

    line = membrane.readline()
    while(line):
        if(line[0]=='v'):
            # This is a vertex line
            # Now i just want the vertices with certain values 
            x = float(line.split(" ")[1])
            if(x>1.3 and x<2.3):
                Output.write(line[2:])
            line = membrane.readline()
            continue
        break


    membrane.close()
    Output.close()

    return 


def Tube_growth_check(folder,Kb):
    # Now its easier cause i have the data
    Data = np.loadtxt(folder+"Vertex_tube_KB_{:.6f}.txt".format(Kb))
    # print(Data[:,0])
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.axis('equal')
    ax.set(xlim =(1.0, 2.5),ylim = (-0.5,0.5),zlim =(-0.5,0.5))
    ax.scatter(Data[:,0],Data[:,1],Data[:,2],color='purple')
    plt.show()
    plt.xlim(-0.1,0.1)
    plt.ylim(-0.1,0.1)
    plt.scatter(Data[:,1],Data[:,2])
    plt.show()


    return


def Tube_growth_radius(folder,Kb):
    Data = np.loadtxt(folder+"Vertex_tube_KB_{:.6f}.txt".format(Kb))

    y = Data[:,1]
    z = Data[:,2]
    # Now i need to find the radius of this circle.
    avg_y=np.mean(y)
    avg_z=np.mean(z)

    r= np.mean(np.sqrt((np.array(y)-avg_y)**2+(np.array(z)-avg_z)**2))
    err=np.std(np.sqrt((np.array(y)-avg_y)**2+(np.array(z)-avg_z)**2))

    print("The radius is {} with an error of {}".format(r,err))

    return r

folder_path_growth = "../Results/Mem3DG_Bead_pulling_oct_growth_arcsim/"


Strengths =[1.0, 2.0, 3.0, 4.0, 5.0, 9.0, 10.0, 12.0, 15.0, 16.0, 20.0, 22.0, 25.0, 28.0, 32.0, 36.0]
for strength in Strengths:
    Tube_growth_data(folder_path_growth,5.0)


# Tube_growth_check(folder_path_growth,5.0)

# Tube_growth_radius(folder_path_growth,5.0)




# main()
# cmap = plt.colormaps['viridis']
# intermediate_points(3.5)
# intermediate_points(4.0)
# intermediate_points(4.5)
# intermediate_points(5.0)
# intermediate_points(5.5)
# intermediate_points(6.0)

# plt.legend()
# plt.xlabel(r"$\Delta L$")
# plt.ylabel(r"Force")
# plt.savefig(base+"Pulling_force_plot.jpg",bbox_inches = 'tight')


