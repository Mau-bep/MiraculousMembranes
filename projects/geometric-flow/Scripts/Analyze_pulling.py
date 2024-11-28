import numpy as np 
import matplotlib.pyplot as plt 
import os


Strengths = [0.1, 0.25,0.5,0.75,1.0,1.25,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,5.75,5.0,5.25,5.5,5.75,6.0]
# Learn to do it with linspace
KB = 1.0
KA = 100000
Nsim = 1
Init_cond = 2

strength = 0.7

base ="../Results/Mem3DG_Bead_pulling_relaxation_arcsim/"

def measure_force():
    dLs = []
    Force = []
    Nsims = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 50, 55, 60, 65, 70,75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170,175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270 ]
    Nsims = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70,75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170,175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270 ]
    
    for Nsim in Nsims:
        folder = "../Results/Mem3DG_Bead_pulling_relaxation_arcsim/nu_1.000_radius_0.200_KA_100000.000_KB_{:.6f}_strength_{:.6f}_Init_cond_4_Nsim_{}/".format(KB,strength,Nsim)
        filename_bead = folder + "Bead_0_data.txt"
        f = open(filename_bead)
        Lines = f.readlines()
        line_first_bead = Lines[-1]
        line_first_bead = line_first_bead.split(" ")
        F_bead = np.array([line_first_bead[3], line_first_bead[4], line_first_bead[5] ],dtype=float)
        f.close()

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
        F_tot = F_bead
        F_tot = np.sqrt(np.sum(F_tot*F_tot))
        # print("THen this value turns into {}".format(F_tot))
        
        dLs.append(dL)
        
        Force.append(F_tot)
    plt.xlim()
    plt.scatter(dLs[1:],Force[1:],c="purple")
    plt.xlabel(r"$\Delta L$",fontsize=15) 
    plt.ylabel(r"Force",fontsize=15)
    plt.savefig(base+"Mem_resistance_pulling_9.png",bbox_inches = 'tight')
    plt.clf()
    plt.scatter(range(1,len(dLs)),dLs[1:],c="purple")
    plt.xlabel(r"Frame",fontsize=15)
    plt.ylabel(r"$\Delta L$",fontsize=15)
    plt.savefig(base+"Extensions_9.png",bbox_inches = 'tight')
    plt.clf()
    




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
            if(x>1.4 and x<2.3):
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
    plt.clf()
    plt.close()
    # plt.xlim(-0.1,0.1)
    # plt.ylim(-0.1,0.1)
    
    plt.title("Kb = {:.2f}".format(Kb))
    
    
    plt.scatter(Data[:,1],Data[:,2])
    plt.axis('equal')
    plt.savefig(folder+"Radius_plot_{:.1f}.png".format(Kb),bbox_inches = 'tight')
    plt.show()


    return


def Tube_growth_radius(folder,Kb):
    Data = np.loadtxt(folder+"Vertex_tube_KB_{:.6f}.txt".format(Kb))

    y = Data[:,1]
    z = Data[:,2]
    # Now i need to find the radius of this circle.
    avg_y=np.mean(y)
    avg_z=np.mean(z)
    plt.title("Kb = {}".format(Kb))
    r= np.mean(np.sqrt((np.array(y)-avg_y)**2+(np.array(z)-avg_z)**2))
    err=np.std(np.sqrt((np.array(y)-avg_y)**2+(np.array(z)-avg_z)**2))

    print("The radius is {} with an error of {}".format(r,err))

    return r

folder_path_growth = "../Results/Mem3DG_Bead_pulling_oct_growth_arcsim/"


def fit():
    Strengths =[1.0, 2.0, 3.0, 4.0, 5.0, 9.0, 10.0, 12.0, 15.0, 16.0, 20.0, 22.0, 25.0, 28.0, 32.0, 36.0]
    radius = []
    for strength in Strengths:
        # Tube_growth_data(folder_path_growth,strength)


        # Tube_growth_check(folder_path_growth, strength)

        r = Tube_growth_radius(folder_path_growth, strength)
        radius.append(r)


    plt.clf()
    plt.xlabel("Kb",fontsize=15.0)
    plt.ylabel("Tube radius",fontsize=15.0)
    plt.scatter(Strengths,radius,color="black")
    plt.savefig(folder_path_growth+"NoFit_radius_curve.png",bbox_inches='tight')
    plt.show()


    # ok so i want to do a linear fit to the thing

    x = np.log(Strengths)
    y = np.log(radius)


    p = np.polyfit(x,y,1)

    x_fit = np.linspace(np.log(1.0),np.log(36))
    y_fit = p[1]+x_fit*p[0]
    plt.plot(x_fit,y_fit,ls='dashed')
    plt.scatter(x,y,color='black')
    print("The values for the fit are m = {} and c = {}".format(p[0],p[1]))
    plt.show()

    y_fit = np.exp(y_fit)
    x_fit = np.exp(x_fit)

    plt.plot(x_fit,y_fit,ls='dashed',color='magenta')
    plt.scatter(Strengths,radius,color='black')
    plt.savefig(folder_path_growth+"Fit_radius_curve.png",bbox_inches='tight')
    plt.xlabel("Kb",fontsize=15.0)
    plt.ylabel("Tube radius",fontsize=15.0)
    plt.show()





# measure_force()



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



def Tube_growth_radius_2(folder):
    Data = np.loadtxt(folder+"Tube_radius.txt")

    y = Data[:,1]
    z = Data[:,3]
    # Now i need to find the radius of this circle.
    # avg_y=np.mean(y)
    # avg_z=np.mean(z)
    # plt.title("Kb = {}".format(Kb))
    # r= np.mean(np.sqrt((np.array(y)-avg_y)**2+(np.array(z)-avg_z)**2))
    # err=np.std(np.sqrt((np.array(y)-avg_y)**2+(np.array(z)-avg_z)**2))

    # print("The radius is {} with an error of {}".format(r,err))

    return [y,z]


def fit2():
    # Strengths =[ 6.0, 10.0, 14.0, 18.0, 22.0, 26.0, 30.0, 34.0, 38.0, 42.0, 46.0, 50.0]
    radius = []
    
    folder_path_growth = "../Results/Mem3DG_Bead_pulling_radius_growth_arcsim/"
    [Strengths,radius] = Tube_growth_radius_2(folder_path_growth)
        # radius.append(r)



    plt.clf()
    plt.xlabel("Kb",fontsize=15.0)
    plt.ylabel("Tube radius",fontsize=15.0)
    plt.scatter(Strengths,radius,color="black")
    plt.savefig(folder_path_growth+"NoFit_radius_curve.png",bbox_inches='tight')
    plt.show()


    # ok so i want to do a linear fit to the thing

    x = np.log(Strengths)
    y = np.log(radius)


    p = np.polyfit(x,y,1)

    x_fit = np.linspace(np.log(1.0),np.log(50.0))
    y_fit = p[1]+x_fit*p[0]


    
    plt.plot(x_fit,y_fit,ls='dashed')
    plt.scatter(x,y,color='black')
    print("The values for the fit are m = {} and c = {}".format(p[0],p[1]))
    plt.show()

    y_fit = np.exp(y_fit)
    x_fit = np.exp(x_fit)

    y_fit2 = np.sqrt( x_fit/(4*0.05))/100.0

    plt.plot(x_fit,y_fit,ls='dashed',color='magenta')
    plt.plot(x_fit,y_fit2,ls='dashed',color='purple')
    plt.scatter(Strengths,radius,color='black')
    plt.savefig(folder_path_growth+"Fit_radius_curve.png",bbox_inches='tight')
    plt.xlabel("Kb",fontsize=15.0)
    plt.ylabel("Tube radius",fontsize=15.0)
    plt.show()



fit2()