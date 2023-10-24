import numpy as np 
import os 
import sys 


# bead_file = open("Bead_data.xyz")

# timestep= bead_file.readline()



def create_bead(Center,radius,frame):

    New_obj=open("bead_{}.obj".format(frame),'+a')

    New_obj.write("##This is a mesh file of the {} frame".format(frame))

    # Now i need to read the sphere obj and write those lines

    Sphere_obj=open("Simple_sphere_bead.obj",'r')
    line=Sphere_obj.readline()
    line=Sphere_obj.readline()
    while(line):
        # Now i need to write this in the other
        Data = line.split(' ')
        Pos = np.array(Data[1:],dtype=float)*radius+Center

        New_obj.write("v {} {} {}\n".format(Pos[0],Pos[1],Pos[2]))


        line=Sphere_obj.readline()




    Sphere_obj.close()
    New_obj.close()

    return 





radius=3.0
Center=np.array([0,0,0])




