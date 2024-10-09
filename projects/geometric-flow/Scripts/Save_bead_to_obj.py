import numpy as np 
import os 
import sys 


# bead_file = open("Bead_data.xyz")

# timestep= bead_file.readline()

strg = float(sys.argv[1])
Init_cond = int(sys.argv[2])
Nsim = int(sys.argv[3])
KA = 100000



def create_bead(Center,radius,frame):

    prev_rad=1.0
    bead_dir=pre_folder+"/Bead_traj/"
    os.makedirs(bead_dir,exist_ok=True)

    New_obj=open(bead_dir+"bead_{}.obj".format(frame),'+w')

    New_obj.write("##This is a mesh file of the {} frame\n".format(frame))

    # Now i need to read the sphere obj and write those lines

    # Sphere_obj=open("../../../input/basic_sphere.obj",'r')

    Sphere_obj=open("../../../input/sphere.obj",'r')
    line=Sphere_obj.readline()
    # print(line)
    line=Sphere_obj.readline()
    # print(line)
    while(line):
        if(line[0]=='#'):
            line=Sphere_obj.readline()
            continue

        # Now i need to write this in the other
        if(line[0]=='v'):
            Data = line.split(' ')
            
            if(len(Data)<4):
                break
            # print(np.array(Data[1:],dtype=float))
            # print(Center)
            Pos = np.array(Data[1:],dtype=float)*radius/prev_rad+Center
            New_obj.write("v {} {} {}\n".format(Pos[0],Pos[1],Pos[2]))

        if(line[0]=='f'):
            New_obj.write(line)
                    
        

        line=Sphere_obj.readline()




    Sphere_obj.close()
    New_obj.close()

    return 




dz=np.array([0.1,0.1,0.1])


radius=1.1
Center=np.array([0,0,0])


pre_folder='../Results/Mem3DG_Bead_Reciprocal_arcsim_up_oct/nu_1.000_radius_1.000_KA_{:.3f}_KB_1.000000_strength_{:.6f}_Init_cond_{}_Nsim_{}/'.format(KA,strg,Init_cond,Nsim)

Bead_orig_data=open(pre_folder+'Bead_0_data.txt')
line_bead=Bead_orig_data.readline()
print(line_bead)
line_bead=Bead_orig_data.readline()
i=0
while(line_bead):
    # I have the data here
    # print(line_bead.split(' '))
    Center=np.array(line_bead.split(' ')[:-4],dtype=float)
    create_bead(Center,radius,i)
    i+=1

    line_bead=Bead_orig_data.readline()
    
# for i in range(100):


#     create_bead(Center+dz*i,radius,i)


