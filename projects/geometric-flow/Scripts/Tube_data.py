import numpy as np
import os 
import matplotlib.pyplot as plt 

import sys




Parent_folder="../Results/Mem3DG_Bead_Pulling_rc/"





def get_radius(file):
    file_1=Parent_folder+file
    file = open(file_1)

    line = file.readline()
    Y_pos=[]
    Z_pos=[]
    X_pos=[]
    while(line):
        split_line=line.split(" ")
        if(split_line[0]=="vn"):
            break 
        
        if(line[0]=="v"):
            X_pos.append(float(split_line[1]))
            Y_pos.append(float(split_line[2]))
            Z_pos.append(float(split_line[3]))
        # print(line)
        line=file.readline()
        

    # plt.scatter(Z_pos,Y_pos)
    # plt.show()
    avg_y=np.mean(Y_pos)
    avg_z=np.mean(Z_pos)

    r= np.mean(np.sqrt((np.array(Y_pos)-avg_y)**2+(np.array(Z_pos)-avg_z)**2))
    err=np.std(np.sqrt((np.array(Y_pos)-avg_y)**2+(np.array(Z_pos)-avg_z)**2))
    print(r)
    return [r,err] 


# filenames=["1_TUBE.obj","2_TUBE.obj","3_TUBE.obj","4_TUBE.obj","5_TUBE.obj","6_TUBE.obj"]
rs=[]
KBs=[0.01,0.04,0.09,0.16,0.25,0.36,0.49,0.64,0.81,1.0]
errs=[]
for i in range(1,11):
    filename = "{}_TUBE.obj".format(i)
    rs.append(get_radius(filename)[0])
    errs.append(get_radius(filename)[1])

rs=np.array(rs)
# plt.scatter(KBs,rs)
plt.errorbar(KBs,rs,yerr=errs,capsize=10,fmt='none',ecolor="purple")
plt.scatter(KBs,rs,color='black')

plt.xlabel("KB")
plt.ylabel("Microtubule radius")
plt.show()