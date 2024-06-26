import numpy as np
import os 
import matplotlib.pyplot as plt 

import sys



# nu=float(sys.argv[1])
# KB=float(sys.argv[2])
# strg=float(sys.argv[3])

nu=1.0
KB=0.005
def main(strg,Init_cond,Nsim):
    

    pre_folder='../Results/Mem3DG_Bead_Reciprocal/nu_1.000_c0_0.000_KA_10.000_KB_0.005000_strength_{:.6f}_Init_cond_{}_Nsim_{}/'.format(strg,Init_cond,Nsim)
                # ../Results/Mem3DG_Bead_Reciprocal/nu_1.000_c0_0.000_KA_10.000_KB_0.005000_strength_0.100000_Init_cond_1_Nsim_1/
    # dir='nu_{:.3f}_c0_0.000_KA_10.000_KB_{:.6f}_Inter_{:.6f}'.format(nu,KB,strg)
    # folder_path= pre_folder+dir+'/'
    folder_path = pre_folder
    Bead_path=folder_path+"Bead_data.txt"
    print(Bead_path)
    file_xyz = open(folder_path+'Bead_data.xyz','w')



    Bead_file = open(Bead_path)
    line = Bead_file.readline()
    r = 1.0

    counter=0
    frame=0
    line=Bead_file.readline()
    while line:
        frame=counter*1000
        splitted_line=line.split(' ')
        if(len(splitted_line)<3):
            break
        Positions=np.array( [splitted_line[0],splitted_line[1],splitted_line[2]],dtype=float)

        file_xyz.write('1\n')
        file_xyz.write('{}\n'.format(frame))
        file_xyz.write('B {} {} {}\n'.format(Positions[0],Positions[1],Positions[2]))

        line=Bead_file.readline()
        counter+=1

    # OK so from here we have the data



    Bead_file.close()
    file_xyz.close()


# for strg in [0.002, 0.6]:
for strg in [0.1,1.0,0.01, 0.001,10.0]:
    main(strg,1,1)
    main(strg,2,1)

        