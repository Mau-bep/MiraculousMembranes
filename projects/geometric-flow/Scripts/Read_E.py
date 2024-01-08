import numpy as np
import os 
import matplotlib.pyplot as plt 

import sys



nu=float(sys.argv[1])
# KB=float(sys.argv[2])
# strg=float(sys.argv[3])

# nu=1.0
KB=0.1
Initial_conds=int(sys.argv[2])

Nsim=1

# def main_bead(strg):


#     pre_folder='../Results/Mem3DG_Cell_Shape/'
#     dir='nu_{:.3f}_c0_0.000_KA_10.000_KB_{:.6f}_init_cond_{}_Nsim_{}'.format(nu,KB,Ini_cond,Nsim)
#         #   nu_1.000_c0_0.000_KA_10.000_KB_0.005000_init_cond_1_Nsim_2
#     # dir_2='nu_{:.3f}_c0_0.000_KA_10.000_KB_{:.6f}_init_cond_{}_Nsim_{}'.format(nu,KB,Ini,Nsim)
    
#     folder_path= pre_folder+dir+'/'

#     # Bead_path=folder_path+"Bead_data.txt"
#     Output_path=folder_path+"Output_data.txt"
#     # file_xyz = open(folder_path+'Bead_data.xyz','w')

#     # Bead_file = open(Bead_path)
#     Output_file= open(Output_path)
#     # line = Output_file.readline()
#     r = 1.0

#     counter=0
#     frame=0
#     line=Output_file.readline()

#     targ_vol=0.0
#     targ_area=0.0


#     # Sim_data << V_bar<<" "<< A_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " " << E_Ben <<" " << E_bead << " "<< grad_norm<<" " << backtrackstep<< " "<< Bead_1.Pos.x << " "<< Bead_1.Pos.y << " "<< Bead_1.Pos.z <<" \n";

#     time_evol=[]
#     Volumes=[]
#     Areas=[]
#     E_vol=[]
#     E_sur=[]
#     E_ben=[]
#     E_bead=[]
#     grad_norm=[]
#     backtrackstep=[]
#     line=Output_file.readline()
#     print(line)
#     while line:
#         # frame=counter*1000
#         splitted_line=line.split(' ')
#         if(len(splitted_line)<13):
#             break
#         time_evol.append(float(splitted_line[2]))
#         Volumes.append(float(splitted_line[3]))
#         Areas.append(float(splitted_line[4]))
#         E_vol.append(float(splitted_line[5]))
#         E_sur.append(float(splitted_line[6]))
#         E_ben.append(float(splitted_line[7]))
#         E_bead.append(float(splitted_line[8]))
#         grad_norm.append(float(splitted_line[9]))
#         backtrackstep.append(float(splitted_line[10]))

#         line=Output_file.readline()

#     V_bar=float(splitted_line[0])
#     A_bar=float(splitted_line[1])


#     # Ok so we have the data and its time to plot !

#     os.makedirs(pre_folder+"Imgs/",exist_ok=True)
    

#     plt.plot(time_evol,Volumes,color='purple')
#     plt.axhline(V_bar,color='black',ls='dashed') 
#     plt.title(r'Volume evolution',fontsize=22)
#     plt.xlabel(r't',fontsize=20)
#     plt.ylabel(r'V',fontsize=20)
#     plt.yscale('log')
#     plt.savefig(pre_folder+"/Imgs/Volume_str_{}.jpg".format(strg),bbox_inches='tight')
#     plt.clf()



#     plt.plot(time_evol,Areas,color='purple')
#     plt.axhline(A_bar,color='black',ls='dashed') 
#     plt.title(r'Area evolution',fontsize=22)
#     plt.xlabel(r't',fontsize=20)
#     plt.ylabel(r'A',fontsize=20)
#     plt.yscale('log')
#     plt.savefig(pre_folder+"/Imgs/Area_str_{}.jpg".format(strg),bbox_inches='tight')
#     plt.clf()


#     E_vol=np.array(E_vol)
#     E_sur=np.array(E_sur)
#     E_ben=np.array(E_ben)
#     E_bead=np.array(E_bead)

#     plt.plot(time_evol,E_vol,color='cyan',label='volume')
#     plt.plot(time_evol,E_sur,color='green',label='surface',ls='dashed')
#     plt.plot(time_evol,E_ben,color='purple',label='bending')
#     plt.plot(time_evol,E_bead,color='red',label='bead')
#     plt.plot(time_evol,E_vol+E_sur+E_ben+E_bead,color='black',label='total')
#     # plt.axhline(V_bar,color='black',ls='dashed') 
#     plt.legend()
#     plt.title(r'Energy evolution',fontsize=22)
#     plt.xlabel(r't',fontsize=20)
#     plt.ylabel(r'E',fontsize=20)
#     # plt.yscale('log')
#     plt.savefig(pre_folder+"/Imgs/Energy_str_{}.jpg".format(strg),bbox_inches='tight')
#     plt.clf()



#     plt.scatter(time_evol,backtrackstep,color='purple')
#     # plt.axhline(A_bar,color='black',ls='dashed') 
#     plt.title(r'timestep evolution',fontsize=22)
#     plt.xlabel(r't',fontsize=20)
#     plt.ylabel(r'dt',fontsize=20)
#     plt.yscale('log')
#     plt.savefig(pre_folder+"/Imgs/Timestep_str_{}.jpg".format(strg),bbox_inches='tight')
#     plt.clf()


#     plt.scatter(time_evol,grad_norm,color='purple')
#     # plt.axhline(A_bar,color='black',ls='dashed') 
#     plt.title(r'gradient norm evolution',fontsize=22)
#     plt.xlabel(r't',fontsize=20)
#     plt.ylabel(r'D E',fontsize=20)
#     plt.yscale('log')
#     plt.savefig(pre_folder+"/Imgs/Gradient_norm_str_{}.jpg".format(strg),bbox_inches='tight')
#     plt.clf()







#     Output_file.close()




def main_shape(Ini_cond):


    pre_folder='../Results/Mem3DG_Cell_Shape/'
    dir='nu_{:.3f}_c0_0.000_KA_10.000_KB_{:.6f}_init_cond_{}_Nsim_{}'.format(nu,KB,Ini_cond,Nsim)
        #   nu_1.000_c0_0.000_KA_10.000_KB_0.005000_init_cond_1_Nsim_2
    # dir_2='nu_{:.3f}_c0_0.000_KA_10.000_KB_{:.6f}_init_cond_{}_Nsim_{}'.format(nu,KB,Ini,Nsim)
    
    folder_path= pre_folder+dir+'/'

    # Bead_path=folder_path+"Bead_data.txt"
    Output_path=folder_path+"Output_data.txt"
    # file_xyz = open(folder_path+'Bead_data.xyz','w')

    # Bead_file = open(Bead_path)
    Output_file= open(Output_path)
    # line = Output_file.readline()
    r = 1.0

    counter=0
    frame=0
    line=Output_file.readline()

    targ_vol=0.0
    targ_area=0.0


    # Sim_data << V_bar<<" "<< A_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " " << E_Ben <<" " << E_bead << " "<< grad_norm<<" " << backtrackstep<< " "<< Bead_1.Pos.x << " "<< Bead_1.Pos.y << " "<< Bead_1.Pos.z <<" \n";

    time_evol=[]
    Volumes=[]
    Bar_Areas=[]
    Areas=[]
    E_vol=[]
    E_sur=[]
    E_ben=[]
    E_bead=[]
    grad_norm=[]
    backtrackstep=[]
    line=Output_file.readline()
    counter=0
    while line and counter<500000:
        # frame=counter*1000
        splitted_line=line.split(' ')
        # print(len(splitted_line))
        if(len(splitted_line)<11):
            break
        if(len(splitted_line)>13):
            continue
        # Bar_Areas.append(float(splitted_line[1]))
        # time_evol.append(float(splitted_line[2]))
        # Volumes.append(float(splitted_line[3]))
        # Areas.append(float(splitted_line[4]))
        # E_vol.append(float(splitted_line[5]))
        # E_sur.append(float(splitted_line[6]))
        # E_ben.append(float(splitted_line[7]))
        # # E_bead.append(float(splitted_line[8]))
        # grad_norm.append(float(splitted_line[8]))
        # backtrackstep.append(float(splitted_line[9]))
        


        # V_bar=float(splitted_line[0])
        # A_bar=float(splitted_line[1])
        
        
        Bar_Area=float(splitted_line[1])
        time=float(splitted_line[2])
        Volume=float(splitted_line[3])
        Area=float(splitted_line[4])
        E_vol=float(splitted_line[5])
        E_sur=float(splitted_line[6])
        E_ben=float(splitted_line[7])
        # E_bead.append(float(splitted_line[8]))
        

        
        V_bar=float(splitted_line[0])
        A_bar=float(splitted_line[1])
        
        
        
        line=Output_file.readline()
        counter+=1


    


    # Ok so we have the data and its time to plot !

    os.makedirs(pre_folder+"Imgs/",exist_ok=True)
    os.makedirs(pre_folder+"Data/",exist_ok=True)


    E_data = open(pre_folder+"Data/Final_energies.txt","a+")

    E_data.write("{} {} {} {} {} {} {} {} \n".format(nu,Initial_conds,Nsim,E_vol,E_sur,E_ben,Area,Volume))
    # E_data.write("{} {} {} {} {} {} {} {} \n".format(nu,Initial_conds,Nsim,E_vol[-1],E_sur[-1],E_ben[-1],Areas[-1],Volumes[-1]))

    E_data.close()


    # plt.plot(time_evol,Volumes,color='purple')
    # plt.axhline(V_bar,color='black',ls='dashed') 
    # plt.title(r'Volume evolution',fontsize=22)
    # plt.xlabel(r't',fontsize=20)
    # plt.ylabel(r'V',fontsize=20)
    # plt.yscale('log')
    # plt.savefig(pre_folder+"/Imgs/Volume_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    # plt.clf()


    # plt.plot(time_evol,Bar_Areas,color='black',ls='dashed')
    # plt.plot(time_evol,Areas,color='purple')
    # # plt.axhline(A_bar,color='black',ls='dashed') 
    # plt.title(r'Area evolution',fontsize=22)
    # plt.xlabel(r't',fontsize=20)
    # plt.ylabel(r'A',fontsize=20)
    # plt.yscale('log')
    # plt.savefig(pre_folder+"/Imgs/Area_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    # plt.clf()


    # E_vol=np.array(E_vol)
    # E_sur=np.array(E_sur)
    # E_ben=np.array(E_ben)
    # # E_bead=np.array(E_bead)

    # plt.plot(time_evol,E_vol,color='cyan',label='volume')
    # plt.plot(time_evol,E_sur,color='green',label='surface',ls='dashed')
    # plt.plot(time_evol,E_ben,color='purple',label='bending')
    # # plt.plot(time_evol,E_bead,color='red',label='bead')
    # plt.plot(time_evol,E_vol+E_sur+E_ben,color='black',label='total')
    # # plt.axhline(V_bar,color='black',ls='dashed') 
    # plt.legend()
    # plt.ylim(0,0.25)
    # plt.title(r'Energy evolution',fontsize=22)
    # plt.xlabel(r't',fontsize=20)
    # plt.ylabel(r'E',fontsize=20)
    # # plt.yscale('log')
    # plt.savefig(pre_folder+"/Imgs/Energy_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    # plt.clf()



    # plt.scatter(time_evol,backtrackstep,color='purple')
    # # plt.axhline(A_bar,color='black',ls='dashed') 
    # plt.title(r'timestep evolution',fontsize=22)
    # plt.xlabel(r't',fontsize=20)
    # plt.ylabel(r'dt',fontsize=20)
    # plt.yscale('log')
    # plt.savefig(pre_folder+"/Imgs/Timestep_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    # plt.clf()


    # plt.scatter(time_evol,grad_norm,color='purple')
    # # plt.axhline(A_bar,color='black',ls='dashed') 
    # plt.title(r'gradient norm evolution',fontsize=22)
    # plt.xlabel(r't',fontsize=20)
    # plt.ylabel(r'D E',fontsize=20)
    # plt.yscale('log')
    # plt.savefig(pre_folder+"/Imgs/Gradient_norm_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    # plt.clf()



    # Volumes=np.array(Volumes)
    # Areas=np.array(Areas)


    # Output_file.close()

    return 
    # return [time_evol,3*Volumes/(4*np.pi)/ ( Areas/(4*np.pi))**(3/2)]

# for strg in [0.002, 0.6]:
# for strg in [0.00001,0.0001,0.05, 0.0005, 0.1 ,1.0 ,5.0,20.0]:
    # main(strg)

# data_1=main_shape(1)
# data_2=main_shape(2)
# data_3=main_shape(3)


main_shape(Initial_conds)


# OK soo

# pre_folder='../Results/Mem3DG_Cell_Shape/'
# plt.axhline(nu,color='black',ls='dashed')
# plt.plot(data_1[0],data_1[1],color=(255/256,51/256,109/256),label='Prolate')
# plt.plot(data_2[0],data_2[1],color=(123/256,40/256,125/256),label='Oblate')
# plt.plot(data_3[0],data_3[1],color=(255/256,166/256,48/256),label='Stomatocyte')

# plt.axhline(data_1[1][-1],ls='dashed',color=(255/256,51/256,109/256))
# plt.axhline(data_2[1][-1],ls='dashed',color=(123/256,40/256,125/256))
# plt.axhline(data_3[1][-1],ls='dashed',color=(255/256,166/256,48/256))
# plt.legend()
# # plt.xlim(0,np.min([data_1[0][-1],data_2[0][-1]]))
# plt.savefig(pre_folder+'Imgs/Reduced_vol_evol_nu_{}.jpg'.format(nu))




