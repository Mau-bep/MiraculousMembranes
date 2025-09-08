import numpy as np
import os 
import matplotlib.pyplot as plt 

import sys



# nu=float(sys.argv[1])
nu = 1.0
# KB=float(sys.argv[2])
# strg=float(sys.argv[3])

# nu=1.0
# Initial_conds=int(sys.argv[2])
Initial_conds= 1
# KB=float(sys.argv[3])
KB=0.01
# strength = float(sys.argv[3])
# Nsim=int(sys.argv[4])

# KA = float(sys.argv[5])



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
    pre_folder='../Results/Mem3DG_Cell_Shape_KB_evol/'
    # dir='nu_{:.3f}_c0_0.000_KA_{:.3f}_KB_{:.6f}_init_cond_{}_Nsim_{}'.format(nu,KA,KB,Ini_cond,Nsim)
    # dir = 'Particle_wrapping_on_plane_phase_space_sept/Surface_tension_15.0000_Bending_10.0000_Edge_reg_1.0000_Bead_radius_0.5000_Frenkel_Normal_nopush_str_140.0000_Switch_Newton_Nsim_1'
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
    while line and counter<400000:
        # frame=counter*1000
        splitted_line=line.split(' ')
        # print(len(splitted_line))
        if(len(splitted_line)<11):
            break
        if(len(splitted_line)>13):
            line=Output_file.readline()
            line=Output_file.readline()
            continue
        Bar_Areas.append(float(splitted_line[1]))
        time_evol.append(float(splitted_line[2]))
        Volumes.append(float(splitted_line[0]))
        Areas.append(float(splitted_line[1]))
        E_vol.append(float(splitted_line[5]))
        E_sur.append(float(splitted_line[3]))
        E_ben.append(float(splitted_line[4]))
        # E_bead.append(float(splitted_line[8]))
        grad_norm.append(float(splitted_line[8]))
        backtrackstep.append(float(splitted_line[9]))
        


        # V_bar=float(splitted_line[0])
        # A_bar=float(splitted_line[1])
        

        # i can do a safety check 
        iftruecontinue=False
        if(splitted_line[1]=='3.558910.568273'):
            print(len(splitted_line[1].split('.')))
        for k in range(8):
            if(len(splitted_line[k].split('.'))>2):
                iftruecontinue=True
                # line=Output_file.readline()
                # continue
        if(iftruecontinue):
            line=Output_file.readline()
            continue
        Bar_Area=float(splitted_line[1])
        time=float(splitted_line[2])
        Volume=float(splitted_line[3])
        Area=float(splitted_line[4])
        E_vol1=float(splitted_line[5])
        E_sur1=float(splitted_line[6])
        E_ben1=float(splitted_line[7])
        # E_bead.append(float(splitted_line[8]))
        V_bar=float(splitted_line[0])
        A_bar=float(splitted_line[1])
        
        
        
        line=Output_file.readline()
        counter+=1


    


    # Ok so we have the data and its time to plot !

    os.makedirs(pre_folder+"Imgs/",exist_ok=True)
    os.makedirs(pre_folder+"Data/",exist_ok=True)


    E_data = open(pre_folder+"Data/Final_energies_KB_{}_Nsim_{}.txt".format(KB,Nsim),"a+")

    E_data.write("{} {} {} {} {} {} {} {} \n".format(nu,Initial_conds,Nsim,E_vol1,E_sur1,E_ben1,Area,Volume))
    # E_data.write("{} {} {} {} {} {} {} {} \n".format(nu,Initial_conds,Nsim,E_vol[-1],E_sur[-1],E_ben[-1],Areas[-1],Volumes[-1]))

    E_data.close()


    plt.plot(time_evol,Volumes,color='purple')
    plt.axhline(V_bar,color='black',ls='dashed') 
    plt.title(r'Volume evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'V',fontsize=20)
    plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Volume_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    plt.clf()


    plt.plot(time_evol,Bar_Areas,color='black',ls='dashed')
    plt.plot(time_evol,Areas,color='purple')
    # plt.axhline(A_bar,color='black',ls='dashed') 
    plt.title(r'Area evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'A',fontsize=20)
    plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Area_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    plt.clf()


    E_vol=np.array(E_vol)
    E_sur=np.array(E_sur)
    E_ben=np.array(E_ben)
    # E_bead=np.array(E_bead)

    plt.plot(time_evol,E_vol,color='cyan',label='volume')
    plt.plot(time_evol,E_sur,color='green',label='surface',ls='dashed')
    plt.plot(time_evol,E_ben,color='purple',label='bending')
    # plt.plot(time_evol,E_bead,color='red',label='bead')
    plt.plot(time_evol,E_vol+E_sur+E_ben,color='black',label='total')
    # plt.axhline(V_bar,color='black',ls='dashed') 
    plt.legend()
    # plt.ylim(0,0.25)
    plt.title(r'Energy evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'E',fontsize=20)
    # plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Energy_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    plt.clf()



    plt.scatter(time_evol,backtrackstep,color='purple')
    # plt.axhline(A_bar,color='black',ls='dashed') 
    plt.title(r'timestep evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'dt',fontsize=20)
    plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Timestep_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    plt.clf()


    plt.scatter(time_evol,grad_norm,color='purple')
    # plt.axhline(A_bar,color='black',ls='dashed') 
    plt.title(r'gradient norm evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'D E',fontsize=20)
    plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Gradient_norm_nu_{}_init_cond_{}.jpg".format(nu,Ini_cond),bbox_inches='tight')
    plt.clf()



    Volumes=np.array(Volumes)
    Areas=np.array(Areas)


    Output_file.close()

    return 







def main_bead(Ini_cond,strength):

    nu = 1.0
    KB=0.01
    
    Nsim=100

    KA = 500

    pre_folder='../Results/Mem3DG_Bead_Reciprocal/'
    # pre_folder='../Results/Mem3DG_Cell_Shape_KB_evol/'
    dir='nu_{:.3f}_c0_0.000_KA_{:.3F}_KB_{:.6f}_strength_{:.6f}_Init_cond_{}_Nsim_{}'.format(nu,KA,KB,strength,Ini_cond,Nsim)
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
    while line and counter<400000:
        # frame=counter*1000
        splitted_line=line.split(' ')
        # print(len(splitted_line))
        if(len(splitted_line)<11):
            break
        if(len(splitted_line)>13):
            line=Output_file.readline()
            line=Output_file.readline()
            continue
        Bar_Areas.append(float(splitted_line[1]))
        time_evol.append(float(splitted_line[2]))
        Volumes.append(float(splitted_line[3]))
        Areas.append(float(splitted_line[4]))
        E_vol.append(float(splitted_line[5]))
        E_sur.append(float(splitted_line[6]))
        E_ben.append(float(splitted_line[7]))
        E_bead.append(float(splitted_line[8]))
        grad_norm.append(float(splitted_line[9]))
        backtrackstep.append(float(splitted_line[10]))
        


        # V_bar=float(splitted_line[0])
        # A_bar=float(splitted_line[1])
        

        # i can do a safety check 
        iftruecontinue=False
        if(splitted_line[1]=='3.558910.568273'):
            print(len(splitted_line[1].split('.')))
        for k in range(8):
            if(len(splitted_line[k].split('.'))>2):
                iftruecontinue=True
                # line=Output_file.readline()
                # continue
        if(iftruecontinue):
            line=Output_file.readline()
            continue
        Bar_Area=float(splitted_line[1])
        time=float(splitted_line[2])
        Volume=float(splitted_line[3])
        Area=float(splitted_line[4])
        E_vol1=float(splitted_line[5])
        E_sur1=float(splitted_line[6])
        E_ben1=float(splitted_line[7])
        E_bead1=float(splitted_line[8])
        V_bar=float(splitted_line[0])
        A_bar=float(splitted_line[1])
        
        
        
        line=Output_file.readline()
        counter+=1


    


    # Ok so we have the data and its time to plot !

    os.makedirs(pre_folder+"Imgs/",exist_ok=True)
    os.makedirs(pre_folder+"Data/",exist_ok=True)


    E_data = open(pre_folder+"Data/Final_energies_strength_{}_Nsim_{}.txt".format(strength,Nsim),"a+")

    E_data.write("{} {} {} {} {} {} {} {} {} \n".format(nu,Initial_conds,Nsim,E_vol1,E_sur1,E_ben1,E_bead,Area,Volume))
    # E_data.write("{} {} {} {} {} {} {} {} {} \n".format(nu,Initial_conds,Nsim,E_vol[-1],E_sur[-1],E_ben[-1],Areas[-1],Volumes[-1]))

    E_data.close()


    plt.plot(time_evol,Volumes,color='purple')
    plt.axhline(V_bar,color='black',ls='dashed') 
    plt.title(r'Volume evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'V',fontsize=20)
    plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Volume_nu_{}_init_cond_{}_strength_{}_Nsim_{}.jpg".format(nu,Ini_cond,strength,Nsim),bbox_inches='tight')
    plt.clf()


    plt.plot(time_evol,Bar_Areas,color='black',ls='dashed')
    plt.plot(time_evol,Areas,color='purple')
    # plt.axhline(A_bar,color='black',ls='dashed') 
    plt.title(r'Area evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'A',fontsize=20)
    plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Area_nu_{}_init_cond_{}_strength_{}_Nsim_{}.jpg".format(nu,Ini_cond,strength,Nsim),bbox_inches='tight')
    plt.clf()


    E_vol=np.array(E_vol)
    E_sur=np.array(E_sur)
    E_ben=np.array(E_ben)
    E_bead=np.array(E_bead)

    plt.plot(time_evol,E_vol,color='cyan',label='volume')
    plt.plot(time_evol,E_sur,color='green',label='surface',ls='dashed')
    plt.plot(time_evol,E_ben,color='purple',label='bending')
    plt.plot(time_evol,E_bead,color='red',label='bead')
    plt.plot(time_evol,E_vol+E_sur+E_ben+E_bead,color='black',label='total')
    # plt.axhline(V_bar,color='black',ls='dashed') 
    plt.legend()
    # plt.ylim(0,0.25)
    plt.title(r'Energy evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'E',fontsize=20)
    # plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Energy_nu_{}_init_cond_{}_strength_{}_Nsim_{}.jpg".format(nu,Ini_cond,strength,Nsim),bbox_inches='tight')
    plt.clf()



    plt.scatter(time_evol,backtrackstep,color='purple')
    # plt.axhline(A_bar,color='black',ls='dashed') 
    plt.title(r'timestep evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'dt',fontsize=20)
    plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Timestep_nu_{}_init_cond_{}_strength_{}_Nsim_{}.jpg".format(nu,Ini_cond,strength,Nsim),bbox_inches='tight')
    plt.clf()


    plt.scatter(time_evol,grad_norm,color='purple')
    # plt.axhline(A_bar,color='black',ls='dashed') 
    plt.title(r'gradient norm evolution',fontsize=22)
    plt.xlabel(r't',fontsize=20)
    plt.ylabel(r'D E',fontsize=20)
    plt.yscale('log')
    plt.savefig(pre_folder+"/Imgs/Gradient_norm_nu_{}_init_cond_{}_strength_{}_Nsim_{}.jpg".format(nu,Ini_cond,strength,Nsim),bbox_inches='tight')
    plt.clf()



    Volumes=np.array(Volumes)
    Areas=np.array(Areas)


    Output_file.close()

    return 
    # return [time_evol,3*Volumes/(4*np.pi)/ ( Areas/(4*np.pi))**(3/2)]

# for strg in [0.002, 0.6]:
# for strg in [0.00001,0.0001,0.05, 0.0005, 0.1 ,1.0 ,5.0,20.0]:
    # main(strg)

# data_1=main_shape(1)
# data_2=main_shape(2)
# data_3=main_shape(3)


# main_shape(Initial_conds)

# for strength_val in [0.001,0.1,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]:
#     main_bead(Initial_conds,strength_val)
# # OK soo





def bead_attachment_Energy(Ini_cond):

    pre_folder='../Results/Mem3DG_Bead_Reciprocal/'
    # pre_folder='../Results/Mem3DG_Cell_Shape_KB_evol/'
    dir='nu_{:.3f}_c0_0.000_KA_{:.3F}_KB_{:.6f}_strength_{:.6f}_Init_cond_{}_Nsim_{}'.format(nu,KA,KB,strength,Ini_cond,Nsim)
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
    while line:
        # frame=counter*1000
        splitted_line=line.split(' ')
        # print(len(splitted_line))
        if(len(splitted_line)<11):
            break
        if(len(splitted_line)>13):
            line=Output_file.readline()
            line=Output_file.readline()
            continue
        Bar_Areas.append(float(splitted_line[1]))
        time_evol.append(float(splitted_line[2]))
        Volumes.append(float(splitted_line[3]))
        Areas.append(float(splitted_line[4]))
        E_vol.append(float(splitted_line[5]))
        E_sur.append(float(splitted_line[6]))
        E_ben.append(float(splitted_line[7]))
        E_bead.append(float(splitted_line[8]))
        grad_norm.append(float(splitted_line[9]))
        backtrackstep.append(float(splitted_line[10]))
        


        # V_bar=float(splitted_line[0])
        # A_bar=float(splitted_line[1])
        

        # i can do a safety check 
        iftruecontinue=False
        for k in range(8):
            if(len(splitted_line[k].split('.'))>2):
                iftruecontinue=True
                # line=Output_file.readline()
                # continue
        if(iftruecontinue):
            line=Output_file.readline()
            continue
        
        Bar_Area=float(splitted_line[1])
        time=float(splitted_line[2])
        Volume=float(splitted_line[3])
        Area=float(splitted_line[4])
        E_vol1=float(splitted_line[5])
        E_sur1=float(splitted_line[6])
        E_ben1=float(splitted_line[7])
        E_bead1=float(splitted_line[8])
        V_bar=float(splitted_line[0])
        A_bar=float(splitted_line[1])
        
        
        
        line=Output_file.readline()
        counter+=1


    


    # Ok so we have the data and its time to plot !

    os.makedirs(pre_folder+"Imgs/",exist_ok=True)
    os.makedirs(pre_folder+"Data/",exist_ok=True)


    E_data = open(pre_folder+"Data/Final_bead_energies_Nsim_{}.txt".format(Nsim),"a+")

    E_data.write("{} {} {} {} {} {} {} {} \n".format(nu,Initial_conds,Nsim,KA,strength,E_bead1,Area,Volume))
    # E_data.write("{} {} {} {} {} {} {} {} {} \n".format(nu,Initial_conds,Nsim,E_vol[-1],E_sur[-1],E_ben[-1],Areas[-1],Volumes[-1]))

    E_data.close()

    Output_file.close()

    return E_bead[-1]




def test_tension():
    KAs = [10, 100, 200, 500, 1000]
    Energies=[]
    strengths = [0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0] 
    pre_folder='../Results/Mem3DG_Bead_Reciprocal/'
        # pre_folder='../Results/Mem3DG_Cell_Shape_KB_evol/'
    Nsim=100
    for i in range(len(KAs)):
        KA=KAs[i]
        En=[]
        for j in range(len(strengths)):
            strength=strengths[j]
            En.append(bead_attachment_Energy(1))
        Energies.append(En)
        plt.plot(strengths,En,label='KA = {}'.format(KA))
    plt.legend()
    plt.xlabel("Interaction strength")
    plt.ylabel("Binding Energy")
    plt.savefig(pre_folder+"/Imgs/Bead_Energy_differen_KA.png",bbox_inches='tight')

    print(Energies)





def main_bead_energy(Nsim,KB,strength):

    nu = 1.0


    KA = 100000
    rc = 2.0
    Ini_cond = 1
    pre_folder='../Results/Mem3DG_Bead_Pulling_rc_august_arcsim/'
    pre_folder = '../Results/Particle_wrapping_on_plane_phase_space_sept/'
    pre_folder = '../Results/Wrapping_single_TEST_meet/'
    # pre_folder='../Results/Mem3DG_Cell_Shape_KB_evol/'
    # dir='nu_{:.3f}_rc_{:.3f}_KA_{:.3F}_KB_{:.6f}_strength_{:.6f}_Init_cond_{}_Nsim_{}'.format(nu,rc,KA,KB,strength,Ini_cond,Nsim)
    
    dir = 'Surface_tension_15.0000_Bending_10.0000_Edge_reg_1.0000_Bead_radius_0.5000_Frenkel_Normal_nopush_str_140.0000_Switch_Newton_Nsim_1'
    # dir = 'Surface_tension_10.0000_Bending_20.0000_Edge_reg_1.0000_Bead_radius_0.3000_Frenkel_Normal_nopush_str_325.0000_Switch_Newton_Nsim_1'
        #   nu_1.000_c0_0.000_KA_10.000_KB_0.005000_init_cond_1_Nsim_2
    # dir_2='nu_{:.3f}_c0_0.000_KA_10.000_KB_{:.6f}_init_cond_{}_Nsim_{}'.format(nu,KB,Ini,Nsim)
    # folder_path = '../Results/Particle_wrapping_on_plane_Sept/Surface_tension_10.0000_Bending_20.0000_Edge_reg_10.0000_Bead_radius_0.3000_Frenkel_Normal_nopush_str_500.0000_Switch_Newton_Switch_t_10000_Nsim_14/'
    
    dir = 'Bending_10.0000_Surface_tension_15.0000_Edge_reg_1.0000_Bead_radius_0.3000_Frenkel_Normal_nopush_str_300.0000_Switch_Newton_Nsim_1'
    dir2 = 'Bending_10.0000_Surface_tension_15.0000_Edge_reg_1.0000_Bead_radius_0.3000_Frenkel_Normal_nopush_str_300.0000_Nsim_1'

    folder_path= pre_folder+dir+'/'
    # folder_path = '../Results/Particle_wrapping_on_plane_Sept/Surface_tension_10.0000_Bending_20.0000_Edge_reg_10.0000_Bead_radius_0.3000_Frenkel_Normal_nopush_str_500.0000_Switch_Newton_Switch_t_10000_Nsim_14/'
    
    # Bead_path=folder_path+"Bead_data.txt"
    Output_path=folder_path+"Output_data.txt"
    Output_path2 = pre_folder + dir2 + '/Output_data.txt'

    Output_file= open(Output_path)
    Output_file2 = open(Output_path2)

    line = Output_file.readline()
    line2 = Output_file2.readline()

    line = Output_file.readline()
    line2 = Output_file2.readline()

    E_bead = []
    times = []
    E_tot = []
    
    E_bead2 = []
    times2 = []
    E_tot2 = []

    while(line and line2):



        data = line.split(' ')
        times.append(float(data[2]))
        E_bead.append(float(data[6]))
        E_tot.append(float(data[7]))
        line = Output_file.readline()

        data = line2.split(' ')
        times2.append(float(data[2]))
        E_bead2.append(float(data[6]))
        E_tot2.append(float(data[7]))
        line2 = Output_file2.readline()



    # plt.yscale('log')
    plt.plot(E_tot,color='black',label = 'With Newton')
    plt.plot(E_tot2,color='purple',label = 'No Newton')
    plt.axvline(40,color='black',ls = 'dashed')
    # plt.xlabel("Time")
    plt.legend()
    plt.ylabel("Total energy")
    plt.show()


    Bead_data_path = folder_path+ "Bead_data_SS.txt"

    Bead_data_file = open(Bead_data_path)

    line = Bead_data_file.readline()

    line = Bead_data_file.readline()
    Mag_f  = []
    X_bead = []
    while(line):

        data = line.split(' ')
        Force = np.array( [float(data[3]), float(data[4]), float(data[5])] )
        Mag_f.append( np.sqrt(np.sum(Force*Force)))
        X_bead.append(float(data[0]))

        line = Bead_data_file.readline()

    # print(Mag_f)
    plt.xlabel("Bead position")
    plt.ylabel("Force magnitude")
    plt.plot(X_bead,Mag_f,color='purple')
    plt.show()
















    plt.xlabel('timestep')
    plt.ylabel('Bead energy')
    plt.scatter(times,E_bead,color='black')
    plt.show()
    differences = []
    # x_useful = []
    # for i in range(len(E_bead)):
    #     differences.append(abs(E_bead[i-1]-E_bead[i]))
    #     if(differences[i]<1e-3):
    #         # print(i)
    #         x_useful.append(i)
    # plt.plot(differences,color='purple')
    # plt.yscale('log')
    # plt.show()

    # plt.axvline(x_useful[0],color='purple')
    # plt.plot(E_bead)

    # plt.show()








def Get_energies_data(Interaction_strengths):
    # This function creates a file that gets the following
    # Volume_E Area_E Bending_E Interaction_E Interaction_strength
    parent_folder= "../Results/Mem3DG_Bead_Reciprocal_finemesh_varKB/"

    E_bead_list = []
    E_area_list = []
    E_vol_list = []
    E_bend_list = []

    counter=0

    Storedata_dir = parent_folder + "Wrapping_data_energies.txt"
    Storedata = open(Storedata_dir)
        
    for strength in Interaction_strengths:
        E_bend_list.append([])
        E_bead_list.append([])
        E_area_list.append([])
        E_vol_list.append([])

        
        Sim_folder = "nu_1.000_radius_1.000_curvadap_0.00_minrel_0.1000_KA_100000.000_KB_1.000000_strength_{:.6f}_Init_cond_2_Nsim_1/".format(strength)
        Output_filename = parent_folder + Sim_folder + "Output_data.txt"

        Output_file = open(Output_filename)

        line = Output_file.readline()

        line = Output_filename.readline()

        while(line):


            data = line.split(' ')
            if(not(len(data)>10)):
               break


            E_vol = data[5]
            E_area = data[6]
            E_bend = data[7]
            E_bead = data[8]


            line = Output_filename.readline()

        Storedata.write("{} {} {} {} {}".format(strength,E_vol,E_area,E_bend,E_bead))
        E_vol_list[counter].append(E_vol)
        E_area_list[counter].append(E_area)
        E_bead_list[counter].append(E_bead)
        E_bend_list[counter].append(E_bend)
    
        counter+=1

    Storedata.close()








    return 




# strengths = [10,50,100,150]

# Get_energies_data(strengths)





main_bead_energy(20,1.0,2.0)

