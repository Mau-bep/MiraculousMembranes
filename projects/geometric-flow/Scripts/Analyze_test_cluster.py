import numpy as np
import matplotlib.pyplot as plt 
import os
import sys
import time


nu=float(sys.argv[1])
KB=float(sys.argv[2])
prefolder='./../Results/Tests_cil_regular/'


def Area_analisis(prefolder,nu,KB):
    dir=prefolder+'nu_{0:.3f}_c0_0.000_KA_10.000_KB_{1:.6f}/'.format(nu,KB)

    grad_theory_comp_t=[]
    grad_finite_comp_t=[]
    tot_grad_theory_t=[]
    tot_grad_finite_t=[]

    for step in range(0,100000,1000):
        prov_grad_theory=[]
        prov_grad_finite=[]
        filename=dir+'Area_Gradient_evaluation_{}.txt'.format(step)
        with open(filename,'r') as file:
            line=file.readline()
            while(line):
                
                Arr=np.array(line.split(' '))
                if(line[0]=='A'):
                    line=file.readline()
                    continue
                if(len(Arr)>2):
                    prov_grad_finite.append(float(Arr[3]))
                    line=file.readline()
                    Arr=np.array(line.split(' '))
                    prov_grad_theory.append(float(Arr[3]))
                else:
                    tot_grad_theory_t.append(float(Arr[0])/2)
                    tot_grad_finite_t.append(float(Arr[1][:-1]))
                    

                line=file.readline()    
            
        grad_finite_comp_t.append(prov_grad_finite)
        grad_theory_comp_t.append(prov_grad_theory)

    # I can do the plots here
    plt.plot(tot_grad_finite_t,color='purple',label='Finite diff')
    plt.plot(tot_grad_theory_t,color='pink',ls='dashed',label='Theory')
    plt.xlabel('timestep')
    plt.yscale('log')
    plt.title('Area E gradient')
    plt.ylabel('Gradient norm')
    plt.legend()
    plt.savefig(prefolder+'Area_grad_norm_evolution_full_nu_{}_KB_{}.jpg'.format(nu,KB),bbox_inches='tight')
    # plt.show()
    plt.clf()
    return 

def Volume_analisis(prefolder,nu,KB):
    dir=prefolder+'nu_{0:.3f}_c0_0.000_KA_10.000_KB_{1:.6f}/'.format(nu,KB)

    grad_theory_comp_t=[]
    grad_finite_comp_t=[]
    tot_grad_theory_t=[]
    tot_grad_finite_t=[]

    for step in range(0,100000,1000):
        prov_grad_theory=[]
        prov_grad_finite=[]
        filename=dir+'Vol_Gradient_evaluation_{}.txt'.format(step)
        with open(filename,'r') as file:
            line=file.readline()
            while(line):
                
                Arr=np.array(line.split(' '))
                if(line[0]=='V'):
                    line=file.readline()
                    continue
                if(len(Arr)>2):
                    line=file.readline()
                    continue
                    prov_grad_finite.append(float(Arr[3]))
                    line=file.readline()
                    Arr=np.array(line.split(' '))
                    prov_grad_theory.append(float(Arr[3]))
                else:
                    tot_grad_theory_t.append(float(Arr[0]))
                    tot_grad_finite_t.append(float(Arr[1][:-1]))
                    

                line=file.readline()    
            
        grad_finite_comp_t.append(prov_grad_finite)
        grad_theory_comp_t.append(prov_grad_theory)

    # I can do the plots here
    plt.plot(tot_grad_finite_t,color='purple',label='Finite diff')
    plt.plot(tot_grad_theory_t,color='pink',ls='dashed',label='Theory')
    plt.xlabel('timestep')
    plt.yscale('log')
    plt.ylabel('Gradient norm')
    plt.title('Volume E gradient')
    plt.legend()
    plt.savefig(prefolder+'Volume_grad_norm_evolution_full_nu_{}_KB_{}.jpg'.format(nu,KB),bbox_inches='tight')
    # plt.show()
    plt.clf()
    return 



def Bending_analisis(prefolder,nu,KB):
    dir=prefolder+'nu_{0:.3f}_c0_0.000_KA_10.000_KB_{1:.6f}/'.format(nu,KB)

    grad_theory_comp_t=[]
    grad_finite_comp_t=[]
    tot_grad_theory_t=[]
    tot_grad_finite_t=[]

    for step in range(0,100000,1000):
        prov_grad_theory=[]
        prov_grad_finite=[]
        filename=dir+'Bending_Gradient_evaluation_{}.txt'.format(step)
        with open(filename,'r') as file:
            line=file.readline()
            while(line):
                
                Arr=np.array(line.split(' '))
                # print(Arr)
                if(line[0]=='B'):
                    line=file.readline()
                    continue
                if(len(Arr)>3):
                    line=file.readline()
                    continue
                    prov_grad_finite.append(float(Arr[3]))
                    line=file.readline()
                    Arr=np.array(line.split(' '))
                    prov_grad_theory.append(float(Arr[3]))
                else:
                    # print('Ever here? \n')
                    tot_grad_theory_t.append(float(Arr[0]))
                    tot_grad_finite_t.append(float(Arr[1]))
                    

                line=file.readline()    
            
        grad_finite_comp_t.append(prov_grad_finite)
        grad_theory_comp_t.append(prov_grad_theory)

    # I can do the plots here
    plt.plot(tot_grad_finite_t,color='purple',label='Finite diff')
    plt.plot(tot_grad_theory_t,color='pink',ls='dashed',label='Theory')
    plt.xlabel('timestep')
    plt.title('Bending E gradient')
    plt.yscale('log')
    plt.ylabel('Gradient norm')
    plt.legend()
    plt.savefig(prefolder+'Bending_grad_norm_evolution_full_nu_{}_KB_{}.jpg'.format(nu,KB),bbox_inches='tight')
    # plt.show()
    plt.clf()
    return 


Area_analisis(prefolder,nu,KB)
print('Finished area grad\n')
Volume_analisis(prefolder,nu,KB)
print('finished volume grad\n')
Bending_analisis(prefolder,nu,KB)
print('Finished bending grad\n')

time.sleep(1)
