import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino'], 'size':14})
rc('text', usetex=True)


KB=0.01
Data_file = open("../Results/Mem3DG_Cell_Shape_KB_evol/Data/Final_energies_KB_{}_Nsim_11_thisdataisgood.txt".format(KB))

line=Data_file.readline()
nu_1=[]
nu_2=[]
nu_3=[]
nus_target=[[],[],[]]
nus=[[],[],[]]
Ebs=[[],[],[]]
Eb_1=[]
Eb_2=[]
Eb_3=[]
Eas=[[],[],[]]
while(line):
    # print(line)

    splitted_line=line.split(' ')
    Init_cond = int(splitted_line[1])
    Volume=float(splitted_line[7])
    Area=float(splitted_line[6])
    nus[Init_cond-1].append( 3*Volume/(4*np.pi*pow( Area/(4*np.pi) ,1.5 ))  )
    nus_target[Init_cond-1].append( float(splitted_line[0]))

    Eas[Init_cond-1].append(float(splitted_line[4]))
    Ebs[Init_cond-1].append(float(splitted_line[5]))

    line=Data_file.readline()


# print(len(nus[0]))
# print(len(nus[1]))
# print(len(nus[2]))

# print(len(Ebs[0]))
# print(len(Ebs[1]))
# print(len(Ebs[2]))



Data_file_2 = open("../Results/Mem3DG_Cell_Shape_KB_evol/Data/Final_energies_KB_{}_Nsim_11_second.txt".format(KB))

line=Data_file_2.readline()



while(line):
    # print(line)

    splitted_line=line.split(' ')
    Init_cond = int(splitted_line[1])
    Volume=float(splitted_line[7])
    Area=float(splitted_line[6])
    nus[Init_cond-1].append( 3*Volume/(4*np.pi*pow( Area/(4*np.pi) ,1.5 ))  )
    nus_target[Init_cond-1].append( float(splitted_line[0]))
    Eas[Init_cond-1].append(float(splitted_line[4]))
    Ebs[Init_cond-1].append(float(splitted_line[5]))

    line=Data_file_2.readline()








plt.scatter(nus[0],Ebs[0],color='pink')
plt.scatter(nus[1],Ebs[1],color='purple')
plt.scatter(nus[2],Ebs[2],color='blue')
plt.ylabel(r'$E_b$')
plt.xlabel(r'$\nu$')
plt.axvline(0.59,color='black',ls='dashed')
plt.axvline(0.65,color='black',ls='dashed')
plt.show()


plt.scatter(nus_target[0],nus[0],color='pink')
plt.scatter(nus_target[1],nus[1],color='purple')
plt.scatter(nus_target[2],nus[2],color='blue')
diag=np.linspace(0,1,10)
plt.plot(diag,diag,color='black',ls='dashed')
plt.ylabel(r'$\nu$')
plt.xlabel(r'$\bar{\nu}$')
plt.show()

plt.scatter(nus[0],Eas[0],color='pink')
plt.scatter(nus[1],Eas[1],color='purple')
plt.scatter(nus[2],Eas[2],color='blue')
plt.yscale('log')
plt.ylabel(r'$E_a$')
plt.xlabel(r'$\nu$')
plt.axvline(0.59,color='black',ls='dashed')
plt.axvline(0.65,color='black',ls='dashed')
plt.show()