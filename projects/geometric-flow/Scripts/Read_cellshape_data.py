import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino'], 'size':14})
rc('text', usetex=True)


Data_file = open("../Results/Mem3DG_Cell_Shape/Data/Final_energies.txt")

line=Data_file.readline()
nu_1=[]
nu_2=[]
nu_3=[]
nus=[[],[],[]]
Ebs=[[],[],[]]
Eb_1=[]
Eb_2=[]
Eb_3=[]

while(line):
    # print(line)

    splitted_line=line.split(' ')
    Init_cond = int(splitted_line[1])
    nus[Init_cond-1].append(float(splitted_line[0]))
    Ebs[Init_cond-1].append(float(splitted_line[5])+ float(splitted_line[4])+float(splitted_line[3]))

    line=Data_file.readline()


# print(len(nus[0]))
# print(len(nus[1]))
# print(len(nus[2]))

# print(len(Ebs[0]))
# print(len(Ebs[1]))
# print(len(Ebs[2]))


plt.scatter(nus[0],Ebs[0],color='pink')
plt.scatter(nus[1],Ebs[1],color='purple')
plt.scatter(nus[2],Ebs[2],color='blue')
plt.axvline(0.59,color='black',ls='dashed')
plt.axvline(0.65,color='black',ls='dashed')
plt.show()

