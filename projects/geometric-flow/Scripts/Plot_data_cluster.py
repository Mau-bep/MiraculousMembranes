import numpy as np
import matplotlib.pyplot as plt
import sys
# plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
# plt.rc('font',**{'family':'serif','serif':['Palatino'], 'size':14})
# plt.rc('text', usetex=True)
KA=0.1
KB=0.00001

# dir='nu_1.000_c0_0.100_KA_10.000_KB_0.000010'

nu=float(sys.argv[1])
KB=float(sys.argv[2])
# dir='nu_1.000_c0_0.000_KA_10.000_KB_0.001000'
dir='nu_{:.3f}_c0_0.000_KA_10.000_KB_{:.6f}'.format(nu,KB)


# Data=np.loadtxt('./Curv_adap_0.10Min_rel_length_0.50/nu_1.00_c0_0.00_KA_1.00_KB_2.00/Output_data.txt',delimiter=" ")
# f=open('./Mem3DG_IMG/Curv_adap_0.10Min_rel_length_0.50/'+dir+'/Output_data.txt')
folder_path='./Mem3DG_IMG_correct/Curv_adap_0.10Min_rel_length_0.50/'
#folder_path='./Tests_cil_regular_finite/'
file_path=folder_path+dir+'/Output_data.txt'

t=[]
dts=[]
V_evol=[]
A_evol=[]
E_sur=[]
E_vol=[]
E_bend=[]
counter=0
A_bar=[]
norm_grad=[]
# counter=0
with open(file_path, "r") as file:
    # Read the first line


    line = file.readline()
    V_bar=float(line.split(' ')[0])
    # A_bar=float(line.split(' ')[1])
    # Loop until the end of the file
    while line:
        counter+=1
        
            # strip() removes the newline character at the end
        Data = file.readline()
        

        # print(Data.split(' ')[2:])
        Arr=np.array(Data.split(' ')[2:-1],dtype=float)
        #print(Arr)
        # print(Arr)
        if(len(Arr)<8):
            break
        # if(counter%100==0):
        # print(counter)
        A_bar.append(float(Data.split(' ')[1]))
        t.append(Arr[0])
        
        V_evol.append(Arr[1])
        A_evol.append(Arr[2])
        E_vol.append(Arr[3])
        E_sur.append(Arr[4])
        E_bend.append(Arr[5])
        norm_grad.append(Arr[6])
        dts.append(Arr[7])

t=np.array(t)
A_bar=np.array(A_bar)
V_evol=np.array(V_evol)
A_evol=np.array(A_evol)
E_vol=np.array(E_vol)
E_sur=np.array(E_sur)
E_bend=np.array(E_bend)


plt.xlabel('t')
plt.title('Volume Evolution')
plt.scatter(t[::100],V_evol[::100],color='red',label='V')
plt.axhline(V_bar,color='black',ls='dashed',label='Target V')
plt.savefig(folder_path+'Volume_evolution_'+dir+'.jpg',bbox_inches='tight')
# plt.show()
plt.clf()
plt.title('Area Evolution')
plt.xlabel('t')
plt.scatter(t[::100],A_evol[::100],color='purple',label='A')
# plt.axhline(A_bar)
plt.plot(t[::100],A_bar[::100],color='black',ls='dashed',label='Target A')
plt.legend()
plt.savefig(folder_path+'Area_evolution_'+dir+'.jpg',bbox_inches='tight')
# plt.show()
plt.clf()


plt.xlabel('t')
plt.title('Energies Evolution')
plt.plot(t[::100],E_vol[::100]+E_sur[::100]+E_bend[::100],color='black',ls='dashed',label='total E')
# plt.xlim(50)
plt.scatter(t[::100],E_vol[::100],color='red',label='E_Vol')
plt.scatter(t[::100],E_sur[::100],color='purple',label='E_Sur')
plt.scatter(t[::100],E_bend[::100],color='pink',label='E_Bend')
# plt.axvline(100,color='pink')
plt.legend()
plt.savefig(folder_path+'Energies_evolution_'+dir+'.jpg',bbox_inches='tight')
# # plt.show()
# print(dts)
plt.clf()
dts=np.array(dts)
# plt.plot(t)
plt.title('Time evolution')
# plt.ylim(1e-2,1e-15)
# plt.axhline(1e-7)
# # print(dts)
# # plt.hist(dts,bins=100)
# # print(np.mean(dts[-100,0]))

# plt.yscale('log')
plt.plot(dts[::100])
plt.yscale('log')
# plt.plot(dts[::100],color='black')
plt.savefig(folder_path+'timestep_evolution_'+dir+'.jpg',bbox_inches='tight')
plt.xlabel('t')
plt.show()
plt.clf()

norm_grad=np.array(norm_grad)

plt.scatter(t[::100],norm_grad[::100],color='green',label='Gradient norm')
plt.legend()
plt.yscale('log')
plt.savefig(folder_path+'gradient_evolution_'+dir+'.jpg',bbox_inches='tight')
plt.show()
# print(np.unique(np.array(dts)[-10000:]))


# I can do the reduce volume evolution hehe

# A_evol=np.array(A_evol)
# V_evol=np.array(V_evol)

plt.clf()
nu_evol= 3*V_evol/(4*np.pi)/ ( A_evol/(4*np.pi))**(3/2)
plt.title("Reduced volume evolution")
plt.axhline(nu,color='black',ls='dashed')
plt.plot(t[::100],nu_evol[::100],color='cyan')
plt.savefig(folder_path+'Reduced_volume_evol_'+dir+'.jpg',bbox_inches='tight')
# plt.show()
