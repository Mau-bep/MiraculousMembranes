import numpy as np 
import matplotlib.pyplot as plt 


nus=np.linspace(0.01,1.0,100)
V=1

A=4*np.pi*(3*V/(4*np.pi*nus))**(2.0/3.0)


# plt.plot(nus,A)
# plt.yscale('log')
# plt.show()

K=10


x=np.linspace(0.0,10.0,100)
r1=5
r2=8
first_seg=np.ones_like(x)*(x<r1)*K
second_seg=(K/(r1-r2))*(x-r2)*(x>r1)*(x<r2)

force = first_seg+second_seg

plt.plot(x,-1*force)

# plt.show()
plt.clf()
# plt.plot(x,first_seg)
# plt.plot(x,second_seg)

U_first=-K*x*(x<r1)+-K*(r1-r2)*(x<r1)/2+K*r1*(x<r1)
U_second= -1*(K/(r1-r2))*(x>r1)*(x<r2)*(x-r2)*(x-r2)/2

U = U_first+U_second

plt.plot(x,-1*U)
# plt.show()

plt.clf()


data =np.loadtxt("../Results/Mem3DG_Bead_Reciprocal_arcsim/Some_radial_dist.txt")


plt.hist(data,bins=100)
plt.xlim(0,10)
plt.show()





