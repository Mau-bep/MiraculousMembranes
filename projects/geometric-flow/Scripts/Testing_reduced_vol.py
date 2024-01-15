import numpy as np 
import matplotlib.pyplot as plt 


nus=np.linspace(0.01,1.0,100)
V=1

A=4*np.pi*(3*V/(4*np.pi*nus))**(2.0/3.0)


plt.plot(nus,A)
plt.yscale('log')
plt.show()
