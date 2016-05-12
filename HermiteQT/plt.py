

import numpy as np 
import matplotlib.pyplot as plt 

plt.subplot(211)

dat = np.genfromtxt('traj.dat')

for i in range(dat.shape[-1]-1):
    plt.plot(dat[:,0],dat[:,i+1])


plt.subplot(212)

dat = np.genfromtxt('xAve.dat')

plt.plot(dat[:,0],dat[:,1],'-')
plt.plot(dat[:,0],dat[:,2],'--')


plt.show() 
