##!/usr/bin/python

import numpy as np
import pylab as plt
import seaborn as sns 

sns.set_context('poster')

plt.subplot(1,1,1)
data = np.genfromtxt(fname='energy.dat') 
#data = np.loadtxt('traj.dat')
#for x in range(1,data.shape[-1]):
plt.plot(data[:,0],data[:,1],label='kinetic')
plt.plot(data[:,0],data[:,2],label='potential')
plt.plot(data[:,0],data[:,3],label='quantum potential')
plt.plot(data[:,0],data[:,4],label='total')

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
plt.xlabel('time')
plt.ylabel('$x_i$')
#plt.title('traj')
plt.legend() 


plt.show() 

