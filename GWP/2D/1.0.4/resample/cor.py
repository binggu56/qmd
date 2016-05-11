##!/usr/bin/python

import numpy as np
import pylab as plt
import seaborn as sns 

sns.set_context('poster')

plt.subplot(1,1,1)
data = np.genfromtxt(fname='cor.dat') 

for x in range(1,data.shape[-1]):
    plt.plot(data[:,0],data[:,x],'k--',lw=2)

plt.xlabel('$time$')
plt.ylabel('$C(t)$')
#plt.title('traj')

#plt.subplot(2,1,2)
#data = np.genfromtxt(fname='../spo/1.0.2/cor.dat') 
#
##for x in range(1,3):
#plt.plot(data[:,0],data[:,1],label='$\Re(C(t))$',lw=1)
#plt.plot(data[:,0],data[:,2],label='$\Im(C(t))$',lw=1)
#z = np.sqrt(data[:,1]**2+data[:,2]**2)
#plt.plot(data[:,0],z,label='$|C(t)|$',lw=1)

plt.legend() 
plt.show() 

