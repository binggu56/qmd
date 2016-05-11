##!/usr/bin/python

import numpy as np
import pylab as plt
import seaborn as sns 

sns.set_context('poster')

#plt.subplot(1,1,1)
data = np.genfromtxt(fname='wft.dat') 
data1 = np.genfromtxt(fname='../../spo/1.0.3/wft3.dat') 
#data0 = np.genfromtxt('../spo_1d/t100')
plt.plot(data[:,0],data[:,1],linewidth=2, label='$|\psi(x,t)|^2$, GWP')
plt.plot(data1[:,0],data1[:,1],'--',linewidth=2,label='QM')
#plt.scatter(data0[:,0],data0[:,1],label='QM')
#plt.plot(data1[:,0],data1[:,1],'k-.',linewidth=2, label='t0')


#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#plt.xlim(0.5,2.4)
plt.xlabel('$x$')
#plt.ylabel('$|\psi(x,t)|^2$')
#plt.title('traj')

#plt.subplot(2,1,2)
#data = np.genfromtxt(fname='c.dat') 
#data = np.loadtxt('traj.dat')
#for x in range(1,10):
#    plt.plot(data[:,0],data[:,x])
#plt.xlabel('time')

plt.savefig('wft.pdf')
plt.legend()
plt.show() 

