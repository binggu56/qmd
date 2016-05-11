##!/usr/bin/python

import numpy as np
import pylab as plt

plt.subplot(1,1,1)
data = np.genfromtxt(fname='wf.dat') 
data1 = np.genfromtxt(fname='wf0.dat') 
#data0 = np.genfromtxt('../spo_1d/t300')
plt.plot(data[:,0],data[:,1],'r--',linewidth=2)
#plt.plot(data0[:,0],data0[:,1],'k-',linewidth=2)
plt.plot(data1[:,0],data1[:,1],'k-.',linewidth=2)


#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
plt.xlim(-8,8)
plt.xlabel('x')
plt.ylabel('$\psi(x,t)$')
#plt.title('traj')

#plt.subplot(2,1,2)
#data = np.genfromtxt(fname='c.dat') 
#data = np.loadtxt('traj.dat')
#for x in range(1,10):
#    plt.plot(data[:,0],data[:,x])
#plt.xlabel('time')



plt.show() 

