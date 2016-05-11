##!/usr/bin/python

import numpy as np
import pylab as plt
import seaborn as sns 

sns.set_context('poster')

plt.subplot(1,1,1)
data = np.genfromtxt(fname='cor.dat') 

#for x in range(1,data.shape[-1]):
    
plt.plot(data[:,0],data[:,1],lw=2,label='Re')
plt.plot(data[:,0],data[:,2],lw=2,label='Im')
#plt.plot(data[:,0],data[:,3],lw=2,label='$|C(t)|$')

#dat = np.genfromtxt(fname='../spo/1.0.2/corr')
#plt.plot(dat[:,0],dat[:,1],'--',label='Re, QM')
#plt.plot(dat[:,0],dat[:,2],'--',label='Im, QM')

#x = np.linspace(0,4,100) 
#y = -np.sin(x)
#plt.plot(x,y,lw=2,label='sin(x)')
plt.xlabel('$Time$')
plt.ylabel('$C(t)$')
#plt.title('traj')

#plt.subplot(2,1,2)
dat = np.genfromtxt(fname='/home/bing/gwp/spo_2d/1.0.0/cor1') 

#for x in range(1,3):
plt.plot(dat[:,0],dat[:,1],'--',label='$\Re(C(t))$',lw=2)
plt.plot(dat[:,0],dat[:,2],'--',label='$\Im(C(t))$',lw=2)
#z = np.sqrt(data[:,1]**2+data[:,2]**2)
#plt.plot(data[:,0],z,label='$|C(t)|$',lw=1)
plt.ylim(-1,1) 
plt.legend()
plt.savefig('cor.pdf')
plt.show() 

