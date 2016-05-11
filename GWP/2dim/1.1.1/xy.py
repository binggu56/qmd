##!/usr/bin/python

import numpy as np
import pylab as plt
import seaborn as sns 

sns.set_context('poster')

plt.subplot(1,1,1)
dat = np.genfromtxt(fname='xy0.dat') 
dat1 = np.genfromtxt(fname='xyt.dat') 

#for x in range(1,data.shape[-1]):
    
plt.plot(dat[:,0],dat[:,1],'ko')
plt.plot(dat1[:,0],dat1[:,1],'ro')
#plt.plot(data[:,0],data[:,3],lw=2,label='sym')
#plt.plot(data[:,0],data[:,5],lw=2,label='asym')
#plt.plot(data[:,0],data[:,3],lw=2,label='$|C(t)|$')

#dat = np.genfromtxt(fname='cor_backup.dat')
#plt.plot(dat[:,0]*2.0,dat[:,1],'--',label='Re, QM')
#plt.plot(dat[:,0]*2.0,dat[:,2],'--',label='Im, QM')

#x = np.linspace(0,4,100) 
#y = -np.sin(x)
#plt.plot(x,y,lw=2,label='sin(x)')
plt.xlabel('$Time$')
plt.ylabel('$C(t)$')
#plt.title('traj')

#plt.subplot(2,1,2)
#dat = np.genfromtxt(fname='/home/bing/gwp/spo_2d/1.0.3/cor1') 

#for x in range(1,3):
#plt.plot(dat[:,0],dat[:,1],'k--',label='$\Re(C(t)),~ QM$',lw=2)
#plt.plot(dat[:,0],dat[:,2],'r--',label='$\Im(C(t)),~ QM$',lw=2)
#z = np.sqrt(data[:,1]**2+data[:,2]**2)
#plt.plot(data[:,0],z,label='$|C(t)|$',lw=1)
plt.ylim(-4,1) 
plt.legend()
#plt.xlim(0,6)
plt.savefig('cor.pdf')
plt.show() 

