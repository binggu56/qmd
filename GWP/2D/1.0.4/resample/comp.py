##!/usr/bin/python

import numpy as np
import pylab as plt

data = np.genfromtxt(fname='t100/wf.dat')
data1 = np.genfromtxt(fname='t300/wf.dat')
data2 = np.genfromtxt(fname='t500/wf.dat')
data3 = np.genfromtxt(fname='t600/wf.dat')
data00 = np.genfromtxt('../spo_1d/t100')
data01 = np.genfromtxt('../spo_1d/t300')
data02 = np.genfromtxt('../spo_1d/t500')
data03 = np.genfromtxt('../spo_1d/t600')

plt.subplot(2,2,1)
plt.xlim(0.5,2.5)

plt.title('t = 100 a.u.')
plt.plot(data[:,0],data[:,1],'r--',linewidth=2,label='LQF')
plt.plot(data00[:,0],data00[:,1],'k-',linewidth=2, label='Exact')
plt.xlabel('x')
plt.ylabel('$\psi^*\psi$')

plt.subplot(2,2,2)
plt.title('t = 300 a.u.')
plt.xlim(0.5,2.5)
plt.plot(data1[:,0],data1[:,1],'r--',linewidth=2)
plt.plot(data01[:,0],data01[:,1],'k-',linewidth=2)
plt.xlabel('x')
plt.ylabel('$\psi^*\psi$')

plt.subplot(2,2,3)
plt.title('t = 500 a.u.')
plt.xlim(0.5,2.5)
plt.plot(data2[:,0],data2[:,1],'r--',linewidth=2)
plt.plot(data02[:,0],data02[:,1],'k-',linewidth=2)
plt.xlabel('x')
plt.ylabel('$\psi^*\psi$')

plt.subplot(2,2,4)
plt.title('t = 600 a.u.')
plt.xlim(0.5,2.5)
plt.plot(data3[:,0],data3[:,1],'r--',linewidth=2)
plt.plot(data03[:,0],data03[:,1],'k-',linewidth=2)
plt.xlabel('x')
plt.ylabel('$\psi^*\psi$')

plt.savefig('wft.pdf')
plt.show() 



