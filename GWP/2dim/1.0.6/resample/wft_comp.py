##!/usr/bin/python

import numpy as np
import pylab as plt

data = np.genfromtxt(fname='dat1/wft.dat')
data1 = np.genfromtxt(fname='dat2/wft.dat')
data2 = np.genfromtxt(fname='dat3/wft.dat')
data3 = np.genfromtxt(fname='dat4/wft.dat')
data4 = np.genfromtxt(fname='dat5/wft.dat')
data5 = np.genfromtxt(fname='dat6/wft.dat')
data6 = np.genfromtxt(fname='dat7/wft.dat')
data7 = np.genfromtxt(fname='dat8/wft.dat')

plt.subplot(111)

plt.xlim(0.5,2.5)

plt.plot(data[:,0],data[:,1],'r--',linewidth=2,label='LQF')
plt.plot(data1[:,0],data1[:,1],'r--',linewidth=2)
plt.plot(data2[:,0],data2[:,1],'r--',linewidth=2)
plt.plot(data3[:,0],data3[:,1],linewidth=2)
plt.plot(data4[:,0],data4[:,1],linewidth=2)
plt.plot(data5[:,0],data5[:,1],linewidth=2)
plt.plot(data6[:,0],data6[:,1],linewidth=2)
plt.plot(data7[:,0],data7[:,1],linewidth=2)
#plt.plot(data03[:,0],data03[:,1],'k-',linewidth=2)



plt.xlim(0.5,2.5)
plt.xlabel('x')
plt.ylabel('$\psi^*\psi$')

plt.savefig('wft.pdf')
plt.show() 



