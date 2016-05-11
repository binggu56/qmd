##!/usr/bin/python

import numpy as np
import pylab as plt
import seaborn as sns 

sns.set_context('poster')

#with open("traj.dat") as f:
#    data = f.read()
#
#    data = data.split('\n')
#
#    x = [row.split(' ')[0] for row in data]
#    y = [row.split(' ')[1] for row in data]
#
#    fig = plt.figure()
#
#    ax1 = fig.add_subplot(111)
#
#    ax1.set_title("Plot title...")    
#    ax1.set_xlabel('your x label..')
#    ax1.set_ylabel('your y label...')
#
#    ax1.plot(x,y, c='r', label='the data')
#
#    leg = ax1.legend()
#fig = plt.figure()
f, (ax1, ax2) = plt.subplots(2, sharex=True)
#f.subplots_adjust(hspace=0.1)
#plt.subplot(211)
ax1.set_ylim(0,4)
data = np.genfromtxt(fname='q.dat') 
#data = np.loadtxt('traj.dat')
for x in range(1,data.shape[1]):
    ax1.plot(data[:,0],data[:,x], linewidth=1)

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#plt.xlabel('time')
ax1.set_ylabel('position [bohr]')
#plt.title('traj')

#plt.subplot(212)
data = np.genfromtxt(fname='c.dat') 
#data = np.loadtxt('traj.dat')
for x in range(1,16):
    ax2.plot(data[:,0],data[:,x], linewidth=1)
ax2.set_xlabel('time [a.u.]')
ax2.set_ylabel('$|c_i|$')
#plt.ylim(-0.2,5)

#plt.subplot(2,2,3)
#data = np.genfromtxt(fname='norm')
#plt.plot(data[:,0],data[:,1],'r-',linewidth=2)
#plt.ylim(0,2)

#plt.subplot(2,2,4)
#data = np.genfromtxt(fname='wf.dat')
#data1 = np.genfromtxt(fname='wf0.dat')
#data0 = np.genfromtxt('../spo_1d/t500')
#plt.plot(data[:,0],data[:,1],'r--',linewidth=2)
#plt.plot(data0[:,0],data0[:,1],'k-',linewidth=2)
#plt.plot(data1[:,0],data1[:,1],'k-.',linewidth=2)
#plt.title('t=100')

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#plt.xlim(0.8,2.1)
#plt.xlabel('x')
#plt.ylabel('$\psi^*\psi$')

plt.savefig('traj.pdf')
plt.show() 



