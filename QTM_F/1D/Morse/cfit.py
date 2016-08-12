#!/usr/bin/python

import numpy as np
import pylab as plt
import seaborn as sns

sns.set_context("poster")
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
plt.subplot(221)
#plt.ylim(-8,8)
data = np.genfromtxt(fname='cfit.dat') 
#data = np.loadtxt('traj.dat')
for x in range(1,6):
    plt.plot(data[:,0],data[:,x],lw=1)

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
plt.xlabel('time')
#plt.ylabel('position')
#plt.title('traj')
plt.ylim(-5,5)

plt.subplot(222)
#data = np.genfromtxt(fname='../2.0.2/traj.dat') 
#data = np.loadtxt('traj.dat')
for x in range(1,6):
    plt.plot(data[:,0],data[:,x+5],lw=1)
#plt.xlabel('time') 
#ax2.yaxis.tick_right()
#ax2.yaxis.set_ticks_position('both')
plt.ylim(-5,5)

plt.subplot(223)
#data = np.genfromtxt(fname='norm')
#for x in range(1,5):
plt.plot(data[:,0],data[:,1+10],lw=1,label='$x$')
plt.plot(data[:,0],data[:,2+10],lw=1,label='$x^2$')
plt.plot(data[:,0],data[:,3+10],lw=1,label='$x^3$')
plt.plot(data[:,0],data[:,4+10],lw=1,label='$1$')
plt.legend()
#plt.ylabel('Norm')
plt.ylim(-2,2)

plt.subplot(224)
#data = np.genfromtxt(fname='norm')
for x in range(1,5):
    plt.plot(data[:,0],data[:,x+14],lw=1)
plt.ylim(-2,2)

#plt.legend()
plt.savefig('cfit.pdf')

plt.show() 



