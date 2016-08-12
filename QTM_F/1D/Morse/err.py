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
plt.subplot(111)
#plt.ylim(-8,8)
dat = np.genfromtxt(fname='err.dat') 
plt.plot(dat[:,0],dat[:,1],lw=1,label='Error p')
plt.plot(dat[:,0],dat[:,2],lw=1,label='Error r')

#plt.ylim(0,5)
#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
plt.xlabel('Time')
#plt.ylabel('position')
#plt.title('traj')

#plt.subplot(212)
#data = np.genfromtxt(fname='../2.0.2/traj.dat') 
#data = np.loadtxt('traj.dat')


plt.legend()

plt.savefig('err.pdf')

plt.show() 



