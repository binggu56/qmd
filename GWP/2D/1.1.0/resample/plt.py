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
plt.subplot(2,1,1)
data = np.genfromtxt(fname='q.dat') 
#data = np.loadtxt('traj.dat')
for x in range(1,data.shape[-1]):
    plt.plot(data[:,0],data[:,x])

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
plt.xlabel('time')
plt.ylabel('$x_i$')
plt.title('traj')

plt.subplot(2,1,2)
data = np.genfromtxt(fname='c.dat') 
#data = np.loadtxt('traj.dat')
for x in range(1,10):
    plt.plot(data[:,0],data[:,x])
plt.xlabel('time')



plt.show() 

