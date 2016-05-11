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
data = np.genfromtxt(fname='wft.dat') 
dat = np.genfromtxt(fname='/home/bing/spo/spo_1d/wft') 
#data = np.loadtxt('traj.dat')
plt.plot(data[:,0],data[:,1],'r-',lw=2)
plt.plot(dat[:,0],dat[:,1],'k-',lw=2)

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')


plt.show() 

