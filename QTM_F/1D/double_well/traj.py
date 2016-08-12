#!/usr/bin/python

import numpy as np
import pylab as plt

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
def plt_traj():
    plt.subplot(111)
    #plt.ylim(-8,8)
    data = np.genfromtxt(fname='traj.dat') 
    #data = np.loadtxt('traj.dat')
    for x in range(1,20):
        plt.plot(data[:,0],data[:,x],lw=1)
    
    #plt.figure(1) 
    #plt.plot(x,y1,'-')
    #plt.plot(x,y2,'g-')
    plt.xlabel('time')
    #plt.ylabel('position')
    #plt.title('traj')
    
    #plt.subplot(212)
    #data = np.genfromtxt(fname='../2.0.2/traj.dat') 
    #data = np.loadtxt('traj.dat')
def plt_r(ax):
    
    x, r, p = np.genfromtxt('r.out', unpack=True, usecols=(0,1,2))
    #ax.set_xlim(-4,4)
    #ax.set_ylim(-10,10)
    ax.plot(x,r,'o')
    ax.plot(x,p,'.')
    
def plt_ave(ax):
     
    data = np.genfromtxt(fname='xAve.out') 
    ax.plot(data[:,0],data[:,1],lw=2)  
    
def plt_MSE(ax):
    
    data = np.genfromtxt(fname='rMSE.out')
    ax.plot(data[:,0], data[:,1],lw=2)
    ax.plot(data[:,0], data[:,2],'--',lw=2)


fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1) 


#plt_ave(ax1)
#plt_MSE(ax2) 
#plt_r(ax1)

plt.show()

plt_traj() 
 



