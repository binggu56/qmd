# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 16:11:08 2016

@author: bing
"""

import numpy as np 
import matplotlib.pyplot as plt 

import scipy 



def plot(ax1):

    x,psi,psi_approx = np.genfromtxt(fname='fit.out',unpack=True)
    
    #ax1.plot(data[:,0],data[:,2],linewidth=2,label='Potential')
    ax1.plot(x,psi,linewidth=2,label='Psi')
    ax1.plot(x,psi_approx,'--',linewidth=2,label='Fit')
      
    #ax1.legend() 
	#pl.legend(bbox_to_anchor=(0.5, 0.38, 0.42, .302), loc=3,ncol=1, mode="expand", borderaxespad=0.)
    #ax1.set_yscale('log')
    ax1.legend() 

def plot_p(ax):
    
    x,u = np.genfromtxt(fname='u.out',unpack=True)
    ax.plot(x,u,'.', markersize=12,label='U(x)')
    #ax.plot(x,du,'.', markersize=12,label='dU(x)')
    #ax.plot(x,-0.5 * (x/2.0 - 0),'--', lw=2,label='Exact')
    #ax.plot(x,-0.5 * (x**2/4.0 - 0.5),'o', lw=2,label='Exact')
    ax.set_xlim(-4,4)
    ax.set_ylim(0,2)
    
def plot_cor(ax):
    
    time,cor = np.genfromtxt(fname='cor.out',unpack=True,dtype=np.complex128)
    ax.plot(time, cor.imag,'-', lw=2,label='$Im C(t)$')
    

    #ax.legend(loc=0) 
def plt_traj():
    plt.subplot(111)
    #plt.ylim(-8,8)
    data = np.genfromtxt(fname='traj.out') 
    #data = np.loadtxt('traj.dat')
    for x in range(1,20):
        plt.plot(data[:,0],data[:,x],lw=1)
    
    #plt.figure(1) 
    #plt.plot(x,y1,'-')
    #plt.plot(x,y2,'g-')
    plt.xlabel('time')    
    
#plt.savefig('energy.pdf')
fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1) 
	
plot(ax1)
plot_cor(ax2)
plt.show()
plt_traj()
plt.show() 
