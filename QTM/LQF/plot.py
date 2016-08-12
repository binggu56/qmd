# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:46:52 2016

@author: bing
"""

import numpy as np 

import matplotlib.pyplot as plt 


dat = np.genfromtxt('traj.dat')

#for i in range(1,10):
#    plt.plot(dat[:,0],dat[:,i])


def plt_e():

    plt.figure()     
    dat = np.genfromtxt('en.out')

    for i in range(1,5):
        plt.plot(dat[:,0],dat[:,i])

def plt_rMSE():

    plt.figure()     
    dat = np.genfromtxt('rMSE.out')

    plt.plot(dat[:,0],dat[:,1],label='rMSE')        



plt_e()
plt.legend() 
plt.show() 