# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:46:52 2016

@author: bing
"""

import numpy as np 

import matplotlib.pyplot as plt 


dat = np.genfromtxt('traj.dat')

for i in range(1,4):
    plt.plot(dat[:,0],dat[:,i])

plt.show() 