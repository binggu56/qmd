

import numpy as np 
import matplotlib.pyplot as plt 

#plt.subplot(211)

#dat = np.genfromtxt('traj.dat')

#for i in range(dat.shape[-1]-1):
#    plt.plot(dat[:,0],dat[:,i+1])

def xAve(ax):

    dat = np.genfromtxt('coherent_nonlinear_n8m2/xAve.dat')
    ax.plot(dat[:,0],dat[:,1],'-',label='N8M2')

    dat = np.genfromtxt('coherent_nonlinear2/xAve.dat')
    ax.plot(dat[:,0],dat[:,1],'-',label='N16M3')
    
    dat = np.genfromtxt('/home/bingg/pyQMD/SPO/spo_2d/1.0.1/xave.dat')     
    ax.plot(dat[:,0],dat[:,1],'--',label='Exact')



fig, ax = plt.subplots()

xAve(ax) 
plt.legend()
plt.show() 
