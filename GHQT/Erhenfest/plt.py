

import numpy as np 
import matplotlib.pyplot as plt 

#plt.subplot(211)
#
#dat = np.genfromtxt('traj.dat')
#
#for i in range(dat.shape[-1]-1):
#    plt.plot(dat[:,0],dat[:,i+1])


plt.subplot(111)

dat = np.genfromtxt('xAve.dat')

#plt.plot(dat[:,0],dat[:,1],'-',label='<x>')
plt.plot(dat[:,0],dat[:,2],'--',label='GWP Center')

dat = np.genfromtxt('../../SPO/spo_1d/xAve.out')
plt.plot(dat[:,0],dat[:,1],'-',label='QM')

plt.legend() 
plt.show() 


def plot_corr():
    plt.subplot(111)
    
    dat = np.genfromtxt('corr.out')
    
    
    plt.plot(dat[:,0],dat[:,1],'--',label='GHQT')
    plt.plot(dat[:,0],dat[:,2],'-')
    
    dat = np.genfromtxt('../../SPO/spo_1d/corr.out')
    plt.plot(dat[:,0],dat[:,1],'-',label='QM')
    
    plt.legend() 
    plt.show() 
    
def plot_wft():   
    plt.subplot(111)

    dat = np.genfromtxt('wft.dat')
    #plt.plot(dat[:,0],dat[:,1],'-',label='wft')
    #plt.plot(dat[:,0],dat[:,2],'-',label='wft')
    plt.plot(dat[:,0],(dat[:,1]**2+dat[:,2]**2),'b--',lw=2,label='GH')
    
    dat = np.genfromtxt('../../SPO/spo_1d/wft_t16.out')
    #plt.plot(dat[:,0],dat[:,1],'--',label='QM')
    #plt.plot(dat[:,0],dat[:,2],'--',label='QM')
    
    plt.plot(dat[:,0],dat[:,3] ** 2 ,'k-',lw=3,label='QM')
    
    
    plt.legend()
    plt.show() 

plot_corr()
plot_wft()