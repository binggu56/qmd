

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

plt.plot(dat[:,0],dat[:,1],'--',label='<x>')
#plt.plot(dat[:,0],dat[:,2],'--',label='GWP Center')

#dat = np.genfromtxt('../../SPO/spo_1d/xAve.out')
#plt.plot(dat[:,0],dat[:,1],'-',label='QM')

plt.legend() 
plt.show() 


def plot_corr():
    plt.subplot(111)
    
    dat = np.genfromtxt('corr.out')
    
    
    plt.plot(dat[:,0],dat[:,1],'--',label='GHQT')
    plt.plot(dat[:,0],dat[:,2],'-')
    
#    dat = np.genfromtxt('../../SPO/spo_1d/corr.out')
#    plt.plot(dat[:,0],dat[:,1],'-',label='QM')
    
    plt.legend() 
    plt.show() 

def plot_coeff():
    plt.subplot(111)
    
    dat = np.genfromtxt('coeff.dat',dtype=np.complex128)
    
    for i in range(1,dat.shape[-1]):
        plt.plot(dat[:,0],dat[:,i])

    
    plt.legend() 
    plt.show() 


def plot_wft():   
    plt.subplot(111)

    dat = np.genfromtxt('wft.dat')
    #plt.plot(dat[:,0],dat[:,1],'--',label='$Re \psi$')
    #plt.plot(dat[:,0],dat[:,2],'-',label='wft')
    plt.plot(dat[:,0],(dat[:,1]**2+dat[:,2]**2),'r--',lw=2,label='GH')
    
    #x = np.linspace(-6,6)
    #t = 0.4 
    #plt.plot(x, (1./np.pi)**0.25 * np.exp(-x**2/2.) * np.cos(t / 2.0) , label='Exact')
    dat = np.genfromtxt('../../SPO/spo_1d/dwell/wft.dat')
    plt.plot(dat[:,0],dat[:,1],'k-',label='QM')
    #plt.plot(dat[:,0],dat[:,2],'--',label='QM')
    
    #plt.plot(dat[:,0],dat[:,3] ** 2 ,'k-',lw=3,label='QM')
    
    
    plt.legend()
    plt.show() 

plot_corr()
plot_wft()
plot_coeff() 