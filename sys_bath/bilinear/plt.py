

import numpy as np 
import matplotlib.pyplot as plt 

#plt.subplot(211)

#dat = np.genfromtxt('traj.dat')

#for i in range(dat.shape[-1]-1):
#    plt.plot(dat[:,0],dat[:,i+1])

def xAve(ax):

    dat = np.genfromtxt('dat1/xAve.dat')
    ax.plot(dat[:,0],dat[:,1],'-',lw=1,label='Nfit=2')

#   dat = np.genfromtxt('dat10/xAve.dat')
#   ax.plot(dat[:,0],dat[:,1],'-',label='linear')
    
    #dat = np.genfromtxt('dat11/xAve.dat')
#ax.plot(dat[:,0],dat[:,1],'-',label='cubic')
    
    dat = np.genfromtxt('/home/bingg/pyQMD/SPO/spo_2d/1.0.2/xave.dat')     
    ax.plot(dat[:,0],dat[:,1],'-',lw=2,label='Exact')

def cy(ax):
    dat = np.genfromtxt('dat4/c.dat',dtype=np.complex128)
    x = dat[:,0]
    y = dat[:,1].real
    p = np.polyfit(x,y,1)
    print(p)
    ax.plot(x,dat[:,1].real,'.',label='real')
    ax.plot(x,np.poly1d(p)(x),'.',lw=2)
    #ax.plot(dat[:,0],dat[:,1].imag,'.',label='imag')
    

fig, ax = plt.subplots()

xAve(ax) 
#cy(ax)
plt.legend()
plt.show() 
