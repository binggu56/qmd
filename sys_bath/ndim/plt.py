

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl 

#plt.subplot(211)

#dat = np.genfromtxt('traj.dat')

#for i in range(dat.shape[-1]-1):
#    plt.plot(dat[:,0],dat[:,i+1])
def set_style():
    plt.style.use(['classic'])
    #print(plt.style.available)
#    font = {'family' : 'Times',
#        'weight' : 'bold',
#        'size'   : 22}

    mpl.rc('font',**{'family':'serif','sans-serif':['Times']})
    mpl.rc('text', usetex=True)

    
    mpl.rc('xtick', labelsize=16) 
    mpl.rc('ytick', labelsize=16)    
    
    
    
    SIZE = 16
    plt.rc('font', size=SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=SIZE)  # fontsize of the x any y labels
    plt.rc('xtick', labelsize=SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SIZE)  # legend fontsize
    plt.rc('figure', titlesize=SIZE)  # # size of the figure title
    
def xAve(ax):

    dat = np.genfromtxt(fname=fname+'/xAve.dat')
    ax.plot(dat[:,0],dat[:,1],'--',lw=2,label='11')

#dat = np.genfromtxt('../1.1.2/dat4/xAve.dat')
#   ax.plot(dat[:,0],dat[:,1],'-',label='nfit=3')
    
#   dat = np.genfromtxt('dat11/xAve.dat')
#   ax.plot(dat[:,0],dat[:,1],'-',label='cubic')
    
    dat = np.genfromtxt('/home/bingg/pyQMD/SPO/spo_2d/1.0.4/xave.dat')     
    ax.plot(dat[:,0],dat[:,1],'k',lw=2,label='Exact')
    ax.set_ylabel('$< x > $')
    ax.set_yticks([-0.4,0,0.4])
    #plt.savefig('xAve.pdf')

def cy(ax):
    dat = np.genfromtxt('dat2/c.dat',dtype=np.complex128)
    x = dat[:,0]
    y = dat[:,1].real
    p = np.polyfit(x,y,1)
    print(p)
    ax.plot(x,dat[:,1].real,'.',label='real')
    ax.plot(x,np.poly1d(p)(x),'.',lw=2)
    
    ax.plot(dat[:,0],dat[:,-1].real,'.',label='imag')
    
def den(ax):
    dat = np.genfromtxt(fname=fname+'/den.dat',dtype=np.complex128)
    #ax.plot(dat[:,0], dat[:,1].real)
    ax.plot(dat[:,0], dat[:,1].imag,'--',lw=2)
    
#    dat = np.genfromtxt(fname='../1.1.2/dat4/den.dat',dtype=np.complex128)
#ax.plot(dat[:,0], dat[:,1].real)
#   ax.plot(dat[:,0], dat[:,1].imag,'--')

    # exact results

    dat = np.genfromtxt('/home/bingg/pyQMD/SPO/spo_2d/1.0.4/den.dat')
    #ax.plot(dat[:,0], dat[:,1],label='real, QM')
    ax.plot(dat[:,0], dat[:,2],'k',lw=2)

    # results without coupling
    x = np.linspace(0,20,200)
    ax.plot(x,np.exp(-0.5) / 2.0 * np.sin(x),'r--',linewidth=1)
    ax.set_ylabel(r'$\rho_{01}$')
    ax.set_yticks([-0.4,0.0,0.4])
    #plt.savefig('den.pdf')
    ax.set_xlim(0,20) 

set_style()

fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
fig.subplots_adjust(hspace=0)

fname = '.'

xAve(ax1) 
#cy(ax)
den(ax2)
#plt.legend()
plt.xlabel('Time [a.u.]')
plt.savefig('nonlinear.pdf')
plt.show() 
