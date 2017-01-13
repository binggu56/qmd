# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 19:57:45 2016

Solve TDSE for a tight-binding model with non-resonant light 

@author: bingg
"""

# -*- coding: utf-8 -*-

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib  import cm
import matplotlib as mpl 

#import seaborn as sns 


norm = mpl.colors.Normalize(vmin=0, vmax=1, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.Blues)

unit_length = 2.0


def set_style():
    plt.style.use(['grayscale'])
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

def set_size(fig):
    fig.set_size_inches(8, 6)
    plt.tight_layout()




def Hamiltonian(t, Nsites):
    N = Nsites
    ta = 1.0 
    tb = 0.7 
    e1 = -2.0 
    e2 = 2.0  
    
    efactor = 1.0 
    
    fock0 = np.zeros((2*N,2*N))
    
    for i in range(N):
        fock0[2*i, 2*i] = e1
        fock0[2*i+1,2*i+1] = e2
    
    for i in range(N):
        fock0[2*i,2*i+1] = -ta
        fock0[2*i+1, 2*i] = -ta
        
    for i in range(N-1):
        fock0[2*i+2, 2*i+1]  = -tb
        fock0[2*i+1, 2*i+2] = -tb  
    
    fock1 = np.zeros((2*N,2*N))
    for i in range(N):
        fock1[2*i,2*i] = (i+1-float(N+1)/2.0) * unit_length                                                                       
        fock1[2*i+1,2*i+1] = fock1[2*i, 2*i] 
    
    fock = fock0 + efactor*field(t) * fock1
   
    return fock 
   
def diag(fock):
    eigs, eigvecs = np.linalg.eigh(fock)
    
    print(eigs[-3],eigvecs[0,-3],eigvecs[19,-3])
    
    xmin = -2
    xmax = 2*Nsite+2 
    dx = 1.0/(xmax-xmin)
    
    set_style()
    fig ,ax = plt.subplots()
    for j in range(2*N):
        for i in range(2*N):
 
            color = mapper.to_rgba(np.abs(eigvecs[i,j])**2)
            ax.axhline(y = eigs[j], xmin=float(i-xmin)*dx, xmax=float(i+1-xmin)*dx, linewidth=2,color=color)
    
    ax.set_xlim(xmin,xmax)


    ECP_L, ECP_R = ECP(N,t)
    plt.axhspan(ymin=-8,ymax=ECP_L,xmin=0.0,xmax=1.6*dx,color='k')
    plt.axhspan(ymin=-8,ymax=ECP_R,xmin=1.0-1.5*dx,xmax=1.0,linewidth=3,color='k')
    
    Z = np.linspace(0,1,5)
    mapper.set_array(Z)
    clb = plt.colorbar(mapper)
    clb.set_ticks([0,0.5,1])    
    #plt.plot(range(2*N),eigs,'.') 
    ax.set_xlabel('Sites')
    ax.set_ylabel('Energy')
    ax.set_ylim(-16,16)
    

    set_size(fig)
    plt.show()         

def field(t):
    
    sigma = 100.0 
    hbarinv = 1.51924888335207073622
    A = 0.5 # strength of the field
    T0 = 500.0 # center of time range for the field 
    phase = 0.0 
    freq = 0.10

    fieldtype = 'monochromatic' 
    
    if fieldtype == 'pulse': 
        phi = A * np.exp(-(t-T0)**2/(2*sigma**2))/(freq * hbarinv)    \
              * ( (t-T0)/(sigma**2) * np.sin( freq *(t-T0) * hbarinv + phase) \
              - freq * hbarinv * np.cos( freq * (t-T0) * hbarinv + phase) )
    elif fieldtype == 'monochromatic':
        omega = 1.0
        A = 4.0 
        phi = A * np.sin(omega * t)
        
    return phi

def ECP(N,t,bias=1.2):
    efactor = 1.0 
    E = field(t)
    left = bias - efactor * E * (float(N+1)/2.0) * unit_length
    right = -bias + efactor * E * (float(N+1)/2.0) * unit_length
    
    return left, right 
    

def basis_trans():
    """
    transform a matrix written in adiabatic eigenstates to site basis and vice versa  
    """
#def SOD(X, dt = 0.001):
#    """
#    second-order propagation of d/dt X = H(t) X 
#    """
#    
#    for kt in range(Nt):
#        k1 = dt * Hamiltonian(t)
#        k2 = dt * Hamiltonian()

def TDSE(Nt ,dt=0.001):
    
    Norbs = Nsites * 2 # two orbs per site 
    #psi = np.zeros((Norbs, Norbs), dtype=np.complex128)
    
    t = 0.0 
    
    # psi(t=0) is the eigenstates of H(t=0)
    #psi = diagH(0.0)
    
    # initial density matrix (full valence band, empty conduction band)
    rho = np.zeros((Norbs,Norbs), dtype = np.complex128)
    
    #rho = basis_trans(rho) # from eigenbasis to AOs 
    
    I = np.identity(Norbs, dtype=np.complex128) # propagator, at t0 is identity matrix  

    U = I 
    psi = np.array([1.0,0.0,0.0,0.0])
    
    f = open('population.dat','w')
    
    
    for k in range(Nt):

        t += dt 
        
        H = Hamiltonian(t, Nsites)
        U = np.linalg.solve(I + 1j * H * dt/2.0, np.dot(I - 1j * H * dt/2.0, U)) 
        f.write('{} {} {} \n'.format(t,*np.abs(np.dot(U,psi)[0:2])**2))

    f.close() 
        
    return U 
        
    


Nsites = 2 
Nt = 4000; dt = 0.005
 
U = TDSE(Nt,dt)
print(U[:,0])

dat = np.genfromtxt('population.dat')
plt.plot(dat[:,0], dat[:,1])
plt.plot(dat[:,0], dat[:,2])

plt.show() 
#T = np.linspace(300,500,400)
##for t in T:
##    print(t,field(t))
#plt.plot(T,field(T))
#Hamiltonian(377.20)
#Hamiltonian(458.90)
#Hamiltonian(417.794486216)
#Hamiltonian(438.345864662)

plt.show()
