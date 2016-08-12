# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 09:35:39 2016

@author: bing
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:42:22 2016

@author: bing
"""
import numpy as np
import scipy  
import numba 
import sys 




bohr_angstrom = 0.52917721092
hartree_wavenumber = 219474.63 

#hartree_wavenumber = scipy.constants.value(u'hartree-inverse meter relationship') / 1e2 


Vmin = -24.2288 

b = np.array([-6.631e-02, 1.346e-01, -3.300e-02, 6e0, -1.4e01, -1.193e02, 2.290e02, \
            1.110e03, -1.850e03, -3.5e03, 6.0e03])



@numba.autojit
def derivs(x):
    """
    Morse potential     
    """ 
    
    PES = 'pH2' 
    
    if PES == 'Morse':
        
        a, x0 = 1.02, 1.4 
        De = 0.176 / 100.0 
    
        d = (1.0-np.exp(-a*x))
        
        v0 = De*d**2
            
        dv = 2. * De * d * a * np.exp(-a*x)
        
    elif PES == 'HO':
        
        v0 = x**2/2.0 
        dv = x 

        #ddv = 2.0 * De * (-d*np.exp(-a*((x-x0)))*a**2 + (np.exp(-a*(x-x0)))**2*a**2)

    elif PES == 'pH2':
        
        dx = 1e-4
        
        v0 = np.zeros(Ntraj)
        dv = np.zeros(Ntraj)
        
        for i in range(Ntraj):
            v0[i] = vpot(x[i])
            dv[i] = ( vpot(x[i] + dx) - v0[i])/dx
        
    return v0,dv

@numba.autojit
def qpot(x,p,r,w):

    """
    Linear Quantum Force : direct polynomial fitting of derivative-log density (amplitude)    
    curve_fit : randomly choose M points and do a nonlinear least-square fitting to a 
            predefined functional form      
    """
    
    #tau = (max(xdata) - min(xdata))/(max(x) - min(x))
    #if tau > 0.6:
    #    pass 
    #else: 
    #    print('Data points are not sampled well.'
    
    Nb = 6
    S = np.zeros((Nb,Nb))
    
    for j in range(Nb):
        for k in range(Nb):
            S[j,k] = np.dot(x**(j+k), w)  
    
    bp = np.zeros(Nb)
    br = np.zeros(Nb)
    
    for n in range(Nb):
        bp[n] = np.dot(x**n * p, w)
        br[n] = np.dot(x**n * r, w)
        
        
    cp = np.linalg.solve(S,bp)
    cr = np.linalg.solve(S,br)  

    #unit = np.identity(Nb)
    #r_approx = cr[0] * unit + cr[1] * x + cr[2] * x**2 + cr[3] * x**3 
    #p_approx = cp[0] * unit + cp[1] * x + cp[2] * x**2 + cp[3] * x**3

    dr = cr[1] + 2. * cr[2] * x + 3. * cr[3] * x**2 + 4.0 * cr[4] * x**3 + 5.0 * cr[5] * x**4 
    dp = cp[1] + 2. * cp[2] * x + 3. * cp[3] * x**2 + 4.0 * cp[4] * x**3 + 5.0 * cp[5] * x**4
    
    ddr = 2. * cr[2] + 6. * cr[3] * x + 12.0 * cr[4] * x**2 + 20.0 * cr[5] * x**3
    ddp = 2. * cp[2] + 6. * cp[3] * x + 12.0 * cp[4] * x**2 + 20.0 * cp[5] * x**3
    
    fr =  -1./2./am * (2. * r * dp + ddp)
    fq = 1./2./am * (2. * r * dr + ddr)  
    
    Eu = -1./2./am * np.dot(r**2 + dr,w)
        
    return Eu,fq,fr 


@numba.autojit     
def sym(V):

    n = V.shape[-1] 
    
    for i in range(n):
        for j in range(i):
            V[j,i] = V[i,j] 
    return V 

@numba.autojit 
def vpot(r):
    
    re = 3.47005
    De = 24.2288 
	
    r = r * bohr_angstrom 
    
    beta_inf = np.log(2.0 * De / u_LR(re)) 
    
    s = 0.0        
    for j in range(11):
        s += b[j] * y_ref(r,1)**j
    
      
    beta = y_ref(r,6) * beta_inf + (1.0 - y_ref(r,6))  * s  
    
    vpot = De * (1.0 - u_LR(r)/u_LR(re) * np.exp(- beta * y_eq(r,6)))**2
    
    vpot = vpot + Vmin 
    
    vpot = vpot / hartree_wavenumber 
    
    return vpot 
	
@numba.autojit 
def y_eq(r,n):
    
    re = 3.47005
     
    y_eq = (r**n - re**n)/(r**n + re**n) 
 
    return y_eq 
    
@numba.autojit    
def y_ref(r,n):
    
    r_ref = 4.60
     
    z = (r**n - r_ref**n)/(r**n + r_ref**n)     

    return z

@numba.autojit      
def u_LR(r):
    
    C6 = 5.820364e04
    C8 = 2.87052154e05 
    C10 = 1.80757343e06 
    
    z = damp(r,6) * C6/r**6 + damp(r,8) * C8/r**8 + damp(r,10) * C10 / r**10 
      
    return z
    


@numba.autojit    	
def damp(r,n):
    
   den = 1.10 
		 
   z = (1.0 - np.exp(-3.30 * den * r / n - 0.423 * (den * r)**2/np.sqrt(float(n))))**(n-1) 
   
   return z 


# initialization    
Ntraj = 4000 
a0 = 0.25  
x0 = 9.0  


x = np.random.randn(Ntraj) 

#x = np.zeros(Ntraj)  
#for k in range(Ntraj):
#    x[k] = np.random.randn() 
#    while x[k] > 3.0:
#        x[k] = np.random.randn()
    
        
x = x / np.sqrt(2.0 * a0) + x0 
print('x range',np.min(x),np.max(x))
p = np.zeros(Ntraj)
r = - a0 * (x-x0) 

w = np.array([1./Ntraj]*Ntraj)
am = 1837.0 




f_MSE = open('rMSE.out','w')
nout = 20       # number of trajectories to print 
fmt =  ' {}' * (nout+1)  + '\n'  
Eu = 0.  

Ndim = 1        # dimensionality of the system  
fric_cons = 0.0008      # friction constant  


def evolve(x,p,r): 
    
    Nt = 80000
    dt = 0.5
    dt2 = dt/2.0 
    t = 0.0 
    
    f = open('traj.dat','w')
    fe = open('en.out','w')
    
    v0, dv = derivs(x)
    Eu,fq,fr = qpot(x,p,r,w)
    
    for k in range(Nt):
        t = t + dt 
    
        p += (- dv + fq) * dt2 - fric_cons * p * dt2   
        r += fr * dt2
        
        x +=  p*dt/am
    
        # force field 
        Eu, fq, fr = qpot(x,p,r,w)
        if Eu < 0:
            print('Error: U = {} should not be negative. \n'.format(Eu))
            #sys.exit()
            
        v0, dv = derivs(x)
    
        p += (- dv + fq) * dt2 - fric_cons * p * dt2 
        r += fr * dt2 
           
        f.write(fmt.format(t,*x[0:nout]))
        Ek = np.dot(p*p,w)/2./am  * hartree_wavenumber
        Ev = np.dot(v0,w) * hartree_wavenumber
        Eu = Eu * hartree_wavenumber
        Etot = Ek + Ev + Eu
        
        fe.write('{} {} {} {} {} \n'.format(t,Ek,Ev,Eu,Etot))
        
        if k == Nt-1:
            print('The total energy = {} cm-1. \n'.format(Etot))
    
    fe.close()
    f.close() 

evolve(x,p,r)


#a, x0, De = 1.02, 1.4, 0.176/100 
#print('The well depth = {} cm-1. \n'.format(De * hartree_wavenumber))
#
#omega  = a * np.sqrt(2. * De / am )
#E0 = omega/2. - omega**2/16./De
#dE = (Etot-E0) * hartree_wavenumber 
#print('Exact ground-state energy = {} Hartree. \nEnergy deviation = {} cm-1. \n'.format(E0,dE))
#    


    
