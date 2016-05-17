# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:42:22 2016

@author: bing
"""
import numpy as np 
import numba 
from math import * 



@numba.autojit
def derivs(x):
    
    g = 0.4 
    
    v0 = x**2/2.0 + g * x**4 / 4.0 
    dv = x + g * x**3 
    
    return v0,dv 

@numba.autojit
def qpot(x,w):
    c = np.array([-0.5,0.0])    
    
    S = np.zeros((Ndim+1,Ndim+1))
    
    for i in range(Ntraj):

        f = np.array([x[i],1.0])
        
        for m in range(Ndim+1):
            for n in range(m+1):
                
                S[m,n] += w[i] * f[m] * f[n]
    S = sym(S)
    
    c = np.linalg.solve(S,c)
    
    u = np.zeros(Ntraj)
    fq = np.zeros(Ntraj)
    
    for i in range(Ntraj):
        
        r = c[0]*x[i] + c[1]
        
        # calculate quantum potential
        u[i] = - r**2/(2.0*am) - c[1]/(2.*am)

        fq[i] = r*c[0]/am

    return fq 

def sym(V):

    n = V.shape[-1] 
    
    for i in range(n):
        for j in range(i):
            V[j,i] = V[i,j] 
    return V 
    
Ntraj = 1000 
a = 1.0
x0 = -1.0  
x = np.random.randn(Ntraj)
x = x / sqrt(2.0 * a) + x0 

p = np.zeros(Ntraj)
w = np.array([1./Ntraj]*Ntraj)
am = 1.0 

Nt = 3000
dt = 0.004 
dt2 = dt/2.0 
t = 0.0 

f = open('traj.dat','w')
nout = 10       # number of trajectories to print 
fmt =  ' {}' * (nout+1)  + '\n'   

Ndim = 1        # dimensionality of the system  
fr = 1.6       # friction constant  

v0, dv = derivs(x)
fq = qpot(x,w)

for k in range(Nt):
    t = t + dt 

    p += (- dv + fq - fr * p)*dt2 
    
    x = x + p*dt/am

    # force field 
    fq = qpot(x,w)
    v0, dv = derivs(x)

    p += (- dv + fq - fr * p)*dt2 
           
    f.write(fmt.format(t,*x[0:nout]))

f.close() 



    