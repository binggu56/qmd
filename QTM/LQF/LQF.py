# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:42:22 2016

@author: bing
"""
import numpy as np 
import numba 

Ndim = 1 

@numba.autojit
def derivs(x):
    v0 = x**2/2.0 
    dv = x
    return v0,dv 

@numba.autojit
def qpot(x,w):
    c = np.zeros(Ndim+1)
    c[0], c[1] =  -0.5, 0.0     
    
    S = np.zeros((Ndim+1,Ndim+1))
    
    for i in range(Ntraj):
        f = np.array([x[i],1.0])
        for m in range(Ndim+1):
            for n in range(Ndim+1):
                S[m,n] = w[i]*f[m]*f[n] + S[m,n]
    
    c = np.linalg.solve(S,c)
    
    u = np.zeros(Ntraj)
    fq = np.zeros(Ntraj)
    
    for i in range(Ntraj):
        r = c[0]*x[i] + c[1]
        
        # calculate quantum potential
        u[i] = - r**2/(2.0*am) - c[1]/(2.*am)

        fq[i] = r*c[0]/am

    return fq 

    
Ntraj = 1000 
x = np.random.randn(Ntraj)
p = np.zeros(Ntraj)
w = np.array([1./Ntraj]*Ntraj)
am = 1.0 

Nt = 1000
dt = 0.005 
t = 0.0 

f = open('traj.dat','w')
fmt =  ' {}' * 4  + '\n'    
print fmt

for k in range(Nt):
    t = t + dt 
    x = x + p*dt/am

    fq = qpot(x,w)
    
    for j in range(Ntraj):
        v0, dv = derivs(x[j])
        p[j] = p[j] - dv*dt + fq[j]*dt 
    
        
    f.write(fmt.format(t,x[0], x[1], x[2],x[3]))

f.close() 



    