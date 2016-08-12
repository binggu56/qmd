# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:42:22 2016

@author: bing
"""
import numpy as np 
import numba 
from math import * 
import sys 

from scipy.optimize import curve_fit

@numba.autojit
def derivs(x):
    
    g = 0.
    
    v0 = x**2/2.0 + g * x**4 / 4.0 
    dv = x + g * x**3 
    
    return v0,dv 
    
def func(x,a,b,c,d):
    return a * x**3 + b * x**2 + c*x + d   

def df(x,a,b,c,d):
    
    return b * 2. * x + c + 3. * a * x**2, 2. * b + 6. * a * x
      

def pade(x,a,b,b1):
    return (a + b * x) / (1. + b1 * x)

def df_pade(x,a,b,b1):
    p = a + b * x     
    dp = b
    ddp = 0.
    q = 1. + b1 * x 
    dq = b1
    df = dp/q - p * dq / q**2
    ddf = ddp/q - 2. * dp * dq / q**2 + p * dq**2 / q**3 
    return df, ddf 
    


@numba.autojit
def qpot(x,p,r,w):

    """
    Linear Quantum Force : direct polynomial fitting of derivative-log density (amplitude)    
    curve_fit : randomly choose M points and do a nonlinear least-square fitting to a 
            predefined functional form      
    """
    
    M = Ntraj         # number of points for curve fitting
    
    # randomly choose M points 
    #xdata = np.zeros(M)
    #pdata = np.zeros(M)
    #rdata = np.zeros(M)   # fitting data 
    
    #for i in range(M):        
        #j = np.random.rand() * Ntraj
        #xdata[i], pdata[i], rdata[i] = x[j], p[j], r[j]

    # choose the first M trajectories 
    xdata, pdata, rdata = x[0:M], p[0:M], r[0:M]
    
    #tau = (max(xdata) - min(xdata))/(max(x) - min(x))
    #if tau > 0.6:
    #    pass 
    #else: 
    #    print('Data points are not sampled well.'
    
    method = 'pade' 
        
    if method == 'poly':    
        # curve_fit of p 
        #p0 = [0.,0.,-1.0,0.0]
        popt, pcov = curve_fit(func, xdata, pdata, method='lm')
        # p_approx = func(x,*popt)
        
        a,b,c,d = popt 
        dp, ddp = df(x,a,b,c,d)
          
        # curve_fit of r with poly 
        popt, pcov = curve_fit(func, xdata, rdata)
    
        a,b,c,d = popt 
        r_approx = func(x,a,b,c,d)
        dr, ddr = df(x,a,b,c,d)
        
        MSE = sum((r-r_approx)**2)/M
        
        if MSE > 10: 
            print('Fitting fails. Mean Square Error = {} \n'.format(MSE))
            sys.exit() 
        else:
            f_MSE.write('{} {} \n'.format(t, MSE))
            
    elif method == 'pade':
        
        popt, pcov = curve_fit(pade, xdata, pdata, method='lm')
        # p_approx = func(x,*popt)
   
        
        a,b,b1 = popt 
        dp, ddp = df_pade(x,a,b,b1)
          
        # curve_fit of r with poly 
        popt, pcov = curve_fit(pade, xdata, rdata, method='lm')
       
        
        a,b,b1 = popt 
        r_approx = pade(x,a,b,b1)
        dr, ddr = df_pade(x,a,b,b1)
        
        MSE = sum((r-r_approx)**2)/M
        if MSE > 10: 
            print('Fitting fails. Mean Square Error = {} \n'.format(MSE))
            sys.exit() 
        else:
            f_MSE.write('{} {} \n'.format(t, MSE))
    
    fr = -1./2./am * (2. * dp * r_approx + ddp)
    fq = 1./2./am * (2. * r_approx * dr + ddr)  
    
    Eu = -1./2./am * np.dot(r**2 + dr,w)
        
    return Eu,fq,fr 


    
def sym(V):

    n = V.shape[-1] 
    
    for i in range(n):
        for j in range(i):
            V[j,i] = V[i,j] 
    return V 






# initialization    
Ntraj = 1000 
a = 1.0
x0 = -1.0  
x = np.random.randn(Ntraj)
x = x / np.sqrt(2.0 * a) + x0 

p = np.zeros(Ntraj)
r = -a * (x-x0)

w = np.array([1./Ntraj]*Ntraj)
am = 1.0 

Nt = 1500
dt = 0.004 
dt2 = dt/2.0 
t = 0.0 


f = open('traj.dat','w')
fe = open('en.out','w')
f_MSE = open('rMSE.out','w')
nout = 20       # number of trajectories to print 
fmt =  ' {}' * (nout+1)  + '\n'  
Eu = 0.  

Ndim = 1        # dimensionality of the system  
fric_cons = 1.6       # friction constant  

v0, dv = derivs(x)
Eu,fq,fr = qpot(x,p,r,w)

for k in range(Nt):
    t = t + dt 

    p += (- dv + fq - fric_cons * p)*dt2   
    r += fr * dt2
    
    x +=  p*dt/am

    # force field 
    Eu, fq, fr = qpot(x,p,r,w)
    v0, dv = derivs(x)

    p += (- dv + fq - fric_cons * p)*dt2 
    r += fr * dt2 
       
    f.write(fmt.format(t,*x[0:nout]))
    Ek = np.dot(p*p,w)/2./am 
    Ev = np.dot(v0,w)
    fe.write('{} {} {} {} {} \n'.format(t,Ek,Ev,Eu,Ek+Ev+Eu))

fe.close()
f.close() 



    