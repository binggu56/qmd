# -*- coding: utf-8 -*-

# Frozen Gaussian and quantum trajectories 
# The GWPs is used to approximate the amplitude instead of the wavefunction 
# Quantum force is computed from the approximated amplitude 

# A problem is that the quantum force obtained from the approximate amplitude is not accurate at the tail 



#import urllib2 as ur
#import re, time
#import os
#import pandas
import matplotlib.pyplot as plt 
#from matplotlib import animation 
import numpy as np 
import numba 
#from vmat import evolve  

#import matplotlib.dates as mdates
#import datetime as dt
#from matplotlib.figure import Figure

#import matplotlib as mpl


import sys 


    
def integrate(x,px,aj,y,py,ak):
    dp = py - px 
    dq = y - x
    
    return (aj*ak)**0.25 * np.sqrt(2./(aj+ak)) * np.exp(    \
            -0.5 * aj*ak/(aj+ak) * (dp**2/aj/ak + dq**2  \
            + 2.0*1j* (self.p/aj + py/ak) *dq)
            )
    
def integrate_x(x,px,aj,y,py,beta):
    
    d = integrate(x,px,aj,y,py,beta)

    return d*(-1j*(px-py) + aj*x+y*beta)/(aj+beta)
    
def integrate_x2(x,px,aj,y,py,beta):
    
    z = -1j*(px-py) + aj*x + y*beta
    
    return z*z*integrate(x,px,aj,y,py,beta)
    
#    def value(z):
#        return np.sqrt(np.sqrt(self.alpha/np.pi))*np.exp(-self.alpha*(z-self.x)**2/2.0)
    
#def set_x(self,z):
#        
#    self.x = z 
#    
#def plot(self,ax,xmin=-4,xmax=4,Np=200):
#    
#    x = np.linspace(-4,4,Np)
#    
#    ax.plot(x,np.abs(self.value(x))**2)
    



#w = np.array([1.0/float(Ntraj) for i in range(Ntraj)])

        
cut = 1e-6
#xmin = x0-np.sqrt(-np.log(cut/np.sqrt(self.alpha/np.pi))/self.alpha)      
#self.xmax = x0+np.sqrt(-np.log(cut/np.sqrt(self.alpha/np.pi))/self.alpha)

#pow = 4 
#self.xmin = self.x - np.sqrt(2.0*pow/self.alpha)
#self.xmax = self.x + np.sqrt(2.0*pow/self.alpha)

#dx = (self.xmax-self.xmin)/(Ntraj-1)

#self.grid = np.linspace(self.xmin,self.xmax,Ntraj)


    
            
def overlap(x):
    
    S = np.zeros((Ntraj,Ntraj),dtype=np.complex128)
    
    for j in range(Ntraj):
        
        aj,qj,pj = alpha, x[j], 0.0  
        
        for k in range(Ntraj):
            
            ak, qk, pk = alpha, x[k], 0.0  
            
            dq = qk - qj 
            dp = pk - pj                 
           
            S[j,k] = (aj*ak)**0.25 * np.sqrt(2./(aj+ak)) * np.exp(    \
            -0.5 * aj*ak/(aj+ak) * (dp**2/aj/ak + dq**2  \
            + 2.0*1j* (pj/aj + pk/ak) *dq)   )
                
    return S
    
def projection(x):
        
    S = overlap(x)
    
    b= np.array([integrate(x[i],0.0,alpha,x0,0.0,a0) for i in range(Ntraj)])
    
    try:
        c = np.linalg.solve(S, b)
    except:
        print("Error: ill-conditioned overlap matrix of initial basis.")
        sys.exit()
    
    return c 


    
#    def get_p(self, beta, smooth=True):
#        """
#        compute p = grad S at grid points  
#        """
#        g = self.basis 
#        
#        dq = 0.005
#        
#        p = np.zeros(Ntraj) 
#        dp = np.zeros(Ntraj) 
#        
#        #print 'x = ',self.grid 
#        #print 'c = ',self.c 
#        
#        if smooth:
#                       
#            for j in range(Ntraj):
#                
#                dz, z = 0.0, 0.0 
#                dz1, z1 = 0.0, 0.0 
#                dz0, z0 = 0.0, 0.0 
#                 
#                
#                q = g[j].x 
#                
#                q1 = q + dq 
#                q0 = q - dq 
#                
#                for k in range(Ntraj):
#                    
#                    
#                    qk, ak = g[k].x, g[k].alpha 
#                    
#                    y = q - qk 
#                    
#                    alfa = (ak*beta)/(ak+beta) 
#                    
#                    an = (ak/np.pi)**0.25 * np.sqrt(beta/(ak+beta))
#                        
#                    z = z + self.c[k] * np.exp(-0.5*alfa*y**2) *an
#                             
#                    dz = dz - self.c[k] * alfa * y * np.exp(-0.5*alfa*y**2) * an
#                    
#                    #ddz = ddz + self.c[k] * (- alfa + alfa**2 * y**2) * np.exp(-0.5*alfa*y**2) * an  
#                           
#                        
#                p[j] = (dz/z).imag
#                #dp.append((-dz/z**2+ddz/z).imag)
#
#                for k in range(Ntraj):
#                    
#                    
#                    qk, ak = g[k].x, g[k].alpha 
#                    
#                    y0, y1 = q0 - qk, q1 - qk  
#                    
#                    alfa = (ak*beta)/(ak+beta) 
#                    
#                    an = (ak/np.pi)**0.25 * np.sqrt(beta/(ak+beta))
#                        
#                    z0 = z0 + self.c[k] * np.exp(-0.5*alfa*y0**2) *an
#                    z1 = z1 + self.c[k] * np.exp(-0.5*alfa*y1**2) *an
#                    
#                    dz0 = dz0 - self.c[k] * alfa * y0 * np.exp(-0.5*alfa*y0**2) * an
#                    dz1 = dz1 - self.c[k] * alfa * y1 * np.exp(-0.5*alfa*y1**2) * an
#                    
#                    
#                dp[j] = (((dz1/z1).imag-(dz0/z0).imag)/2./dq)
#
#        #print 'p = ',p 
#        
#        return p, dp 
#    
#    def norm(self):
#        
#        c = self.c 
#        S = self.overlap()
#        
#        z = np.vdot(c,S.dot(c))
#        #for j in range(Ntraj):
#        #    for k in range(Ntraj):
#        #        z += np.conj(c[j]) * S[j,k] * c[k]
#        return z.real
#    
#    def overlap(self):
#        
#        g = self.basis 
#        
#        S = np.zeros((Ntraj,Ntraj),dtype=np.complex128)
#        
#        for j in range(Ntraj):
#            
#            aj,qj,pj = g[j].alpha,g[j].x,g[j].p 
#            
#            for k in range(Ntraj):
#                
#                ak, qk, pk = g[k].alpha, g[k].x, g[k].p 
#                
#                dq = qk - qj 
#                dp = pk - pj                 
#               
#                S[j,k] = (aj*ak)**0.25 * np.sqrt(2./(aj+ak)) * np.exp(    \
#                -0.5 * aj*ak/(aj+ak) * (dp**2/aj/ak + dq**2  \
#                + 2.0*1j* (pj/aj + pk/ak) *dq)   )
#                    
#        return S
 
    


class Tdse(Psi):

    def __init__(self,c,g, **kwargs):
    
        Psi.__init__(self,c,g)
        
        try:
            self.dt = kwargs['dt']
            self.Nt = kwargs['Nt']
            self.am = kwargs['am']
            self.modelName = kwargs['modelName']
            self.beta = kwargs['beta']
            
        except:
            print("ERROR: missing args for Tdse.")
        


    def solve(self):
        
        gx = np.array([basis.x for basis in self.basis]) 
        gp = [x.p for x in self.basis]
        alpha = [x.alpha for x in self.basis]
        
        am = self.am 
        
        S = self.overlap() 
        
        p, dp = self.get_p(beta = self.beta)
        
        D = dmat(am,gx,gp,p,dp,alpha,S)
         
        V = vmat(gx,gp,alpha,S) 
        K = kmat(am,gx,gp,alpha,S)
        
        c = self.c 
        
        #enk = np.vdot(c,K.dot(c))
        #env = np.vdot(c,V.dot(c))               
        
        H = (K+V)-1j*D
    
        b = np.dot(H,c)
        
        b = -1j*b
        
        #print 'ham = ', H
        
        S = self.overlap()
        
        #print 'overlap',S
     
        try:
            dc = np.linalg.solve(S, b)
            
        except:
            print("Error: ill-conditioned overlap matrix for dc")
            print(S) 
            sys.exit()
        
        return dc 
    

    
    def propagate_x(self):
        
        """
        update quantum trajectories using momentum obtained from current wavefunction
        also construct new GWP basis 
        """
        dt = self.dt 
        am = self.am 
        g = self.basis 
        
        p, dp = self.get_p(beta = self.beta)
        
        
        self.grid = self.grid + p*dt/am 
        
        #g = [gwp(self.grid[i],0.0) for i in range(Ntraj)]
        #print 'p',p
        #print 'dp',dp 
        #print '\n'
        
        for k in range(Ntraj):
            g[k].set_x(self.grid[k])
            # g[k].alpha = g[k].alpha - dp[k]/am * g[k].alpha * dt 
        
       
    
        
    def evolve(self, integrator='SOD'):
        
        """
        SOD   : second-order difference 
        Euler : first-order integration 
        """
        
        dt = self.dt 
        Nt = self.Nt
        am = self.am 
        g = self.basis 
        
        
        f = open('x.dat','w')
        f1 = open('cor.dat','w')
        f2 = open('energy.dat','w')
        
        t = 0.
        
#        if integrator == 'Euler':
#            for i in range(Nt):
#            
#                t = t + dt
#            
#                dc = self.solve()
#                self.c = self.c + dc*dt 
#                self.propagate_x()
#                
#                
#        
#        
#                
#                f.write(str(t) + ' ' + ' '.join([str(self.grid[i]) for i in range(Ntraj)]) +'\n')
#                
#        elif integrator == 'SOD':
            
        cold = self.c 
        dc= self.solve()
                    
        
        self.c = self.c + dc*dt 
        
        for i in xrange(Nt):
        
            t = t + dt
            
            #self.propagate_x()
            #p, dp = self.get_p(beta = self.beta)
            
            for k in range(Ntraj):
        
                g[k].x += p[k]*dt/am 
                #g[k].alpha = g[k].alpha - dp[k]/am * g[k].alpha * dt 

            
            dc = self.solve()
            
            
            cnew = cold + dc*2.0*dt
            cold = self.c 
            self.c = cnew
            
            #self.c = self.c/np.sqrt(self.norm())
            
            if i% int(Nt/3) == 0:
            #if False:
            
                
                print '\n ------- step {:2d} --------\n'.format(i)
                
                print 'norm {:6.2f} \n'.format(self.norm()) 
                #print 'energy Ek = {:6.2f}, Ev = {:6.2f} '.format(enk.real,env.real) 
                
            #f.write(str(t) + ' ' + ' '.join([str(self.grid[i]) for i in range(Ntraj)]) +'\n')
            
            z = self.corr()
            
            f1.write('{:12.6f} {:12.6f} {:12.6f} \n'.format(t, z.real, z.imag))
            #f2.write('{:12.6f} {:12.6f} {:12.6f}'.format(enk,env,enk+env))
                
                
                    
        f.close()
        f1.close() 
        f2.close() 
        

    
    def corr(self):
        
        c = self.c
        S = self.overlap()
        
        z = c.dot(S.dot(c))
        
        return z


def quantum_force(x,c):
    """
    compute quantum potential from the approximate amplitude     
    """
    dz,z = 0.0, 0.0 
    u = np.zeros(len(x))
    du = np.zeros(len(x))
    
    for j in range(len(x)):
    
        dz, z = 0.0, 0.0 
             
        q = x[j] 
    
        for k in range(len(x)):
        
        
            qk, ak = x[k], alpha 
            
            y = q - qk 
            
            an = (ak/np.pi)**0.25
                
            z = z + c[k] * np.exp(-0.5*ak*y**2) * an
                     
            dz = dz + c[k] * (ak**2 * y**2 - ak) * np.exp(-0.5*ak*y**2) * an
            
            dddz += c[k] * (-ak**2 * y * (ak * y**2 - 3.0)) * np.exp(-0.5*ak*y**2) * an
            
            #ddz = ddz + self.c[k] * (- alfa + alfa**2 * y**2) * np.exp(-0.5*alfa*y**2) * an  
                   
                
        u[j] = dz/z * (-0.5)
        du[j] = (-0.5) * (dddz/z - dz/z**2)
    
    return du 


# global params

x0 = 0.0 
a0 = 1.0 
p0 = 0.0 
Ntraj = 16 

Ntraj = 20
alpha = 8.0  # basis width 

x = np.random.randn(Ntraj) / np.sqrt(2.0 * a0) + x0 

print('initial wavefunction x0 = {:6.2f}, p0={:6.2f}, a0 = {:6.2f} \n'.format(x0,p0,a0))
        
print('particle mass'.format(am))  
print 'number of basis {:4d} \n'.format(Ntraj)
print('initial grid \n') 
print(x)

print 'time interval = {:6.2f} time steps = {:4d} \n'.format(dt,Nt)

c =  projection(x)


print('c(0) =  ',c) 
print('trajectory interval '.format(x[1]-x[0]))
print('variation of basis'.format(np.sqrt(1./alpha))

if start.grid[1]-start.grid[0] < np.sqrt(1./g[0].alpha):
    print("WARNING: grid spacing is too small.")
    

# expanded initial wavefunction 

print('\n initialization succeed ... \n')

# propagate quantum trajectories 
for kt in range(Nt):
    
    x += p * dt/m 

    # obtain quantum force from approximated wavefunction 
    du = quantum_force(x,c)
    dv = dv(x)
    p += (-dv - du) * dt  
    
    # update the amplitude 
    dp = linear_regression(p)
    dc = update_amplitude(x,c,dp)
    c += dc * dt 
    
  

# output final data       
Np = 200 
x = np.linspace(-5,5,Np)
     


with open('pot.dat','w') as f:
    for i in xrange(Np):
        f.write('{:12.6f} {:12.6f} \n'.format(x[i], dv(x[i])[0]))
        
with open('wft.dat','w') as f:
    for i in xrange(Np):
        f.write('{:12.6f} {:12.6f} {:12.6f} \n'.format(x[i],np.abs(wf.value(x[i]))**2, \
                np.abs(se.value(x[i]))**2 ))

# final wavefunction 
ax.plot(x,np.abs(se.value(x))**2,'k--')

#plt.legend(loc=0)
plt.show()  

@numba.jit 
def dv(x):
        
    pot = 'Double_well'
    
    if pot == 'Morse':
    
        a, x0, De = 1.02, 1.4, 0.176
    
        d = (1.0-np.exp(-a*(x-x0)))
    
        v0 = De*d**2
    
        ddv = 2.0 * De * (-d*np.exp(-a*((x-x0)))*a**2 + (np.exp(-a*(x-x0)))**2*a**2)
        
    elif pot == 'Double_well':
        
        eta = 1.3544 
        v0 = x**4/16.0/eta - x**2/2.0 
        ddv = 3./4./eta * x**2 - 1.0 
    
    elif pot == 'Harmonic':
        
        v0 = x**2/2.0
        ddv = 1.0 
    
    else:
        print "ERROR: there is no such potential."
            
    return v0,ddv
    
    

@numba.jit
def vmat(x,p,alpha,S):
    """
    local harmonic approximation (LHA) 
    v(x) = v0 + v'(x-x0) + v''(x-x0)**2/2
    x0 = (alpha*x1+beta*x2)/(alpha+beta)
    """
    Ntraj = len(x)
    l = np.zeros((Ntraj,Ntraj),dtype=np.complex128)
     
    for j in range(Ntraj):
        
        
        aj,qj,pj = alpha[j], x[j], p[j] 
    
        for k in range(Ntraj):
        
            ak,qk,pk = alpha[k], x[k], p[k] 
   
            x0 = (aj*qj+ak*qk)/(aj+ak)
        
            v0,ddv = dv(x0)
        
            d0 = v0 + ddv/2.0*((pj-pk)**2/(aj+ak)**2 + 1.0/(aj+ak))
        
            l[j,k] = d0*S[j,k]
    
       

    return l

@numba.jit 
def kmat(am,x,p,alpha,S):
    
    Ntraj = len(x)
    
    l = np.zeros((Ntraj,Ntraj),dtype=np.complex128)
     
    for j in range(Ntraj):
                 
        aj,qj,pj = alpha[j], x[j], p[j] 
        
        for k in range(Ntraj):
        
            ak,qk,pk = alpha[k], x[k], p[k]
            
            p0 = (aj*pk + ak*pj)/(aj+ak)
            d0 = 0.5/am * ( (p0+1j*aj*ak/(aj+ak)*(qj-qk))**2 + aj*ak/(aj+ak) )
            
            l[j,k] = d0*S[j,k]
        
    
    return l

@numba.jit           
def dmat(am,x,gp,p,dp,alpha,S):
    """
    time-derivative matrix D = <gj|i d/dt | gk> = <gj| -dp/(2m)*ak*(x-qk)^2 + ak*pk/m * (x-qk) - dp/(4m) | gk>
    the trajectory momentum p is different from the momentum for basis (gp)    
    """
    Ntraj = len(x)
      
 
    l = np.zeros((Ntraj,Ntraj),dtype=np.complex128)
    
    
    for j in range(Ntraj): 
        
        aj,qj,pj = alpha[j], x[j], gp[j] 
        
        for k in range(Ntraj):
            
            ak,qk,pk = alpha[k], x[k], gp[k]
        
            d2 = dp[k]/2./am * ak * 1./(aj+ak)**2 * ( 2j*aj*(pj-pk) * (qj-qk) - aj**2 * (qj-qk)**2 + (pj-pk)**2 - (aj+ak) )
            
            d1 = ak*p[k]/am * (-(1j*(pj-pk)-aj*(qj-qk)))/(aj+ak) 
            
            d0 =  - dp[k]/4.0/am
            
            l[j,k] = (d0+d1+d2)*S[j,k]

     
    return l 
        


d = dict({'dt': 0.01,  'Nt' : 300, 'am' : 1.0, 'modelName':'Double_well', 'beta' : 1.0})
np.set_printoptions(precision=4,threshold=5)

