

import numpy as np 

Ntraj = 2000 
Nt= 100 
dt = 0.01 
am = 1.0 
x0 = -1.0 

x = np.random.randn(Ntraj)+x0 
p = np.zeros(Ntraj)
w = np.array([1./Ntraj]*Ntraj)

def derivs(x):
    c = 0.0
    v = x**2/2. + c * x**4/4.0 
    dv = x + c * x**3  
    return v,dv 

def LQF(x,w,xAve,var):
    
    r = -(x-xAve)/var/2.0 

    dr = -1.0/2.0/var
    uAve =  (np.dot(r**2,w))/2./am 
    
    du = -1./am * (r*dr)

    return uAve,du 

def qpot(x,w,xAve,xVar):

    alpha =  1.0 / xVar  
    z = x - xAve 

    g0 = (alpha/2./np.pi)**0.5 * np.exp(-alpha*z**2/2.)
    g1 = z*g0

    s00 = np.sqrt(alpha/np.pi)/2.0
    s11 = 1./4./np.sqrt(alpha*np.pi) 

    b0 = np.dot(g0,w) 
    b1 = np.dot(g1,w) 

    # c0 = b0/s00 
    c0 = 1.0 
    c1 = b1/s11
    
    print 'coeff',c0,c1

    P = c0*(-alpha*(x-xAve)) + c1*(-alpha + alpha**2*(x-xAve)**2)
    dP = - alpha*c0 + c1 * alpha**2 * 2.0*(x-xAve) 
    ddP = 2.0 * c1 * alpha**2

    Q = c0 + c1*(x-xAve)     
    dQ = c1 
    # ddQ = 0.0, ignore P*ddQ/Q
    r = 0.5 * P/Q
    dr = 0.5 *(-P*dQ/Q**2 + dP/Q)
    ddr = 0.5 * (-(dP*dQ/Q**2 - P*dQ*dQ/Q**3) + (ddP/Q - dP*dQ/Q**2))

    du = -(2.0*r*dr + ddr)/2.0/am 

    uAve = np.dot(r**2,w)/2./am 

    return uAve,du  
    

f1 = open('traj.dat', 'w')
f2 = open('xAve.dat', 'w')
f3 = open('energy.dat', 'w')
f4 = open('coeff.dat', 'w')

t = 0.0
dt2 = dt/2.0 

fmt = ' {} '*11 + '\n'

xAve = np.dot(x,w)
xSqdAve = np.dot(x*x,w) 
xVar = (xSqdAve - xAve**2) 

v,dv = derivs(x)
uAve,du = qpot(x,w,xAve,xVar)

for k in range(Nt):
    t += dt 

    

    p = p - dv*dt2 - du*dt2 
    
    x += p * dt / am 
    
    # compute observables 
    xAve = np.dot(x,w)
    xSqdAve = np.dot(x*x,w) 
    xVar = (xSqdAve - xAve**2) 
    
    # force fields
    v,dv = derivs(x)
    uAve,du = qpot(x,w,xAve,xVar)

    p = p - dv*dt2 - du*dt2 

    vAve = np.dot(v,w)
    kAve = np.dot(p*p/2./am, w) 
    
    # overlap with a GWP ~ (xAve, xVar)
    #uAve = overlap(x,w,xAve,xVar) 
    #print('overlap with GWP {} \n'.format(ov))

    f3.write('{} {} {} {} \n'.format(t,kAve,vAve,uAve))
    f1.write(fmt.format(t,*x[0:10]))
    f2.write(' {} {} {} \n'.format(t,xAve,xVar))
    #f4.write('{} {} {} \n'.format(t,c0,c1))

f1.close()
f2.close() 


