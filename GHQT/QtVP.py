# Time-dependent quantum dynamics with TD Gauss-Hermite basis 
# 1D anharmonic system 
# V',V'' follows DFVP 

import numpy as np
import math
import sys 

Ntraj = 16000 
print('Number of trajectories = {} \n'.format(Ntraj))

dt = 0.004
am = 1.0 
x0 = -0.0 
alpha = 1.0+0j 
a = alpha.real 
b = alpha.imag

g = 0.0 # anharmonic constant in potential 

Nt = input('Time interval = {} \n How many time steps? '.format(dt))
x = np.random.randn(Ntraj)/np.sqrt(2.0*alpha.real) + x0 
p = np.zeros(Ntraj)
w = np.array([1./Ntraj]*Ntraj)

def derivs(x,xAve):
    """
    Local Harmonic Approximation V 
    """
    
    V1, V2 = Hessian_VP(xAve)

    v = V0 + V1*(x-xAve) + V2 * (x-xAve)**2/2.
    dv = V1 + V2 * (x-xAve) 

    return v,dv

def Hessian_VP(a,x,xAve,c,w):
    """
    1. Classical trajectory, V',V'' at position xAve 
    2. Variational principle 

    """
    #V0 =  x**2/2.0 + g*x**4/4. 
    #V1 =  x + g * x**3
    #V2 =  1.0 + 3.0 * g * x**2

    M1 = M1mat(a) 
    M2 = M2mat(a) 
    M3 = M3mat(a) 
    M4 = M4mat(a) 

    Vm0 = Vmat(a,x,xAve,w) 
    Vm1 = V1mat(a,x,xAve,w) 
    Vm2 = V2mat(a,x,xAve,w)
    
    #print(np.vdot(c,Vm0.dot(c)))
    #print(np.vdot(c,Vm1.dot(c)))
    #print(np.vdot(c,Vm2.dot(c)))


    W1 = np.vdot(c,(np.dot(M1,Vm0) - Vm1).dot(c))
    W2 = np.vdot(c,(np.dot(M2,Vm0) - Vm2).dot(c))

    S11 = np.vdot(c,(np.dot(M1,M1) - M2).dot(c))
    S12 = np.vdot(c,(np.dot(M1,M2) - M3).dot(c))
    S21 = np.vdot(c,(np.dot(M2,M1) - M3).dot(c))
    S22 = np.vdot(c,(np.dot(M2,M2) - M4).dot(c))

    #print('\nmatrix\n')
    #print(M1,M2)

    #print(S11,S12,S21,S22)
    print('W1,W2 = ',W1,W2)
    print('Det S = ',S11*S22-S12*S21)

    V1 = ((W1 * S22 - W2 * S12)/(S11*S22 - S12*S21)).real 
    V2 = 2.0 * ((W2 * S11 - W1 * S21)/(S11*S22 - S12*S21)).real  
    
    return V1,V2

def LQF(x,w,xAve,var):
    
    a = 1. / 2. / var 
    r = - a * (x-xAve) 
    dr = - a 

    uAve =  (np.dot(r**2,w))/2./am 
    
    du = -1./am * (r*dr)

    return uAve,du 

def qpot(x,w,xAve,xVar):

    alpha =  1.0 / xVar / 2.0 
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

def expand(alpha,x,xAve,pAve,w,c): 

    z = x - xAve 
    # V = < m | DeltaV | n > 
    
    # hermite polynomails 
    #H0 = 1. 
    #H1 = z * np.sqrt(2.) 
    #H2 = (4. * z**2 - 2.) / 4. / np.sqrt(2.) 
    a = alpha.real 

    V = Vmat(a,x,xAve,w)
    K = Kmat(alpha,x,pAve)
    D = Dmat(alpha,x,xAve,pAve,c,w)

    f6.write(' {} {} '.format(t,np.vdot(c,(K+V).dot(c)).real))

    #V[0,2] = np.dot(DeltaV * H0 * H2, w)
    #V[1,2] = np.dot(DeltaV * H1 * H2, w) 

    #V = symmetrize(V)

    #V[1,1] = np.dot(DeltaV * H1 * H1, w) 
    #V[2,2] = np.dot(DeltaV * H2 * H2, w) 

    #print 'iH + D', 1j*(K+V) + D 

    dc = -(1j*(K+V) + D).dot(c)

    return dc 

def Vmat(a,x,xAve,w):

    z = (x - xAve)*np.sqrt(a)
    #H0 = 1.0 
    #H1 = z * np.sqrt(2.)

    V = Potential(x)
    #V0 = Potential(xAve) 
    #V = V - V0 

    #v = V0 * np.identidy(Nb) 
    
    Vm = np.zeros((Nb,Nb))
    
    H = Hermite(z) # Hermite polynomails 

    for i in range(Nb):
        for j in range(i+1):
            Vm[j,i] =  np.dot(V * H[i] * H[j], w)

    Vm = Sym(Vm) 

    #Vm[0,0] = np.dot(V,w) 
    #Vm[0,1] = np.dot(V * H[] * H1, w) 
    #Vm[1,1] = np.dot(V * H1 * H1, w)    
    #Vm[1,0] = Vm[0,1] 

    return Vm 

def V1mat(a,x,xAve,w):

    z = (x - xAve)*np.sqrt(a)

    V = Potential(x)
    Vm = np.zeros((Nb,Nb))
    
    H = Hermite(z) # Hermite polynomails 

    for i in range(Nb):
        for j in range(i+1):
            Vm[j,i] =  np.dot(V * (x-xAve) * H[i] * H[j], w)

    Vm = Sym(Vm) 
    
    return Vm 

def V2mat(a,x,xAve,w):

    z = (x - xAve)*np.sqrt(a)
    V = Potential(x)
    Vm = np.zeros((Nb,Nb))
    
    H = Hermite(z) # Hermite polynomails 

    for i in range(Nb):
        for j in range(i+1):
            Vm[j,i] =  np.dot(V * (x-xAve)**2 *  H[i] * H[j], w)

    Vm = Sym(Vm) 

    return Vm 

def M1mat(a):
    
    M1 = np.zeros((Nb,Nb)) 

    for m in range(Nb-1):
        M1[m,m+1] = np.sqrt(float(m+1)/2.0/a)

    M1 = Sym(M1) 

    return M1 

def M2mat(a):
    
    M2 = np.zeros((Nb,Nb)) 

    for m in range(Nb):
        M2[m,m] = (float(m) + 0.5)/a 
    
    if Nb > 1: 
        for m in range(Nb-2):
            M2[m,m+2] = np.sqrt(float((m+1)*(m+2)))/2.0/a 

    M2 = Sym(M2)

    return M2 

def M3mat(a):
    
    M3 = np.zeros((Nb,Nb)) 

    for m in range(Nb-1):
        M3[m,m+1] = 3.0 * (float(m+1)/2./a)**1.5 

    if Nb > 2:
        for m in range(Nb-3):
            M3[m,m+3] = np.sqrt(float((m+1)*(m+2)*(m+3))) / (2.0*a)**1.5 
    
    M3 = Sym(M3) 

    return M3 

def M4mat(a):
    
    M4 = np.zeros((Nb,Nb))

    for m in range(Nb):
        M4[m,m] =  float(3.0 * m**2 + 3.0 * (m+1)**2) / (2.*a)**2
    
    if Nb > 1: 
        for m in range(Nb-2):
            M4[m,m+2] = (4.0*m + 6.0) * np.sqrt(float((m+1)*(m+2))) / (2.*a)**2
            
    if Nb > 3: 
        for m in range(Nb-4):
            M4[m,m+4] = np.sqrt(float((m+1)*(m+2)*(m+3)*(m+4))) / (2.0*a)**2

    M4 = Sym(M4) 

    if Nb > 1:
        if not M4[0,1] == M4[1,0]: 
            print(M4) 
            print('\n ERROR: Not symmetric matrix M4.\n')
            sys.exit() 
    return M4


def Hermite(x):

    cons = np.array([1. / np.sqrt(float(2**n) * float(math.factorial(n))) for n in range(Nb)])
    
    H = [] 
    H.append(1.0) 
    H.append( x * 2.0 ) 
    if Nb > 2:
        for n in range(2,Nb):
            Hn = 2.0 * x * H[n-1] - 2.0*(n-1) * H[n-2]
            H.append(Hn)
    
    for n in xrange(Nb):
        H[n] = H[n]*cons[n] 

    return H

#    if n == 0: 
#        H.append(1.)  
#    elif n == 1: 
#        return 2. * x * cons 
#    elif n == 2: 
#        return (4. * x**2 - 2.) * cons   
#    elif n == 3: 
#        return (8.0 * x**3 - 12.0 * x) * cons 
#    elif n == 4:
#        return (16.0 * x**4 - 48.0 * x**2 + 12.0) * cons 
#    elif n == 5:
#        return (32.0*x**5 - 160.0*x**3 + 120.0*x) * cons 
#    elif n == 6: 
#        return ()

def Potential(x):
    
    return  x**2/2.0 + g * x**4 / 4.0

def Kmat(alpha,x,pAve):

    K = np.zeros((Nb,Nb),dtype=complex)

    ar = alpha.real 

    for j in range(Nb): 
        K[j,j] = np.abs(alpha)**2 / ar * (2. * j + 1.)/2. +  pAve**2 
    
    for j in range(1,Nb):
        K[j-1,j] = -1j*np.conj(alpha) * pAve * np.sqrt(2. * j / ar)
        K[j,j-1] = np.conj(K[j-1,j])

    if Nb > 2: 
        for j in range(2,Nb):
            K[j-2,j] = - np.sqrt(float((j-1)*j)) * np.conj(alpha)**2 / 2. / ar  
            K[j,j-2] = np.conj(K[j-2,j])
    

    #K[0,0] = np.abs(alpha)**2/alpha.real / 2. + pAve**2
    #K[1,1] = np.abs(alpha)**2/alpha.real * 3.0 / 2. + pAve**2 

    #K[0,1] = -1j*np.conj(alpha) * pAve * np.sqrt(2.*j/alpha.real)
    #K[1,0] = np.conj(K[0,1])
    K = K / (2.*am) 

    return K 

def Dmat(alpha,x,xAve,pAve,c,w):

    D = np.zeros((Nb,Nb),dtype=complex)
    
    a = alpha.real 

    V1, V2  = Hessian_VP(a,x,xAve,c,w) 

    # time derivative 
    dq = pAve / am
    dp = - V1
    da = (-1.0j/am * alpha**2 + 1.j * V2) 
    
    for k in range(Nb):
        D[k,k] = - 1j*da.imag/2./a * (float(k) + 0.5) - 1j * pAve**2/am    
    
    for k in range(1,Nb):
        D[k-1,k] = np.sqrt(float(k)/2./a) * ( - np.conj(alpha) * dq + 1j * dp) 
        D[k,k-1] = np.sqrt(float(k)/2./a) * ( alpha * dq + 1j * dp )
    
    if Nb > 2: 
        for k in range(2,Nb):
            D[k-2,k] = np.conj(da)/2./a * np.sqrt(float(k*(k-1)))/2.0 
            D[k,k-2] = - da/2./a * np.sqrt(float(k*(k-1)))/2.0
            


    #D[0,0] = - 1j * pAve**2 / am 
    #D[1,1] = - 1j * pAve**2 / am 
    #D[0,1] = - (np.conj(alpha)*pAve/am + 1j * V1) / np.sqrt(2.*alpha.real)
    #D[1,0] = (alpha * pAve / am - 1j * V1) / np.sqrt(2.*alpha.real)

    return D 


def Sym(V):
    n = V.shape[-1]
    
    for i in range(n):
        for j in range(i):
            V[i,j] = V[j,i] 
    return V 
    


def SaveWf(alpha, pAve, S,c,xAve,xVar,fname='wft.dat'):
    """
    save wavefunction to file 
    """
    f = open(fname,'w')
    
    x = np.linspace(-6,6,200) 

    a = alpha.real 
    z = (x - xAve) * np.sqrt(a) 

    #print('GWP width parameter {}, Real alpha {}'.format(a,alpha.real))

    phi0 = (a/np.pi)**0.25 * np.exp( - alpha * (x-xAve)**2/2.0 \
            + 1j*pAve*(x-xAve))
    
    H = Hermite(z) 

    basis = [] 
    for i in range(Nb):
        basis.append(H[i]*phi0)
    #phi1 = Hermite(z,1) * phi0 
    #phi2 = Hermite(z,2) * phi0 
    #phi2 = (4. * z*z - 2.) / 4. / np.sqrt(2.) * phi0 

    for i in xrange(len(x)):
        wf = 0.+0.j 
        for j in xrange(Nb):
            wf += c[j]*basis[j][i] 
        f.write('{} {} {} \n'.format(x[i], wf.real,wf.imag))

    f.close() 

    
f1 = open('traj.dat', 'w')
f2 = open('xAve.dat', 'w')
f3 = open('energy.dat', 'w')
f4 = open('coeff.dat', 'w')
f5 = open('norm.dat', 'w')
f6 = open('energy.dat', 'w')

#def corr(c0,c) 

t = 0.0
dt2 = dt/2.0 
S = 0. # phase in GWP 


xAve = np.dot(x,w)
xSqdAve = np.dot(x*x,w) 
xVar = (xSqdAve - xAve**2) 


pAve = np.dot(p,w) 
print('\n Initial Time \n ')
print('Mean position = {}, Variance = {} '.format(xAve,xVar))
print('Initial momentum {}'.format(pAve))

#ak = Hessian(xAve)
#Nb = 10 # number of basis function 
Nb = input('Please enter number of basis function \n ')

c = np.zeros(Nb,dtype=complex) # initial expansion coeffs
c[0] = 1.0

SaveWf(alpha, pAve, S, c, xAve,xVar,fname='wf0.dat') # save initial wavefunction 

# LHA
V1,V2 = Hessian_VP(a,x,xAve,c,w) 
print('Before time V1 = {}, V2 = {}'.format(V1,V2))
dv =  V2* (x-xAve) + V1 
uAve,du = LQF(x,w,xAve,xVar)

V = Potential(x) 
vAve = np.dot(V,w) 
print(' Quantum potential = {} \n '.format(uAve))
print(' Potential = {} \n '.format(vAve))
print(' Total energy = {} \n '.format(uAve+vAve))


# format for output data 
fmt = ' {} '*11 + '\n'
fmtC = ' {} '*(Nb+1) + '\n'

# update c for one timestep 
cold = c 
dc = expand(alpha,x,xAve,pAve,w,c)
c = c + dc*dt 

for k in range(Nt):
    
    t += dt 

    # leap-frog alg for {x,p}
    p = p - dv*dt2 - du*dt2 
    
    x += p * dt / am 
    
    # compute observables 
    xAve = np.dot(x,w)
    xSqdAve = np.dot(x*x,w) 
    xVar = (xSqdAve - xAve**2) 
    
    # force fields
    V1,V2 = Hessian_VP(a,x,xAve,c,w)
    print('c = ',c)
    print('V1 = {}, V2 = {}'.format(V1,V2))
    dv =  V2 * (x-xAve) + V1 
    uAve, du = LQF(x,w,xAve,xVar)

    p = p - dv*dt2 - du*dt2 

    # average momentum, update phase parameters 
    pAve = np.dot(p,w) 

    #alpha += (-1.0j/am * alpha**2 + 1.j * ak) * dt 
    a = 1.0 / 2.0 / xVar 
    b += (- (a**2 - b**2) / am + V2)*dt 
    alpha = a + 1j*b 

    # S is the real part of the complex phase term, imaginary part is absorbed 
    # into the normalization constant N 
    #S += ( pAve**2/2./am - V0 - a/2./am ) * dt 
    #S += (pAve**2/2./am - V0 ) * dt 

    # update c, second-order difference 
    dc = expand(alpha,x,xAve,pAve,w,c)
    cnew = cold + 2.0*dc*dt 
    cold = c 
    c = cnew 

    f5.write( '{} {} \n'.format(t,np.vdot(cold,cold)))
    

    #vAve = np.dot(v,w)
    kAve = np.dot(p*p/2./am, w) 
    
    # overlap with a GWP ~ (xAve, xVar)
    #uAve = overlap(x,w,xAve,xVar) 
    #print('overlap with GWP {} \n'.format(ov))
    
    # save data 
    SaveWf(alpha, pAve,S, cold, xAve,xVar)

    f3.write('{} {} {} \n'.format(t,kAve,uAve))
    f1.write(fmt.format(t,*x[0:10]))
    f2.write(' {} {} {} {} \n'.format(t,xAve,pAve,xVar))
    f4.write(fmtC.format(t,*c[0:Nb]))

f1.close()
f2.close() 


