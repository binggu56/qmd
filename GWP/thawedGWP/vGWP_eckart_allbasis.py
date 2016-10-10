"""
orthogonal matching-pursuit algorithm for quantum dynamics 

thawed Guassian wavepackets

width parameters treated within TDVP 

@date: June 25, 2016 
@author: Bing Gu 

""" 
import numpy as np 
import numba
import sys 
import math 

@numba.autojit 
def gint(aj,x,px,ak,y,py):
    """
    Gaussian integrals     
    """
    dp = py - px 
    dq = y - x 
    
    return (aj*ak)**0.25 * np.sqrt(2./(aj+ak)) * np.exp(    \
            -0.5 * aj*ak/(aj+ak) * (dp**2/aj/ak + dq**2  \
            + 2.0*1j* (px/aj + py/ak) *dq) )

@numba.autojit 
def iproj(a,x,p,a0,x0,p0):
    """
    the intial wavefunction is expanded into a linear combination of GWPs
    """

    N = len(x)
    b = np.zeros(N,dtype=np.complex128)
    
    # overlap matrix 
    S = overlap(a,x,p)
    
    # projection onto each basis 
    for k in range(N):
        
        ak, qk, pk = a[k], x[k], p[k] 
        
        b[k] = gint(ak, qk, pk, a0, x0, p0)

    c = np.linalg.solve(S,b)
        
    return c


def overlap_FGA(x,p): 
    """
    construct overlap matrix from GWPs defined by {x,p}
    """
    N = len(x)
    S = np.zeros((N,N),dtype=np.complex128)
    
    for j in range(N):
        
        aj, qj, pj = a, x[j], p[j] 
        
        for k in range(N):
            
            ak, qk, pk = a, x[k], p[k]  
            
            dq = qk - qj 
            dp = pk - pj                 
           
            S[j,k] = (aj*ak)**0.25 * np.sqrt(2./(aj+ak)) * np.exp(    \
            -0.5 * aj*ak/(aj+ak) * (dp**2/aj/ak + dq**2  \
            + 2.0*1j* (pj/aj + pk/ak) *dq)   )
                
    return S

def overlap(a,x,p): 
    """
    construct overlap matrix from GWPs defined by {a,x,p}
    """
    N = len(x)
    S = np.zeros((N,N),dtype=np.complex128)
    
    for j in range(N):
        
        aj, qj, pj = a[j], x[j], p[j] 
        
        for k in range(N):
            
            ak, qk, pk = a[k], x[k], p[k]  
            
            dq = qk - qj 
            dp = pk - pj                 
           
            S[j,k] = (aj*ak)**0.25 * np.sqrt(2./(aj+ak)) * np.exp(    \
            -0.5 * aj*ak/(aj+ak) * (dp**2/aj/ak + dq**2  \
            + 2.0*1j* (pj/aj + pk/ak) *dq)   )
                
    return S


def omp(aold, xold, pold, cold):
    
    """
    orthogonal matching pursuit 
    target function is a linear combination of GWPs with parameters {aold, xold, pold, cold}    
    """
    
    global a, gx, gp
    
    N = len(gx)
    
    tmp = np.zeros(N, dtype=complex)
    c0 = np.zeros(N, dtype=complex)

    
    x = [] 
    c = [] 
    p = [] 
    index = [] 
    res = 1.0 
    
    # projection of the target function onto each basis 
    for k in range(N):        
        c0[k] = iproj(gx[k],gp[k],a,aold,xold,pold,cold)
    
    tmp = c0 
            
    i = np.argmax(np.abs(tmp))
    index.append(i)
    x.append(gx[i])
    p.append(gp[i])
    c.append(c0[i])
        
    nmax = 40 # maximum number of basis defined by user 
        
    res = 1.0 - np.abs(c[0])**2
    
    b = c  # initial orthogonal expansion coefficients 
    
    while res > 1e-6 and len(x) < nmax: 
        
        tmp = c0 
        
        for i in range(len(x)):
            d = [gint(gx[k],gp[k],a, x[i], p[i], a) for k in range(N)]
            tmp = tmp - b[i] * np.array(d)  
        
        i = np.argmax(np.abs(tmp))
    
        index.append(i)
        x.append(gx[i])
        p.append(gp[i])
        c.append(c0[i])
        
        # orthogonalization of obtained basis so far 
        
        S = overlap(x,p)
    
        b = np.linalg.solve(S,c)
 
        res = 1.0 + np.vdot(b,S.dot(b)).real - np.vdot(b,c).real * 2.0  
        
        print('Residual Error = {}'.format(res))
     
    
    return x,p,b, index 

def qpot(xgrid, x,c):
    """
    compute quantum potential from the approximate amplitude     
    """
    dz,z = 0.0, 0.0 
    u = np.zeros(len(xgrid))
    du = np.zeros(len(xgrid))
    
    delta = 0.0
    
    
    for j in range(len(xgrid)):
    
        dz, z = 0.0, 0.0
        ddz = 0.0 
        dddz = 0.0 
             
        q = xgrid[j] 
    
        for k in range(len(x)):
        
        
            qk, ak = x[k], a 
            
            y = q - qk 
            
            an = (ak/np.pi)**0.25
                
            z += c[k] * np.exp(-0.5*ak*y**2) * an
            
            dz += c[k] * (-ak * y) * c[k] * np.exp(-0.5*ak*y**2) * an
                     
            ddz += c[k] * (ak**2 * y**2 - ak) * np.exp(-0.5*ak*y**2) * an
            
            #ddz = ddz + self.c[k] * (- alfa + alfa**2 * y**2) * np.exp(-0.5*alfa*y**2) * an  
            dddz += c[k] * (-ak**2 * y * (ak * y**2 - 3.0)) * np.exp(-0.5*ak*y**2) * an
                
        u[j] = dz/z * (-0.5)
        du[j] = (-0.5) * (dddz/(z+delta) - dz*ddz/(z**2+delta))
    
    return u,du 
    

def get_p(xgrid, a,x,p,c, smooth=True):
    """
    compute p = grad S at grid points  
    """
    
    Ntraj = len(x)
    
    phase_p = np.zeros(len(xgrid)) 

    delta = 1e-6 # soft parameter
    beta = 128.0 
    
    if smooth:
                   
        for j in range(len(xgrid)):
            
            dz, z = 0.0, 0.0 
                     
            q = xgrid[j] 
            
            for k in range(Ntraj):
                
                
                qk, ak = x[k], a[k] 
                
                y = q - qk 
                
                alfa = (ak*beta)/(ak+beta) 
                
                an = (ak/np.pi)**0.25 * np.sqrt(beta/(ak+beta))
                    
                z = z + c[k] * np.exp(-0.5*alfa*y**2) *an
                         
                dz = dz - c[k] * alfa * y * np.exp(-0.5*alfa*y**2) * an
                
                #ddz = ddz + self.c[k] * (- alfa + alfa**2 * y**2) * np.exp(-0.5*alfa*y**2) * an  
                       
                    
            phase_p[j] = (dz/z).imag
            #dp.append((-dz/z**2+ddz/z).imag)

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
    else: 
        
        for j in range(len(xgrid)):
            
            dz, z = 0.0+0j, 0.0+0j  
                     
            q = xgrid[j] 
            
            for k in range(Ntraj):
                
                pk, qk, ak = p[k], x[k], a 
                
                y = q - qk 
                
                g = (ak/np.pi)**0.25 * np.exp(- 0.5 * ak * y**2 + 1j * pk * y ) 
                    
                z += c[k] * g
                         
                dz += c[k] * (- ak * y + 1j*pk ) * g
                
                #ddz = ddz + self.c[k] * (- alfa + alfa**2 * y**2) * np.exp(-0.5*alfa*y**2) * an  
                       
                    
            phase_p[j] = (dz/(z+delta)).imag
  
    return phase_p 
    
def solve(a,x,p,c):
        
    
    S = overlap(a,x,p) 
    
    phase_p = get_p(x,a,x,p,c) # phase momentum dS/dx 
    
    D1 = d1mat(a,x,p,phase_p,S)
    D2 = d2mat(a,x,p,c,phase_p,S)
    
    print('D matrix \n')
    print(D1)
    print(D2)
     
    V = vmat(a,x,p,S) 
    
    print('potential matrix \n')
    print(V)
    
    V1 = v1mat(a,x,p,c,S)  
    
    print('potential matrix 1 \n')
    print(V1)
    
    K = kmat(a,x,p,S)
    
    print('kinetic energy \n')
    print(K)
    
    K1 = k1mat(a,x,p,c,S) 
    print('kinetic energy 1 \n')
    print(K1)    
    
    J = Jmat(a,x,p,c,S)
    
    print('J matrix \n')
    print(J)
    
    G = Gmat(a,x,p,c,S)
    
    print('G matrix \n')
    print(G)
    print('condition number of G',cnum(G))
    
    # computes energy expectation value 
    #enk = np.vdot(c,K.dot(c))
    #env = np.vdot(c,V.dot(c))      
    
    
    N = len(x)

    b1 = - np.dot(1j * (K+V) + D1,c)
    b2 = - np.dot(1j * (K1+V1) + D2, c)
    
    b  = np.zeros(2*N,dtype=np.complex128)

    b[0:N] = b1
    b[N:2*N] = b2 
    
    print(b)
    
    M = blkMatrix(S,np.conj(J).T,J,G)    
    #print 'overlap',S
 
    print('\n condition number of M',cnum(M))
    try:
        dz = np.linalg.solve(M, b)
        
        print('TDSE optimized EOM for z \n')
        print(dz)
        
    except:
        print("Error: Equation of motion cannot be solved, possibaly ill-conditioned overlap matrix for dc")
        print(S) 
        sys.exit()
    
    return dz



@numba.jit 
def dv(x):
        
    pot = 'Eckart'
    
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
        
    elif pot == 'Eckart':
        
        D = 16.0
        #a = 1.3624d0 
        a = 1.0
        v0 = D/np.cosh(a*x)**2
        dv = -2.0 * a * D / np.cosh(a*x)**2 * np.tanh(a*x) 
        ddv = 6.0*D*a**2*np.sinh(a*x)**2/np.cosh(a*x)**4-2.*D*a**2/np.cosh(a*x)**2

    
    else:
        print("ERROR: there is no such potential.")
            
    return v0,dv,ddv
    
    

@numba.jit
def vmat(a,x,p,S):
    """
    local harmonic approximation (LHA) 
    v(x) = v0 + v'(x-x0) + v''(x-x0)**2/2
    x0 = (alpha*x1+beta*x2)/(alpha+beta)
    """
    N = len(x)
    l = np.zeros((N,N),dtype=np.complex128)
     
    for j in range(N):
        
        aj,qj,pj = a[j], x[j], p[j] 
    
        for k in range(N):
        
            ak,qk,pk = a[k], x[k], p[k] 
   
            x0 = (aj*qj+ak*qk)/(aj+ak)
        
            v0,v1,ddv = dv(x0)
        
            tmp2 = 1.0/(aj+ak)

            d0 = v0 + ddv/2.0 * tmp2 
        
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
def kmat(a,x,p,S):
    
    Ntraj = len(x)
    
    l = np.zeros((Ntraj,Ntraj),dtype=np.complex128)
     
    for j in range(Ntraj):
                 
        aj,qj,pj = a[j], x[j], p[j] 
        
        for k in range(Ntraj):
        
            ak,qk,pk = a[k], x[k], p[k]
            
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
    
@numba.jit           
def d1mat(a,x,gp,p,S):
    """
    time-derivative matrix D = <gj| dq/dt * (x-qk) | gk> 
    the trajectory momentum p is different from the momentum for basis (gp)    
    """
    N = len(x)
      
 
    l = np.zeros((N,N),dtype=np.complex128)
    
    
    for j in range(N): 
        
        aj,qj,pj = a[j], x[j], gp[j] 
        
        for k in range(N):
            
            ak,qk,pk = a[k], x[k], gp[k]
        
            #d2 = dp[k]/2./am * ak * 1./(aj+ak)**2 * ( 2j*aj*(pj-pk) * (qj-qk) - aj**2 * (qj-qk)**2 + (pj-pk)**2 - (aj+ak) )
            
            d1 = ak * p[k]/am * (-(1j*(pj-pk)-aj*(qj-qk)))/(aj+ak) 
            
            l[j,k] = d1*S[j,k]
     
    return l 

@numba.jit           
def d2mat(a,x,gp,c,p,S):
    """
    time-derivative matrix D = <d/daj gj| dq/dt * (x-qk) | gk> = conjugate(cj) * pk/am <gj | (1/4/aj - (x-qj)^2/2) * ak * (x-qk) | gk>
    the trajectory momentum p is different from the momentum for basis (gp)    
    """
    N = len(x)
       
    l = np.zeros((N,N),dtype=np.complex128)
    
    for j in range(N): 
        
        aj,qj,pj = a[j], x[j], gp[j] 
        
        for k in range(N):
            
            ak,qk,pk = a[k], x[k], gp[k]

            tmp1 = aj*(qj-qk)/(aj+ak) 
            
            dq = qj - qk 

            tmp3 = ((aj**2 - 2.0 * ak**2 + aj * ak * (-1.0 + ak * (qj - qk)**2)) * dq)/(aj + ak)**3
        
            d1 = p[k]/am * np.conj(c[j]) * ( ak/4.0/aj * tmp1 - ak/2.0 * tmp3) 
            
            l[j,k] = d1*S[j,k]
     
    return l 

def Jmat(a,x,p,c,S):
    """
        matrix element < d/daj gj | gk> = conj(cj) * <gj| 1/4/aj - (x-qj)^2/2 | gk> 
    """

    N = len(x)
       
    l = np.zeros((N,N),dtype=np.complex128)
    
    for j in range(N): 
        
        aj,qj,pj = a[j], x[j], p[j] 
        
        for k in range(N):
            
            ak,qk,pk = a[k], x[k], p[k]

            dq = qj - qk 

            tmp2 = 1./(aj+ak) + aj**2*dq**2/(aj+ak)**2 

            d1 = np.conj(c[j]) * ( 1.0/4.0/aj  - tmp2/2.0) 
            
            l[j,k] = d1*S[j,k]
     
    return l 
    
    
def Gmat(a,x,p,c,S):
    """
        matrix element < d/daj gj | d/dak gk> = conj(cj) * ck * <gj| (1/4/aj - (x-qj)^2/2) * (1/4/ak - (x-qk)^2/2) | gk> 
    """
    N = len(x)
       
    l = np.zeros((N,N),dtype=np.complex128)
    
    for j in range(N): 
        
        aj,qj,pj = a[j], x[j], p[j] 
        
        for k in range(N):
            
            ak,qk,pk = a[k], x[k], p[k]

            dq = qj - qk 

            tmp22 = (1.0/((aj + ak)**4)) * (-3.0* aj * ak * (-2.0 + ak * dq**2) + \
                    ak**2 * (3.0 + ak * dq**2) + aj**2 * (3. - 3. * ak * dq**2 + ak**2 * dq**4) + aj**3 * dq**2)

            tmp2j = 1./(aj+ak) + aj**2*dq**2/(aj+ak)**2
            tmp2k = 1./(aj+ak) + ak**2*dq**2/(aj+ak)**2
            
            d = np.conj(c[j]) * c[k] * (1.0/16.0/aj/ak - 1.0/8.0/aj * tmp2k - 1.0/8.0/ak * tmp2j + 1./4.0*tmp22) 

            l[j,k] = d * S[j,k] 
    
    return l 


def norm(a,x,p,c):
    
    S = overlap(a,x,p)
    
    z = np.vdot(c,S.dot(c))
    
    return z.real

def corr(a,x,p,c):
    
    S = overlap(a,x,p)
    
    z = c.dot(S.dot(c))
    
    return z

def blkMatrix(A,B,C,D):
    """
    put four matrices into blocks to form larger matrix P =  A B 
                                                             C D 
    """
    N = len(A)
    l = np.zeros((2*N,2*N),dtype=np.complex128)

    for i in range(N):
        for j in range(N):
            
            l[i,j] = A[i,j] 
            l[i,j+N] = B[i,j]              
            l[i+N,j] = C[i,j] 
            l[i+N,j+N] = D[i,j]  
    
    return l 

def k1mat(a,x,p,c,S):
    """
        matrix element conj(cj) * < d/daj gj | K | gk> = conj(cj) * <gj| (1/4/aj - (x-qj)^2/2) * (-1/2m) * (ak^2*(x-qk)^2 - ak) | gk> 
    """
    N = len(x)
       
    l = np.zeros((N,N),dtype=np.complex128)
    
    for j in range(N): 
        
        aj,qj,pj = a[j], x[j], p[j] 
        
        for k in range(N):
            
            ak,qk,pk = a[k], x[k], p[k]

            dq = qj - qk 

            tmp22 = (1.0/((aj + ak)**4)) * (-3.0* aj * ak * (-2.0 + ak * dq**2) + \
                    ak**2 * (3.0 + ak * dq**2) + aj**2 * (3. - 3. * ak * dq**2 + ak**2 * dq**4) + aj**3 * dq**2)

            tmp2j = 1./(aj+ak) + aj**2*dq**2/(aj+ak)**2
            tmp2k = 1./(aj+ak) + ak**2*dq**2/(aj+ak)**2
            
            d = -1.0/2.0/am * ( 1.0/4.0/aj*ak**2 * tmp2k - ak/4.0/aj - ak**2/2.0 * tmp22 + ak * tmp2j) 
            
            l[j,k] = d * S[j,k] 
    
    return l
    
def v1mat(a,x,p,c,S):
    """
        matrix element conj(cj) * < d/daj gj | V | gk> = conj(cj) * <gj| (1/4/aj - (x-qj)^2/2) * ( V0 + V'(x-x0) + V''/2 * (x-x0)^2) | gk> 
    """
    N = len(x)
       
    l = np.zeros((N,N),dtype=np.complex128)
    
    for j in range(N): 
        
        aj,qj,pj = a[j], x[j], p[j] 
        
        for k in range(N):
            
            ak,qk,pk = a[k], x[k], p[k]
            
            x0 = (aj*qj + ak*qk)/(aj+ak) 

            v0,v1,ddv = dv(x0)

            dq = qj - qk 

            tmp2j = 1./(aj+ak) + aj**2*dq**2/(aj+ak)**2
            tmp2ave = 1.0/(aj+ak) 
            tmp2j1ave = (- 2.0 * ak * dq)/(aj + ak)**2
            tmp2j2ave = ( 3.0 * aj + ak * (3.0 + ak * dq**2) )/(aj + ak)**3 
            
            d = 1.0/4.0/aj * (v0 + ddv/2.0 * tmp2ave) - v0 * tmp2j - v1/2.0*tmp2j1ave - ddv/4.0 * tmp2j2ave
            
            l[j,k] = d * S[j,k] 
    
    return l


def cnum(S):
    """
    check condition number of a matrix S 
    """
    d = np.linalg.cond(S)
    
    if d > 1e6: 
        sys.exit('Singular matrix')        
        #print('number of basis = {} \n'.format(len(x)))
    return d
##################### MAIN CODE ###################

### define the initial wavefunction, which usually is a Gaussian 

# initial wavefunction is expanded as a linear combination of GWPs
nold = 1
a0 = 2.0
x0 = -3.0
p0 = 3.0
am = 1.0

### define a over-complete dictionary for matching-pursuit 

nb = 8     # number of elements in the dictionary 
#gx = np.random.randn(N) / np.sqrt(2.0 * a0) + x0 
cut = 1e-5
xmin = x0 - np.sqrt(-np.log(cut/np.sqrt(a0/np.pi))/a0)      
xmax = x0 + np.sqrt(-np.log(cut/np.sqrt(a0/np.pi))/a0)

x = np.linspace(xmin,xmax,nb)
print('configuration space = {}, {} \n'.format(xmin,xmax))

p = np.zeros(nb)
#gp = np.random.randn(N)
#gp += p0
a = np.array([4.0+0j] * nb) # basis of GWPs width 

### OMP, generate a set of GWPs 
#x, p, c, index = omp(a0, x, p, c)

print('Start the propagation ...')  

# check singularity of overlap matrix

print('number of basis = {} \n'.format(nb))

print('basis sets ... \n')
print('width parameters ',a,'\n')


### INITIAL PARAMETERS 
Nt, dt = 1, 0.001
t = 0.0 
dt2 = dt/2.0 

print('time steps = {}, time interval = {}'.format(Nt,dt))

f_traj = open('traj.dat', 'w')
f_cor = open('cor.dat', 'w')

nout = 20 
fmt = ' {} ' * (nout+1) + ' \n' 


c = iproj(a,x,p,a0,x0,p0)
print('initial expansion coefficients ...\n')
print(c)

#cold = c 
#dc = solve(x,p,c)
#c += dc*dt 
#
#for k in range(Nt):
#
#    t = t + dt
#    
#    P_all = get_p(gx,x,p,cold)
#    
#    gx += P_all * dt / am 
#    x = [gx[i] for i in index]
#    p = [gp[i] for i in index]
#
#    dc = solve(x,p,c)
#    
#    
#    cnew = cold + dc*2.0*dt
#    cold = c 
#    c = cnew
    
for k in range(Nt):
    
    t += dt 
    
    P_all = get_p(x,a,x,p,c)
    x += P_all * dt / am 

    
    # update c     
    dz = solve(a,x,p,c)
    
    dc = dz[0:nb]
    da = dz[nb:2*nb]
    
    c += dc * dt
    a += da * dt 
    
    
    
    print('width',a,'\n')
    
    #print('norm = {} \n'.format(norm(x,p,c)))
    
    # renormalization 
    S = overlap(a,x,p)

    anm = np.vdot(c,S.dot(c)).real    
    c = c/np.sqrt(anm)

    #f_cor.write(' {} {} \n'.format(t*2,corr(a,x,p,c)))
    #f_traj.write(fmt.format(t,*gx[0:nout]))
print('============================')
print('*   propagation complete   *')
print('============================.')    

# close files 
f_traj.close() 
f_cor.close() 



# output data 
xgrid = np.linspace(-6,3,200)
#phase_p = get_p(xgrid, x,p,c,beta)
#f = open('u.out','w')
#for i in range(len(xgrid)):
#    f.write('{} {} \n'.format(xgrid[i], phase_p[i]))
#f.close() 




f = open('fit.out', 'w')

psi = (a0/np.pi)**0.25 * np.exp(-a0*(xgrid - x0)**2/2.0 + 1j * p0 * (xgrid-x0))
#        + cold[1] * np.exp(-a0*(xgrid - xold[1])**2/2.0 + 1j * pold[1] * (xgrid-xold[1])))

#psi_approx = c[0] * (a/np.pi)**0.25 * np.exp(-a*(xgrid - x[0])**2/2.0 + 1j * p[0] * (xgrid-x[0]))
#if len(x) > 1:
#    for j in range(1,len(x)):
#        psi_approx += c[j] * (a/np.pi)**0.25 * np.exp(-a*(xgrid - x[j])**2/2.0 + 1j * p[j] * (xgrid - x[j]))

#for i in range(len(xgrid)):
#    f.write('{} {} {} \n'.format(xgrid[i], np.abs(psi[i]), np.abs(psi_approx[i])))

f.close()
