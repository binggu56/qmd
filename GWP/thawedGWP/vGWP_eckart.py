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
def gint(x,px,aj,y,py,ak):
    """
    Gaussian integrals     
    """
    dp = py - px 
    dq = y - x 
    
    return (aj*ak)**0.25 * np.sqrt(2./(aj+ak)) * np.exp(    \
            -0.5 * aj*ak/(aj+ak) * (dp**2/aj/ak + dq**2  \
            + 2.0*1j* (px/aj + py/ak) *dq) )

@numba.autojit 
def iproj(x,px,ax,aold,xold,pold,cold):
    """
    compute the projection of intial wavefunction on basis k 
    the intial wavefunction is assumed to be a linear combination of GWPs
    """
    z = 0.+0j 
    for j in range(len(xold)):
        d = gint(x, px, a, xold[j], pold[j], aold)
        z += d * cold[j]
    return z 


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
    

def get_p(xgrid, x,p,c, smooth=False):
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
                
                
                qk, ak = x[k], a 
                
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
    
def solve(x,p,c):
        
    
    S = overlap(x,p) 
    
    P = get_p(x,x,p,c) # phase momentum dS/dx 
    
    D = dmat_FGA(x,p,P,S)
     
    V = vmat_FGA(x,p,S) 
    K = kmat_FGA(x,p,S)
    
    #enk = np.vdot(c,K.dot(c))
    #env = np.vdot(c,V.dot(c))               
    
    H = (K+V)-1j*D

    b = np.dot(H,c)
    
    b = -1j*b
    
    #print 'overlap',S
 
    try:
        dc = np.linalg.solve(S, b)
        
    except:
        print("Error: ill-conditioned overlap matrix for dc")
        print(S) 
        sys.exit()
    
    return dc 



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
        ddv = 6.0*D*a**2*np.sinh(a*x)**2/np.cosh(a*x)**4-2.*D*a**2/np.cosh(a*x)**2

    
    else:
        print("ERROR: there is no such potential.")
            
    return v0,ddv
    
    

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
        
            v0,ddv = dv(x0)
        
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
    
    
    for j in range(Ntraj): 
        
        aj,qj,pj = a[j], x[j], gp[j] 
        
        for k in range(Ntraj):
            
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

            tmp3 = ((aj**2 - 2.0 * ak**2 + aj * ak * (-1.0 + ak * (qj - qk)**2)) * (xj - xk))/(aj + ak)**3
        
            d1 = p[k]/am * np.conj(c[j]) * ( ak/4.0/aj * tmp1 - ak/2.0 * tmp3) 
            
            l[j,k] = d1*S[j,k]
     
    return l 

def Jmat(a,x,p,c,S):
    """
        matrix element < d/daj gj | gk> = conj(cj) * <gj| 1/4/aj - (x-qj)^2/2 | gk> 
    """
    
    tbd 

def norm(x,p,c):
    
    S = overlap(x,p)
    
    z = np.vdot(c,S.dot(c))
    
    return z.real

def corr(x,p,c):
    
    S = overlap(x,p)
    
    z = c.dot(S.dot(c))
    
    return z

def blkMatrix(A,B,C,D):
    """
    put four matrices into blocks to form larger matrix P =  A B 
                                                             C D 
    """
    N = A.size[-1] 
    l = np.zeros(2*N,2*N)

    for i in range(N):
        for j in range(N):
            
            l[i,j] = A[i,j] 
            l[i,j+N] = B[i,j]  

    for i in range(N+1,2*N):
        for j in range(N):
            
            l[i,j] = C[i,j] 
            l[i,j+N] = D[i,j]  
    
    return l 

##################### MAIN CODE ###################

# define the initial wavefunction, which usually is a Gaussian 
# initial wavefunction is written as a linear combination of GWPs
nold = 1
a0 = 2.0
x0 = -3.0
p0 = 3.0
am = 1.0
 
x = [x0]
p = [p0]

#s = gint(xold[0],0.,a0,xold[1],0.,a0)
#cold = [1./np.sqrt(2. + 2. * s), 1./np.sqrt(2. + 2. * s)]
c = [1.0] 


### define a over-complete dictionary for matching-pursuit 

N = 400     # number of elements in the dictionary 
#gx = np.random.randn(N) / np.sqrt(2.0 * a0) + x0 
cut = 1e-5
xmin = x0 - np.sqrt(-np.log(cut/np.sqrt(a0/np.pi))/a0)      
xmax = x0 + np.sqrt(-np.log(cut/np.sqrt(a0/np.pi))/a0)

gx = np.linspace(xmin,xmax,N)
print('configuration space = {}, {} \n'.format(xmin,xmax))

gp = np.zeros(N)
#gp = np.random.randn(N)
#gp += p0
a = 8.0 # basis of GWPs width 

# OMP, generate a set of GWPs 
x, p, c, index = omp(a0, x, p, c)

# solve TDSE 

# check singularity of overlap matrix 
     
nb = len(x) # the number of basis function 
print('number of basis = {} \n'.format(nb))


Nt, dt = 600, 0.001
t = 0.0 
dt2 = dt/2.0 


f_traj = open('traj.out', 'w')
f_cor = open('cor.out', 'w')

nout = 20 
fmt = ' {} ' * (nout+1) + ' \n' 


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
    
    P_all = get_p(gx,x,p,c)
    gx += P_all * dt / am 

    
    # update c     
    dc = solve(x,p,c)
    
    c += dc * dt 

    
    # update basis position and momentum 
    x = [gx[i] for i in index]
    p = [gp[i] for i in index]
    
    #print('norm = {} \n'.format(norm(x,p,c)))

    S = overlap(x,p)

    # renormalization 
    anm = np.vdot(c,S.dot(c)).real    
    c = c/np.sqrt(anm)
    
    cnum = np.linalg.cond(S)
    if cnum > 1e5: 
        print('OMP at {} timestep'.format(k))
        #sys.exit('Singular matrix')
        x,p,c,index = omp(a, x,p,c)
        
        print('number of basis = {} \n'.format(len(x)))
     

    f_cor.write(' {} {} \n'.format(t*2,corr(x,p,c)))
    f_traj.write(fmt.format(t,*gx[0:nout]))
    
print('Mission Complete.')    

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

psi_approx = c[0] * (a/np.pi)**0.25 * np.exp(-a*(xgrid - x[0])**2/2.0 + 1j * p[0] * (xgrid-x[0]))
if len(x) > 1:
    for j in range(1,len(x)):
        psi_approx += c[j] * (a/np.pi)**0.25 * np.exp(-a*(xgrid - x[j])**2/2.0 + 1j * p[j] * (xgrid - x[j]))

for i in range(len(xgrid)):
    f.write('{} {} {} \n'.format(xgrid[i], np.abs(psi[i]), np.abs(psi_approx[i])))

f.close()
