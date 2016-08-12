"""
orthogonal matching-pursuit algorithm for quantum dynamics 

@date: June 25, 2016 
@author: Bing Gu 

""" 
import numpy as np 
import numba

#define a over-complete dictionary 

N = 400 # number of elements in dictionary  
gx = np.random.randn(N) * 2.0
gp = np.random.randn(N)

@numba.autojit 
def gint(x,px,aj,y,py,ak):
    dp = py - px 
    dq = y - x 
    
    return (aj*ak)**0.25 * np.sqrt(2./(aj+ak)) * np.exp(    \
            -0.5 * aj*ak/(aj+ak) * (dp**2/aj/ak + dq**2  \
            + 2.0*1j* (px/aj + py/ak) *dq) )

@numba.autojit 
def iproj(x,px,ax):
    """
    compute the projection of intial wavefunction on basis k 
    the intial wavefunction is assumed to be a linear combination of GWPs
    """
    z = 0.+0j 
    for j in range(nold):
        d = gint(gx[k],gp[k],a, xold[j], pold[j],a0)
        z += d * cold[j]
    return z 


def overlap(x,p,a = 1.): 
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

# initial wavefunction is written as a linear combination of GWPs
nold = 2
a0 = 2.0 
xold = [-1.0, 1.0]
pold = [0.0, 0.0]
s = gint(xold[0],0.,a0,xold[1],0.,a0)
cold = [1./np.sqrt(2. + 2. * s), 1./np.sqrt(2. + 2. * s)]


a = 1.0 

tmp = np.zeros(N,dtype=complex)
c0 = tmp
 
x = [] 
c = [] 
p = [] 
res = 1.0 


# projection of the target function onto each basis 
for k in range(N):        
    c0[k] = iproj(gx[k],gp[k],a)

tmp = c0 
        
i = np.argmax(np.abs(tmp))
x.append(gx[i])
p.append(gp[i])
c.append(c0[i])
    
nmax = N # maximum number of basis defined by user 
    
res = 1.0 - np.abs(c[0])**2

b = c  # initial orthogonal expansion coefficients 

while res > 1e-4 and len(x) < nmax: 
    
    tmp = c0 
    
    for i in range(len(x)):
        d = [gint(gx[k],gp[k],a, x[i], p[i], a) for k in range(N)]
        tmp = tmp - b[i] * np.array(d)  
    
    i = np.argmax(np.abs(tmp))

    x.append(gx[i])
    p.append(gp[i])
    c.append(c0[i])
    
    # orthogonalization of obtained basis so far 
    
    S = overlap(x,p,a)

    b = np.linalg.solve(S,c)
    
    res = 1.0 + np.vdot(b,S.dot(b)).real - np.vdot(b,c).real * 2.0  


    
nb = len(x) # the number of basis function 
print('number of basis = {} \n residue = {} \n'.format(nb,res))



f = open('fit.out', 'w')
xgrid = np.linspace(-4,4,200)
psi = (a0/np.pi)**0.25 * (cold[0] * np.exp(-a0*(xgrid - xold[0])**2/2.0 + 1j * pold[0] * (xgrid-xold[0])) \
        + cold[1] * np.exp(-a0*(xgrid - xold[1])**2/2.0 + 1j * pold[1] * (xgrid-xold[1])))

psi_approx = b[0] * (a/np.pi)**0.25 * np.exp(-a*(xgrid - x[0])**2/2.0 + 1j * p[0] * (xgrid-x[0]))
if nb > 0:
    for j in range(1,nb):
        psi_approx += b[j] * (a/np.pi)**0.25 * np.exp(-a*(xgrid - x[j])**2/2.0 + 1j * p[j] * (xgrid - x[j]))

for i in range(len(xgrid)):
    f.write('{} {} {} \n'.format(xgrid[i], np.abs(psi[i]), np.abs(psi_approx[i])))

f.close()
