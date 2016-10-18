import numpy as np 
import matplotlib.pyplot as plt 

#use parameters

pi=3.141592653589793238e0
auev=27.2113845e0
hbarinv=1.51924888335207073622e0
hbar=0.65822e0
qe = 1.6022e-19
permi0 = 8.854187817e-12
mgauleg1=8

boltz  = 1.3806506e-23
ev2j   = 1.60210e-19
dsmall = 1.0e-8

freqau=1.378999779e6
aufs=2.418884326505e-2
cunity=1.0+0j
czero=0.0+0j 
eye= 0.0+1.0j

dpt = [-0.96028986e0, -0.79666648e0, -0.52553241e0, -0.18343464e0,  0.18343464e0, 0.52553241e0,  \
    0.79666648e0, 0.96028986e0]
dw  = [ 0.10122854e0, 0.22238103e0, 0.31370665e0, 0.36268378e0, 0.36268378e0,  0.31370665e0,\
        0.22238103e0, 0.10122854e0]

#use blkMatrix
#use dataType

norbs = 20
nsite = 10
npade = 20
onsite1 = 2.0
onsite2 = -2.0
hop1 = -1.0
hop2 = -0.7 
alpha=4.1
gammaL=0.05
gammaR=0.05

lc = unit_length = 2.0 

telec=300.0e0
fieldtyp=0
tmax=1500.0e0 
dt=0.05

al=0.1
ar=0.1
nener = 2000
engyl = -10.0
engyr = 10.0 

nLead = 2
miu0 = 0.0 

beta = ev2j / (boltz * telec) 
im = 0.0+1.0j



#type(dBlkSequence), intent(in) :: Lambda(nLead)

fermi = np.zeros(nLead)
curr = np.zeros((nLead, nLead))

downrange=0.0
uprange=0.0

AL = np.zeros((norbs,norbs))
AR = np.zeros((norbs,norbs))
AL[0,0] = gammaL/2.0
AL[1,0] = gammaL/2.0
AL[1,1] = gammaL/2.0
AL[0,1] = gammaL/2.0

AR[0,0] = gammaR/2.0
AR[0,1] = gammaR/2.0
AR[1,0] = gammaR/2.0
AR[1,1] = gammaR/2.0

Lambda = [AL, AR]

def field(t):
    
    sigma = 100.0 
    hbarinv = 1.51924888335207073622
    A = 0.5 # strength of the field
    T0 = 500.0 # center of time range for the field 
    phase = 0.0 
    freq = 0.10
    phi = A * np.exp(-(t-T0)**2/(2*sigma**2))/(freq * hbarinv)    \
              * ( (t-T0)/(sigma**2) * np.sin( freq *(t-T0) * hbarinv + phase) \
              - freq * hbarinv * np.cos( freq * (t-T0) * hbarinv + phase) )

    return phi

def ECP(t,bias=3e0):

    global nsite, unit_length
    
    efactor = 1.0e0
    E = field(t)
    left = bias - efactor * E * (float(nsite+1)/2.0) * unit_length
    right = -bias + efactor * E * (float(nsite+1)/2.0) * unit_length
    
    return left, right 
    
def fock(t):
    
    global onsite1, onsite2, hop1, hop2, norbs, nsite

    unit_length = 1.22e0 
    
    efactor = 1.0 
    
    fock0 = np.zeros((norbs,norbs))
    
    for i in range(nsite):
        fock0[2*i, 2*i] = onsite1
        fock0[2*i+1,2*i+1] = onsite2
    
    for i in range(nsite):
        fock0[2*i,2*i+1] = hop1
        fock0[2*i+1, 2*i] = hop1
        
    for i in range(nsite-1):
        fock0[2*i+2, 2*i+1]  = hop2
        fock0[2*i+1, 2*i+2] = hop2
    
    fock1 = np.zeros((norbs,norbs))
    for i in range(nsite):
        fock1[2*i,2*i] = (i+1-float(nsite+1)/2.0) * unit_length                                                                       
        fock1[2*i+1,2*i+1] = fock1[2*i, 2*i] 
    
    fock = fock0 + efactor*field(t) * fock1
    
    return fock


#fock0 = fock(438.345864662)

t = 458.90 # minimum of field
#t = 438.345864662 # maximum of field

fock0 = fock(t)


ampl, ampr = ECP(t,bias=1.0)
amp = np.array([ampl, ampr])

class Lead:
    def __init__(self, nStart, nSigma):
        self.nStart = nStart 
        self.nSigma = nSigma 

Lead = [Lead(0,2), Lead(norbs-2,2)]
print('\n ========= enter  Landauer_m.f90 =========== \n')

def landauer(nLead, miu0, beta, amp, norbs, fock0, Lambda, Lead):

#read(5,nml=land)
    tlamda = np.zeros((norbs,norbs)) # self energy matrix 
    for iLead in range(nLead):
        i = Lead[iLead].nStart
        j = Lead[iLead].nSigma
        tlamda[i:i+j, i:i+j] = tlamda[i:i+j, i:i+j] + Lambda[iLead][0:j,0:j]
    
    
    for i in range(norbs):
        tlamda[i, i] += 2.0e-7 # add a small imaginary part
    
    f = open('transmission.dat', 'w')
    
    curr = np.zeros((nLead, nLead))
    
    for iLead in range(nLead):
        for jLead in range(iLead+1, nLead):

            print('from  {} to {} \n'.format(jLead,iLead))
            
            enerstar= miu0 + min(amp[iLead],amp[jLead]) - downrange   #initial point of energy
            enerend = miu0 + max(amp[iLead],amp[jLead]) + uprange     #end point of energy
    
            dener = (enerend-enerstar)/float(nener)     # integration step
            
            print(' integration range = {} {} \n '.format(enerstar, enerend))
            print(' integration step = {} \n'.format(dener))
            
            for k in range(nener):
                
                energy = enerstar + float(k)*dener
                
                Id = np.identity(norbs, dtype=np.complex128)
                
                gr = energy * Id - (fock0 - im * tlamda)

        #call zgetrf(norbs, norbs, greenr, norbs, ipiv, info)
        #call zgetri(norbs, greenr, norbs, ipiv, ctmp, norbs*norbs, info)
                gr = np.linalg.inv(gr)
                
                tmp1 = (energy-miu0-amp[iLead])*beta
                tmp2 = (energy-miu0-amp[jLead])*beta
                
                fermi = np.zeros(nLead)
                fermi[iLead] = Fermi(tmp1)
                fermi[jLead] = Fermi(tmp2)
    
    
                GammaL = np.zeros((norbs,norbs),dtype=np.complex128)
                i = Lead[iLead].nStart
                j = Lead[iLead].nSigma
                GammaL[i:i+j, i:i+j] += Lambda[iLead][0:j,0:j] + 0e0j 
#                # Lambda : linewidth function 
#                #call zgemm('c', 'n', norbs, norbs, norbs, cunity, greenr, norbs, ctmp0, norbs, czero, ctmp2, norbs)
#                ctmp0 = czero
                GammaR = np.zeros((norbs,norbs),dtype=np.complex128)
                i = Lead[jLead].nStart
                j = Lead[jLead].nSigma
                GammaR[i:i+j, i:i+j] += Lambda[jLead][0:j,0:j] + 0.0j 
                #call zgemm('n', 'n', norbs, norbs, norbs, cunity, ctmp0, norbs, ctmp2, norbs, czero, ctmp1, norbs)
                #call zgemm('n', 'n', norbs, norbs, norbs, cunity, greenr, norbs, ctmp1, norbs, czero, ctmp, norbs)
                #GammaL, GammaR = Lambda                 
                T = np.dot(gr, np.dot(GammaR, np.dot(np.conjugate(np.transpose(gr)), GammaL))) 
                                
                tmp = trace(T).real/np.pi**2
                
                tmp2 = 2.0*np.pi * tmp * (fermi[iLead]-fermi[jLead]) * dener
            
                f.write('{} {:14.7e} {:14.7e} \n'.format(energy, 4.0*np.pi**2*tmp, tmp2))
            
                curr[iLead,jLead] += tmp2 
                
            
    curr = 2e0 * curr * 1.6022e5 * hbarinv # to nA, 2.d0 comes from spin 
    
    for iLead in range(nLead):
        for jLead in range(iLead + 1, nLead):
            print('current from {} lead to {} lead is {}'.format(jLead,iLead,curr[iLead,jLead]))
    
    print('\n ========= leave  Landauer_m.f90 =========== \n')

def Fermi(x):
    """
    Fermi distribution function 1/(1+exp(x))
    """
    if x > 100:
        return 0e0 
    else:
        return 1e0/(1e0+np.exp(x))

def trace(A):
    """
    compute trace of a matrix 
    """
    tmp = 0.0+0.0j
    for i in range(A.shape[0]):
        tmp += A[i,i]
    
    return tmp 

def plt_trans(ax):
    dat = np.genfromtxt('transmission.dat')
    ax.plot(dat[:,0],dat[:,1])
    ax.set_yscale('log')
    plt.show()
    
print('beta = {} eV^-1'.format(beta))
landauer(nLead, miu0, beta, amp, norbs, fock0, Lambda, Lead)

fig,ax = plt.subplots()
plt_trans(ax) 