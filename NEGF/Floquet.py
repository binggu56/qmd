# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 14:54:11 2017

Floquet theory for periodic Hamiltonian 

@author: bingg
"""

import numpy as np 
import matplotlib.pyplot as plt 

import sys 

def delta(i,j):
    if i==j: 
        return 1
    else:
        return 0 
        
def ham(Norbs):
    """
    Norbs: number of atomic orbitals 

    """

    
    h = np.zeros((Norbs,Norbs))
    # diagonal elements     
    #for i in range(Norbs):
    #    h[i,i] = E0 

    h[0,0] = E0 
    h[1,1] = E1 
    
    # off-diagonal elements     
    for i in range(Norbs-1):
        h[i,i+1] = t 
    for i in range(1,Norbs):
        h[i+1,i] = t
    
    return h 

def Floquet_ham(Norbs, Nt):
    """
    construct a Floquet hamiltonian of the size Norbs * Nt  
    """
    # frequency of the radiation  

    
    NF = Norbs * Nt 
    F = np.zeros((NF,NF))
    
#    for n in range(Nt):
#        for k in range(Norbs):
#        # we need to map the index i to double-index (n,k) n : time Fourier component 
#        # k : atomic bassi index
#        
#        # a relationship for this is Norbs * n + k ---> i
#            i = Norbs * n + k 
#            for m in range(Nt):
#                for l in range(Norbs):
#                    j = Norbs * m + l 
#                    
#                    F[i,j] =  
    
    # for a two-state model 
    
    # starting point of Fourier component exp(i n\omega t)
    N0 = -(Nt-1)/2 
    for n in range(Nt):
        for m in range(Nt):
            F[n * Norbs, m * Norbs] = t * (N0 + n) * omega * delta(n,m)
            F[n * Norbs + 1, m * Norbs + 1] = (E1 + t * (N0+n) * omega) * delta(n,m)
            F[n * Norbs, m * Norbs + 1] = t * delta(n,m+1)
            F[n * Norbs + 1, m * Norbs] = t * delta(n,m-1)
            
    # after construct the Floquet Hamiltonian, the eigenvalues are computed 
    eigvals, eigvecs = np.linalg.eigh(F)
    
    # specify a range to choose the quasienergies, choose [-hbar omega/2, hbar * omega/2]
    eigvals_subset = np.zeros(Norbs)
    eigvecs_subset = np.zeros((NF , Norbs))
    
    for i in range(NF):
        j = 0
        if  -omega/2.0 <= eigvals[i] <= omega/2.0:
            eigvals_subset[j] = eigvals[i]
            eigvecs_subset[:,j] = eigvecs[:,i]
            j += 1 
    if j != Norbs: 
        print("Error: nuber of Floquet states is not equal to the number of orbitals.")
        sys.exit() 
        
        
    # now we have a complete linear independent set of solutions for the TD problem 
    # to compute the coefficients before each Floquet state if we start with |alpha>
    # at time 0, this is done by solving a matrix equation CG = 1 
    G = np.zeros((Norbs,Norbs))
    for i in range(Norbs):
        for j in range(Norbs):
            tmp = 0.0 
            for m in range(Nt):
                tmp += eigvecs_subset[m * Norbs + j, i]
            G[i,j] = tmp
            
    C = np.linalg.inv(G)
    
    
    return F 
    
Nt = 9 # has to be odd integer 
Norbs = 2 

omega = 1.0
 
E0 = 0.0
E1 = 2.0
t = 1.0

quasiE0 = (1.0-np.sqrt(5.0))/2.0 + 1.0
quasiE1 = (1.0+np.sqrt(5.0))/2.0 - 2.0

print('Quasienergies')
print(quasiE0, quasiE1)

F = Floquet_ham(Norbs, Nt)


eigvals.sort()



x = range(Norbs * Nt)
plt.plot(x , eigvals,'-o',markersize=8)
plt.axhline(y = quasiE0,lw=1)
plt.axhline(y = quasiE1,lw=1)

plt.show()




                    