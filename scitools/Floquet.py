import numpy as np
from scipy import linalg
import sys


def quasiE(H0, H1, Nt, omega):
    """
    Construct the Floquet hamiltonian of size Norbs * Nt and Diagonalize

    INPUT: (in order)
        H0 : system Hamiltonian
        H1    : drive H1 * cos(omega * t)
        Nt    : number of Fourier components
        omega: drive frequency

    OUTPUT:
        eigvals_subset: Quasienergies in First BZ
        eigvecs_subset: Floquet modes
    """
    # diagonalize system H
    #U, sysEigvals = Hamilton(Norbs)

    # transform basis to eigenstats of system H
    #M = basisTransform(U)

    Norbs = H0.shape[-1]


    # dimensionality of the Floquet matrix
    NF = Norbs * Nt
    F = np.zeros((NF,NF), dtype=np.complex128)

    N0 = -(Nt-1)/2 # starting point for Fourier companent of time exp(i n w t)

    # construc the Floquet H for a general tight-binding Hamiltonian
    for n in range(Nt):
        for m in range(Nt):

#            # atomic basis index
#            for k in range(Norbs):
#                for l in range(Norbs):
#
#                # map the index i to double-index (n,k) n : time Fourier component
#                # with relationship for this :  Norbs * n + k = i
#
#                    i = Norbs * n + k
#                    j = Norbs * m + l
            # fill a block of the Floquet Matrix
            istart = Norbs * n
            iend = Norbs * (n+1)
            jstart = Norbs * m
            jend = Norbs * (m+1)

            F[istart:iend, jstart:jend] = HamiltonFT(H0, H1, n-m) + (n + N0) \
                     * omega * delta(n, m) * np.eye(Norbs)


    # for a two-state model

#    for n in range(Nt):
#        for m in range(Nt):
#            F[n * Norbs, m * Norbs] = (N0 + n) * omega * delta(n,m)
#            F[n * Norbs + 1, m * Norbs + 1] = (onsite1 + (N0+n) * omega) * delta(n,m)
#            F[n * Norbs, m * Norbs + 1] = t * delta(n,m+1)
#            F[n * Norbs + 1, m * Norbs] = t * delta(n,m-1)
    #print('\n Floquet matrix \n', F)

    # compute the eigenvalues of the Floquet Hamiltonian,
    eigvals, eigvecs = linalg.eigh(F)

    #print('Floquet quasienergies', eigvals)

    # specify a range to choose the quasienergies, choose the first BZ
    # [-hbar omega/2, hbar * omega/2]
    eigvals_subset = np.zeros(Norbs)
    eigvecs_subset = np.zeros((NF , Norbs), dtype=np.complex128)


    # check if the Floquet states is complete
    j = 0
    for i in range(NF):
        if  eigvals[i] < omega/2.0 and eigvals[i] > -omega/2.0:
            eigvals_subset[j] = eigvals[i]
            eigvecs_subset[:,j] = eigvecs[:,i]
            j += 1
    if j != Norbs:
        print("Error: Number of Floquet states {} is not equal to \
              the number of orbitals {} in the first BZ. \n".format(j, Norbs))
        sys.exit()


#    # now we have a complete linear independent set of solutions for the time-dependent problem
#    # to compute the coefficients before each Floquet state if we start with |alpha>
#    # At time t = 0, constuct the overlap matrix between Floquet modes and system eigenstates
#    # G[j,i] = < system eigenstate j | Floquet state i >
#    G = np.zeros((Norbs,Norbs))
#    for i in range(Norbs):
#        for j in range(Norbs):
#            tmp = 0.0
#            for m in range(Nt):
#                tmp += eigvecs_subset[m * Norbs + j, i]
#            G[j,i] = tmp
#
#
#    # to plot G on site basis, transform it to site-basis representation
#    Gsite = U.dot(G)

    return eigvals_subset, eigvecs_subset

def delta(i,j):
    if i==j:
        return 1
    else:
        return 0

def HamiltonFT(H0, H1, n):
    """
    Fourier transform of the full Hamiltonian, required to construct the
    Floquet Hamiltonian

    This is particulary for drive with cos(omega *t)

    INPUT
        n : Fourier component index
        H0: system Hamiltonian
        H1 : perturbation Hamiltonian
    """
    Norbs = H0.shape[-1]

    if n == 0:

        return H0

    elif n == 1 or n == -1:
        return 0.5 * H1

    else:
        return np.zeros((Norbs,Norbs))




