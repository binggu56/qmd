# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 10:20:48 2016

@author: bing
"""

import numpy as np 

DIM = 3 

def V(r):
    rSqd = r.dot(r)
    return 0.5 * rSqd 

def ensureCapacity(N):

    maxN = 0       # remember the size of the array

    if N < maxN:          # no need to expand array
        return                # do nothing

    oldMaxN = maxN            # remember the old capacity
    
    if maxN > 0:
        maxN *= 2             # double capacity
    else:
        maxN = 1
        
    if N > maxN - 1:      # if this is not sufficient
        maxN = N + 1     # increase it so it is sufficient

    # allocate new storage
    rNew = np.zeros(maxN, DIM) 
    newAlive = [] 
    for n in range(maxN):
        if n < oldMaxN:      # copy old values into new arrays
            rNew[n][0:DIM] = r[n][0:DIM]
            newAlive[n] = alive[n] 

    #delete [] r;               // release old memory
    #r = rNew;                  // point r to the new memory
    #delete [] alive;
    #alive = newAlive;


# We need to measure the energy, its variance, and the wave function of the ground state.

# observables
rMax = 4                # max value of r to measure psi
NPSI = 100              # number of bins for wave function
psi = np.zeros(NPSI)              # wave function histogram

def zeroAccumulators():
    global ESum, ESqdSum, psi 

    ESum = 0 
    ESqdSum = 0 
    
    psi = np.zeros(NPSI)
    return 



N = N_T                   # set N to target number specified by user
for n in range(N):
    ensureCapacity(n)
    for d in range(DIM): 
        r[n,d] = np.random.uniform() - 0.5    
    alive[n] = True

zeroAccumulators()
E_T = 0                   # initial guess for the ground state energy



def oneMonteCarloStep(n):

    # Diffusive step
    for d in range(DIM):
        r[n,d] += np.random.normal() * np.sqrt(dt);

    x = r[n,0:DIM]
    
    # Branching step
    q = np.exp(- dt * (V(x) - E_T))
    
    survivors = int(q)
    if q - survivors > np.random.uniform():
        ++survivors

    # append survivors-1 copies of the walker to the end of the array
    for n in range(survivors - 1):
        ensureCapacity(N);
        for (int d = 0; d < DIM; d++)
            r[N][d] = r[n][d];
        alive[N] = true;
        ++N;
    }

    // if survivors is zero, then kill the walker
    if (survivors == 0)
        alive[n] = false;
}

def oneTimeStep() {

    // DMC step for each walker
    int N_0 = N;
    for (int n = 0; n < N_0; n++)
        oneMonteCarloStep(n);

    // remove all dead walkers from the arrays
    int newN = 0;
    for (int n = 0; n < N; n++)
    if (alive[n]) {
        if (n != newN) {
            for (int d = 0; d < DIM; d++)
                r[newN][d] = r[n][d];
            alive[newN] = true;
        }
        ++newN;
    }
    N = newN;

    // adjust E_T
    E_T += log(N_T / double(N)) / 10;

    // measure energy, wave function
    ESum += E_T;
    ESqdSum += E_T * E_T;
    for (int n = 0; n < N; n++) {
        double rSqd = 0;
        for (int d = 0; d < DIM; d++)
            rSqd = r[n][d] * r[n][d];
        int i = int(sqrt(rSqd) / rMax * NPSI);
        if (i < NPSI)
            psi[i] += 1;
    }
}


int main() {

    cout << " Diffusion Monte Carlo for the 3-D Harmonic Oscillator\n"
         << " -----------------------------------------------------\n";
    cout << " Enter desired target number of walkers: ";
    cin >> N_T;
    cout << " Enter time step dt: ";
    cin >> dt;
    cout << " Enter total number of time steps: ";
    int timeSteps;
    cin >> timeSteps;

    initialize();

    // do 20% of timeSteps as thermalization steps
    int thermSteps = int(0.2 * timeSteps);
    for (int i = 0; i < thermSteps; i++)
        oneTimeStep();

    // production steps
    zeroAccumulators();
    for (int i = 0; i < timeSteps; i++) {
        oneTimeStep();
    }

    // compute averages
    double EAve = ESum / timeSteps;
    double EVar = ESqdSum / timeSteps - EAve * EAve;
    cout << " <E> = " << EAve << " +/- " << sqrt(EVar / timeSteps) << endl;
    cout << " <E^2> - <E>^2 = " << EVar << endl;
    double psiNorm = 0, psiExactNorm = 0;
    double dr = rMax / NPSI;
    for (int i = 0; i < NPSI; i++) {
        double r = i * dr;
        psiNorm += pow(r, DIM-1) * psi[i] * psi[i];
        psiExactNorm += pow(r, DIM-1) * exp(- r * r);
    }
    psiNorm = sqrt(psiNorm);
    psiExactNorm = sqrt(psiExactNorm);
    ofstream file("psi.data");
    for (int i = 0; i < NPSI; i++) {
        double r = i * dr;
        file << r << '\t' << pow(r, DIM-1) * psi[i] / psiNorm << '\t'
             << pow(r, DIM-1) * exp(- r * r / 2) / psiExactNorm << '\n';
    }
    file.close();
}

