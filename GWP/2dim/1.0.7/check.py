import pylab as plt
import numpy as np 

dat = np.genfromtxt('qforce.dat') 
for i in range(dat.shape[-1]):
    plt.plot(dat[:,0],dat[:,i]) 

plt.show()
