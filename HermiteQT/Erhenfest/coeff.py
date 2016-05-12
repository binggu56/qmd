

import numpy as np 
import matplotlib.pyplot as plt 

plt.subplot(111)

dat = np.genfromtxt('coeff.dat',dtype=np.complex128)

#plt.plot(dat[:,0].real,dat[:,1].real,'-',label='c0')
#plt.plot(dat[:,0].real,dat[:,1].imag,'-',label='c0')
for i in range(1,dat.shape[-1]):
    plt.plot(dat[:,0].real,np.abs(dat[:,i]),'-',label='c'+str(i))

#plt.plot(dat[:,0].real,np.abs(dat[:,2])**2,'--',label='c1')
#plt.plot(dat[:,0].real,np.abs(dat[:,3])**2,'--',label='c1')

plt.ylim(-2,2)

#plt.legend()
plt.show() 
