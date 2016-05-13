

import numpy as np 
import matplotlib.pyplot as plt 

plt.subplot(111)

dat = np.genfromtxt('energy.dat')

plt.plot(dat[:,0],dat[:,1],'-',label='K')
plt.plot(dat[:,0],dat[:,2],'--',label='V')
plt.plot(dat[:,0],dat[:,3],'--',label='U')

en = dat[:,1]+dat[:,2]+dat[:,3] 

plt.plot(dat[:,0],en,'--',label='E')


plt.legend()
plt.show() 
