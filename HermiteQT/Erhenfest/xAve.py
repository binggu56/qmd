

import numpy as np 
import matplotlib.pyplot as plt 

plt.subplot(111)

dat = np.genfromtxt('xAve.dat')
plt.plot(dat[:,0],dat[:,1],'-',label='wft')
plt.plot(dat[:,0],dat[:,3],'-',label='wft')

dat = np.genfromtxt('/Users/bing/lab/spo/spo_1d/xAve.dat')
plt.plot(dat[:,0],dat[:,1],'--',label='QM')
plt.plot(dat[:,0],dat[:,2],'--',label='QM')


plt.legend()
plt.show() 
