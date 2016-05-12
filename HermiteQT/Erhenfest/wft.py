

import numpy as np 
import matplotlib.pyplot as plt 

plt.subplot(111)

dat = np.genfromtxt('wft.dat')
#plt.plot(dat[:,0],dat[:,1],'-',label='wft')
#plt.plot(dat[:,0],dat[:,2],'-',label='wft')
plt.plot(dat[:,0],(dat[:,1]**2+dat[:,2]**2),'b--',lw=2)

dat = np.genfromtxt('/Users/bing/lab/spo/spo_1d/wft')
#plt.plot(dat[:,0],dat[:,1],'--',label='QM')
#plt.plot(dat[:,0],dat[:,2],'--',label='QM')
plt.plot(dat[:,0],dat[:,3],'k-',lw=3,label='QM')


plt.legend()
plt.show() 
