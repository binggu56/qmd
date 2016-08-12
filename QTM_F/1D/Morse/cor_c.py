##!/usr/bin/python

import numpy as np
import pylab as plt
import matplotlib as mpl
import seaborn as sns

sns.set_context("poster",font_scale=1.5)
sns.set_style({'font.family':'Times New Roman'})

mpl.rcParams['lines.linewidth'] = 2

data = np.genfromtxt(fname='cor_c.dat') 

ncols = data.shape[1]

#for x in range(1,ncols):
#plt.plot(data[:,0],data[:,1],linewidth=2,label='$\Re(C_{xx})$')
fig, (ax1,ax2) = plt.subplots(1,2)
for i in range(1,5):
    ax1.plot(data[:,0],data[:,i],linewidth=2,label='c for DOF 1')
ax1.legend(loc=3)

for i in range(1,5):
    ax2.plot(data[:,0],data[:,i+5],linewidth=2,label='c for DOF 2')

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#plt.xlim(0,10)
plt.legend(loc=3)
plt.xlabel('Time [a.u.]')
#plt.ylabel('Positions')
plt.savefig('cor.pdf')
plt.show() 

