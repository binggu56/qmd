##!/usr/bin/python

import numpy as np
import pylab as plt
import matplotlib as mpl
import seaborn as sns

sns.set_context("poster",font_scale=1.5)
sns.set_style({'font.family':'Times New Roman'})

mpl.rcParams['lines.linewidth'] = 2

data = np.genfromtxt(fname='cor.dat') 

ncols = data.shape[1]

#for x in range(1,ncols):
#plt.plot(data[:,0],data[:,1],linewidth=2,label='$\Re(C_{xx})$')
plt.plot(data[:,0],data[:,2],linewidth=2,label='$\Im(C_{11})$')
plt.plot(data[:,0],data[:,4],linewidth=2,label='$\Im(C_{22})$')
plt.plot(data[:,0],data[:,6],linewidth=2,label='$\Im(C_{33})$')
plt.plot(data[:,0],data[:,8],linewidth=2,label='$\Im(C_{44})$')
plt.plot(data[:,0],data[:,10],linewidth=2,label='$\Im(C_{12})$')

#plt.plot(data[:,0],data[:,3],linewidth=2,label='$\Re(C_{yy})$')
#plt.plot(data[:,0],data[:,4],linewidth=2,label='$\Im(C_{yy})$')

    

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
plt.xlim(0,40)
plt.legend(loc=3)
plt.xlabel('Time [a.u.]')
#plt.ylabel('Positions')
plt.savefig('cor.pdf')
plt.show() 

