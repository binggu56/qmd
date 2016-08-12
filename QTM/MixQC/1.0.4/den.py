##!/usr/bin/python

import numpy as np
import pylab as pl

#with open("traj.dat") as f:
#    data = f.read()
#
#    data = data.split('\n')
#
#    x = [row.split(' ')[0] for row in data]
#    y = [row.split(' ')[1] for row in data]
#
#    fig = plt.figure()
#
#    ax1 = fig.add_subplot(111)
#
#    ax1.set_title("Plot title...")    
#    ax1.set_xlabel('your x label..')
#    ax1.set_ylabel('your y label...')
#
#    ax1.plot(x,y, c='r', label='the data')
#
#    leg = ax1.legend()
#fig = plt.figure()
data = np.genfromtxt(fname='den.dat')
#dat = np.genfromtxt(fname='../2d12/en.dat')
#data = np.loadtxt('traj.dat')
#for x in range(1,10):
pl.subplot(111)
#pl.title('two-steps fitting alg')
#pl.ylabel('Energy [hartree]')
#pl.plot(data[:,0],data[:,1],'b--',linewidth=2,label='QT')
#pl.plot(dat[:,0],dat[:,2],'r-',linewidth=2)

dat = np.genfromtxt(fname='../spo_2d/1.0.0/den.dat') 
pl.plot(dat[:,0],dat[:,1],'k',linewidth=3,label='QM')
pl.plot(dat[:,0],dat[:,2],'k',linewidth=3)
#pl.plot(data[:,0],data[:,4],'k-',linewidth=2,label='Energy')
#pl.legend(bbox_to_anchor=(0.5, 0.38, 0.42, .302), loc=3,ncol=1, mode="expand", borderaxespad=0.)

dat1 = np.genfromtxt(fname='dom4/den.dat') 
pl.plot(dat1[:,0],dat1[:,1],'r',linewidth=2,label='L = 4')
pl.plot(dat1[:,0],dat1[:,2],'r',linewidth=2,label='L = 4')

dat2 = np.genfromtxt(fname='dom8/den.dat') 
pl.plot(dat2[:,0],dat2[:,1],'g',linewidth=2,label='L = 8')
pl.plot(dat2[:,0],dat2[:,2],'g',linewidth=2)

dat3 = np.genfromtxt(fname='dom16/den.dat') 
pl.plot(dat3[:,0],dat3[:,1],'b',linewidth=2,label='L = 16')
pl.plot(dat3[:,0],dat3[:,2],'b',linewidth=2)
#pl.subplot(212)
#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
pl.xlim(0,5)
#pl.xlabel('time [a.u.]')
##pl.ylabel('Energy [hartree]')
#pl.plot(data[:,0],data[:,1],'r--',linewidth=2,label='Kinetic')
#pl.plot(dat[:,0],dat[:,1],'k-',linewidth=2)
#pl.yscale('log')
pl.legend(loc=3)
#pl.title()

pl.savefig('energy.pdf')
pl.show() 

