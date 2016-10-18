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
data = np.genfromtxt(fname='xy.dat')
#dat = np.genfromtxt(fname='../2d12/en.dat')
#data = np.loadtxt('traj.dat')
#for x in range(1,10):
pl.subplot(111)
#pl.xlim(0,4)
pl.title('two-steps fitting alg')
pl.ylabel('Energy [hartree]')
pl.plot(data[:,1],data[:,0],'o',label='')
#pl.legend(bbox_to_anchor=(0.5, 0.38, 0.42, .302), loc=3,ncol=1, mode="expand", borderaxespad=0.)

#pl.subplot(212)
#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#pl.xlim(0,4)
#pl.xlabel('time [a.u.]')
#pl.ylabel('Energy [hartree]')
#pl.plot(dat[:,0],dat[:,1],'k-',linewidth=2)
#pl.yscale('log')
#pl.legend()
#pl.title()

#pl.savefig('energy.pdf')
pl.show() 
