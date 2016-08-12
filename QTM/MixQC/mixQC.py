# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 16:07:54 2016

@author: bing
"""

import numpy as np
import pylab as pl
import matplotlib as mpl 

pl.style.use('ggplot')

font = {'family' : 'Times New Roman',
        'weight' : 'normal', 
        'size'   : '18'}

mpl.rc('font', **font)  # pass in the font dict as kwargs
#import seaborn as sns 
#sns.set_context('poster')
#sns.set_style("whitegrid")

data = np.genfromtxt(fname='/home/bing/qt/1.0.7/den.dat')
#dat = np.genfromtxt(fname='../2d12/en.dat')
#data = np.loadtxt('traj.dat')
#for x in range(1,10):
pl.figure(figsize=(14,9))
pl.subplot(111)

#pl.ylabel('Energy [hartree]')
pl.plot(data[:,0],data[:,1],'b-.',linewidth=2,label='L=8,$\partial_y U =0 $ ')
pl.plot(data[:,0],data[:,2],'b-.',linewidth=2)

dat = np.genfromtxt(fname='/home/bing/qt/spo_2d/1.0.3/den.dat') 
pl.plot(dat[:,0],dat[:,1],'k',linewidth=3,label='QM w/o $U_c$')
pl.plot(dat[:,0],dat[:,2],'k',linewidth=3)
#pl.plot(data[:,0],data[:,4],'k-',linewidth=2,label='Energy')
#pl.legend(bbox_to_anchor=(0.5, 0.38, 0.42, .302), loc=3,ncol=1, mode="expand", borderaxespad=0.)

#dat = np.genfromtxt(fname='dom2/den.dat') 
#pl.plot(dat[:,0],dat[:,1],'c',linewidth=2,label='L = 2 ')
#pl.plot(dat[:,0],dat[:,2],'c',linewidth=2)
#
#dat = np.genfromtxt(fname='dom16/den.dat') 
#pl.plot(dat[:,0],dat[:,1],'r--',linewidth=3,label='L = 16 ')
#pl.plot(dat[:,0],dat[:,2],'r--',linewidth=3)
#
#
dat1 = np.genfromtxt(fname='/home/bing/qt/1.0.6/den.dat') 
pl.plot(dat1[:,0],dat1[:,1],'r--',linewidth=2,label='L = 8, $\partial_y \Lambda =0 $' )
pl.plot(dat1[:,0],dat1[:,2],'r--',linewidth=2)

dat2 = np.genfromtxt(fname='/home/bing/qt/1.0.5/dom8/den.dat') 
pl.plot(dat2[:,0],dat2[:,1],'g',linewidth=2,label='L = 8')
pl.plot(dat2[:,0],dat2[:,2],'g',linewidth=2)

#pl.subplot(212)
#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
pl.xlim(0,4)
#pl.xlabel('time [a.u.]')
##pl.ylabel('Energy [hartree#')
#pl.plot(data[:,0],data[:,1]#'r--',linewidth=2,label='Kinetic')
#pl.plot(dat[:,0],dat[:,1],'k-',linewidth=2)
#pl.yscale('log')
pl.legend(loc=2)
pl.xlabel('Time [a.u.]')
pl.ylabel('Density overlap')
#pl.title()

pl.savefig('den.pdf')
pl.title('QC dynamics with 8 domains')
pl.show() 