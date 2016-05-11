##!/usr/bin/python

import numpy as np
import pylab as plt
import seaborn as sns 

sns.set_context('poster',font_scale=1.5)
sns.set_style({'font.family':'Times New Roman'})

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
plt.subplot(1,1,1)
data = np.genfromtxt(fname='phase.dat') 
#data = np.loadtxt('traj.dat')
nb = int(data.shape[-1]/2)
for x in range(0,nb):
    plt.plot(data[:,x],data[:,x+nb])

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
plt.ylabel('Momentum')
plt.xlabel('Position')
#plt.title('traj')
#plt.xlim(-16,16) 
#plt.ylim(-16,16) 

plt.savefig('phase.pdf')
plt.show() 

