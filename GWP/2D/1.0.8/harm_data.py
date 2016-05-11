##!/usr/bin/python

import numpy as np
import pylab as plt 
import subprocess 
import seaborn as sns 

#sns.set_context('poster')

plt.subplot(1,1,1)
t,dat0,dat1,dat2,dat3,dat4,dat5 = np.genfromtxt(fname='cor.dat',usecols=(0,1,2,3,4,5,6),unpack=True) 

#for x in range(1,data.shape[-1]):
    
#plt.plot(data[:,0],data[:,1],lw=2,label='Re')
#plt.plot(data[:,0],data[:,2],lw=2,label='Im')
#plt.plot(data[:,0],data[:,3],lw=2,label='$|C(t)|$')

#dat = np.genfromtxt(fname='../spo/1.0.2/corr')
#plt.plot(dat[:,0],dat[:,1],'--',label='Re, QM')
#plt.plot(dat[:,0],dat[:,2],'--',label='Im, QM')

#x = np.linspace(0,4,100) 
#y = -np.sin(x)
#plt.plot(x,y,lw=2,label='sin(x)')
#plt.xlabel('$Time$')
#plt.ylabel('$C(t)$')
#dat0 = ' '.join(map(str, corr.tolist()))
#dat1 = ' '.join(map(str, cori.tolist()))
dat = ''
str1 = ''
str2 = '' 
for i in range(0,len(dat0)):
    dat = dat + str(dat0[i])+'+'+str(dat1[i])+'i ' 
#    dat = dat+str(dat0[i])+' '  

for i in range(0,len(dat0)):
    str1 = str1 + str(dat2[i])+'+'+str(dat3[i])+'i ' 

for i in range(0,len(dat0)):
    str2 = str2 + str(dat4[i])+'+'+str(dat5[i])+'i ' 

#plt.subplot(2,1,2)
#dat = np.genfromtxt(fname='/home/bing/gwp/spo_2d/1.0.0/cor1') 
f = open('harm_cor.dat', 'w')
f.write(str(dat))

f1 = open('harm_cor_sys.dat', 'w')
f1.write(str(str1))

f2 = open('harm_cor_asym.dat', 'w')
f2.write(str(str2))

dt = t[1]-t[0]
#cmd = 'harminv -t'+str(dt)+'0-3 < harm_cor.dat'
print('harminv -t 0.002 0-1 < harm_cor.dat')
#subprocess.Popen("harminv -t 0.004 0-1 < harm_cor.dat",shell=True)

print('harminv -t 0.002 0-1 < harm_cor_sys.dat')
#subprocess.Popen("harminv -t 0.004 0-1 < harm_cor_sys.dat",shell=True)

print('harminv -t 0.002 0-1 < harm_cor_asym.dat')
#subprocess.Popen("harminv -t 0.004 0-1 < harm_cor_asym.dat",shell=True)

#f2 = open('harm_cori.dat','w')
#f2.write(str(dat1))
#np.savetxt('harm_cor.dat',str(dat),fmt="%s")
#for x in range(1,3):
#plt.plot(dat[:,0],dat[:,1],'--',label='$\Re(C(t))$',lw=2)
#plt.plot(dat[:,0],dat[:,2],'--',label='$\Im(C(t))$',lw=2)
#z = np.sqrt(data[:,1]**2+data[:,2]**2)
#plt.plot(data[:,0],z,label='$|C(t)|$',lw=1)
#plt.ylim(-0.2,0.2) 
#plt.legend()
#plt.xlim(0,36)
#plt.savefig('cor.pdf')
#plt.show() 

