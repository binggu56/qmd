import numpy as np
import pylab as plt
import seaborn as sns 

# make the fig pretty 

sns.set_context("poster",font_scale=1.5)
sns.set_style({'font.family':'Times New Roman'})


from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

minorLocator = MultipleLocator(1) 

en,c1,c2,c3,c4 = np.loadtxt('fft.dat',unpack=True)

# make plots 

fig, ax  = plt.subplots(1,1) 
#ax.xaxis.set_minor_locator(minorLocator)

ax.plot(en,c1,linewidth=2,label='$\mathcal{F}[C_{11}(t)](\omega)$')
ax.plot(en,c2,linewidth=2,label='$\mathcal{F}[C_{22}(t)](\omega)$')
ax.plot(en,c3,linewidth=2,label='$\mathcal{F}[C_{33}(t)](\omega)$')
ax.plot(en,c4,linewidth=2,label='$\mathcal{F}[C_{44}(t)](\omega)$')
#plt.plot(dat[:,0],dat[:,2],linewidth=2,label='$\widetilde{X^2(t)}(\omega)$')
#plt.plot(dat[:,0],dat[:,5],linewidth=2,label='$\widetilde{X(t)}(\omega)$')
#plt.plot(dat[:,0],dat[:,6],'--',linewidth=2,label='$\widetilde{XX(t)}(\omega)$')
plt.legend()
#plt.xlim(0,3)
plt.xlabel('$\omega$')
plt.savefig('fft.pdf')
plt.show()
