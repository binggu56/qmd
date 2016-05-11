import numpy as np
import pylab as plt

c1 = 16.
c2 = 2e6
d = 0.16
c4 = 0.00001

x = np.linspace(0.00001,0.001,200)
y1 = d/np.exp(c1*x)/x
y2 = 3.0*d/x
y3 = np.exp(-x)/np.tan(x)/3
y4 = np.exp(-c4*x)/np.tan(x)/3

plt.plot(x,y1,'r-',linewidth=2)
plt.plot(x,y2,'r--',linewidth=2)
plt.plot(x,y3,'k--',linewidth=2)
plt.plot(x,y4,'g--',linewidth=2)

plt.savefig('force.pdf')
plt.show()
