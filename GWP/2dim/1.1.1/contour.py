
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
#import seaborn as sns 
#sns.set_context("paper",font_scale=1.5)
#sns.set_style({'font.family':'Times New Roman'})

#matplotlib.rcParams['xtick.direction'] = 'out'
#matplotlib.rcParams['ytick.direction'] = 'out'

matplotlib.rcParams.update({'font.size': 20})
font = {'family' : 'Times New Roman', 'weight' : 'normal', 'size'   : 18}                
matplotlib.rc('font', **font)

delta = 0.02
xmin = -4.0 
xmax = 4.0 
ymin = -3.0 
ymax = 3.0 
X = np.arange(xmin, xmax, delta)
Y = np.arange(ymin, ymax, delta)
x, y = np.meshgrid(X, Y)
#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
a = 16 
#Z1 = np.exp(-a*(X-1)**2 -a*(Y+1)**2) 
#Z2 = np.exp(-a*(X+1)**2 -a*(Y-1)**2) 
#z = (1.0-(x-np.sqrt(2.0)*y)**2/3.0)**2/8.0 + (np.sqrt(2.0)*x-y)**2/6.0-1.0/8.0 

a = 1.0 
b = 4.0 
c = 4.0
z = y**2*(a*y**2-b)+c*(x-y)**2/2.0+b**2/4.0/a 
cmap = cm.get_cmap('ocean')

# Create a simple contour plot with labels using default colors.  The
# inline argument to clabel will control whether the labels are draw
# over the line segments of the contour, removing the lines beneath
# the label
plt.figure()
levels = [0.0, 1, 2, 4.0, 8.0,12.0,16.0,20]
CS = plt.contour(x, y, z, levels,cmap=cmap)
#sns.kdeplot(z)
#plt.contour(X, Y, Z2,cmap=cmap)
plt.clabel(CS, inline=1, fontsize=9)
#cb = plt.colorbar(CS)
#cb.set_label('')

#plt.title('Simplest default with labels')
dat = np.genfromtxt(fname='xyt.dat')
plt.plot(dat[:,0],dat[:,1],'ko',markersize=4)
plt.xlabel('x [Bohr]')
plt.ylabel('y [Bohr]')
plt.xlim(-4,4)
plt.ylim(-3.6,3)

plt.savefig('contour.pdf')
plt.show()
