import numpy as np
import cmath 
import pylab as plt 

Emin = 0.01e0
Emax = 4e0 
nk = 400 # number of points of FT  
de = (Emax-Emin)/(nk-1) 

def fft_sin(nt,dt,data,alfa=0e0):  
  """
  fourier transform with sin(t)
  INPUT: 
    nt: time steps 
    dt: time interval 
    alfa: damping coefficient 

  OUTPUT:
    out: array(Nt), freqency 
  
  """ 
  out = [0e0 for i in range(nk)]
  
  for i in range(0,nk):
      en = Emin + i*de 
      t = 0e0 
      for j in range(0,nt): 
        t = t + dt   
        out[i] = (data[j]*np.exp(-1j*en*t)).real*np.exp(-alfa*t**2)*dt+out[i] 
        
  for i in range(nk): 
    out[i] = out[i]/np.pi  
  
  return out 

# data analysis 
time, real, imag = np.loadtxt('cor.dat',usecols=(0,1,2),unpack=True)

nt = len(time)
dt = time[1]-time[0]
alfa = -np.log(10e-4)/time[-1]**2
print(alfa)
cor = real+1j*imag 
out0 = []
out0 = fft_sin(nt,dt,cor,alfa) 
#out1 = fft_sin(nt,dt,cyy,alfa)
#out2 = fft_sin(nt,dt,cxy,alfa) 

en = [0e0 for i in range(nk)]
for i in range(nk):
  en[i] = Emin + (i-1)*de


# save data 
DataOut = np.column_stack((en,out0)) 
np.savetxt('fft.dat', DataOut, fmt=('%14.7e'))


# make figure 
plt.subplot(111)
plt.plot(en,out0,lw=2)

plt.savefig('fft.pdf')
plt.show()

