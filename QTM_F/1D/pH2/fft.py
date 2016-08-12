import numpy as np

IM = complex(0e0,1e0)
Emin = 0e0
Emax = 3e0 
nk = 400 
de = (Emax-Emin)/(nk-1) 

def fft_sin(t,data,alfa=0e0):  
  """
  fourier transform with sin(t)
  INPUT: 
    Nt: time steps 
    dt: time interval 
    alfa: damping coefficient 

  OUTPUT:
    out: array(Nt), freqency 
  
  """ 
  Nt = len(t) 
  out = [0e0 for i in range(nk)]  
  dt = t[1]-t[0]
  
  for i in range(0,nk):
      en = Emin + i*de 
      for j in range(0,Nt): 
        out[i] = 2e0*data[j]*np.sin(en*t[j])*np.exp(-alfa*t[j]**2)*dt+out[i] 
        
  for i in range(nk): 
    out[i] = out[i]/(2.0*np.pi)  
  
  return out 

# data analysis 
t, c1, c2, c3, c4 = np.loadtxt('cor.dat',usecols=(0,2,4,6,8),unpack=True)

Nt = len(t)
alfa = -np.log(1e-2)/t[-1]**2
print('scaling factor = ',alfa)

out1 = []
out2 = []
out3 = []
out4 = []

out1 = fft_sin(t,c1,alfa) 
out2 = fft_sin(t,c2,alfa) 
out3 = fft_sin(t,c3,alfa) 
out4 = fft_sin(t,c4,alfa) 

en = [0e0 for i in range(nk)]

for i in range(nk):
  en[i] = Emin + (i-1)*de

DataOut = np.column_stack((en,out1,out2,out3,out4)) 
np.savetxt('fft.dat', DataOut, fmt=('%14.7e'))

