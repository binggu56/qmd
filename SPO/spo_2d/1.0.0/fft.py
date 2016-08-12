import numpy as np

PI = 3.1415926e0 
IM = complex(0e0,1e0)
Emin = 0.5e0
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
        out[i] = data[j]*np.sin(en*t)*np.exp(-alfa*t**2)*dt+out[i] 
        
  for i in range(nk): 
    out[i] = out[i]/np.pi  
  
  return out 

# data analysis 
time, cxx, cyy, cxy = np.loadtxt('cor.dat',delimiter=',',usecols=(0,2,4,6),unpack=True)

nt = len(time)
dt = time[1]-time[0]
alfa = -np.log(10e-4)/time[-1]**2
print(alfa)
out0 = []
out1 = []
out2 = []
out0 = fft_sin(nt,dt,cxx,alfa) 
out1 = fft_sin(nt,dt,cyy,alfa)
out2 = fft_sin(nt,dt,cxy,alfa) 

en = [0e0 for i in range(nk)]
for i in range(nk):
  en[i] = Emin + (i-1)*de

DataOut = np.column_stack((en,out0,out1,out2)) 
np.savetxt('fft.dat', DataOut, fmt=('%14.7e'))

