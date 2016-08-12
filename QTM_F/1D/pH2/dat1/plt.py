import numpy as np 
import matplotlib.pyplot as plt 

def plt_pes():
 
	data = np.genfromtxt('pes.out')

	plt.figure()  
	plt.plot(data[:,0], data[:,1],'b-',lw=3) 

	plt.ylim(-0.001,0.01)
	plt.xlim(2,10) 

def plt_en(ax1,ax2):

	data = np.genfromtxt(fname='energy.out')

	ax1.plot(data[:,0],data[:,2],linewidth=2,label='Potential')
	ax1.plot(data[:,0],data[:,3],linewidth=2,label='Quantum Potential')
	ax1.plot(data[:,0],data[:,4],linewidth=2,label='Energy')
	#pl.legend(bbox_to_anchor=(0.5, 0.38, 0.42, .302), loc=3,ncol=1, mode="expand", borderaxespad=0.)
	
	
	ax2.plot(data[:,0],data[:,1],linewidth=2,label='K')


#plt.legend()
	ax2.set_ylabel('Energy [hartree]')


#plt.savefig('energy.pdf')
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1) 
	
plt_en(ax1, ax2) 
#plt_pes()
plt.show() 


