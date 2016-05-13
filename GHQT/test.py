import numpy as np 

def Sym(V):
    n = V.shape[-1]

    for i in range(n):
        for j in range(i):
            V[i,j] = V[j,i]
    
    return V


a = np.zeros((2,2))

a[0,1] = 1 
a[1,1] = 3 
print(a) 
a = Sym(a) 
print(a) 
