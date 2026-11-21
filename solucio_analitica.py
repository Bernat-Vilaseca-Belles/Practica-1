import numpy as np
import matplotlib as plt
from constants import T_c, sigma, V_0, k, t_a, T_0

def sol_analít(dz, N, t, error):
    vector=T_c*np.ones(N)
    E=np.ones(N)
    n=0
    b_1=True #No convergent
    while b_1:
        for i in range(N):
            vector[i]+=T_0*4*((1-np.exp(-((2*n+1)**2)*(np.pi**2)*t))
                                        /(((2*n+1)**3)*(np.pi**3)))*np.sin((2*n+1)*np.pi*(i*dz))
            if n>1:
                E[i]=abs(((1-np.exp(-((2*n+1)**2)*(np.pi**2)*t))
                                 /(((2*n+1)**3)*(np.pi**3)))*np.sin(((2*n)+1)*np.pi*(i*dz)))/vector[i]
        if max(E)<=error:
            b_1=False
        n+=1
    return vector