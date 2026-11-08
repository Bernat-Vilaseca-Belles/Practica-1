import numpy as np
import matplotlib as plt
from constants import T_c, sigma, V_0, k, t_a
N=101
z_tilde=np.linspace (0,1,N) #creem els punst on avaluarem la funció temperatura
n=500 #número fins al que trunquem el sumatori
beta=(T_c*sigma*V_0**2)/(k)
def solucio_analitica(z,t,n):
    T=T_c*np.ones_like(z)
    for i in range(n):
        T+=(k/(sigma*V_0**2))*(4*(1-np.exp(-(2*n+1)**2*np.pi**2*t)))*(np.sin((2*n+1)*np.pi*z))/(np.pi**3*(2*n+1)**3)
    return T    
T_analitica=solucio_analitica(z_tilde,t_a,n)