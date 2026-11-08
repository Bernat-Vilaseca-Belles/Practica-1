import numpy as np
import matplotlib as plt
from constants import T_c, sigma, V_0, k, t_a
N=101
z_tilde=np.linspace (0,1,N) #creem els punst on avaluarem la funció temperatura
beta=(T_c*sigma*V_0**2)/(k)
epsilon=0.0000000001 #definim un epslon per assegurar-nos que per un n gran, sumem mensy que epsilon
n_max=1000
def solucio_analitica(z,t,epsilon,n_max):
    T=T_c*np.ones_like(z)
    n=0
    menor_que_epsilon=False
    while not menor_que_epsilon and n<n_max:
        T+=(k/(sigma*V_0**2))*(4*(1-np.exp(-(2*n+1)**2*np.pi**2*t)))*(np.sin((2*n+1)*np.pi*z))/(np.pi**3*(2*n+1)**3)
        if np.max(np.abs((k/(sigma*V_0**2))*(4*(1-np.exp(-(2*n+1)**2*np.pi**2*t)))*(np.sin((2*n+1)*np.pi*z))/(np.pi**3*(2*n+1)**3)))<epsilon:
            menor_que_epsilon=True
            print(f"La sèrie ha convergit a {n+1} termes (per a epsilon = {epsilon}).")
        n+=1
    return T    

T_analitica=solucio_analitica(z_tilde,t_a,epsilon,n_max)