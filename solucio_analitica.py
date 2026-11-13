import numpy as np
import matplotlib as plt
from constants import T_c, sigma, V_0, k, t_a, T_0
# N=101
# z_tilde=np.linspace (0,1,N) #creem els punst on avaluarem la funció temperatura
# beta=(2*k)/(T_c*sigma*V_0**2)
# epsilon=0.0001 #definim un epslon per assegurar-nos que per un n gran, sumem menys que epsilon
# n_max=1000
# def solucio_analitica(z,t,epsilon,n_max):
#     T=T_c*np.ones_like(z)
#     n=0
#     menor_que_epsilon=False
#     while not menor_que_epsilon and n<n_max:
#         T+=((sigma*V_0**2)/(2*k))*(4*(1-np.exp(-(2*n+1)**2*np.pi**2*t)))*(np.sin((2*n+1)*np.pi*z))/(np.pi**3*(2*n+1)**3)
#         if np.max(np.abs(((sigma*V_0**2)/(2*k))*(4*(1-np.exp(-(2*n+1)**2*np.pi**2*t)))*(np.sin((2*n+1)*np.pi*z))/(np.pi**3*(2*n+1)**3)))<epsilon:
#             menor_que_epsilon=True
#             print(f"La sèrie ha convergit a {n+1} termes (per a epsilon = {epsilon}).")
#         n+=1
#     return T    

# T_analitica=solucio_analitica(z_tilde,t_a,epsilon,n_max)

def sol_analít(dz, N, t, error):
    vector=T_c*np.ones(N)
    E=np.ones(N)
    n=0
    b_1=True #No convergent
    while b_1:
        for i in range(N):
            vector[i]+=T_0*(4*((1-np.exp(-((2*n+1)**2)*(np.pi**2)*t))
                                        /(((2*n+1)**3)*(np.pi**3)))*np.sin((2*n+1)*np.pi*(i*dz)))
            if n>1:
                E[i]=abs(T_0*(4*((1-np.exp(-((2*n+1)**2)*(np.pi**2)*t))
                                 /(((2*n+1)**3)*(np.pi**3)))*np.sin((2*n+1)*np.pi*(i*dz))))/vector[i]
        if max(E)<=error:
            b_1=False
        n+=1
    return vector