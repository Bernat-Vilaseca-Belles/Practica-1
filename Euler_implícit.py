import numpy as np
import matplotlib.pyplot as plt
import solucio_analitica as sa
import constants as const
import math
from gauss_seidel import gauss_seidel as gs


#Càlcul de les variables necessàries.:
tau = (36.5+273.15)/const.T_0 #Temperatura inicial normalitzada. tau=T(t=0)/T_0, T_0= k/sigma*(V_0**2).

#Establim les dades de la discratització:
N=101
delta_z = 1/(N-1) #valor de la discretització espacial normalitzada. delta_z'= h/(L_0*N). L_0=h ==> delta_z'=1/N.
delta_t= ((delta_z**2), 0.5*(delta_z**2)) #Conjunt de valors de les diferents discretitzacions temporals.
gamma= [1,0.55] #coficient
t_a=0.025
p=np.ones(2)
for i in range(2): p[i]=math.floor(t_a/delta_t[i]) #Límit fins a on es realitza la simulació.

foxtrot=np.ones(2)
for o in range(2):
    for j in range(int(p[o])):
        phi=tau*np.ones(N)
        v=phi[1:N-1]
        charlie=np.ones(N)
        for k in range(N): charlie[k]=phi[k]
        u=gs(gamma[j], delta_z, v, N, 0.0001)
        for k in range(N-2): v[k]=u[k]
    foxtrot[o]=np.append(phi[0], v, phi[-1])


x=np.linspace(0, const.L_0, N) #Definim l'eix x de la gràfica

Temp_1, Temp_2= foxtrot
Temp_real=np.ones(N) #Llista per calcular la temperatura a partir de la funció analítica.
for i in range(N):
    Temp_real[i]=sa.sol_analít(x[i]/const.L_0,t_a, 0.00000000001)-273.15 #Temperatura a cada punt del mallat segons la solució analítica

E_1, E_2=0, 0 #Errors de cada solució
for i in range(N):
    E_1+=abs(Temp_1[i]- Temp_real[i]) #Es calcula la suma dels errors de cada punt del mallat de la solució 1.
    E_2+=abs(Temp_2[i]- Temp_real[i]) #Suma errors solució 2.
E_1=E_1/N #Es normalitza l'error 1
E_2=E_2/N #Normalització de l'error 2

