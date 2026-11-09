import numpy as np
import math
import constants as const
import matplotlib.pyplot as plt
import solucio_analitica as sa
"""
En aquest codi s'ha calculat respecte a les funcions normalitzades. 
A raó d'això, cal tenir en compte els valors dels coeficients de normalització.
"""

#Establim una funció que ens doni la temperatura al llarg de tota la discretització de z a partir de la solució de l'equació diferencial i fins a un temps determinat.
def calcul_temperatura(T0, dz, N, p, g): #Temps inicil, interval entre dos punts de z, nombre de punts, temps màxim i gamma (que es dedueix en el document de calculs d'euler explícit.)
    phi = T0*np.ones(N) #Temperatura normalitzada a cada iteració.
    charlie= T0*np.ones(N) #Eina per trobar la temperatura a la partició de temps n+1. Nom escollit a l'atzar.
    for i in range(int(p)):
        for j in range(1, N-1):
            charlie[j] = phi[j] + g*(phi[j+1]- (2*phi[j])+ phi[j-1] + (dz**2)) #Es calcula la temperatura al temps n+1
        for j in range(1,N-1): phi[j] =charlie[j]  #S'assigna la nova temperatura a la llista phi.
    T=np.ones(N) #Llista que serà la temperatura resultant en ºC
    for j in range(N):
        T[j]=phi[j]*const.T_0-273.15 #Temperatura resultant en cada punt
    return T

#Càlcul de les variables necessàries.:
tau = (36.5+273.15)/const.T_0 #Temperatura inicial normalitzada. tau=T(t=0)/T_0, T_0= k/sigma*(V_0**2).

#Establim les dades de la discratització:
N=101
delta_z = 1/(N-1) #valor de la discretització espacial normalitzada. delta_z'= h/(L_0*N). L_0=h ==> delta_z'=1/N.
delta_t= (0.51*(delta_z**2), 0.49*(delta_z**2), 0.25*(delta_z**2)) #Conjunt de valors de les diferents discretitzacions temporals.
gamma= [0.51,0.49,0.25] #coficient
t_a=0.025
p=np.ones(3)
for i in range(3): p[i]=math.floor(t_a/delta_t[i]) #Límit fins a on es realitza la simulació.


#Apliquem el càlcul de la temperatura a les posicions de la discretització de z normalitzada
Temp_1=calcul_temperatura(tau, delta_z, N, p[0], gamma[0]) #Temperatura per a la primera discretització
Temp_2=calcul_temperatura(tau, delta_z, N, p[1], gamma[1]) #Temperatura per a la segona discretització
Temp_3=calcul_temperatura(tau, delta_z, N, p[2], gamma[2]) #Temperatura per a la tercera discretització

x=np.linspace(0, const.L_0, N)

Temp_real=np.ones(N) #Llista per calcular la temperatura a partir de la funció analítica.
for i in range(N):
    Temp_real[i]=sa.sol_analít(x[i]/const.L_0,t_a, 0.00000000001)-273.15 #Temperatura a cada punt del mallat segons la solució analítica

E_1, E_2, E_3=0, 0, 0 #Errors de cada solució
for i in range(N):
    E_1+=abs(Temp_1[i]- Temp_real[i]) #Es calcula la suma dels errors de cada punt del mallat de la solució 1.
    E_2+=abs(Temp_2[i]- Temp_real[i]) #Suma errors solució 2.
    E_3+=abs(Temp_3[i]- Temp_real[i]) #Suma errors solució 3.
E_1=E_1/N #Es normalitza l'error 1
E_2=E_2/N #Normalització de l'error 2
E_3=E_3/N #Normalització de l'error 3

plt.figure(figsize=(8,6))
plt.scatter(x=x, y=Temp_1, marker='o', color="purple", s=2, label=fr"$\Delta \hat t=0.51 (\Delta \hat z)^2$. Error={E_1}")
plt.scatter(x=x, y=Temp_2, marker='o', color="red", s=2, label=fr"$\Delta \hat t=0.49 (\Delta \hat z)^2$. Error={E_2}")
plt.scatter(x=x, y=Temp_3, marker='o', color="blue", s=2, label=fr"$\Delta \hat t=0.25 (\Delta \hat z)^2$. Error={E_3}")
plt.scatter(x, Temp_real, color='k', s=2, label="Funció real")
plt.xlabel(r"$z$ (m)")    #(o $x$, depenent de com es vulgui anomenar)
plt.ylabel(r"$T$ (ºC)")
plt.title(r"Temperatures a $t_a=0.025$ amb Euler explícit.")
plt.legend(loc="best")
plt.yscale('log')
plt.show()