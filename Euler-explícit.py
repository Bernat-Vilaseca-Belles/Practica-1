import numpy as np
import math
import constants as const
import matplotlib.pyplot as plt

# fitxer = open("Resultats ER.txt", "w") #Fitxer on es guardaran els resultats de les simulacions.

#En aquest codi s'ha calculat respecte a les funcions normalitzades. A raó d'això, cal tenir en compte els valors dels coeficients de normalització. En aquest cas, s'ha utilitzat la funció de normalització de la distribució normal, que és 1/(sigma*sqrt(2*pi)), on sigma és la desviació estàndard. Això significa que els valors dels coeficients de normalització i els valors de les constants.

#Càlcul de les variables necessàries.:
tau = (36.5+273)/const.T_0 #Temperatura inicial normalitzada. tau=T(t=0)/T_0, T_0= k/sigma*(V_0**2).

#Comencem establint les dades de la discratització:
N=101
delta_z = 1/(N-1) #valor de la discretització espacial normalitzada. delta_z'= h/(L_0*N). L_0=h ==> delta_z'=1/N.
delta_t= (0.51*(delta_z**2), 0.49*(delta_z**2), 0.25*(delta_z**2)) #Conjunt de valors de les diferents discretitzacions temporals.
gamma= [0.51,0.49,0.25] #coficient
t_a=0.025
p=np.ones(3)
for i in range(3): p[i]=math.floor(t_a/delta_t[i]) #Límit fins a on es realitza la simulació.


#Elements per imprimir els resultats:
foxtrot = np.ones(3) #Llista de temps que ha aguantat el teixit en les condicions correctes per a cada discretització. Nom escollit a l'atzar.

def calcul_temperatura(T0, dz, N, p, g):
    phi = T0*np.ones(N)
    charlie= T0*np.ones(N)
    for i in range(int(p)):
        for j in range(1, N-1):
            charlie[j] = phi[j] + g*(phi[j+1]- (2*phi[j])+ phi[j-1] + (dz**2))
        for j in range(1,N-1): phi[j] =charlie[j] 
    T=np.ones(N)
    for j in range(N):
        T[j]=phi[j]*const.T_0-273
    return T

Temp_1=calcul_temperatura(tau, delta_z, N, p[0], gamma[0]) #Temperatura per a la primera discretització
Temp_2=calcul_temperatura(tau, delta_z, N, p[1], gamma[1]) #Temperatura per a la segona discretització
Temp_3=calcul_temperatura(tau, delta_z, N, p[2], gamma[2]) #Temperatura per a la tercera discretització

x=np.ones(N)

for i in range(N):
    x[i]= i*delta_z*const.L_0

print(len(x))
print(len(Temp_1))
plt.figure(figsize=(8,6))
plt.errorbar(x=x, y=Temp_1, xerrr=None, yerr=None, fmt='o', color="purple", markersize=2, ecolor="gray", elinewidth=0.8, capsize=1.5, label=r"$\Delta \hat t=0.51 \Delta (\hat z)^2$")
plt.errorbar(x=x, y=Temp_2, xerrr=None, yerr=None, fmt='o', color="red", markersize=2, ecolor="gray", elinewidth=0.8, capsize=1.5, label=r"$\Delta \hat t=0.49 \Delta (\hat z)^2$")
plt.errorbar(x=x, y=Temp_3, xerrr=None, yerr=None, fmt='o', color="blue", markersize=2, ecolor="gray", elinewidth=0.8, capsize=1.5, label=r"$\Delta \hat t=0.25 \Delta (\hat z)^2$")
plt.xlabel(r"$z$ (m)")    #(o $x$, depenent de com es vulgui anomenar)
plt.ylabel(r"$T$ (ºC)")
plt.title(r"Temperatures a $t_a=0.025$ amb Euler explícit.")
plt.legend(loc="best")
plt.show()