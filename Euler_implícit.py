import numpy as np
import matplotlib.pyplot as plt
import solucio_analitica as sa
import constants as const
import math
from gauss_seidel import gauss_seidel as gs

#Definició del vector sortida (S'aplica més endevant)
def calcul_b(vector, dz, g, tau):
    b=np.ones(N-2) #Definició del vector corresponent de Ax=b
    for i in range(N-2): b[i]=g*(dz**2)+vector[i] #Assignem els valors de b
    b[0]+=g*tau #Com que és el valor de j=2 i no està definida la operació per a j=1, cal afegir-hi el valor esctremal
    b[N-3]+=g*tau #Afegim el valor extremal de la temperatura per a j=N
    return b

#Càlcul de les variables necessàries.:
tau = (36.5+273.15)/const.T_0 #Temperatura inicial normalitzada. tau=T(t=0)/T_0, T_0= k/sigma*(V_0**2).

#Establim les dades de la discratització:
N=101
delta_z = 1/(N-1) #valor de la discretització espacial normalitzada. delta_z'= h/(L_0*N). L_0=h ==> delta_z'=1/N.
delta_t= ((delta_z**2), 0.5*(delta_z**2)) #Conjunt de valors de les diferents discretitzacions temporals.
gamma= [1,0.5] #coficient
t_a=0.025
p=np.ones(2)
for i in range(2): p[i]=math.floor(t_a/delta_t[i]) #Límit fins a on es realitza la simulació.

#Definim la Matriu A
A=[]
for k in range(len(gamma)):
    Bravo=[np.zeros(N-2) for _ in range(N-2)] #Definició de la matriu.
    for i in range(N-2):
        for j in range(N-2):
            if j==i:
                Bravo[i][j]=1+(2*gamma[k]) #Definim els valors de A a la diagonal.
            if ((j-1)==i) or ((j+1)==i):
                Bravo[i][j]=-gamma[k] #Definim els valors de A al voltant de la diagonal. Els altres valors són 0
    A.append(Bravo)







foxtrot=[[],[]] #Llista formada per dues llistes que contindran les temperatures en les diferents posicions. Nom aleatori.
for o in range(2):
    phi=tau*np.ones(N) # Assignem la temperatura normalitzada com a llista.
    v=phi[1:N-1] #Agafem tots els valors de T que no siguin els valors 
    for j in range(int(p[o])):
        u, total=gs(A[o], gamma[o], delta_z, v, calcul_b(v, delta_z, gamma[o], tau), tau, N, 0.001) #Apliquem el mètode de gauss-seidel
        for k in range(N-2): v[k]=u[k] #Reescribim els valors de v.
    foxtrot[o]=np.append(np.array([phi[0]]), np.array(v)) #Posem el valor de la temperatura a z=0 al començament de la llista.
    foxtrot[o]=np.append(np.array(foxtrot[o]), np.array([phi[-1]])) #Posem el valor de la temperatua a z=h al final de la llista.
    for i in range(len(foxtrot[o])):
        foxtrot[o][i]=(foxtrot[o][i]*const.T_0)-273.15 #Canviem la llista a temperatures no normalitzades.

x=np.linspace(0, const.L_0, N) #Definim l'eix x de la gràfica

Temp_1, Temp_2= foxtrot #Assignem els valors de cada temperatura a una llista corresponent.
Temp_real=sa.sol_analít(delta_z, N,t_a, 0.00001) #Llista per calcular la temperatura a partir de la funció analítica.
for i in range(N):
    Temp_real[i]=Temp_real[i]-273.15 #Temperatura a cada punt del mallat segons la solució analítica

E_1, E_2=np.ones(N), np.ones(N) #Errors de cada solució
for i in range(N):
    E_1[i]=abs(Temp_1[i]- Temp_real[i])/Temp_real[i] #Es calcula la suma dels errors de cada punt del mallat de la solució 1.
    E_2[i]=abs(Temp_2[i]- Temp_real[i])/Temp_real[i] #Suma errors solució 2.
# E_1=E_1/N #Es normalitza l'error 1
# E_2=E_2/N #Normalització de l'error 2

#Gràfic.
plt.figure(figsize=(16,6))
plt.subplot(1,2,1)
plt.scatter(x, Temp_1, color='g', s=2, label=rf"Discretització $\Delta \tilde t = (\Delta \tilde z)^2$")
plt.scatter(x, Temp_2, color='b', s=2, label=rf"Discretització $\Delta \tilde t = 0.5 \, (\Delta \tilde z)^2$")
plt.scatter(x, Temp_real, color='k', s=3, label=rf"Funció real")
plt.xlabel(r"$\Delta z$ (m)")
plt.ylabel(r"$T$ ($^\circ $C)")
plt.title(r"Mètode d'Euler implícit.")
plt.legend(loc='best')

plt.subplot(1,2,2)
plt.scatter(x, E_1, color='g', s=2, label=rf"Error de $\Delta \tilde t = (\Delta \tilde z)^2$")
plt.scatter(x, E_2, color='b', s=2, label=rf"Error de $\Delta \tilde t = 0.5(\Delta \tilde z)^2$")
plt.xlabel(r"$\Delta z$ (m)")
plt.ylabel(r"Error normalitzat")
plt.title(r"Errors del mètode d'Euler implícit.")
plt.legend(loc='best')
plt.show()