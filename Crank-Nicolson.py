import numpy as np
from constants import T_c,T_0
import matplotlib.pyplot as plt
import solucio_analitica as sa
import constants as const
import math
from gauss_seidel import gauss_seidel as gs

ee = np.genfromtxt('Temperatura_ee.csv', delimiter=';', filling_values=np.nan)
ei = np.genfromtxt('Temperatura_ei.csv', delimiter=';', filling_values=np.nan)
ee_T1=ee[0,:]
ee_T2=ee[1,:]
ee_T3=ee[2,:]

ei_T1=ei[0,:]
ei_T2=ei[1,:]



#Definim la matriu A
def calcul_A(g):
    Bravo = [np.zeros(N-2) for _ in range(N-2)]
    for i in range(N-2):
        for j in range(N-2):
            if i==j:
                Bravo[i][i]= 2*(1+g)
            if i>0:
                Bravo[i][i-1]=-g
            if i<N-3:
                Bravo[i][i+1]=-g
    return Bravo

#Definim el vector b
def calcul_b(v, g, dt, tau):
    b=np.ones(N-2)
    deca=2*(1-g)
    golf=2*dt
    for i in range(1, N-3):
        b[i]=g*v[i-1]+deca*v[i]+ g*v[i+1] + golf
    b[0] = 2*g*tau+ deca*v[0] + g*v[1] + golf
    b[N-3] = g*v[N-4] + deca*v[N-3] + 2*g*tau + golf 
    return b


#Càlcul de les variables necessàries.:
tau = T_c/T_0 #Temperatura inicial normalitzada. tau=T(t=0)/T_0, T_0= k/sigma*(V_0**2).

#Establim les dades de la discratització:
N=101
delta_z = 1/(N-1) #valor de la discretització espacial normalitzada. delta_z'= h/(L_0*N). L_0=h ==> delta_z'=1/N.
delta_t= ((delta_z**2), 0.5*(delta_z**2)) #Conjunt de valors de les diferents discretitzacions temporals.
gamma= [1,0.5] #coficient
t_a=0.025
p=np.ones(2)
for i in range(2): p[i]=math.floor(t_a/delta_t[i]) #Límit fins a on es realitza la simulació.

A=[calcul_A(gamma[0]), calcul_A(gamma[1])]



foxtrot=[[],[]] #Llista formada per dues llistes que contindran les temperatures en les diferents posicions. Nom aleatori.
for o in range(2):
    phi=tau*np.ones(N) # Assignem la temperatura normalitzada com a llista.
    v=phi[1:N-1] #Agafem tots els valors de T que no siguin els valors 
    for j in range(int(p[o])):
        u, total=gs(A[o], gamma[o], v, calcul_b(v, gamma[o], delta_t[o], tau), N, 0.00001) #Apliquem el mètode de gauss-seidel
        for k in range(N-2): v[k]=u[k] #Reescribim els valors de v.
    foxtrot[o]=np.append(np.array([phi[0]]), np.array(v)) #Posem el valor de la temperatura a z=0 al començament de la llista.
    foxtrot[o]=np.append(np.array(foxtrot[o]), np.array([phi[-1]])) #Posem el valor de la temperatua a z=h al final de la llista.
    for i in range(len(foxtrot[o])):
        foxtrot[o][i]=(foxtrot[o][i]*const.T_0)-273.15 #Canviem la llista a temperatures no normalitzades.

x=np.linspace(0, const.L_0, N) #Definim l'eix x de la gràfica

Temp_1, Temp_2= foxtrot #Assignem els valors de cada temperatura a una llista corresponent.
Temp_real=sa.sol_analít(delta_z, N,t_a, 0.0001) #Llista per calcular la temperatura a partir de la funció analítica.
for i in range(N):
    Temp_real[i]=Temp_real[i]-273.15 #Temperatura a cada punt del mallat segons la solució analítica

E_1, E_2=np.ones(N), np.ones(N) #Errors de cada solució
for i in range(N):
    E_1[i]=abs(Temp_1[i]- Temp_real[i])/Temp_real[i] #Es calcula la suma dels errors de cada punt del mallat de la solució 1.
    E_2[i]=abs(Temp_2[i]- Temp_real[i])/Temp_real[i] #Suma errors solució 2.
# E_1=E_1/N #Es normalitza l'error 1
# E_2=E_2/N #Normalització de l'error 2

#Gràfic
plt.figure(figsize=(16,6))
plt.subplot(1,2,1)
plt.scatter(x, Temp_1, color='g', s=2, label=rf"Discretització $\Delta \tilde t = (\Delta \tilde z)^2$")
plt.scatter(x, Temp_2, color='b', s=2, label=rf"Discretització $\Delta \tilde t = 0.5 \, (\Delta \tilde z)^2$")
plt.scatter(x, Temp_real, color='k', s=3, label=rf"Funció real")
plt.xlabel(r"$\Delta z$ (m)")
plt.ylabel(r"$T$ ($^\circ $C)")
plt.title(r"Mètode de Crack-Nicolson")
plt.legend(loc='best')

plt.subplot(1,2,2)
plt.scatter(x, E_1, color='g', s=2, label=rf"Error de $\Delta \tilde t = (\Delta \tilde z)^2$")
plt.scatter(x, E_2, color='b', s=2, label=rf"Error de $\Delta \tilde t = 0.5(\Delta \tilde z)^2$")
plt.xlabel(r"$\Delta z$ (m)")
plt.ylabel(r"Error normalitzat")
plt.title(r"Errors del mètode de Crack-Nicolson")
plt.legend(loc='best')
plt.show()


"""
Ara farem el gràfic on comparem tots els altres mètodes respecte aquest mètode amb la
discretització de delta t=0.5 (delta z)**2
"""
###Euler explícit
Err_11, Err_12, Err_13=np.ones(N), np.ones(N), np.ones(N)
for i in range(N):
    Err_11[i]=abs(ee_T1[i]- Temp_2[i])/Temp_2[i]
    Err_12[i]=abs(ee_T2[i]- Temp_2[i])/Temp_2[i]
    Err_13[i]=abs(ee_T3[i]- Temp_2[i])/Temp_2[i]

##Euler implícit
Err_21, Err_22=np.ones(N), np.ones(N)
for i in range(N):
    Err_21[i]=abs(ei_T1[i]- Temp_2[i])/Temp_2[i]
    Err_22[i]=abs(ei_T2[i]- Temp_2[i])/Temp_2[i]

plt.figure(figsize=(16,6))
plt.subplot(1,2,1)
plt.scatter(x, ee_T1, marker='o', color='g', s=2, label=rf"Euler explícit amb $\Delta \tilde t = 0.51 \,(\Delta \tilde z)^2$")
plt.scatter(x, ee_T2, marker='o', color='b', s=2, label=rf"Euler explícit amb $\Delta \tilde t = 0.49 \, (\Delta \tilde z)^2$")
plt.scatter(x, ee_T3, marker='o', color="purple", s=2, label=fr"Euler explícit amb $\Delta \hat t=0.25 \,(\Delta \hat z)^2$")
plt.scatter(x, Temp_2, color='k', s=3, label=rf"Crack-Nickolson amb $\Delta \tilde t = 0.5 \, (\Delta \tilde z)^2$")
plt.xlabel(r"$z$ (m)")
plt.ylabel(r"$T$ ($^\circ $C)")
plt.legend(loc="best")
plt.yscale('log')
plt.title(r"Comparació Euler explícit amb Crack-Nicolson")
plt.legend(loc='best')

plt.subplot(1,2,2)
plt.scatter(x, Err_11, color='g', s=2, label=rf"Euler explícit amb $\Delta \tilde t = 0.51 \, (\Delta \tilde z)^2$")
plt.scatter(x, Err_12, color='b', s=2, label=rf"Euler explícit amb $\Delta \tilde t = 0.49 \, (\Delta \tilde z)^2$")
plt.scatter(x, Err_13, color='purple', s=2, label=rf"Euler explícit amb $\Delta \tilde t = 0.25 \,(\Delta \tilde z)^2$")
plt.xlabel(r"$\Delta z$ (m)")
plt.ylabel(r"Error normalitzat")
plt.title(r"Errors del mètode d'Euler explícit comparat amb Crack Nicolson")
plt.yscale('log')
plt.legend(loc='best')
plt.show()

##Gràfics Euler implícit

plt.figure(figsize=(16,6))
plt.subplot(1,2,1)
plt.scatter(x, ei_T1, marker='o', color='r', s=2, label=rf"Euler explícit amb $\Delta \tilde t =(\Delta \tilde z)^2$")
plt.scatter(x, ei_T2, marker='o', color='indigo', s=2, label=rf"Euler explícit amb $\Delta \tilde t = 0.5 \, (\Delta \tilde z)^2$")
plt.scatter(x, Temp_2, color='k', s=3, label=rf"Crack-Nickolson amb $\Delta \tilde t = 0.5 \, (\Delta \tilde z)^2$")
plt.xlabel(r"$z$ (m)")
plt.ylabel(r"$T$ ($^\circ $C)")
plt.legend(loc="best")
plt.title(r"Comparació Euler explícit amb Crack-Nicolson")
plt.legend(loc='best')

plt.subplot(1,2,2)
plt.scatter(x, Err_21, color='r', s=2, label=rf"Euler explícit amb $\Delta \tilde t=  (\Delta \tilde z)^2$")
plt.scatter(x, Err_22, color='indigo', s=2, label=rf"Euler explícit amb $\Delta \tilde t = 0.5 \, (\Delta \tilde z)^2$")
plt.xlabel(r"$\Delta z$ (m)")
plt.ylabel(r"Error normalitzat")
plt.title(r"Errors del mètode d'Euler explícit comparat amb Crack Nicolson")
plt.legend(loc='best')
plt.show()