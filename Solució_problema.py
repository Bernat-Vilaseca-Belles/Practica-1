import numpy as np
from constants import T_c,T_0
import matplotlib.pyplot as plt
import solucio_analitica as sa
import constants as const
import math
from gauss_seidel import gauss_seidel as gs

"""
Resolem el problema amb el mètode de Crank-Nicolson amb gamma=0.5. És la discretització i el mètode que millor funciona, entre els que hem probat.
"""


#Definim el vector b. S'extreu del mètode de Crank-Nicolson.
def calcul_b(v, g, dt, tau): #Calculem la b per a cada iteració.
    b=np.ones(N-2)
    deca=2*(1-g) #Ho definim per estalviar-nos operacions
    golf=2*dt #Ho definim per estalviar-nos operacions.
    for i in range(1, N-3): #Agafem els valors que no són dels extrems. En els altres hem de tenir en compte els valors que cauen fora de la matriu A.
        b[i]=g*v[i-1]+deca*v[i]+ g*v[i+1] + golf
    b[0] = 2*g*tau+ deca*v[0] + g*v[1] + golf
    b[N-3] = g*v[N-4] + deca*v[N-3] + 2*g*tau + golf 
    return b


#Càlcul de les variables necessàries.:
tau = T_c/T_0 #Temperatura inicial normalitzada. tau=T(t=0)/T_0, T_0= k/sigma*(V_0**2).
T_1=(50+273.15)/T_0 #Temperatura normalitzada de 50ºC.
T_2=(80+273.15)/T_0 #Temperatura normalitzada de 80ºC

#Establim les dades de la discratització:
N=101
delta_z = 1/(N-1) #valor de la discretització espacial normalitzada. delta_z'= h/(L_0*N). L_0=h ==> delta_z'=1/N.
delta_t= 0.5*(delta_z**2) #Conjunt de valors de les diferents discretitzacions temporals.
gamma= 0.5 #coficient

#Punts
z_1=0.0075/const.L_0 #Punt on comença la regió malalta.
z_2=0.0125/const.L_0 #Punt on s'acaba la regió malalta.
j_1=math.ceil(z_1/delta_z)-1 #Valor de la iteració on comença la regió malalta.
j_2=math.floor(z_2/delta_z)-1 #Valor de la iteració on acaba la regió malalta.

A = [np.zeros(N-2) for _ in range(N-2)] #Creem una matriu de (N-2)x(N-2) plena de zeros.
for i in range(N-2):
    for j in range(N-2):
        if i==j: #Valors de la diagonal.
            A[i][i]= 2*(1+gamma)
        if i>0: #Valors sota de la diagonal. No agafem el valor de la primera fila, perquè no està dins de la matriu. El tenim en compte en el càlcul de la b.
            A[i][i-1]=-gamma
        if i<N-3: #Valors sota de la diagonal. No agafem el valor de l'última fila, perquè no està dins de la matriu. El tenim en compte en el càlcul de la b.
            A[i][i+1]=-gamma


###Implementació del mètode per solucionar el problema.
foxtrot=[] #Llista formada per dues llistes que contindran les temperatures en les diferents posicions. Nom aleatori.
phi=tau*np.ones(N) # Assignem la temperatura normalitzada com a llista.
v=phi[1:N-1] #Agafem tots els valors de T que no siguin els valors
b_1=True #No s'ha superat la temperatura de 80ºc. Boolean 1
b_2=True #Fora de l'interval no s'ha superat la temperatura de 50ªC. Boolean 2
comptador_1=0
comptador_2=0
while b_1 and b_2:
    v, total=gs(A, gamma, v, calcul_b(v, gamma, delta_t, tau), N, 10**(-12)) #Apliquem el mètode de gauss-seidel
    u=v[j_1:j_2] #Llistat de temperatures dins de la regió malalta.
    w=np.append(v[:(j_1-1)], v[(j_2+1):]) #Llistat de temperatures fora de la regió malalta.
    if (max(v)>T_2): #Si la temperatura a quasevol punt és major a 80ºC
        b_1=False
    if (max(w)>T_1): #Si fora de la regió malalta algun punt arriba als 50ºC.
        b_2=False
    if b_1 and (max(u)>T_1): #Comptem quan temps passa la regió malalta per sobre dels 50ºC.
        comptador_1+=1
    comptador_2+=1
print("Temps d'exposició: ", (comptador_1*delta_t*const.t_0))
print("Temps amb la màquina engegada: ", (comptador_2*delta_t*const.t_0))
foxtrot=np.append(np.array([phi[0]]), np.array(v)) #Posem el valor de la temperatura a z=0 al començament de la llista.
foxtrot=np.append(np.array(foxtrot), np.array([phi[-1]])) #Posem el valor de la temperatua a z=h al final de la llista.
for i in range(len(foxtrot)):
    foxtrot[i]=(foxtrot[i]*const.T_0)-273.15 #Canviem la llista a temperatures no normalitzades.

Temperatura=foxtrot.copy()
x=np.linspace(0, const.L_0, N) #Definim l'eix x de la gràfica

#recuadres. Fem els rectangles per dibuixar les tres regions.
X_1=(0, x[j_1], x[j_1], 0)
X_2=(x[j_1], x[j_2], x[j_2], x[j_1])
X_3=(x[j_2], const.h, const.h, x[j_2])
Y_1=(min(Temperatura), min(Temperatura), max(Temperatura), max(Temperatura))
Y_2=(min(Temperatura), min(Temperatura), max(Temperatura), max(Temperatura))
Y_3=(min(Temperatura), min(Temperatura), max(Temperatura), max(Temperatura))


plt.figure(figsize=(16,6))
plt.scatter(x, Temperatura, color='purple', s=3.5, label=rf"Temperatura")
plt.vlines(x[j_1], ymin= min(Temperatura), ymax=max(Temperatura), colors='k') #Delimitem la primera frontera entre la regió sana i la malalta.
plt.vlines(x[j_2], ymin= min(Temperatura), ymax=max(Temperatura), colors='k') #Delimitem la segona frontera entre la regió sana i la malalta.
plt.hlines(y=50, xmin=0, xmax=const.h, linewidths=0.8, colors='k')
plt.fill(X_1, Y_1, color='skyblue', alpha=0.5, label="zona sana") #Emplenem el rectangle de la regió sana.
plt.fill(X_2, Y_2, color='red', alpha=0.3, label="zona malalta") #Emplenem el recuadre de la regió malalta.
plt.fill(X_3, Y_3, color='skyblue', alpha=0.5) #Emplenem el recuadre de l'altra regió sana.
plt.xlabel(r"$z$ (m)")
plt.ylabel(r"$T$ ($^\circ $C)")
plt.title(r"Solució del problema")
plt.legend(loc='best')
plt.show()
