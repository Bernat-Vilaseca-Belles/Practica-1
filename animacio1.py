import numpy as np
from constants import T_c,T_0
import matplotlib.pyplot as plt
import constants as const
import math
from gauss_seidel import gauss_seidel as gs


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

#definim la funció color:
def color_T(T_a,T_b,T_1,T_2,T_3):
    d_1=T_2-T_1
    d_2=T_3-T_2
    T=(T_b+T_a)/2
    col=[0,0,0]
    if T<=T_2:
        col[0]=135+(120*T/d_1)
        col[1]=206-(206*T/d_1)
        col[2]=235-(235*T/d_1)
    if (T>T_2) and (T<=T_3):
        col[0]=255-(255*T/d_2)
    print(col)
    return (col[0]/255.0,col[1]/255.0,col[2]/255.0)

#Càlcul de les variables necessàries.:
tau = T_c/T_0 #Temperatura inicial normalitzada. tau=T(t=0)/T_0, T_0= k/sigma*(V_0**2).
T_1=(50+273.15)/T_0
T_2=(80+273.15)/T_0

#Establim les dades de la discratització:
N=101
delta_z = 1/(N-1) #valor de la discretització espacial normalitzada. delta_z'= h/(L_0*N). L_0=h ==> delta_z'=1/N.
delta_t= 0.5*(delta_z**2) #Conjunt de valors de les diferents discretitzacions temporals.
gamma= 0.5 #coficient

#Punts
z_1=0.0075/const.L_0
z_2=0.0125/const.L_0
j_1=math.ceil(z_1/delta_z)-1
j_2=math.floor(z_2/delta_z)-1

A = [np.zeros(N-2) for _ in range(N-2)]
for i in range(N-2):
    for j in range(N-2):
        if i==j:
            A[i][i]= 2*(1+gamma)
        if i>0:
            A[i][i-1]=-gamma
        if i<N-3:
            A[i][i+1]=-gamma


###Implementació del mètode per solucionar el problema.
foxtrot=[] #Llista formada per dues llistes que contindran les temperatures en les diferents posicions. Nom aleatori.
phi=tau*np.ones(N) # Assignem la temperatura normalitzada com a llista.
v=phi[1:N-1] #Agafem tots els valors de T que no siguin els valors
b_1=True #No s'ha superat la temperatura de 80ºc. Boolean 1
b_2=True #Fora de l'interval no s'ha superat la temperatura de 50ªC. Boolean 2
comptador=0

plt.figure(figsize=(16,9))
# creem els requadres on pintarem l'animació
X, Y=[], []

for i in range(N-1):
    london=(i*delta_z*const.L_0, (i+1)*delta_z*const.L_0, (i+1)*delta_z*const.L_0, i*delta_z*const.L_0) #En x
    paris=(0,0,1,1) #Eix y
    X.append(london)
    Y.append(paris)

x=np.linspace(0, const.L_0, N) #Definim l'eix x de la gràfica
while b_1 and b_2:
    v, total=gs(A, gamma, v, calcul_b(v, gamma, delta_t, tau), N, 0.0000000001) #Apliquem el mètode de gauss-seidel
    # for k in range(N-2): v[k]=u[k] #Reescribim els valors de v.
    u=v[j_1:j_2]
    w=np.append(v[:(j_1-1)], v[(j_2+1):])
    if (max(v)>T_2):
        b_1=False
    if (max(w)>T_1):
        b_2=False
    if b_1 and (max(u)>T_1):
        comptador+=1
    foxtrot=np.append(np.array([phi[0]]), np.array(v)) #Posem el valor de la temperatura a z=0 al començament de la llista.
    foxtrot=np.append(np.array(foxtrot), np.array([phi[-1]])) #Posem el valor de la temperatua a z=h al final de la llista.
    #posem el color a l'animació
    plt.cla()
    plt.vlines(x[j_1], ymin= 0, ymax=1, colors='k')
    plt.vlines(x[j_2], ymin= 0, ymax=1, colors='k')

    for i in range(N-1):
        plt.fill(X[i],Y[i],color=color_T(foxtrot[i],foxtrot[i+1],tau,T_1,T_2))
    plt.xlabel(r"$\Delta z$ (m)")
    plt.ylabel(r"$T$ ($^\circ $C)")
    plt.title(r"Solució del problema")
    plt.legend(loc='best')
    plt.pause(delta_t*const.t_0)
# print("Temps d'exposició: ", (comptador*delta_t*const.t_0))


plt.show()



