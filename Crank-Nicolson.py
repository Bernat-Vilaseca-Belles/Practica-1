import numpy as np
from constants import T_c,T_0
import math

#Definim la matriu A
def calcul_A(g):
    A = [np.zeros(N-2)]*(N-2)
    for i in range(N-2):
        for j in range(N-2):
            if i==j:
                A[i][i]= 2*(1+g)
            if i>0 and ((j-1)==i):
                A[i][j]=-g
            if ((j+1)==i):
                A[i][j]=-g

#Definim el vector b
def calcul_b(v, g, dt, tau):
    b1 = 2*(1-g)*v[0] + g*v[1] + 2*dt + g*tau
    bn = g*v[-1] + 2*(1-g)*v[-2] + 2*dt + g*tau
    br = np.zeros(N-4)
    for i in range(N-4):
        br[i]= g*v[i] + 2*(1-g)*v[i+1] + g*v[i+2] + 2*dt

    return [b1,br,bn]


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

