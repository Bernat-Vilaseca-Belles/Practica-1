import numpy as np
from constants import T_c, T_0
N=101
n_max=1000
epsilon=0.0001 #error per  a poder comprovar la onvergència
T_tilde_0=T_c/T_0
def gauss_seidel(T_anterior,Delta_t,Delta_z,n_max, epsilon): #T_anterior és la temperatura en el temps m
    gamma=Delta_t/Delta_z**2
    T_actual=T_anterior.copy() #el primer valor de T_actual (T en el temps m+1) és el valor anterior
    for n in range(n_max):
        T_i_anterior=T_actual.copy() #guardem el valor anterior

        for i in range(1,N-1):
            #Apliquem la fórmula de Guass-seidel
            T_actual[i]=(T_anterior[i] + Delta_t + gamma*(T_actual[i-1]+T_i_anterior[i+1]))/(1+2*gamma)
    #apliquem les condicons de contron
        T_actual[0]=T_tilde_0
        T_actual[-1]=T_tilde_0
        #comprovem si la diferencia entre el terme actual i l'anterior és més petita que epsilon i el bucle pot aturar-se
        error=np.max(np.abs(T_actual-T_i_anterior))
        if error<epsilon:
            break
    return T_actual