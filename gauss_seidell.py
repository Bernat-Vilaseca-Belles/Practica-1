import numpy as np
from constants import T_c, T_0
# N=101
# n_max=1000
# epsilon=0.0001 #error per  a poder comprovar la onvergència
# T_tilde_0=T_c/T_0
# def gauss_seidel(T_anterior,Delta_t,Delta_z,n_max, epsilon): #T_anterior és la temperatura en el temps m
#     gamma=Delta_t/Delta_z**2
#     T_actual=T_anterior.copy() #el primer valor de T_actual (T en el temps m+1) és el valor anterior
#     for n in range(n_max):
#         T_i_anterior=T_actual.copy() #guardem el valor anterior

#         for i in range(1,N-1):
#             #Apliquem la fórmula de Guass-seidel
#             T_actual[i]=(T_anterior[i] + Delta_t + gamma*(T_actual[i-1]+T_i_anterior[i+1]))/(1+2*gamma)
#     #apliquem les condicons de contron
#         T_actual[0]=T_tilde_0
#         T_actual[-1]=T_tilde_0
#         #comprovem si la diferencia entre el terme actual i l'anterior és més petita que epsilon i el bucle pot aturar-se
#         error=np.max(np.abs(T_actual-T_i_anterior))
#         if error<epsilon:
#             break
#     return T_actual

#definim la funció per fer el gauss seidel
def gauss_seidel(g, dz, vector, N, error): #gamma, interval de z normalitzat, llistat de temperatura normalitzada, excloent els extrems.
    A=[np.zeros(N-2)]*(N-2) #creem la matriu A de mida (N-2)x(N-2)
    bravo=np.ones(N-2)
    b=np.ones(N-2)
    for i in range(N-2): b[i]=g*(dz**2)+vector[i] #creem el vector b de components b[i] i mida N-2

    for i in range(N-2):
        for j in range(N-2):
            if j==i:
                A[i][j]=1+(2*g) #els elements de la diagonal principal valen 1+2*gamma
            if ((j-1)==i) or ((j+1)==i):
                A[i][j]=-g #els elements de les diagonals adjacents valen gamma i menys gamma
    b_1=True #No convergència. Boolean 1.
    k=0
    #fem el while per aplicar la formula de Gauss Seidel
    while(b_1):
        x=vector
        echo=np.ones(N-2)
        resta=np.ones(N-2)
        for i in range(N-2): 
            echo[i]=x[i] #guardem el valor de la iteració anteruor
            suma_1, suma_2=0, 0 #comencem amb els dos sumatoris igual a 0
            if i>1:
                for j in range(0, i-2):
                    suma_1+=A[i][j]*x[j]
            for j in range(i+1, N-2):
                suma_2+=A[i][j]*x[j]
            x[i]=(1/A[i][i])*(b[i]-suma_1-suma_2)
            resta[i]=x[i]-echo[i]
        if max(resta)<=error:
            b_1=False
        k+=1
    return x, k