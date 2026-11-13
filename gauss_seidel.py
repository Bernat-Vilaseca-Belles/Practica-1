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


<<<<<<< HEAD
def gauss_seidel(vector, N, error, A, b): #gamma, interval de z normalitzat, llistat de temperatura normalitzada de la iteració anterior excloent els extrems, temperatura inicial, nombre de punts en z, error que estem disposats a assumir.
    '''
    A=[np.zeros(N-2) for _ in range(N-2)] #Definició de la matriu.
    b=np.ones(N-2) #Definició del vector corresponent de Ax=b
    for i in range(N-2): b[i]=g*(dz**2)+vector[i] #Assignem els valors de b
    b[0]+=g*tau #Com que és el valor de j=2 i no està definida la operació per a j=1, cal afegir-hi el valor esctremal
    b[N-3]+=g*tau #Afegim el valor extremal de la temperatura per a j=N
    for i in range(N-2):
        for j in range(N-2):
            if j==i:
                A[i][j]=1+(2*g) #Definim els valors de A a la diagonal.
            if ((j-1)==i) or ((j+1)==i):
                A[i][j]=-g #Definim els valors de A al voltant de la diagonal. Els altres valors són 0
    '''

=======
def gauss_seidel(A, g, vector, b, N, error): #Matriu del sistema, gamma, llistat de temperatura normalitzada de la iteració anterior excloent els extrems, terme b de l'expressió Ax=b, nombre de punts en z, error que estem disposats a assumir.
>>>>>>> refs/remotes/origin/main
    b_1=True #No convergència. Boolean 1.
    k=0 #Comptador.
    while(b_1):
        x=vector.copy() #Fem una còpia del vector que li subministrem.
        echo=np.ones(N-2) #Eina matemàtica per operar les x. Represneta ser la x a la iteració k quan calculem la iteració k+1. Nom aleatori.
        E=np.ones(N-2) #Llista de l'error per dir quin error estem disposats a assumir.
        for i in range(N-2):
            echo[i]=x[i] #Definim els valors de echo com els de l'iteració k abans d'assignar els valors de k+1 a les x.
            suma_1, suma_2=0.0, 0.0 #Sumes de la fòrmula del gauss-seidel.
            if i>0: #Modifiquem només si la i és major a 0. 
                suma_1=A[i][i-1]*x[i-1] #Apliquem la suma, però com que tots els altres valors de la fila de A i j<i-1 són zero, no els hem de tenir en compte
            if i<N-3: #Només tenim en compte els termes que no siguin a l'extrem
                suma_2=A[i][i+1]*x[i+1] #Fem la suma per tots els valors no zero de j>i+1
            x[i]=(1/(A[i][i]))*(b[i]-suma_1-suma_2) #Apliquem la fórmula
            E[i]=abs(x[i]-echo[i])/x[i] #Calculem la convergència.
        if max(E)<=error: #Mirem si convergeix.
            b_1=False #Veiem que ha convergit, per tant canviem el boolean.
        k+=1 #Augmentem el comptador.
    return x, k