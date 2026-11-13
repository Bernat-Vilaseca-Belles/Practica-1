import numpy as np


def cranknicolson(N,g,dt,Tc,ta):
    #Definim la matriu A
    A = [np.zeros(N-2)]*(N-2)
    for i in range(N-2):
        for j in range(N-2):
            if i==j:
                A[i][i]= 2*(1+g)
            if ((j-1)==i) or ((j+1)==i):
                A[i][j]=-g

    #Definim la matriu B
    B = [np.zeros(N-2)]*(N-2)
    for i in range(N-2):
        for j in range(N-2):
            if i==j:
                A[i][i]= 2*(1-g)
            if ((j-1)==i) or ((j+1)==i):
                A[i][j]=-g

    #Definim C
    C = [dt](N-2)

    #Definim el vector temperatura inicial
    T = [Tc](N-2)

    for p in range(ta):
        for j in range(N-2): 
            T[j] = B[i][j]


