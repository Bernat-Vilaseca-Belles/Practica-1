import numpy as np
import math

fitxer = open("Resultats ER.txt", "w") #Fitxer on es guardaran els resultats de les simulacions.

#En aquest codi s'ha calculat respecte a les funcions normalitzades. A raó d'això, cal tenir en compte els valors dels coeficients de normalització. En aquest cas, s'ha utilitzat la funció de normalització de la distribució normal, que és 1/(sigma*sqrt(2*pi)), on sigma és la desviació estàndard. Això significa que els valors dels coeficients de normalització i els valors de les constants.

#Càlcul de les constants de normalització necessàries:
tau = (36.5+273)*0.472*(40**2)/0.56 #Temperatura inicial normalitzada. tau=T(t=0)/T_0, T_0= k/sigma*(V_0**2).
T_M = (80+273)*0.472*(40**2)/0.56 #temperatura 80º normalitzada.
T_m = (50 +273)*0.472*(40**2)/0.56 #temperatura 50º normalitzada.
L_0 =0.02

#Comencem establint les dades de la discratització:
N=101
h=0.02
delta_z = 1/(N-1) #valor de la discretització espacial normalitzada. delta_z'= h/(L_0*N). L_0=h ==> delta_z'=1/N.
delta_t= (0.51*(delta_z**2), 0.49*(delta_z**2), 0.25*(delta_z**2)) #Conjunt de valors de les diferents discretitzacions temporals.
gamma= [(o/(delta_z**2)) for o in delta_t] #coficient
j_1= math.ceil(0.0075/L_0/delta_z) -1 #Punt de la component z on comença la zona malalta. Utilitzem math.ceil per arrodonir a l'alça. Així ens assegurem que estem dins de la regió malalta.
j_2= math.floor(0.0125/L_0/delta_z) -1 #Punt de la component z on acaba la zona malalta. Utilitzem math.floor per arrodonir a la baixa. Així ens assegurem que estem dins de la regió malalta.

# print(j_1)
# print(j_2)
#Elements per imprimir els resultats:
foxtrot = np.ones(3) #Llista de temps que ha aguantat el teixit en les condicions correctes per a cada discretització. Nom escollit a l'atzar.

for p in range(3): #per cada delta_t i gamma.
    #Establim les condicions inicials:
    phi = tau*np.ones(N) #Temperatura normalitzada a t=0 s
    charlie= np.ones(N) #Nom aleatori per al càlcul de la temperatura a la següent iteració.
    b_1, b_2 =True, True #Conjunt de boleans per a les operacions
    comptador_1 =0 #nombre d'iteracions a partir de la qual la zona malalta supera la temperatura de 50ºC.
    comptador_2 =0 #nombre d'iteracions fins que la zona sana supera els 50ºC o algun punt supera els 80ºC.
    while b_1 and b_2:
        for j in range(N-2):
            charlie[j] = phi[j] + gamma[p]*(phi[j+2]- 2*phi[j+1]+ phi[j] + (delta_z**2)) #Càlcul de la següent iteració de la temperatura en funció del temps.
        phi =charlie
        if max(phi)> T_M:
            b_1=False
        daca = phi[j_1:j_2] #Nom aleatori per denominar els punts de la zona malalta.
        eco=np.append(phi[:(j_1-1)], phi[(j_2+1):]) #Nom aleatori per denominar els punts de la zona sana
        if min(daca) > T_m: #Si la temperatura de la zona malalta supera els 50ºC, comencem a comptar quan temps passa.
            comptador_1+=1 
        if max(eco) > T_m: #Si la temperatura de la zona sana supera els 50ºC, s'atura la simulació.
            b_2=False
        comptador_2+=1 
    foxtrot[p]=comptador_2*delta_t[p]*(h**2)*3698*1081/0.56 #Temps no normalitzat que ha aguantat el teixit en la discretització p+1.
print("*****\nResultats de les 3 simulacions:\n")
print("Temps de resistència del teixit per a discretització normalitzada delta_t=0.51 delta_z^2:", foxtrot[0], "\n")
print("Temps de resistència del teixit per a discretització normalitzada delta_t=0.49 delta_z^2:", foxtrot[1], "\n")
print("Temps de resistència del teixit per a discretització normalitzada delta_t=0.25 delta_z^2:", foxtrot[2])

fitxer.write("*****\nResultats de les 3 simulacions:\n")
fitxer.write("Temps de resistència del teixit per a discretització normalitzada delta_t=0.51 delta_z^2: "+ str(foxtrot[0]) + "\n")
fitxer.write("Temps de resistència del teixit per a discretització normalitzada delta_t=0.49 delta_z^2: "+ str(foxtrot[1]) + "\n")
fitxer.write("Temps de resistència del teixit per a discretització normalitzada delta_t=0.25 delta_z^2: "+ str(foxtrot[2]))
fitxer.close()