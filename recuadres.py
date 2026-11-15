import numpy as np
import constants as const

N=101
delta_z = 1/(N-1) #valor de la discretització espacial normalitzada. delta_z'= h/(L_0*N). L_0=h ==> delta_z'=1/N.
delta_t= 0.5*(delta_z**2) #Conjunt de valors de les diferents discretitzacions temporals.
gamma= 0.5 #coficient

#recuadres
X, Y=[], []

for i in range(N-1):
    london=(i*delta_z*const.L_0, (i+1)*delta_z*const.L_0, (i+1)*delta_z*const.L_0, i*delta_z*const.L_0) #En x
    paris=(0,0,1,1) #Eix y
    X.append(london)
    Y.append(paris)