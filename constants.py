from numpy import sqrt
k=0.56 #Conductivitat tèrmica (W/(m*K))
sigma=0.472 #Conductivitat elèctrica ((A^2*s^3)/(kg*m^3))
rho=1081 #Densitat (kg/m^3)
c_v=3698 #Calor específic (J/(kg*K))
V_eff=40/sqrt(2) #Voltatge eficaç ((kg*m^2)/(A*s^3))
V_0=40 #Voltatge de pic a pic ((kg*m^2)/(A*s^3))
alpha=(k/(c_v*rho))
h=0.02 #gruix total de la paret (m)
T_0=(sigma*(V_0**2))/(2*k) #temperatura característica (K)
L_0=h #longitud característica (m)
t_0=L_0**2/alpha #temps característic (s)
t_a=0.025 #temps adimensional
T_c=273.15+36.5 #temperatura a la frontera (K)
