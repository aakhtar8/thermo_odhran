import numpy as np
import scipy as sp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

minL = 0.05 #inlet mass flow (kg/s)
Tsat = 105 #inlet temperature (K)
cpL = 2 #heat capacity (kJ/kg.K)
cpV = 1 #heat capacity (kJ/kg.K)
rhoL = 800 #liquid density (kg/m3)
volumeV = 0.01 #m3
area = 0.1 #m2
gasConstant = 287.058*28 #J/kgK
hVap = 200 #kJ/kg
P = 200000 #Pa
mW = 28 #amu
mS = 1 #kg
cpS = 1.1 #kJ/kgK
u = 1 #heat transfer coefficient between fluid and solid (W/m2.K)
T0 = 80 #initial fluid temp (K)

Tin = 300 #Inlet fluid temperature (K)

hIn = 200 + cpV*(Tin-Tsat) #Inlet enthalpy (kJ/kg)

def temp(H):
    #returns temperature and mass holdup given enthalpy
    if H<=0:
        T = Tsat + (H/cpL)
        #print("1")
    elif H >= hVap:
        T = Tsat + ((H-hVap)/cpV)
        #print("3")
    else:
        T = Tsat
        #print("2")

    return T

def mass(H,T):
    #returns temperature and mass holdup given enthalpy
    if H<=0:
        m = rhoL*volumeV
        #print("1")
    elif H >= hVap:
        m = P*mW*volumeV/(gasConstant*T)
        #print("3")
    else:
        m = volumeV/((1/rhoL)+((gasConstant*Tsat/(P*mW)) - (1/rhoL))*(H/hVap))
        #print("2")
    print(m)
    return m


def my_system(t, y, minL, mS, cpS, hIn, u, area):
    H, Ts = y
    T = temp(H)
    m = mass(H,T)
    # could also insert any variables that are dependent on t here
    DHdt = ((minL * (hIn - H)) - (u * area * (T - Ts))) / m
    DTsdt = ((u * area * (T - Ts))) / (mS * cpS)
    dydt = [DHdt, DTsdt]  # solving for both derivatives manually
    return dydt

# starting condition
H0 = cpL*(T0-Tsat) #fluid starting tempearature (K)
Ts0 = 80 #solid starting temperature (K)
y0 = [H0, Ts0]

# time steps
t = (0, 200) #time range
t_eval = np.linspace(0, 200, 101)  # optional

solution = solve_ivp(my_system, t, y0, args=(minL, mS, cpS, hIn, u, area), method='LSODA', t_eval=t_eval)
#print(solution.y[0])  # T
#print(solution.y[1])  # Ts

x = solution.y[0]
x1 = np.array(x)

Temperatures = []
for i in x1:
    Temperatures.append(temp(i))
#print("NP Temperatures")
#print(Temperatures)

# plot
#plt.plot(solution.t, solution.y[0], 'b', label='H_fluid') #plotting fluid solution t vs T

#Working with water/steam (assuming steam as ideal gas)

minL2 = 0.05  # inlet liquid flow (kg/s)
U2 = 1  # heat transfer coefficient (
A2 = 0.1  # heat transfer area (m2)
V2 = 0.01  # liquid volume in vessel (m^3)
Tin2 = 300  # inlet gas temperature (K)
Cpl2 = 2  # liquid heat capacity (kJ/kg.K)
Cpv2 = 1  # vapour heat capacity (kJ/kg.K)
Cs2 = 1.1  # solid heat capacity (kJ/kg.K)
heatVap2 = 200  # heat of vaporisation (kJ/kg)
T12 = 100  # lower temperature bound for saturation (K
T22 = 110  # upper temperature bound for saturation
P2 = 200000  # pressure (Pa) (for vG)
MW2 = 28  # molecular mass of fluid
R2 = 287.058*28  # J/kg.K #gas constant
vL2 = 1/1000  # liquid specific volume #m3/kg
Ms2 = 1  # mass of solid in bed

def cp(Tn):
    if Tn<T12:
        return Cpl2
    elif Tn>T22:
        return Cpv2
    else:
        cP2 = (Cpl2*(Tn-T12) + heatVap2 + Cpv2*(T22-Tn))/(T22-T12)
        #print("Cp = ",cP)
        return cP2

def mass2(Tn):
    if Tn<T12:
        return V2/vL2
    elif Tn>T22:
        return (P2*MW2*V2)/(R2*T22)
    else:
        m2 = V2/(vL2+((R2*T22/P2*MW2)-vL2)*((Tn-T12)/(T22-T12)))
        #print("m = ",m)
        return m2


# return [DT/dt, DTs/dt] for values of t, using parameters passed in.
def my_system2(t2, y2, minL2, Ms2, Cs2, Tin2, U2, A2):
        Tn, Ts = y2
        Cp2 = cp(Tn)
        m2 = mass2(Tn)
        # could also insert any variables that are dependent on t here
        DTdt2 = (minL2 * Cp2 * (Tin2 - Tn) - U2 * A2 * (Tn - Ts)) / (m2 * Cp2)
        DTsdt2 = ((U2 * A2 * (Tn - Ts))) / (Ms2 * Cs2)
        dydt2 = [DTdt2, DTsdt2] #solving for both derivatives manually
        return dydt2

# starting condition
T02 = 80 #fluid starting tempearature (K)
Ts02 = 80 #solid starting temperature (K)
y02 = [T02, Ts02]

# time steps
t2 = (0, 200) #time range
t_eval2 = np.linspace(0, 200, 101)  # optional

# integrate
solution2 = solve_ivp(my_system2, t2, y02, args=(minL2, Ms2, Cs2, Tin2, U2, A2), method='LSODA', t_eval=t_eval2)
#print(solution.y[0])  # T
#print(solution.y[1])  # Ts
#print(solution2.y[0])



plt.plot(solution.t, Temperatures, 'b', label='Fluid Temperature (dH/dt model)') #plotting fluid solution t vs T
plt.plot(solution2.t, solution2.y[0], 'b', marker = ".", label='Fluid Temperature (dT/dt model)') #plotting fluid solution t vs T
plt.plot(solution.t, solution.y[1], 'g', label='Fluid Temperature (dH/dt model)') #plotting solid solution t vs Ts
plt.plot(solution2.t, solution2.y[1], 'g', marker = ".", label='Solid Temperature (dT/dt model)') #plotting solid solution t vs Ts
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()

