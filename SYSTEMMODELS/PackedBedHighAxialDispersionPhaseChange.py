import numpy as np
import scipy as sp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

minL = 0.05 # inlet liquid flow (kg/s)
U = 1  # heat transfer coefficient (
A = 0.1  # heat transfer area (m2)
V = 0.01  # liquid volume in vessel (m^3)
Tin = 300  # inlet gas temperature (K)
Cpl = 2  # liquid heat capacity (kJ/kg.K)
Cpv = 1  # vapour heat capacity (kJ/kg.K)
Cs = 1.1  # solid heat capacity (kJ/kg.K)
heatVap = 200  # heat of vaporisation (kJ/kg)
T1 = 105  # lower temperature bound for saturation (K
T2 = 110  # upper temperature bound for saturation
P = 200000  # pressure (Pa) (for vG)
MW = 28  # molecular mass of fluid
R = 287.058*28  # J/kg.K #gas constant
vL = 1/800  # liquid specific volume #m3/kg
Ms = 57.75  # mass of solid in bed

def cp(T):
    if T<T1:
        return Cpl
    elif T>T2:
        return Cpv
    else:
        cP = (Cpl*(T-T1) + heatVap + Cpv*(T2-T))/(T2-T1)
        #print("Cp = ",cP)
        return cP

def mass(T):
    if T<T1:
        return V/vL
    elif T>T2:
        return (P*MW*V)/(R*T2)
    else:
        m = V/(vL+((R*T2/P*MW)-vL)*((T-T1)/(T2-T1)))
        #print("m = ",m)
        return m


# return [DT/dt, DTs/dt] for values of t, using parameters passed in.
def my_system(t, y, minL, Ms, Cs, Tin, U, A):
        T, Ts = y
        Cp = cp(T)
        m = mass(T)
        # could also insert any variables that are dependent on t here
        DTdt = (minL * Cp * (Tin - T) - U * A * (T - Ts)) / (m * Cp)
        DTsdt = ((U * A * (T - Ts))) / (Ms * Cs)
        dydt = [DTdt, DTsdt] #solving for both derivatives manually
        return dydt

# starting condition
T0 = 80 #fluid starting tempearature (K)
Ts0 = 80 #solid starting temperature (K)
y0 = [T0, Ts0]

# time steps
t = (0, 7000) #time range
t_eval = np.linspace(0, 7000, 100)  # optional

# integrate
solution = solve_ivp(my_system, t, y0, args=(minL, Ms, Cs, Tin, U, A), method='LSODA', t_eval=t_eval)
#print(solution.y[0])  # T
#print(solution.y[1])  # Ts
#print(solution.y[0])

# generates data of simulation one. data2 is an array where each row is 
# a time step and the first column is the time
# and all the others are the solution
data2 = np.concatenate((np.array([solution.t]),solution.y)).transpose()

np.savetxt("data2.txt",data2)
print("data 2 loaded")
# # plot
plt.plot(solution.t, solution.y[0], 'b', label='T') #plotting fluid solution t vs T
# #plt.plot(solution.t, solution.y[1], 'g', label='Ts') #plotting solid solution t vs Ts
plt.title("Single Tank Fluid Temperature for U = 1 W/m2.K")
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()


