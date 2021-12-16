import numpy as np
import scipy as sp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

m = 5 #mass of fluid in bed
min = 1 #flowrate in
M = 5 #mass of solid bed
Cp = 1 #heat capacity of fluid (kJ/kgK) #air
Cs = 1.1 #heat capacity of solid (kJ/kgK) #quartzite
U = 1 #heat transfer coefficient between fluid and solid #random value
A = 1 #heat transfer area #random value
Tin = 350 #fluid temperature in

# y = [T, Ts]`
# return [DT/dt, DTs/dt] for values of t, using parameters passed in.
def my_system(t, y, min, m, M, Cp, Cs, Tin, U, A):
        T, Ts = y
        # could also insert any variables that are dependent on t here
        DTdt = (min * Cp * (Tin - T) - U * A * (T - Ts)) / (m * Cp)
        DTsdt = ((U * A * (T - Ts))) / (M * Cs)
        dydt = [DTdt, DTsdt] #solving for both derivatives manually
        return dydt

# starting condition
T0 = 293 #fluid starting tempearature (K)
Ts0 = 293 #solid starting temperature (K)
y0 = [T0, Ts0]

# time steps
t = (0, 100) #time range
t_eval = np.linspace(0, 100, 101)  # optional

# integrate
solution = solve_ivp(my_system, t, y0, args=(min, m, M, Cp, Cs, Tin, U, A), method='LSODA', t_eval=t_eval)
#print(solution.y[0])  # T
#print(solution.y[1])  # Ts

# plot
plt.plot(solution.t, solution.y[0], 'b', label='T') #plotting fluid solution t vs T
plt.plot(solution.t, solution.y[1], 'g', label='Ts') #plotting solid solution t vs Ts
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()






