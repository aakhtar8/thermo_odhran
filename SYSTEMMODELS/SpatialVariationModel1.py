import numpy as np
import scipy as sp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math
from FLUIDS.fluid import Fluid
from FLUIDS.mixture import Mixture

nitrogen = Fluid("n2")
oxygen = Fluid("o2")
Kij = np.array([[0, -0.0119], [-0.0119, 0]]) #interaction parameters
air = Mixture("Air",[oxygen, nitrogen], {"basis": "mole", "fraction": [0.20, 0.80]}, interaction=Kij)

n = 3  # number of elements in spatial array
P = 200000  # Pa

t_bubble = air.ideal_Tbubble(P)        # bubble point temp at default press (200000)
t_dew = air.ideal_Tdew(P)
Tsat = (t_bubble+t_dew)/2

cpL = 2  # heat capacity (kJ/kg.K)
cpV = 1  # heat capacity (kJ/kg.K)
rhoL = 800  # liquid density (kg/m3)
volumeV = 0.01 / n  # m3
area = 0.1 / n  # m2
mS = 1 / n  # kg
gasConstant = 287.058 * 28  # J/kgK
hVap = 200  # kJ/kg

mW = 28  # amu
mS = 57.75 / n  # kg
cpS = 1.1  # kJ/kgK
u = 1  # heat transfer coefficient between fluid and solid (W/m2.K)
T0 = 80  # initial fluid temp (K)

Tin = 90  # Inlet fluid temperature (K)
minL = 0.05  # inlet mass flow (kg/s)

# function for enthalpy
hIn = cpL * (Tin - Tsat)  # Inlet enthalpy (kJ/kg)

def temp(H):
    # returns temperature given enthalpy (after processing function)
    if H <= 0:
        T = Tsat + (H / cpL)

    elif H >= hVap:
        T = Tsat + ((H - hVap) / cpV)
    else:
        T = Tsat
    return T


def mass(H, T):
    # returns mass holdup given enthalpy
    if H <= 0:
        m = rhoL * volumeV
    elif H >= hVap:
        m = P * mW * volumeV / (gasConstant * T)
        # print("mass goldup = ", m)
    else:
        m = volumeV / ((1 / rhoL) + ((gasConstant * Tsat / (P * mW)) - (1 / rhoL)) * (H / hVap))
    return m


def dpdh(H, T):
    # returns mass holdup given enthalpy
    if H <= 0:
        m = rhoL * volumeV
        dpdh = 0
    elif H >= hVap:
        m = P * mW * volumeV / (gasConstant * T)
        dpdh = -1 * (m / (cpV * T * volumeV))
        # print("mass goldup = ", m)
    else:
        m = volumeV / ((1 / rhoL) + ((gasConstant * Tsat / (P * mW)) - (1 / rhoL)) * (H / hVap))
        dpdh = -1 * (m ** 2 / volumeV ** 2) * (((gasConstant * Tsat / (P * mW)) - (1 / rhoL)) / hVap)

    return dpdh


def my_system(t, y, n, hIn, Tin, minL, mS, cpS, u, area):
    dydt = np.zeros(3 * n)

    # y = [h_0, Ts_0, m_0, ... h_n, Ts_n, m_n]
    # y[0] = hIn
    # y[1] = Tin
    # y[2] = minL

    M_dot = minL

    i=0
    T, m = temp(y[i]), mass(y[i], temp(y[i]))
    M_T = u * area * (y[i + 1] - T)
    dydt[i] = (M_dot * (hIn - y[i]) + M_T) / m
    dydt[i + 1] = -M_T / (mS * cpS)
    dydt[i + 2] = M_dot * dpdh(y[i], T) * volumeV / n


    for i in range(3, 3 * n, 3):
        T, m = temp(y[i]), mass(y[i], temp(y[i]))
        M_T = u * area * (y[i + 1] - T)

        # [h, T_S, m]
        dydt[i] = (M_dot * (y[i - 3] - y[i]) + M_T) / m
        dydt[i + 1] = -M_T / (mS * cpS)
        dydt[i + 2] = M_dot * dpdh(y[i], T) * volumeV / n

        M_dot -= dydt[i + 2]

    return dydt


# starting condition
H0 = cpL * (T0 - Tsat)  # fluid starting enthalpy (kJ/kg)
Ts0 = 80  # solid starting temperature (K)
m0 = mass(H0, T0)

y0 = np.copy(np.asarray([H0, Ts0, m0] * n))
y0[0] = H0#hIn
y0[1] = Ts0#Tin
y0[2] = m0#minL

# time steps
t_max = 1000
t = (0, t_max)  # time range
t_eval = np.linspace(0, t_max, 1000 * t_max)  # optional

parameters = (n, hIn, Tin, minL, mS, cpS, u, area)
solution = solve_ivp(my_system, t, y0, args=parameters, method='LSODA', t_eval=t_eval)

h_arr = solution.y[0::3]
Ts_arr = solution.y[1::3]
m_arr = solution.y[2::3]

time = solution.t
# dt = 2*(time[2:]-time[:-2])

# #recalculating derivative arrays
# dhdt_arr = (h_arr[:,2:]-h_arr[:,:-2]) / dt
# dTsdt_arr = (Ts_arr[:,2:]-Ts_arr[:,:-2]) / dt
# dmdt_arr = (m_arr[:,2:]-m_arr[:,:-2]) / dt

dydt = np.zeros((3 * n, time.size))
for i in range(time.size):
    dydt[:, i] = my_system(t, solution.y[:, i], n, hIn, Tin, minL, mS, cpS, u, area)

dhdt_arr = dydt[0::3]
dTsdt_arr = dydt[1::3]
dmdt_arr = dydt[2::3]

# creating the fluid array from the solid array
Tf_arr = np.zeros_like(Ts_arr)
Tf_arr = Ts_arr + dTsdt_arr * (mS * cpS) / (u * area)
# Tf_arr[:, 0] = Ts0  # I.C cannot be calculated as there is no derivative at the boundary. Manually entering here
data1 = np.concatenate([np.array([time]), Tf_arr]).transpose()
# data1 = np.stack((data1,data1),axis=2)

np.savetxt("data1.txt",data1)
print("data 1 loaded")



# snapshot = 100000  # snapshot time spacing: 1000=1s
# style = ['-', '--', '-.', ':', '-', '--', '-.', ':','-', '--']
# for i in range(len(style)):
#     plt.plot(Ts_arr[:, i*snapshot], style[i],
#              color='b',
#              label="t="+str(i*snapshot/1000)+"s")

# plt.ylabel(r'T(K)')
# plt.xlabel(r'position(j)')
# plt.legend(loc='lower left')
# plt.show()

# plt.plot(time,Tf_arr[0, :], '-.', color='b', label="position 0")
# plt.plot(time,Tf_arr[1, :], '--', color='b', label="position 1")
# plt.plot(time,Tf_arr[2, :], ':', color='b', label="position 2")
# from plot_together import plot_together
# plot_together(data1[:,[0,1]:], '--', color='b', label="position 0")
# plot_together(data1[:,[0,2],:], '--', color='b', label="position 1")
# plot_together(data1[:,[0,3],:], ':', color='b', label="position 2")
# plt.ylabel(r'T (K)')
# plt.xlabel(r'Time')
# plt.ylim((80, 91))
# plt.grid('both')
# plt.title("T fluid for all tanks (n=3)")
# plt.legend(loc='lower right')
# plt.show()
