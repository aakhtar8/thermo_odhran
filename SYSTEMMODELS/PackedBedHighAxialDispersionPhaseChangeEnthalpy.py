import numpy as np
import scipy as sp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

n=32

minL = 0.05 #inlet mass flow (kg/s)
Tsat = 105 #inlet temperature (K)
cpL = 2 #heat capacity (kJ/kg.K)
cpV = 1 #heat capacity (kJ/kg.K)
rhoL = 800 #liquid density (kg/m3)
volumeV = 0.01 #m3
minitial = 800*0.01
area = 0.1 #m2
gasConstant = 287.058*28 #J/kgK
hVap = 200 #kJ/kg
P = 200000 #Pa
mW = 28 #amu
mS = 57.75 #kg
cpS = 1.1 #kJ/kgK
u = 1 #heat transfer coefficient between fluid and solid (W/m2.K)
T0 = 80 #initial fluid temp (K)

Tin = 90 #Inlet fluid temperature (K)

hIn = cpL*(Tin-Tsat) #Inlet enthalpy (kJ/kg)

def temp(H):
    #returns temperature given enthalpy (after processing function)
    if H<=0:
        T = Tsat + (H/cpL)

    elif H >= hVap:
        T = Tsat + ((H-hVap)/cpV)
    else:
        T = Tsat
    return T

def mass(H,T):
    #returns mass holdup given enthalpy
    if H<=0:
        m = rhoL*volumeV
    elif H >= hVap:
        m = P*mW*volumeV/(gasConstant*T)
        #print("mass goldup = ", m)
    else:
        m = volumeV/((1/rhoL)+((gasConstant*Tsat/(P*mW)) - (1/rhoL))*(H/hVap))
    return m

def dpdh(H,T):
    # returns mass holdup given enthalpy
    if H <= 0:
        m = rhoL * volumeV
        dpdh = 0
    elif H >= hVap:
        m = P * mW * volumeV / (gasConstant * T)
        dpdh = -1*(m/(cpV*T*volumeV))
        # print("mass goldup = ", m)
    else:
        m = volumeV / ((1 / rhoL) + ((gasConstant * Tsat / (P * mW)) - (1 / rhoL)) * (H / hVap))
        dpdh = -1*(m**2/volumeV**2)*(((gasConstant*Tsat/(P*mW))-(1/rhoL))/hVap)

    return dpdh

def dhdtcalc(H,T,Ts,m):

    Dhdt = ((minL * (hIn - H)) - (u * area * (T - Ts))) / m

    return Dhdt

def dtsdtcalc(T,Ts):

    DTsdt = ((u * area * (T - Ts))) / (mS * cpS)

    return DTsdt

def dmdtcalc(dhdt,dpdh):
    dmdt = volumeV * dhdt * dpdh
    return dmdt

def my_system(t, y, minL, mS, cpS, hIn, u, area):
    H, Ts = y
    T = temp(H)
    m = mass(H,T)
    # could also insert any variables that are dependent on t here
    DHdt = dhdtcalc(H,T,Ts,m)
    DTsdt = dtsdtcalc(T,Ts)
    dydt = [DHdt, DTsdt]  # solving for both derivatives manually
    return dydt

# starting condition
H0 = cpL*(T0-Tsat) #fluid starting enthalpy (kJ/kg)
Ts0 = 80 #solid starting temperature (K)
y0 = [H0, Ts0]

# time steps
t = (0, 7000) #time range
t_eval = np.linspace(0, 7000, 1000)  # optional

solution = solve_ivp(my_system, t, y0, args=(minL, mS, cpS, hIn, u, area), method='LSODA', t_eval=t_eval)
#print(dhdtextract)

#converting fluid enthalpy array back to temperatures
x = solution.y[0] #fluid enthalpy array
x1 = np.array(x) #converting x to numpy array for computation below

xs = solution.y[1] #fluid enthalpy array
x2 = np.array(xs) #converting x to numpy array for computation below

#print(len(x1))
Temperatures = []
for i in x1:
    Temperatures.append(temp(i))

Temperatures = np.array(Temperatures)
massarray = []

for h,t in zip(x1,Temperatures):
    massarray.append(mass(h,t))
massarray = np.array(massarray)

dpdharray = []
for h,t in zip(x1,Temperatures):
    dpdharray.append(dpdh(h,t))

dhdtarray = []
for h,T,Ts,m in zip(x1,Temperatures,x2,massarray):
    dhdtarray.append(dhdtcalc(h,T,Ts,m))


dhdtarray = np.array(dhdtarray)
#print(dpdharray)
print("dpdh len: ",len(dpdharray))
print("dhdt len: ", len(dhdtarray))
print("mass array len: ",len(massarray))

dmdt = volumeV*(np.multiply(dhdtarray,dpdharray))
mout = [minL - dmdt for dmdt in dmdt]



#Creating Semi-log arrays for fluid and solid
solidTs = solution.y[1] #solid temperature array
semilogsolidtemperatures = [Tin - solidTs for solidTs in solidTs] #semi log solid temperatures
semilogfluidtemperatures = [Tin - Temperatures for Temperatures in Temperatures] #converting to semi log fluid temperatures
semilogmassflowout = [minL - mout for mout in mout] #INCORRECT
#print(semilogmassflowout)

#slope:
slope = np.log(semilogfluidtemperatures[300]- semilogfluidtemperatures[350])/(350-400)
#print("slope 200 = ", slope)
#print("time constant = ", 1/abs(slope))

#plt.plot(solution.t, semilogsolidtemperatures, 'g', label='Tin - T_solid') #plotting solid solution t vs Ts
#plt.plot(solution.t, semilogfluidtemperatures, 'b', label='Tin - T_Fluid') #plotting fluid solution t vs T
#plt.yscale('log')
# plot
plt.plot(solution.t, Temperatures, 'b', label='T_Fluid') #plotting fluid solution t vs T
#plt.plot(solution.t, solution.y[1], 'g', label='T_solid') #plotting solid solution t vs Ts
#plt.plot(solution.t, solution.y[0], 'm', label='H_fluid') #plotting fluid solution t vs T
#plt.plot(solution.t, massarray, 'c', label='Mass Holdup') #plotting fluid solution t vs T
#plt.plot(solution.t, mout, 'g', label='m_out') #plotting fluid solution t vs T
#plt.plot(solution.t, massarray, 'b', label='mass holdup', scaley = 'log') #plotting fluid solution t vs T
#plt.plot(solution.t, dpdharray, 'm', label='dp/dh') #plotting fluid solution t vs T
#print(dpdharray)
#plt.yscale('log')



#plt.plot(solution.t, semilogmassflowout, 'b', label='Semi-Log m_out') #plotting fluid solution t vs T

plt.legend(loc='best')
plt.xlabel('time')
plt.ylabel('Temperature (K)')
plt.title("Single Tank Model - Liquid Phase Only")
#plt.ylabel('Tin - T')
plt.grid()
plt.show()



