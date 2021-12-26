import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math


'''Computing Time step using Von-Neumann Stability equation'''
'''----------------------------*********************--------------------------------'''
def time_step(dt, dx): #come back to this
    a = rho_cm
    b = G*cf
    t = np.arange(0, 2*np.pi, 0.001)
    val = abs(1 - 2*km*dt/(a*dx**2) + 2*km*dt/(a*dx**2)*np.cos(t) - 1j*(b*dt)*np.sin(t)/(a*dx))
    s =  np.max(val)
    while s > 1:
        dt = dt/1.5
        val = abs(1 - 2*km*dt/(a*dx**2) + 2*km*dt/(a*dx**2)*np.cos(t) - 1j*(b*dt)*np.sin(t)/(a*dx))
        s =  max(val)
    return dt/2
'''----------------------------*********************--------------------------------'''

L = 3 #input length of packed bed
diameter = 0.20 #length to diameter ratio of 5:1
time = 300 #simulation time in seconds 
m_f = 0.09 #mass flow rate of fluid entering bed

T0 = 293 #initial temperature of packed bed 
T_in = 77 #temperature of fluid entering packed bed

epsilon = 0.4 #void fraction

rho_f = 1.225 #fluid density
kf = 0.025 #fluid thermal conductivity
cf = 1070 #fluid cp (was 700)

rho_s = 1800 #solid density
ks = 30 #solid thermal conductivity
cs = 880 #solid heat capacity (cp)

'''Given Data'''
'''*******************************************************************************'''
'''*******************************************************************************'''
Pr_f = 0.73 #prandtl number of fluid
c1 = 0.115 #elouali et al value
c2 = 1 #elouali et al value


flowarea = ((math.pi)/4)*diameter**2
G = m_f/flowarea
print(G)

d = 0.01 #particle diameter
u = m_f/(rho_f*flowarea) #packed bed fluid velocity
mu_f = 18.6 * 10**(-6) #fluid viscosity 
Re_p = rho_f*u*d/mu_f #fluid reynolds number

kf_star = kf*epsilon*(1 + c1*(Re_p*Pr_f)**c2) #effective thermal conductivity of the fluid

km = ks*(1 - epsilon*(ks - kf_star)/(kf_star + epsilon**(1/3)*(ks - kf_star))) #effective thermal conductivity of the mixture
rho_cm = epsilon*rho_f*cf + (1 - epsilon)*rho_s*cs #effective heat capacity of the packed bed
'''*******************************************************************************'''
'''*******************************************************************************'''

N = 100 #number of elements
x = np.linspace(0, L, N+1) #number of nodes
dx = x[1] - x[0] #length of elements

# Computing Time step using von-neumann stability analysis
dt = time_step(1, dx) #von neumann stable time step 
print(dt)

# Calculating constants
a1 = km*dt/(rho_cm*dx**2) - G*cf*dt/(2*rho_cm*dx) #discretised PDE constant a1
b1 = 2*km*dt/(rho_cm*dx**2) #discretised PDE constant b1
c1 = km*dt/(rho_cm*dx**2) + G*cf*dt/(2*rho_cm*dx) #discretised PDE constant c1


t = 0 #setting initial value of t = 0 
T = np.zeros((int(time/dt)+2, N+1)) #setting up Temperature - Position array

n = 0 #setting position equal to 0 (initial position)
# Initial Condition
T[n, :] = T0 #setting T = T0 at t=0 for all positions (initial temperature of entire packed bed)
# Boundary Condition
T[:, 0] =  T_in #setting all times at n=0 to inlet fluid temperature

while t <= time: #running simulation until required simulation time
    for i in range(1,N): #solving for temperature at every increment
        T[n+1,i] = T[n,i] + a1*T[n,i+1] - b1*T[n,i] + c1*T[n,i-1] #using three nodes (i-1,i,i+1) at time n to solve for temperature at n+1, position i
      
    T[n+1,-1]  = T[n+1,-2] #setting final position temperature equal to second last temperature
    T[n+1, 0] =  T_in #resetting position 0 as inlet fluid temperature
    n = n + 1 #incrementing position
    t = t + dt #incrementing time
    
    
'''Plotiing the Temperature at different time steps'''
'''###################################################################################'''
i = 0 #setting i = 0 for plot increment
y = np.linspace(0, 0.1, len(x)) #setting up temperature plot linspace to be same as position linspace
X,Y = np.meshgrid(x,y) #creating rectangular grid from x and y arrays
fig, [ax1, ax2] = plt.subplots(2,1, figsize = (10,9)) #setting up two plots on one figure
while i <= n:
    ax1.plot(x, T[int(i),:], label = "@ t = {} s".format(np.round(int(i)*dt,1))) 
    #incrementing i such that temperature profile prints 
    #every 20% of total time increment
    i = i + 0.2*n

ax1.set_xlabel("Length (m)", fontsize = 15) #setting x axis, plot 1
ax1.set_ylabel("Temperature (K)", fontsize = 15) #setting y axis, plot 1
ax1.legend(loc=(0.60,0), prop={"size":16}) #setting legend, plot 1
ax1.grid(True) #setting grid plot 1


row, cols = X.shape #setting matrix shape
ax2.set_aspect(1) #setting aspect ratio

z = T[n-1,:] #capturing what the temperature distribution looks like at the end of charge/discharge
Z = np.tile(z, (cols, 1)) 
cp = ax2.contourf(X, Y, Z, 500, cmap=cm.coolwarm)
cbar = fig.colorbar(cp, ax = ax2, orientation = "horizontal") #creating colorbar plot1.2
plt.legend([f'Temperature (K) @ {np.round(t,1)} s']) #colorbar legend

plt.show()