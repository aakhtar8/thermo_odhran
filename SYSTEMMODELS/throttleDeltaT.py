# shows the plot of temperture difference across valve
# liqufecaction model using OOP
import sys
[sys.path.append(i) for i in ['.', '..']]
import copy
from UNITS.exchanger import Exchanger
from UNITS.Throttling_Valve import Throttling_Valve
from UNITS.Flash_Seprator import Flash_Sep
from stream import Stream
name = {"CASRN": False, "ids": ["oxygen", "nitrogen"]}
air = Stream("air",Ids=name, Tcs=[154.4, 126.2], Pcs=[5.01E6, 3.356E6], Omegas=[0.022, 0.038], Mws=[32, 28.0134], Kijs=[[0, -0.0119], [-0.0119, 0]],composition=[0.21, 0.79], temp=298, pres=3000000)
hotair = copy.deepcopy(air)
# getting initial estimates for cold air
tv_initlize = Throttling_Valve(hotair, discharge_pres=101325)
coldair = tv_initlize.discharge

n = 1
k = 0       # counter for smooth plotting
valve_feed_temp = []
valve_discharge_temp = []
exchanger_cold_discharge = []
fraction_liquified = []
while n <= 200:
    print("ITERATION: ", n)
    # now starting the loop
    hx1 = Exchanger("hx1", hotair, coldair, approach_temp=1)
    valve_feed_temp.append(hx1.hot_discharge.temp)
    exchanger_cold_discharge.append(hx1.cold_discharge.temp)

    tv1 = Throttling_Valve(hx1.hot_discharge, discharge_pres=101325)
    valve_discharge_temp.append(tv1.discharge.temp)
    
    flash1 = Flash_Sep("flash1", tv1.discharge)
    fraction_liquified.append(flash1.liquid.molar_flow/flash1.feed.molar_flow)
    
    # check whether solution converged or not
    # temperture, flowrate and composition differnce should be within tolrence for coldair and vapor
    converged = False
    if abs(flash1.vapor.molar_flow - coldair.molar_flow) < 0.0001:        # flowrate
        if abs(flash1.vapor.temp - coldair.temp) < 0.0001:        # temperture
            if [x-y for x,y in zip(flash1.vapor.composition, coldair.composition)] < [0.000001]*len(flash1.vapor.composition):  # composition
                k += 1      # running the simulation to get bit of stright line 

    coldair = flash1.vapor
    n += 1
    if k == 5:
        break


# observing results
import matplotlib.pyplot as plt
import numpy as np
valve_feed_temp = np.array(valve_feed_temp)
valve_discharge_temp = np.array(valve_discharge_temp)
exchanger_cold_discharge = np.array(exchanger_cold_discharge)
temp_difference = valve_feed_temp - valve_discharge_temp
for i in range(len(fraction_liquified)):
    if fraction_liquified[i] is None:
        fraction_liquified[i] = 0.0
fraction_liquified = np.array(fraction_liquified)
fraction_liquified *= 100
x = np.arange(1, n)
fig, ax = plt.subplots()
plt.title('Temperature drop Across Valve')
ax.plot(x, temp_difference, color = 'g')
 
# giving labels to the axises
ax.set_xlabel('Iterations', color = 'r')
ax.set_ylabel('throttle temperture drop (K)', color = 'g')
# defining display layout
plt.tight_layout()
 
# show plot
plt.savefig("temp_drop.png")
plt.show()
