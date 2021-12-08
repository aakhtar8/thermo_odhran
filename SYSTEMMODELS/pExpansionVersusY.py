# plots compressor pressure vs liquid yield
# liqufecaction model using OOP
import sys
[sys.path.append(i) for i in ['.', '..']]
import copy
from UNITS.exchanger import Exchanger
from UNITS.Throttling_Valve import Throttling_Valve
from UNITS.Flash_Seprator import Flash_Sep
from stream import Stream
import numpy as np
name = {"CASRN": False, "ids": ["oxygen", "nitrogen"]}
pressures = np.arange(2500000, 101330, -500000)
liquid_produced = []
for pressure in pressures:

    air = Stream("air",Ids=name, Tcs=[154.4, 126.2], Pcs=[5.01E6, 3.356E6], Omegas=[0.022, 0.038], Mws=[32, 28.0134], Kijs=[[0, -0.0119], [-0.0119, 0]],composition=[0.21, 0.79], temp=298, pres=3000000)
    hotair = copy.deepcopy(air)
    # getting initial estimates for cold air
    tv_initlize = Throttling_Valve(hotair, discharge_pres=pressure)
    coldair = tv_initlize.discharge
    n = 1
    k = 0
    while True:
        print("ITERATION: ", n)
        # now starting the loop
        hx1 = Exchanger("hx1", hotair, coldair, approach_temp=1)

        tv1 = Throttling_Valve(hx1.hot_discharge, discharge_pres=pressure)

        flash1 = Flash_Sep("flash1", tv1.discharge)

        # check whether solution converged or not
        # temperture, flowrate and composition differnce should be within tolrence for coldair and vapor
        if abs(flash1.vapor.molar_flow - coldair.molar_flow) < 0.0001:        # flowrate
            if abs(flash1.vapor.temp - coldair.temp) < 0.0001:        # temperture
                if [x-y for x,y in zip(flash1.vapor.composition, coldair.composition)] < [0.000001]*len(flash1.vapor.composition):  # composition
                    k += 1      # running the simulation to get bit of stright line
                    fraction_liquified = flash1.liquid.molar_flow/flash1.feed.molar_flow 

        coldair = flash1.vapor
        n += 1
        if k == 5:
            break
        elif n == 500:         # SOLUTION HAS NOT CONVERGED IN 500 ITERATIONS
            # most likely there is no liquid yield in 500 iterations
            fraction_liquified = flash1.liquid.molar_flow/flash1.feed.molar_flow
            break
    liquid_produced.append(fraction_liquified)
# observing results
import matplotlib.pyplot as plt
import numpy as np
for i in range(len(liquid_produced)):
    if liquid_produced[i] is None:
        fraction_liquified[i] = 0.0
fraction_liquified = np.array(liquid_produced)
fraction_liquified *= 100
plt.plot(pressures/100000, fraction_liquified)
plt.title("Expansion Pressure vs Liquid Yield")
plt.xlabel("Throttling pressure(bar)")
plt.ylabel("Fraction Liquified (%)")
plt.show()
plt.savefig("Pexp_vs_liq.png")