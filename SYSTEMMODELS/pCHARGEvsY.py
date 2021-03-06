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
pressures = np.arange(500000, 20100000, 500000)
liquid_produced = []
for pressure in pressures:

    air = Stream("air",Ids=name, Tcs=[154.4, 126.2], Pcs=[5.01E6, 3.356E6], Omegas=[0.022, 0.038], Mws=[32, 28.0134], Kijs=[[0, -0.0119], [-0.0119, 0]],composition=[0.21, 0.79], temp=298, pres=pressure)
    hotair = copy.deepcopy(air)
    # getting initial estimates for cold air
    tv_initlize = Throttling_Valve(hotair, discharge_pres=101325)
    coldair = tv_initlize.discharge
    n = 1
    k = 0
    # initilizer for the looping issue for convergence
    temps = []
    flows = []
    compos = []
    while True:
        print("ITERATION: ", n)
        # now starting the loop
        hx1 = Exchanger("hx1", hotair, coldair, approach_temp=1)

        tv1 = Throttling_Valve(hx1.hot_discharge, discharge_pres=101325)

        flash1 = Flash_Sep("flash1", tv1.discharge)
        a = [flash1.vapor.molar_flow - coldair.molar_flow, flash1.vapor.temp - coldair.temp, [x-y for x,y in zip(flash1.vapor.composition, coldair.composition)]]
        flows.append(a)
        if n > 10:  # perform the average check after a minimum of 10 iterations
            if (flows[n-1][0] - flows[n-3][0]) < 0.0001:    # flow have alternating convergence
                if (flows[n-1][1] - flows[n-3][1]) < 0.0001:        # tempertures have alternative convergencee
                    if [x-y for x,y in zip(flows[n-1][2], flows[n - 3][2])] < [0.00001]*len(flows[n-1][2]):
                        # all variables have alternatively converged
                        fraction_liquified = flash1.liquid.molar_flow/flash1.feed.molar_flow
                        break

        # check whether solution converged or not
        # temperture, flowrate and composition differnce should be within tolrence for coldair and vapor
        if abs(flash1.vapor.molar_flow - coldair.molar_flow) < 0.01:        # flowrate
            if abs(flash1.vapor.temp - coldair.temp) < 0.01:        # temperture
                if [x-y for x,y in zip(flash1.vapor.composition, coldair.composition)] < [0.001]*len(flash1.vapor.composition):  # composition
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
        liquid_produced[i] = 0.0
fraction_liquified = np.array(liquid_produced)
fraction_liquified *= 100
plt.plot(pressures/100000, fraction_liquified)
plt.title("Compressor Pressure vs Liquid Yield")
plt.xlabel("Compressor pressure(bar)")
plt.ylabel("Fraction Liquified (%)")
plt.show()
plt.savefig("P_vs_liq.png")