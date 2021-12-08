"""
This script is designed to simulate the liquefection process of air
following units have be included in the this simulations
1) compressor
2) Heat Exchanger
3) throttling valve
4) Vapour-Liquid Seprator

method of creating a mixture and using the units is explained on the go as the script progresses
short documentation for each method is also provided in method description

defining and using a set of instances for the code intilization is given at the end to all methods
how everything is being, is also been tesed. 

* For Now it is assumed that air entering the loop is entering in compressed form and it's cooled to 300 K
** If an import is specific to a method, it explicitly called and used inside that method only
"""
# imports
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
from thermo.interaction_parameters import IPDB
import copy
# methods
"""
Heat Exchanger
Mainly HE has three methods
i) delta_h
ii) discharge_temp
iii) approach_exchanger

In order to model exchanger, one needs to call only "approach_exchanger" method. Rest of methods act as the auxilary functions for this
"""
# # enthalpy change function
def delta_h(flasher, conditions: dict, discharge_t):
    """
    this will act as wrapper to FlashVL and properly return the results based of structred input
    Parameter
    flasher: a FlashVL Object for given system
    conditions: dict of the shape {"temp": T, "pres": P, "composition":[xa, xb], "flow rate": mol/s}
    discharge_t: discharge temp in K, for which enthalpy change is to computed
    ----------------------
    RETURNS
    change in enthalpy for given conditions
    """
    h_in = flasher.flash(T=conditions["temp"], P=conditions["pres"], zs=conditions["composition"]).H()
    h_out = flasher.flash(T=discharge_t, P=conditions["pres"], zs=conditions["composition"]).H()
    dh = (h_out - h_in)*conditions["flow rate"]
    return dh

# discharge temp function
def discharge_temp(flasher, conditions: dict, dh):
    """
    Returns the discharge temp against the given enthalpy change (dh) and feed conditions
    PARAMETERS
    flasher: a FlashVL Object for given system
    conditions: dict of the shape {"temp": T, "pres": P, "composition":[xa, xb], "flow rate": mol/s}
    dh: enthalpy change for which discharge temp is to calculated
    -----------------------------
    
    RETURNS
    discharge temperature in K
    """
    from scipy.optimize import fsolve
    
    def func(temp):
        return delta_h(flasher, conditions, temp) - dh
    discharge_temp = fsolve(func, conditions["temp"])
    return discharge_temp[0]


# Exchanger integration model
def approach_exchanger(hotConditions: dict, coldConditions: dict, flasher, approach_temp=5, flow="counter-current"):

    """
    This will calculate the discharge temperture for hot and cold mixtures given feed tempertures and approach_temp
    --------------
    Parameters
    hotConditions: dict containing temp, press and composition for hot fluid
                {"temp": T, "pres": P, "composition":[xa, xb], "flow rate": mol/s}
    coldConditions: dict, temp, press and composition for cold fluid
                {"temp": T, "pres": P, "composition":[xa, xb], "flow rate": mol/s}
    flash: an instanace of FlashVL object for the given set of components
    optional approach_temp: minimum temperature difference allowed for heat transfer
    optional flow: flow arrangement inside exchanger, "co-current" or "counter-current"

    RETURNS
    ----------------------------
    (hotDischargeTemp, coldDischargeTemp, exchanger_duty) tuple with updated tempertures
    * only implmenting it for mixture as this is our use case
    """
    print("="*50)
    print("Exchanger feed:\n")
    print("Hot Fluid: \n", hotConditions)
    print("Cold Fluid: \n", coldConditions)


    # checking input for validation
    if (hotConditions["temp"] - coldConditions["temp"]) <= approach_temp:
        # heat transfer is not possible, ask user to change the inputs and rerun
        print("given conditions cannot satisfy minimum approach temperture")
        print("please update input streams or minimum approach temperture")
        raise ValueError("improper temperture specifications")
    else:
        # tempertures are valid
         
        if flow == "counter-current":
            # assume hot side discharge temp
            hot_dischargeT = coldConditions["temp"] + approach_temp
            # assume cold side discharge temp
            cold_dischargeT = hotConditions["temp"] - approach_temp
        elif flow == "co-current":
            # has to be done iteratively
            pass
        # calculate enthalpy change for hot and cold
        dh_hot = delta_h(flasher, hotConditions, hot_dischargeT)    # w
        dh_cold = delta_h(flasher, coldConditions, cold_dischargeT)  # w
        # compare the two changes and adjust temp so the enthalpy changes are equal to the smaller value
        if abs(dh_hot) < abs(dh_cold):
            # cold stream will not reach minimum approach temp
            cold_dischargeT = discharge_temp(flasher, coldConditions, -dh_hot)
        else:
            # hot stream will not meet approach conditions
            hot_dischargeT = discharge_temp(flasher, hotConditions, -dh_cold)

    hot_discharge = copy.deepcopy(hotConditions)
    hot_discharge["temp"] = hot_dischargeT
    cold_discharge = copy.deepcopy(coldConditions)
    cold_discharge["temp"] = cold_dischargeT
    print("\n----------------------------------------------------")
    print("Exchanger Discharge:\n")
    print("Hot Discharge Temperture: ", hot_dischargeT)
    print("Cold Discharge Temperture: ", cold_dischargeT)
    print("Exchanger Duty: ", abs(dh_hot))
    return (hot_discharge, cold_discharge, abs(dh_hot))

# model for throttling valve
def throttle_valve(feed_conditions: dict, discharge_pres, flasher):
    """
    for given discharge pressure throttle valve will provide discharge temperture of the stream
    """
    print("="*50)
    print("\nValve Feed: \n\n", feed_conditions)
    pt_flash = flasher.flash(P=feed_conditions["pres"], T=feed_conditions["temp"], zs = feed_conditions["composition"])
    enthalpy = pt_flash.H()
    throttle_T = flasher.flash(P=discharge_pres, H=enthalpy, zs = feed_conditions["composition"]).T
    discharge = copy.deepcopy(feed_conditions)
    discharge["temp"] =throttle_T
    discharge["pres"] = discharge_pres
    print("\n----------------------------------------------------")
    print("Valve Discharge: \n\n", discharge)
    return discharge

# mode for vapor-liquid seprator
def vl_split(feed_conditions: dict, flasher):
    """
    this will flash the mixture at feed conditions and return two condition dictionaries
    1) liquid conditions
    2) vapor conditions
    flow rate of feed will be split based on VF ratio obtained as result of flash
    ---------------------------------------
    RETURNS
    two dictionaries representing vapor and liquid phase and gaseous phase
    """
    print("="*50)
    print("\nFlash Feed: \n", feed_conditions)
    pt_flash = flasher.flash(P=feed_conditions["pres"], T=feed_conditions["temp"], zs = feed_conditions["composition"])
    vf = pt_flash.VF
    # if all vapors(vf=1) then no need to compute flash, replace outlet vapor stream with the feed stream
    # reverse case of all lquid flashing.
    if vf == 0.:
        # all the feed is in liquid
        liquid_composition = pt_flash.zs
        gas_composition = None
        liquid_flowrate = feed_conditions["flow rate"]
        gas_flowrate = None 
    elif vf < 1.:
        # feed is the mixtur of liquid mixture
        liquid_composition = pt_flash.liquid0.zs
        gas_composition = pt_flash.gas.zs
        liquid_flowrate = feed_conditions["flow rate"] * (1 - vf)
        gas_flowrate = feed_conditions["flow rate"] * vf
    elif vf == 1.:
        # all the feed is in vapor
        gas_composition = pt_flash.zs
        liquid_composition = None
        liquid_flowrate = None
        gas_flowrate = feed_conditions["flow rate"]
    # spliting up the feed
    liquid = copy.deepcopy(feed_conditions)
    gas = copy.deepcopy(feed_conditions)

    liquid["flow rate"] = liquid_flowrate
    liquid["composition"] = liquid_composition
    gas["flow rate"] = gas_flowrate
    gas["composition"] = gas_composition
    
    # showing and returning results
    print("\n----------------------------------------------------")
    print("\nFlash Liquid: \n", liquid)
    print("Flash vapors: \n", gas)
    return (liquid, gas)

"""
Defining the system and conditions dictionary
"""
# each stream is defined as dictionary
air_hot = {
    "temp": 290.,
    "pres": 3000000.,
    "composition": [0.21, 0.79],
    "flow rate": 1.0
}
air_cold ={
    "temp": 150.,
    "pres": 101325.,
    "composition": [0.21, 0.79],
    "flow rate": 1.0
}

# system is defined with constraints and application here
kijs = [[0, -0.0119], [-0.0119, 0]]         # Binary interaction parameters for nitrogen and oxygen


# creating air object as a mixture of oxygen and nitrogen
constants, properties = ChemicalConstantsPackage.from_IDs(['oxygen', 'nitrogen'])

kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
# using critical properties from reference hx work
eos_kwargs = dict(Tcs=[154.6, 126.2], Pcs=[50E5, 33.56E5], omegas=[0.021, 0.038], kijs=kijs)

gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)

liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)


#air_cold = throttle_valve(air_hot, 101325, flasher)
# checking eqiupment behavior
exchanger_cold_discharge = []
valve_feed_temp = []
valve_discharge_temp = []
liquid_flowrate = []
# convergence loop
n = 1       # counter for iterations
while(n <= 50):
    print("ITERATION: ", n)   
    (hotAir_discharge, coldAir_discharge, hx_duty) = approach_exchanger(air_hot, air_cold, flasher, approach_temp=1)
    exchanger_cold_discharge.append(coldAir_discharge["temp"])
    valve_feed_temp.append(hotAir_discharge["temp"])
    valveOutAir = throttle_valve(hotAir_discharge, 101325, flasher)
    valve_discharge_temp.append(valveOutAir["temp"])
    (liquid_air, vapor_air) = vl_split(valveOutAir, flasher)
    liquid_flowrate.append(liquid_air["flow rate"])
    air_cold = copy.deepcopy(vapor_air)
    n += 1
#if n == 50:
#    print("convergence loop did not converge in specified iterations")
# observing results
import matplotlib.pyplot as plt
import numpy as np
valve_feed_temp = np.array(valve_feed_temp)
valve_discharge_temp = np.array(valve_discharge_temp)
exchanger_cold_discharge = np.array(exchanger_cold_discharge)
temp_difference = valve_feed_temp - valve_discharge_temp
for i in range(len(liquid_flowrate)):
    if liquid_flowrate[i] is None:
        liquid_flowrate[i] = 0.0
fraction_liquified = np.array(liquid_flowrate)*100/air_hot["flow rate"]
x = np.arange(1, n)




plt.plot(x, valve_discharge_temp, label="Discharge Temperture")
plt.plot(x, exchanger_cold_discharge, label="Cold air discharge temperutre")
#plt.plot(x, temp_difference, label="temperture drop across valve")
#
plt.title("Temperture Change across TV")
plt.xlabel("Iterations")
plt.ylabel("Temperature Differrence (K)")
plt.legend()
plt.show()

plt.plot(x, fraction_liquified, label="Fraction Liquified")
plt.title("Fraction Liqeufied")
plt.xlabel("Iterations")
plt.ylabel("Liquid Fraction (%)")
plt.show()






