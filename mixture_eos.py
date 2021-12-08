"""
This file will contain handling of the mixture using thermo library
"""
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
from thermo.interaction_parameters import IPDB
kijs = [[0, -0.0119], [-0.0119, 0]]         # Binary interaction parameters for nitrogen and oxygen

# creating air object as a mixture of oxygen and nitrogen
constants, properties = ChemicalConstantsPackage.from_IDs(['oxygen', 'nitrogen'])
#air = PRMIX(T=280, P=1E6, Tcs=[154.6, 126.2], Pcs=[50E5, 33.56E5], omegas=[0.021, 0.038], zs=[0.21, 0.79], kijs=kijs)
#print(air.H())
#print(type(air))
eos_kwargs = dict(Tcs=[154.6, 126.2], Pcs=[50E5, 33.56E5], omegas=[0.021, 0.038], kijs=kijs)

# provide initial estimates for liquid and gas
gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)

liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
# provide composition going into the flash
zs = [0.21, 0.79]
# perform a PT flash
PT = flasher.flash(T=108, P=1e6, zs=zs)
print(type(PT))
print("vapor to feed ratio ", PT.VF)
print("vapor composition ", PT.gas.zs)
#print("liquid composition ", PT.liquid0.zs)
print("enthalpy is ", PT.H())

# enthalpy change function
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




"""
Heat exchanger section of the code. this will calculate temperture of the mixture based on inputs provided
--------------------------
Included functions
1) approach_exchanger()
"""
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
    (hotDischarge, coldDischarge) tuple with updated tempertures
    * only implmenting it for mixture as this is our use case
    """

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
        # cold stream will achieve minimum discharge temp
        hot_dischargeT = discharge_temp(flasher, hotConditions, -dh_cold)
    else:
        # hot stream will achieve minimum approach temp
        cold_dischargeT = discharge_temp(flasher, coldConditions, -dh_hot)
    return (hot_dischargeT, cold_dischargeT, dh_hot)


print("="*20)
air_hot = {
    "temp": 310,
    "pres": 1013250,
    "composition": [0.21, 0.79],
    "flow rate": 1
}
air_cold = {
    "temp": 290,
    "pres": 1013250,
    "composition": [0.21, 0.79],
    "flow rate": 1
}

print("Rechecking acutal results for enthalpy change and required temp functions")
air_dh = delta_h(flasher, air_cold, 300)
print("air dh is ", air_dh)
air_cold["temp"] = 300
air_temp = discharge_temp(flasher, air_cold, -air_dh)
print("against dh air temperture from the previous discharge is ", air_temp)

print("testing heat exchanger functionality")
hx_results = approach_exchanger(air_hot, air_cold, flasher)
print("heat exchanger results")
print(f"hot stream: discharge temp: {hx_results[0]}\nDuty: {hx_results[2]}")
print(f"cold stream: discharge temp: {hx_results[1]}\nDuty: {hx_results[2]}")


def throttle_valve(feed_conditions: dict, discharge_pres, flasher):
    """
    for given discharge pressure throttle valve will provide discharge temperture of the stream
    """
    pt_flash = flasher.flash(P=feed_conditions["pres"], T=feed_conditions["temp"], zs = feed_conditions["composition"])
    enthalpy = pt_flash.H()
    throttle_T = flasher.flash(P=discharge_pres, H=enthalpy, zs = feed_conditions["composition"]).T
    return throttle_T

throttle_T = throttle_valve(air_cold, 1e5, flasher)
print("1 bar throttling gives ",throttle_T)

print("Capturing the whole scenario")









