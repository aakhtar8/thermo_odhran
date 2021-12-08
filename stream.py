"""
Stream class will contain mixture(and eventually both component and mixture) conditions its present at
Following methods will be the part of stream class
1) Enthalpy
2) PH Flash
3) enthalpy change against specified temperture
4) temperture change required for given enthalpy change
"""
# imports
import copy
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL, FlashPureVLS
from thermo.interaction_parameters import IPDB
from thermo import *


class Stream:
    """
    Stream class to unify and create basis for all the methods
    PARAMETERS
    * for individual component specify the required properties as list with one element e.g. [3E5]
        name: name of the compound or mixture
        Ids: dictionary containing id method (CASRN or names) and list of IDs
        Ids: CAS_RN or names in the form of a list, for CASRN set CASRN=True
        Tcs: list of critical temperatures for all compoenent
        Pcs: list of critical pressures for all components
        Omegas: list of eccentric factors for all components
        Kijs: binary interaction parameters for the components involved, if unspecified then library will be searched
        user_defined_constants: boolean indicating whether chemical correlations, used by the library, are specified by user or defaults are used.
        temp: temperature of the stream in K
        pres: Pressure of the stream in Pa
        flowRate: dictionary with units of flowrate(key: units) and value (key: amount)
        composition: list of molar fraction of each component in mixture, use [1] for individual component
    """
    def __init__(self, name: str, Ids={"CASRN": False, "ids": []}, Tcs:list=[], Pcs:list=[], Omegas:list=[], Mws=[], Kijs:list=[], user_defined_constants: bool=False, temp: float=298, pres: float = 101325, flowRate={"units": "mol/s", "amount": 1.}, composition=[1.]):
        """
        create an instance of component and mixture to perform unit operation simulation
        """
        self.name = name
        self.material = self.define_pr(Ids=Ids, Tcs=Tcs, Pcs=Pcs, Omegas=Omegas, Mws=[], Kijs=Kijs, user_defined_constants=False)
        if len(Tcs) > 1:
            self.mixture = True
        else:
            self.mixture = False
        # units handler does not have a molecular weight functionlity modified    
        self.molar_flow = self.units_handler(flowRate)
        self.temp = temp
        self.pres = pres
        self.composition = composition
    # define function will create EquilibriumState object from thermo library using PENGROBINSON EOS
    def define_pr(self, Ids={"CASRN": False, "ids": []}, Tcs:list=[], Pcs:list=[], Mws=[], Omegas:list=[], Kijs:list=[], user_defined_constants=False, temp: float=298, pres: float = 101325):
        if len(Tcs) == 1:
            # given substance is a compound
            pr_compound = PR(Tc=Tcs[0], Pc=Pcs[0], omega=Omegas[0], T=temp, P=pres)
            constants = ChemicalConstantsPackage(Tcs=Tcs, Pcs=Pcs, omegas=Omegas, MWs=Mws)
            HeatCapacityGases = HeatCapacityGas(IDs = Ids["ids"][0])
            correlations = PropertyCorrelationsPackage(constants, HeatCapacityGases=HeatCapacityGases, skip_missing=True)
            eos_kwargs = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)

            liquid = CEOSLiquid(PR, HeatCapacityGases=HeatCapacityGases, eos_kwargs=eos_kwargs)
            gas = CEOSGas(PR, HeatCapacityGases=HeatCapacityGases, eos_kwargs=eos_kwargs)

            flasher = FlashPureVLS(constants, correlations, gas=gas, liquids=[liquid], solids=[])
            return flasher
        else:
            # use mixture methods
            # for now support is only added for componenet selection from library
            constants, properties = ChemicalConstantsPackage.from_IDs(Ids["ids"])
            # if interaction parameters are specified, use them
            if len(Kijs) == 0:
                # interaction parameters have not been specified by user
                # use library data
                kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
            else:
                kijs = Kijs
            eos_kwargs = dict(Tcs=Tcs, Pcs=Pcs, omegas=Omegas, kijs=kijs)

            gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
            liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
            
            flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
            return flasher


        
        
    # units handler
    def units_handler(self, stream: dict):
        """
        Analyze units of the stream flow and convert them to molar flow in standard (mol/s)
        * year means 365 days of flow since operational days can vary with project
        units allowed
        1) g/s, g/hr, g/day, g/year
        2) kg/s, kg/hr, kg/day, kg/year
        3) m-tonne/s, m-tonne/hr, m-tonne/day, m-tonne/year
        4) mol/s, mol/hr, mol/day, mol/year
        5) kmol/s, kmol/hr, kmol/day, kmol/year
        """

        values = list(stream.values())
        if values[0] == "mol/s":
            # check whether units are already converted to standard or not
            # if its standarized, then flow would have been replaced by molar_flow
            return values[1]
        else:   # units have not been converted to mol/s
            # perform conversion
            given_unit = values[0]
            [quantity, time] = given_unit.split("/")
            # convert time
            if time == "hr":
                values[1] /= 3600
                time = "s"
            elif time == "day":
                values[1] /= (24*3600)
                time = "s"
            elif time == "year":
                values[1] /= (24 * 3600 * 365)
                time = "s"
            elif time == "s":
                # nothing needs to the be done
                values[1] = values[1]
            else:
                # time units could not the infered in this case
                raise ValueError(f"Units type not infered for time units with given units {time}, valid option are s/hr/day/year")
                values[1] = None
                time = None

            # get the molar mass of the material of stream
            if isinstance(stream["material"], Mixture):
                molar_mass = stream["material"].MW()    # mixture method
            else:
                molar_mass = stream["material"].getMw() # Fluid method
            
            # check if its molar flow or mass flow, give solution in case of molar flow
            if quantity == "mol":
                # nothing needs to be done, specify molar_flow attribute
                values[1] = values[1]
            elif quantity == "kmol":
                values[1] = values[1] * 1000
            elif quantity == "g":
                values[1] /= molar_mass
            elif quantity == "kg":
                values[1] /= molar_mass * 1000
            elif quantity == "m-tonne":
                values[1] /= molar_mass * 1000 * 1000
            else:
                # given quantity could not be classifed
                raise ValueError(f"Given units could be parsed as {quantity} valid options are mol/kmol/g/kg/m-tonne")
                values[1] = None
                values[0] = None
            # return the converted value if conversion is sucessful
            if values[1] is not None:
                return values[1]
            else:   # stop the flow of program
                print("unit conversion could not be performed with given units, check logs")
                raise SystemExit
    
    # print stream properties
    def show(self):
        print(f"Stream specification for {self.name}")
        print(f"Temperature: {self.temp}")
        print(f"Pressure: {self.pres}")
        print(f"Molar Flow: {self.molar_flow}")
        print(f"Composition: {self.composition}")
    # enthalpy change function
    def delta_h(self, discharge_t):
        """
        this will act as wrapper to FlashVL and properly return the results based of structred input
        Parameter
        discharge_t: discharge temperture of the stream for which delta_h is to calculated(K)
        ----------------------
        RETURNS
        change in enthalpy for given conditions
        """
        if self.mixture:
            # user mixture methods
            flasher = self.material
            h_in = flasher.flash(T=self.temp, P=self.pres, zs=self.composition).H()
            h_out = flasher.flash(T=discharge_t, P=self.pres, zs=self.composition).H()
            dh = (h_out - h_in)*self.molar_flow
            return dh
        else:
            # use component methods
            # to be implemented latter
            pass
        
    # discharge temp function
    def discharge_temp(self, dh):
        """
        Returns the discharge temp against the given enthalpy change (dh) and feed conditions
        PARAMETERS
        dh: enthalpy change for which discharge temp is to calculated
        -----------------------------

        RETURNS
        discharge temperature in K
        """
        from scipy.optimize import fsolve

        def func(temp):
            return self.delta_h(temp) - dh
        discharge_temp = fsolve(func, self.temp)
        return discharge_temp[0]
        



