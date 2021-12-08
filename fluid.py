import copy
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
from thermo.interaction_parameters import IPDB
from thermo import ChemicalConstantsPackage, PropertyCorrelationsPackage, PRMIX, SRKMIX, CEOSLiquid, CEOSGas, FlashPureVLS
"""
Fluid class will handle creation of mixtures or individual components
"""
class Fluid(self, consituentS, user_defined: boolean):
    """
    Creation of individual components or mixtures
    PARAMETERS
    consituentS: name of the compound/list of compounds
    user_defied: True or False based on whether or not constants and property package have been chosen by user 
    """
    if user_defined:
        # create constants and property object from critical properties provided by user
        constants = ChemicalConstantsPackage(Tcs=[512.5], Pcs=[8084000.0], omegas=[0.559], MWs=[32.04186], CASs=['67-56-1'])
        # this won't work as HeatCapacityGases object has not been created or defined
        correlations = PropertyCorrelationsPackage(constants, HeatCapacityGases=HeatCapacityGases, skip_missing=True)
    else:
        # use the builtin preferences for chemical properites and constants correlations
        constants, properties = ChemicalConstantsPackage.from_IDs(consituentS)

    # handling mixtures and indivudal components
    if len(consituentS) == 1:
        # only one component, use fluid properties

    else:
        # user mixture specific methods
        if user_defined:
            # specify binary interaction parameters for the system
            kijs = [[0, -0.0119], [-0.0119, 0]
        else:
            # use interaction parameters with specified method
            kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        # create mixture object
        eos_kwargs = dict(Tcs=[154.6, 126.2], Pcs=[50E5, 33.56E5], omegas=[0.021, 0.038], kijs=kijs)

        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)



    # get coefficients for the given property from database
    def get_coefficients(self, required_property = "CP_LIQUID"):
        # connect to database
        con = sqlite3.connect("fluid_properties.db")
        pointer = con.cursor()
        if self.is_common:
            # common or commerical name is specified as input
            querry = 'select * FROM ' + required_property + ' where (generic_name= (?))'
        else:
            # formula of the fluid is specified now querring as formula
            querry = 'select * FROM ' + required_property + ' where (formula= (?))'
        
        coefficients = pointer.execute(querry, (self.name,))
        coefficients = coefficients.fetchall()
        # in case the entered componenet is not present in the database inform the user
        # and give option to add the parameters for new property
        
        if len(coefficients) == 0:
            print(f"No record found against {self.name} for {required_property}")
            add_component = input("Do you want to add properties data for new component (y/n)")
            # if user agrees to add the data, add data for the required componenet in respecitve table
            if add_component == 'y':
                print(f"Please enter the relavent coefficents for {required_property} in following format")
                data = input("(name, chemical formula, **coefficents, min_temp, max_temp); please include parenthesis\n")
                cleaned_data = self.formater(data, required_property)
                self.add_todb(cleaned_data)
                # once the data is added, get the required property coffiencients
                coefficients = pointer.execute(querry, (self.name,))
                coefficients = coefficients.fetchall()
                con.close()
            else:
                return coefficients
        # coefficents found in the database, but are in string format
        # convert them to tuple and then return the value exculding component name and formula
        numeric_coefficients = [float(x) for x in coefficients[0][2:]]
        coefficients = coefficients[0][:2] + tuple(numeric_coefficients)
        return coefficients


    # get chemical formula of the species
    def getName(self):
        if self.is_common:
            # get formula from the database
            con = sqlite3.connect("fluid_properties.db")
            pointer = con.cursor()
            querry = 'select [formula] FROM GENERAL_PARAMETERS where (generic_name= (?))'  
            formula = pointer.execute(querry, (self.name,))
            formula = formula.fetchall()[0][0]
            con.close()
        else:
            formula = self.name
        return formula