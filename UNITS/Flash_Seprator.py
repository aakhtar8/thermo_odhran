# flash seprator will seprate the liquid and vapor phase based on the the flash conditions provided
from stream import Stream
import copy
class Flash_Sep:
    def __init__(self,name: str, stream:Stream):
        """
        Flash class will seprate liquid and vapor at given conditions
        PARAMETERS
        name: name of the flash seprator (string)
        stream: an instance of Stream class which needs to be split up
        -----------------------------------
        ATTRIBUTES
        name: name of the equipment
        feed: stream entering the flash separator
        vapor: vapor stream leaving the separator
        liquid: liquid stream leaving the separator
        """
        self.name = name
        self.feed = stream
        (self.liquid, self.vapor) = self.solve()

    def solve(self):
        print("="*50)
        print("Flash Seprator Feed:\n")
        self.feed.show()
        pt_flash = self.feed.material.flash(P=self.feed.pres, T=self.feed.temp, zs = self.feed.composition)
        liquid = copy.deepcopy(self.feed)
        vapor = copy.deepcopy(self.feed)
        vf = pt_flash.gas_beta          # molar fraction of gas in resultant flash
        # calculate liquid and vapor flowrates
        liquid.molar_flow = liquid.molar_flow*(1. - vf)
        vapor.molar_flow = vapor.molar_flow * vf
        if vf == 1.:
            liquid.composition = [None]*len(self.feed.composition)
            vapor.composition = self.feed.composition
        elif vf > 0.:
            liquid.composition = pt_flash.liquid_zs
            # calculate vapor fraction from this liquid fraction and flowrate
            mass_in = [x*self.feed.molar_flow for x in self.feed.composition]       # list of molar flow for each component
            liquid_out = [x*liquid.molar_flow for x in liquid.composition]     # list of moalr flow in liquid
            
            vapor_out = [x - y for x,y in zip(mass_in, liquid_out)]
            vapor.composition = [x/vapor.molar_flow for x in vapor_out]
        else:       # vf = 0
            liquid.composition = self.feed.composition
            vapor.composition = [None] * len(self.feed.composition)
        print("-"*50)
        print("liquid discharge:\n")
        liquid.show()
        print("vapor discharge:\n")
        vapor.show()
        return (liquid, vapor)


