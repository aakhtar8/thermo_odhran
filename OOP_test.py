# testing of the methods and classes during conversion of thermo library to OOP
from stream import Stream
import copy

from UNITS.Throttling_Valve import Throttling_Valve
# create mixture stream
name = {"CASRN": False, "ids": ["oxygen", "nitrogen"]}
air = Stream("air",Ids=name, Tcs=[154.4, 126.2], Pcs=[5.01E6, 3.356E6], Omegas=[0.022, 0.038], Mws=[32, 28.0134], Kijs=[[0, -0.0119], [-0.0119, 0]],composition=[0.21, 0.79], temp=300)
"""
other parameters of stream like pressure, flowrate, composition, interaction parameters
are either irrlevent or default value is used
"""
#air.show()
#print(290 == air.discharge_temp(air.delta_h(290)))      # true means both enthalpy and temperture functions are working correctly
# testing exchanger functions and classes

from UNITS.exchanger import Exchanger

air.temp = 290
hotair = copy.deepcopy(air)
air.temp = 280
coldair = copy.deepcopy(air)
hx1 = Exchanger("hx1", hotair, coldair)

# throttling valve functionality test
from UNITS.Throttling_Valve import Throttling_Valve
tv1 = Throttling_Valve(hx1.hot_discharge, discharge_pres=10000)

# liquid vapor seprator
from UNITS.Flash_Seprator import Flash_Sep
flash1 = Flash_Sep("flash1", coldair)




