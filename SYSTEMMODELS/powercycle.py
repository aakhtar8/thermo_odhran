import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FLUIDS.air import Air
from UNITS.compressor import Compressor
from UNITS.heat_exchanger import HeatExchanger
from tabulate import tabulate
from UNITS.packedbed import PackedBed

liquidairmass = 54 #96kg of liquid air 
mair = 0.09 #discharge flow from liquid air vessel
dischargetime = liquidairmass/mair
print("discharge time: ", dischargetime)

p_lair = 3000000 #pressure of liquid air after pump
t_lair = 77 #temperature of liquid air after pump

pb1 = PackedBed(Air())
pb1.temperatureProfile(1.4,0.2,dischargetime,mair,303,t_lair)







