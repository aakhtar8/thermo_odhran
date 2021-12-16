import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FLUIDS.propylene_glycol import PropyleneGlycol
from FLUIDS.air import Air
from UNITS.compressor import Compressor
from UNITS.heat_exchanger import HeatExchanger
from UNITS.expansion_valve import ExpansionValve
from tabulate import tabulate

air = Air() #setting up air instance
AtmT = 293 #Atmospheric Temperature
AtmP = 101325 #Atmospheric Pressure 
mair = 0.02 #Inlet air mass flowrate (kg/s)

c1 = Compressor(Air()) #setting up air compressor
c1P = 3000000 #outlet pressure from compressor
c1T = c1.calcOutTemperature(AtmT,AtmP,c1P) #outlet temperature from compressor calculated
c1T = 363 #actual outlet temperature with coolant 

aftcT = 303
condT = 85

ev1 = ExpansionValve(Air()) #setting up expansion valve 
exp1T = ev1.calcOutTemperature(condT,c1P,AtmP) #temperature after expansion
print("outT: ", exp1T)

exp1Y = ev1.calcLiquidQuality(condT,c1P,AtmP)
exp1X = ev1.calcGasQuality(condT,c1P,AtmP)
print("exp1Y: ",exp1Y)
print("exp1X: ",exp1X)

mReturn = mair*exp1X
cond1 = HeatExchanger(Air(),Air()) #setting up condenser
returnT = cond1.calcOutTemperature(mair,c1P,aftcT, condT,mReturn, AtmP, exp1T,5,0,1)
print("returnT: ", returnT)

therminolInT = 303
therminolOutT = 363
therminolHX1 = HeatExchanger(Air(),PropyleneGlycol())
mTherminol = therminolHX1.calcMassFlow2(mair,c1P, c1T, aftcT, AtmP, therminolInT, therminolOutT)
print("mTherminol: ", mTherminol)

mLiquidAir = mair*exp1Y #charge rate to vessel once min temperature reached
mLiquidTotal = 54 #120L = 96kg

print(tabulate([
                ['Inlet',AtmT,AtmP,mair, 0.0],
                ['After Compression/aftercooler', aftcT,c1P, mair, 0.0],
                ['After Condenser', condT, c1P, mair, 0.0],
                ['After Expansion Valve', exp1T, AtmP, mair, 0.0],
                ['Return Air', exp1T, 500000, mReturn, 0.0],
                ['Liquid Air', exp1T, 500000, mair-mReturn, exp1Y],
                ['Inlet Therminol (Propylene Glycol)',therminolInT, 200000,mTherminol, 0.0],
                ['Outlet Therminol', therminolOutT, 200000,mTherminol, 0.0],
                ],

                headers = ['Stage', 'Temperature (K)', 'Pressure (Pa)', 'Mass Flow (kg/s)', 'Liquid Fraction'], 
                tablefmt='orgtbl'))
print("Time to charge vessel to 54kg (min): ", (mLiquidTotal/mLiquidAir)/60)

bubbleT = air.Tsb(AtmP)
dewT = air.Tsd(AtmP)
print("B: ", bubbleT)
print("D: ", dewT)










