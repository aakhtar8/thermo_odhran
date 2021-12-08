"""
This file is used to test and counter check the validity and methods performance of the thermo library
PENGROBINSON EOS is used in these test cases
"""
from thermo import *
nitrogen = PR(Tc=154.6, Pc=15460000, omega=0.288, T=200, P=1E6)
print(nitrogen.raw_volumes)
print(nitrogen.phase)
print(nitrogen.more_stable_phase)
# handling the mixtures
eos = PRMIX(T=115.0, P=1E6, Tcs=[126.1, 190.6], Pcs=[33.94E5, 46.04E5], omegas=[0.04, 0.011], zs=[0.5, 0.5], kijs=[[0.0, 0.0289], [0.0289, 0.0]])


