import numpy as np
from FLUIDS.fluid import Fluid
from FLUIDS.mixture import Mixture
import matplotlib.pyplot as plt

water = Fluid("h2o")
nitrogen = Fluid("n2")
oxygen = Fluid("o2")


Kij = np.array([[0, -0.0119], [-0.0119, 0]])
air = Mixture("Air",[oxygen, nitrogen], {"basis": "mole", "fraction": [0.21, 0.79]}, interaction=Kij)

x = range(100000, 1500001, 100000)

t_bubble = []
t_dew = []
for pressure in x:
    t_bubble.append(air.ideal_Tbubble(pres = pressure))      # bubble point temp at default press (101325)
    t_dew.append(air.ideal_Tdew(pres = pressure))             # dew point temp at default press (101325)

print(t_bubble)
print(t_dew)

plt.plot(x, t_bubble, 'b', label = 'T_BUBBLE')
plt.plot(x, t_dew, 'g', label = 'T_DEW')
plt.legend(loc='best')
plt.xlabel('Pressure (Pa)')
plt.ylabel('Temperature (K)')
plt.grid()
plt.show()