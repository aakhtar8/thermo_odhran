from plot_together import plot_together
import numpy as np
import matplotlib.pyplot as plt

data1=np.loadtxt('data1.txt') # loads data from simulation one
data2=np.loadtxt('data2.txt') # loads data from simulation two

#Position versus temperature for a range of times
snapshot = 1000  # snapshot time spacing: 1000=1s
style = ['-', '--', '-.', ':']
for i in range(len(style)):
    plt.plot(data1[:, i*snapshot], style[i],
             color='k',
             label="t="+str(i*snapshot/1000)+"s")
plt.ylabel(r'T (K)')
plt.xlabel(r'Time')
#plt.ylim((80, 91))
plt.grid('both')
plt.title("T_solid for all tanks (n=3) and the single tank (Tin=300)")
plt.legend(loc='lower right')
plt.show()
