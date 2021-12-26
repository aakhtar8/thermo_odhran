from plot_together import plot_together
import numpy as np
import matplotlib.pyplot as plt

data1=np.loadtxt('data1.txt') # loads data from simulation one
data2=np.loadtxt('data2.txt') # loads data from simulation two

#Ttime versus temperature plot for all positions
plt.plot(data1[:,0],data1[:,1], '--', color='b',label="position 1")
plt.plot(data1[:,0],data1[:,2], '--', color='r',label="position 2")
plt.plot(data1[:,0],data1[:,3], '--', color='g',label="position 3")
#plt.plot(data2[:,0],data2[:,2], ':', color='g',label="Single tank")
plt.ylabel(r'T (K)')
plt.xlabel(r'Time')
#plt.ylim((80, 91))
plt.grid('both')
plt.title("T_solid for all tanks (n=3) (Tin=90)")
plt.legend(loc='lower right')
plt.show()
