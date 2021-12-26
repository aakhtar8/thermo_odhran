import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FLUIDS.air import Air


'''
air = Air()
x = np.arange(80, 293,1)
y = [air.cp(_,101325) for _ in x]
plt.plot(x,y)
plt.show()
'''
air = Air()
cp1 = air.cp(293,101325)
print("cv1hi: ", cp1)
