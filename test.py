from scipy.optimize import fsolve
def func(x):
    return (x**3) + (2*x) - 3
xs = fsolve(func, -3)
print(xs)