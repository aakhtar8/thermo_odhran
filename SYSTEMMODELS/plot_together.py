import matplotlib.pyplot as plt
import numpy as np

def plot_together(data, colors, *args, **kwargs):
    '''Data is a np array of dimension n, where n is the
     number of plots that will be ploted togheter. Every 
     n-th dimension have its first column as x and the 
     second column as y.'''

    for i in range(data.shape[2]):
        plt.plot(data[:,0,i], data[:,1,i],color=colors[i],*args,**kwargs)