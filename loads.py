import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

tc = 0.12
Dy = 26 * tc
clmax = 1.61

rho = 1.225 #kg/m^3
S = 14.14 
W = 21000 #N


def load(V):
    q = 0.5 * rho * V**2
    n = q * clmax / (W/S)
    return n 

n_cruise = load(55.6)
n_Vd = load(55.6 * 1.05)
n_climb = load(300/60)

df = pd.read_csv('lift.csv')
lift = df.iloc[:,0]
lift_size = lift.size
x = np.linspace(-9.66/2, 9.66/2, lift_size)

n_loads = np.array([n_cruise, n_Vd, n_climb])
n_max = np.max(n_loads)
lift_max_load = lift * n_max


plt.plot(x, lift, label='Lift')
plt.plot(x, lift_max_load, label='Lift max load')
plt.show()