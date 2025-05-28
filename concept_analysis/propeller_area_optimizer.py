from propeller_power_calculator import Propeller
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from helper import calculate_total_noise

# Optimize areas

r1 = np.linspace(0.5, 2)
r2 = np.linspace(0.5, 2)
total_powers = []
iteration = 0
plot = 0

r1 = 3.33 
r2 = 2.353

a1 = np.pi * (r1/2)**2
a2 = np.pi * (r2/2)**2

mtow = 2700 * 9.81 / 2 / 10
t1 = mtow * a1/(a1+a2)
t2 = mtow * a2/(a1+a2)

j = []
eff = []
ct = []
cp = []
print(t1)
prop = Propeller(1, 785, 50, 1, 4500, 5)
#prop.compute_maximum_RPM(use_maximum_RPM=True)
prop.high_loading_design()
prop.summary_parameters(plot_chord_twist_distributions=False)
plt.plot(prop.total_circulation)
plt.show()

