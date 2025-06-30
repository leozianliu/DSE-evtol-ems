import numpy as np
from propeller_power_calculator import Propeller


d = 0.254
v = 15.87
p = 68.77
rpm = 6519
rho = 1.225
b = 2
t = 3.5

p = Propeller(t, v, d, rpm, b)
p.design_propeller_geometry()
p.change_flight_regime(v, t)
print(p.shaft_power)
p.summary_parameters(plot_chord_twist_distributions=True)