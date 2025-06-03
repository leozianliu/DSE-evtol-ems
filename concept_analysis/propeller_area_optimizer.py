from propeller_power_calculator import Propeller
import numpy as np


r1 = 3.33 
r2 = 2.353

a1 = np.pi * (r1/2)**2
a2 = np.pi * (r2/2)**2

mtow = 2200 * 9.81 / 2 
t1 = mtow * a1/(a1+a2)
t2 = mtow * a2/(a1+a2)

print(t1)
prop = Propeller(t1, 8, 3.33, 4500, 5)
prop.compute_maximum_RPM(use_maximum_RPM=True)
prop.design_propeller_geometry()
prop.summary_parameters(plot_chord_twist_distributions=True)
prop.compute_coefficients(print_coeff=True)
prop.change_flight_regime(40)
#prop.compute_coefficients()
prop.summary_parameters()




