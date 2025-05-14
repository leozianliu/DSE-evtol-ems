import numpy as np
import matplotlib.pyplot as plt
from propeller_cruise_calculator import Propeller



rho = 1.225
A = np.linspace(1, 80)
vel = 55.5555
thrust = np.array([])
eff_real = np.array([]) 
eff_ideal = np.array([])
for i in range(len(A)):

    prop = Propeller(10000, vel, rho, A[i])
    thrust = np.append(thrust, prop.compute_thrust())
    prop.propeller_thrust = prop.compute_thrust()
    eff_ideal = np.append(eff_ideal, prop.compute_efficiency()[0])
    eff_real = np.append(eff_real, prop.compute_efficiency()[1])

plt.plot(A, eff_real)
plt.xlabel('Propeller Area')
plt.grid()
plt.show()