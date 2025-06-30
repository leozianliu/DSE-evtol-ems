import numpy as np
from propeller_power_calculator import Propeller
import matplotlib.pyplot as plt

# d = v * e / p(v)
e = (124 - 2 * 25) * 1000 * 3600 * 2/3 / 2 # Wh

g = 9.81
radius1 = 3.19
radius2 = 2.286
a1 = np.pi * radius1 ** 2
a2 = np.pi * radius2 ** 2

v = np.linspace(20, 200)
prop = Propeller(2444.44, 64.44, radius1, 2000, 5)
prop.design_propeller_geometry()
prop.compute_maximum_RPM(use_maximum_RPM=True)
d = []
for vel in v:
    preq = 0.5 * 1.225 * vel**3 * 11.27 * 0.014 + (2600*g)**2 / (0.5 * 1.225 * vel * 11.27 * np.pi * 0.9 * 13.1) + vel * (0.09 * 0.5 * 1.225 * vel**2 * 3.436449)
    prop.change_flight_regime(vel, preq/vel / 2 * a1/(a1 + a2))
    d.append(vel * e / prop.shaft_power/1000)

plt.plot(v, d, color = 'green')
plt.xlabel('Velocity [m/s]', fontsize=14)
plt.ylabel('Total Range [km]', fontsize=14)
plt.title('Maximum Range Diagram', fontsize=14)
plt.grid()
plt.tight_layout()
plt.show()
