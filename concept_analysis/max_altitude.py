import numpy as np
from propeller_power_calculator import Propeller
import matplotlib.pyplot as plt

from ambiance import Atmosphere
h = np.linspace(0, 40000)

atmo = Atmosphere(h)

def drag(velocity_wing, velocity_fuselage, density):

    cd_fuselage = 0.1444908713 
    cd_wing = 0.035
    wing_surface = 14
    frontal_area = 3.436449
    return 0.5 * density * (velocity_fuselage ** 2 * cd_fuselage * frontal_area + velocity_wing ** 2 * cd_wing * wing_surface) * 1.5


radius1 = 3.19
radius2 = 2.286
a1 = np.pi * radius1 ** 2
a2 = np.pi * radius2 ** 2

t1 = a1 / (a1 + a2)
t2 = a2 / (a1 + a2)

prop1 = Propeller(2444.44, 64.44, radius1, 2000, 5)
prop2 = Propeller(2444.44, 64.44, radius2, 2000, 5)

v = 55.56
prop1.vi = 0

prop1.design_propeller_geometry()
prop2.design_propeller_geometry()

prop1.compute_maximum_RPM(use_maximum_RPM=True)
prop2.compute_maximum_RPM(use_maximum_RPM=True)

max_power_big = 350 
max_power_small = 223 
mtow = 2600 * 9.81

g = 9.81
cl = mtow / (0.5 * 1.225 * 55.56**2 * 11.27)
p1 = []
p2 = []

for i, altitude in enumerate(h):
    prop1.density = atmo.density[i]
    prop2.density = atmo.density[i]

    v = np.sqrt(2 * mtow / (cl * atmo.density[i] * 11.27))    
    
    preq = 0.5 * atmo.density[i] * v**3 * 11.27 * 0.014 + (2600*g)**2 / (0.5 * atmo.density[i] * v * 11.27 * np.pi * 0.9 * 13.1) + v * (0.09 * 0.5 * atmo.density[i] * v**2 * 3.436449)
    prop1.speed_of_sound = atmo.speed_of_sound[i]
    prop1.compute_maximum_RPM(use_maximum_RPM=True)

    prop1.change_flight_regime(v, preq/v)
    p1.append(preq / 2 * t1)
    p2.append(prop1.shaft_power / 2 * t1)


plt.plot(h/1000, np.array(p1)/1000, color='red', label='Inbound Engine Power Required')
plt.plot(h/1000, np.array(p2)/1000, color='green', label='Outbound Engine Power Required')
plt.axhline(y=max_power_big, color='r', linestyle='--', label='Inbound Engine Power Available')
plt.axhline(y=max_power_small, color='g', linestyle='--', label='Outbound Engine Power Available')
plt.legend()
plt.grid()
plt.savefig('max_altitude.png', dpi=800)
#plt.tight_layout()
plt.xlabel('Altitude [km]')
plt.ylabel('Power [kW]')
plt.show()


