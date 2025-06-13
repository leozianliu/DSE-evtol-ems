from propeller_power_calculator import Propeller

import matplotlib.pyplot as plt
import numpy as np

prop1 = Propeller(1500, 120, 3.19, 2000, 5)
prop2 = Propeller(1500, 120, 2.286, 2000, 5)

prop1.design_propeller_geometry()
prop1.compute_maximum_RPM(use_maximum_RPM=True)

prop2.design_propeller_geometry()
prop2.compute_maximum_RPM(use_maximum_RPM=True)

j = []
eff = []
ct = []
cp = []

for v in np.linspace(0, 200):
    prop1.change_flight_regime(v, 1500)
    j.append(prop1.advance_ratio)
    eff.append(prop1.efficiency)
    ct.append(prop1.thrust_coefficient)
    cp.append(prop1.power_coefficient)


plt.plot(j, eff, color='green', label="Efficiency")
plt.plot(j, ct, color='blue', label="Thrust Coefficient")
plt.plot(j, cp, color='red', label="Power Coefficient")
plt.legend()
plt.xlabel("Advance Ratio")
plt.ylabel("Efficiency and Coefficients")
plt.grid()
plt.tight_layout()
plt.show()
# def drag(velocity_wing, velocity_fuselage):

#     cd_fuselage = 0.1444908713
#     cd_wing = 0.035
#     wing_surface = 14
#     frontal_area = 3.436449

#     return 0.5 * 1.225 * (velocity_fuselage ** 2 * cd_fuselage * frontal_area + velocity_wing ** 2 * cd_wing * wing_surface) * 1.5

# required_power = []
# available_power = []
# dragall = []
# vel = np.linspace(0.001, 150)
# clcd = 11.3

# prop1.change_flight_regime(10, 10)
# prop2.change_flight_regime(10, 10)

# for v in vel:
#     prop1.change_flight_regime(v, drag(v + prop1.vi, v) / 2 * (2/3))
#     prop2.change_flight_regime(v, drag(v + prop2.vi, v) / 2 * (1/3))

#     force = (mtow / (0.5 * 1.225 * v**2 * 14)) / clcd * 0.5 * 1.225 * v**2 * 17.436449
#     dragall.append(force)
#     prop1.change_flight_regime(v, force)
#     prop2.change_flight_regime(v, force)
#     required_power.append(v * force/1000)
#     available_power.append( (prop1.shaft_power + prop2.shaft_power)/1000)

# plt.plot(vel, required_power, color='red', label='Power Required')
# plt.plot(vel, available_power, color='blue', label='Propeller Power')
# plt.title('Power Study')
# plt.ylabel('Power [kW]')
# plt.xlabel('Velocity [m/s]')
# plt.legend()
# plt.grid()
# plt.tight_layout()
# plt.show()