import numpy as np
from propeller_power_calculator import Propeller
import matplotlib.pyplot as plt

v = np.linspace(0, 230)


radius1 = 3.19
radius2 = 2.286
a1 = np.pi * radius1 ** 2
a2 = np.pi * radius2 ** 2

prop1 = Propeller(5357.142, 140, radius1, 3555.55, 5)
prop2 = Propeller(5357.142, 140, radius2, 3555.55, 5)

prop1.design_propeller_geometry()
prop2.design_propeller_geometry()


def drag(velocity_wing, velocity_fuselage):

    cd_fuselage = 0.09
    cd_wing = 0.035
    wing_surface = 14
    frontal_area = 3.436449

    return 0.5 * 1.225 * (velocity_fuselage ** 2 * cd_fuselage * frontal_area + velocity_wing ** 2 * cd_wing * wing_surface) * 1.5

prop1.compute_maximum_RPM(use_maximum_RPM=True)
prop1.vi = 0
eff = []
j= []
ct = []
cp = []

for vel in v:
    t = drag(vel + prop1.vi, vel)
    prop1.change_flight_regime(vel, t)
    eff.append(prop1.efficiency)
    j.append(prop1.advance_ratio)
    ct.append(prop1.thrust_coefficient)
    cp.append(prop1.power_coefficient)



plt.plot(j, eff, color='green', label='Efficiency')
#plt.plot(j, ct, color='blue', label='Thrust Coefficient')
#plt.plot(j, cp, color='red', label='Power Coefficient')
plt.grid()
plt.ylim(top=1)
plt.xlabel('Advance Ratio', fontsize=14)
plt.ylabel('Efficiency, Performance Coefficients', fontsize=14)
plt.title('Propeller Performance Diagram', fontsize=14)
plt.legend()
plt.tight_layout()
#plt.savefig('power_envelope.png', dpi=800)
plt.show()
