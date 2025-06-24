import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import pandas as pd
from scipy.interpolate import interp1d

from propeller_power_calculator import Propeller

_polar_df = pd.read_csv("concept_analysis/polar.txt", sep=r'\s+', skiprows=2,
                         names=["AOA","cl","cd","cm"])
cl_interp = interp1d(_polar_df.AOA, _polar_df.cl, kind='cubic', fill_value="extrapolate")
cd_interp = interp1d(_polar_df.AOA, _polar_df.cd, kind='cubic', fill_value="extrapolate")

def get_aero_coefficients(alpha):
    CL = cl_interp(alpha)
    CD = cd_interp(alpha)
    return CL, CD

radius1 = 3.19
radius2 = 2.286
a1 = np.pi * radius1 ** 2
a2 = np.pi * radius2 ** 2


prop1 = Propeller(2444.44, 64.44, radius1, 2000, 5)
prop2 = Propeller(2444.44, 64.44, radius2, 2000, 5)

prop1.design_propeller_geometry()
prop2.design_propeller_geometry()

prop1.compute_maximum_RPM(use_maximum_RPM=True)
prop2.compute_maximum_RPM(use_maximum_RPM=True)

max_power_big = 350 / 2
max_power_small = 223 / 2

max_rpm_big = 2700
max_rpm_small = 9600

velocity = np.linspace(10, 110)
altitude = np.linspace(0.001, 50000)

number_blades = 5
density = 1.225
g = 9.81

p_req1 = []
p_req2 = []

for v in velocity:
    preq = 0.5 * 1.225 * v**3 * 11.27 * 0.014 + (2600*g)**2 / (0.5 * 1.225 * v * 11.27 * np.pi * 0.9 * 13.1) + v * (0.09 * 0.5 * 1.225 * v**2 * 3.436449)
    prop1.change_flight_regime(v, preq/v/2 * a1/(a1 + a2))
    prop2.change_flight_regime(v, preq/v/2 * a2/(a1 + a2))
    p_req1.append(prop1.shaft_power/1000)
    p_req2.append(prop2.shaft_power/1000)


plt.plot(velocity, p_req1, color='red', label='Inbound Engine Power Required')
plt.plot(velocity, p_req2, color='green', label='Outbound Engine Power Required')
plt.axhline(y=max_power_big, color='r', linestyle='--', label='Inbound Engine Power Available')
plt.axhline(y=max_power_small, color='g', linestyle='--', label='Outbound Engine Power Available')
plt.grid()
plt.xlabel('Velocity [m/s]')
plt.ylabel('Power [kW]')
plt.title('Vehicle Power Study')
plt.legend()
plt.tight_layout()
plt.savefig('max_speed.png', dpi=800)
plt.show()
