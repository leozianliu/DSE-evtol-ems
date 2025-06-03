from propeller_power_calculator import Propeller
import numpy as np
import matplotlib.pyplot as plt

cd_fuselage = 0.1444908713
cd_wing = 0.035
cd_total = cd_fuselage + cd_wing
frontal_area = 3.436449
mtow = 2200 * 9.80665
big_diameter = 3.33 
small_diameter = 2.353

lift_over_drag = 10

prop_area_big = np.pi * (big_diameter/2) ** 2
prop_area_small = np.pi * (small_diameter/2) ** 2

hover_thrust_big = mtow / 2  * prop_area_big / (prop_area_big + prop_area_small)
hover_thrust_small = mtow / 2  * prop_area_small / (prop_area_big + prop_area_small)

big_prop_1 = Propeller(hover_thrust_big, 5, big_diameter, 1, 5)
big_prop_2 = Propeller(hover_thrust_big / lift_over_drag, 55.56, big_diameter, 1, 5)

big_prop_1.compute_maximum_RPM(use_maximum_RPM=True)
big_prop_2.compute_maximum_RPM(use_maximum_RPM=True)

big_prop_1.design_propeller_geometry()
big_prop_2.design_propeller_geometry()