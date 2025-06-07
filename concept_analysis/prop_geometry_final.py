import numpy as np

from propeller_power_calculator import Propeller

big_prop = Propeller(2000, 75.555555555555565, 3.19, 2000, 5)
small_prop = Propeller(8222.222222222223, 75.555555555555565, 2.286, 3666.6666666666665, 5)


big_prop.design_propeller_geometry()
small_prop.design_propeller_geometry()

big_prop.compute_maximum_RPM()
small_prop.compute_maximum_RPM()

big_prop.chord *= 2

big_prop.summary_parameters(plot_chord_twist_distributions=True)
