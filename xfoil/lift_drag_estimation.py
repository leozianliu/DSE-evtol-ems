from xfoil import XFoil
from xfoil.test import naca0012
import numpy as np

xf = XFoil()

xf.airfoil = naca0012

xf.Re = 1e6
xf.max_iter = 40
a, cl, cd, cm = xf.aseq(-20, 20, 0.5)[0:4]

import matplotlib.pyplot as plt
plt.plot(a, cl)
plt.show()