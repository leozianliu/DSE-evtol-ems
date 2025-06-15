import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from tube_skin import TubeGeometry
from loads import LoadsCalculator
import aeroloads as aero
from stress_job_skin import StressCalculations

optimized = np.loadtxt('final_thickness.csv', delimiter=',', skiprows=1)
thickness_ideal = optimized[:,0]
span = optimized[:,1]

plt.plot(span, thickness_ideal)

# Get the current Axes instance
ax = plt.gca()

ax.set_xlabel("Spanwise position [m]")
ax.set_ylabel("Thickness [m]")
ax.set_title("Unified Thickness Distribution")

ax.xaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_major_locator(MultipleLocator(0.0001))

ax.grid()
plt.show()



spanwise_change = [0, 1.5, 2.25 , 3.0, 3.75, 4.5, 5.25, 6.2412]
t_change = [0, 0, 0.0002, 0.0002, 0.0004, 0.0003, 0.0008, 0]

tube = TubeGeometry(6.2412 * 2, 0.27, 0.0031, spanwise_change, t_change, taper_ratio=0.5)

num_points = 300

geometry = tube.get_tube_matrix(num_points)


thickness = geometry[0]     # thickness
inertiaX = geometry[1]      # inertia_Ix
inertaiJ = geometry[2]      # inertia_J
radius = geometry[3]        # radius_out
y_vals = geometry[4]        # y_value
inertiaX_skin = geometry[5] # inertiax_skin
area_skin = geometry[6]     # area_skin

plt.plot(y_vals, thickness*1e3,label='Final thickness', color = 'blue')
plt.plot(span, thickness_ideal*1e3,label='Optimized min thickness', color = 'red')

plt.title("Thickness Distribution")
plt.xlabel("Spanwise location y (m)")
plt.ylabel("Thickness of the tube (mm)")
plt.legend()
plt.grid()
plt.show()