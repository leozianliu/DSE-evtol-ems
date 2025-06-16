import numpy as np
import matplotlib.pyplot as plt
from tube_skin import TubeGeometry
from loads import LoadsCalculator
from stress_job_skin import StressCalculations
import aeroloads as aero

spanwise_change = [0, 1.5, 2.25, 3.0, 3.75, 4.5, 5.25, 6.2412]
t_change = [0, 0, 0.0002, 0.0002, 0.0004, 0.0007, 0.0004, 0]

tube = TubeGeometry(6.2412 * 2, 0.27, 0.0036, spanwise_change, t_change, taper_ratio=0.5)

num_points = 300

geometry = tube.get_tube_matrix(num_points)


thickness = geometry[0]     # thickness
inertiaX = geometry[1]      # inertia_Ix
inertaiJ = geometry[2]      # inertia_J
radius = geometry[3]        # radius_out
y_vals = geometry[4]        # y_value
inertiaX_skin = geometry[5] # inertiax_skin
area_skin = geometry[6]     # area_skin

# Find max widths at each interval
width = []

for i in range(len(spanwise_change) - 1):
    start = spanwise_change[i]
    end = spanwise_change[i + 1]

    # Get indices in this interval
    mask = (y_vals >= start) & (y_vals < end)

    # Get corresponding radii and find max
    if np.any(mask):  # make sure there's data in the interval
        max_r = np.max(radius[mask])
    else:
        max_r = np.nan  # or 0, depending on how you want to handle empty bins

    width.append(max_r)
# Convert result to numpy array
width = np.array(width) * 5.88



calculator = LoadsCalculator('horizontal', [935, 481], 300)
calculator.thrust_loads()
calculator.engine_weight_loads()
calculator.weight_loads(0)
calculator.aero_moment()
shear_lift, shear_drag, moment_lift, moment_drag, normal_lift = calculator.aerodynamic_loads(lift= 2.5 * aero.lift_gull_rh, drag = np.zeros(43, dtype=int))

shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()
loads = [shear_x, shear_z, moment_x, moment_z, torque, normal]

stress = StressCalculations(geometry, loads, 360)
normal_max, normal_min, shear_max, shear_min, skin_top, skin_bottom = stress.get_combined_stresses_tube()

# volume = 2 * np.pi * radius * thickness * 6.2412 * 2/num_points
# mass = np.sum(volume * 2780)
# print(mass)

# plt.plot(y_vals, skin_top, label='Stress top skin', color = 'blue')
plt.plot(y_vals, skin_bottom, label='Stress bottom skin', color = 'red')

plt.title("Normal stress skin")
plt.xlabel("Spanwise location y (m)")
plt.ylabel("Max stress (MPa)")
plt.legend()
plt.grid()
plt.show()





