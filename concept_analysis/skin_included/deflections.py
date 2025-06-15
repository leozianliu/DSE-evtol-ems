import numpy as np
import matplotlib.pyplot as plt
from tube_skin import TubeGeometry
from loads import LoadsCalculator
import aeroloads as aero
import math
from scipy.integrate import cumulative_trapezoid, simps

spanwise_change = [0, 1.5, 2.25, 3.0, 3.75, 4.5, 5.25, 6.2412]
t_change = [0, 0, 0.0002, 0.0002, 0.0004, 0.0003, 0.0008, 0]

tube = TubeGeometry(6.2412 * 2, 0.27, 0.0031, spanwise_change, t_change, taper_ratio=0.5)

calculator = LoadsCalculator('horizontal', [935, 481], 300)
calculator.thrust_loads()
calculator.engine_weight_loads()
calculator.weight_loads(0)
calculator.aero_moment()
shear_lift, shear_drag, moment_lift, moment_drag, normal_lift = calculator.aerodynamic_loads(lift= 2.5 * aero.lift_gull_rh, drag = np.zeros(43, dtype=int))

shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()
loads = [shear_x, shear_z, moment_x, moment_z, torque, normal]

# stress = StressCalculations(matrix, loads, 360)
# normal_max, normal_min, shear_max, shear_min, skin_top, skin_bottom = stress.get_combined_stresses_tube()


tube_matrix = tube.get_tube_matrix(300)

# Extract arrays
thickness = tube_matrix[0]
Ix = tube_matrix[1]
J = tube_matrix[2]
radii_out = tube_matrix[3]
y_vals = tube_matrix[4]
inertiax_skin = tube_matrix[5]
area_enc = np.pi * radii_out**2

Ix = Ix + inertiax_skin



def get_twist(torque, area_enc, thickness, radius, G):
    dtwist = torque * 2 * math.pi * radius / (4 * (area_enc ** 2) * G * thickness)
    return dtwist


def integrate_twist(dtwist, half_span, num_points):
    n = len(dtwist)
    dz = half_span / (n - 1)
    dz = np.full(num_points, dz)
    y = np.linspace(0, 6.126, len(dtwist))

    twist = cumulative_trapezoid(dtwist, y)

    return dz, twist


dz, twist = integrate_twist(get_twist(torque, area_enc, thickness, radii_out, 28 * 10 ** 9), 6.126, 300)

# print(thickness, twist)
# y = np.linspace(0, 6.126, len(twist))
#
# plt.plot(y, twist*180/np.pi)
# plt.xlabel("Spanwise location z (m)")
# plt.ylabel("Twist angle θ (deg)")
# plt.title("Twist distribution along span")
# plt.grid(True)
# plt.show()

def get_angle(moment, inertia, E_mod, y_values, half_span):
    n = len(y_values)
    dy = half_span/(n-1)
    y = np.linspace(0, half_span, n)

    moment_over_EI = moment/(inertia*E_mod)
    angle = cumulative_trapezoid(moment_over_EI, y_values)

    return angle

angle = get_angle(moment_x, Ix, 72.4 * 1e9, y_vals, 6.126)

def get_deflection(angle, y_values):
    deflection = cumulative_trapezoid(angle, y_values)
    return deflection

deflection = get_deflection(angle, y_vals[:-1])

plt.plot(y_vals[:-2], deflection)
plt.show()

