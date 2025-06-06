import numpy as np
import matplotlib.pyplot as plt
from tube_geometry import TubeGeometry
from loads import LoadsCalculator
import aeroloads as aero
import math
from scipy.integrate import cumulative_trapezoid, simps

tube = TubeGeometry(12.252, 0.26, 0.005, [2, 4], [0.002, 0.001], [1, 2, 3, 4, 5, 6],
                    [0.026, 0.026, 0.026, 0.026, 0.026, 0.026], taper_ratio=0.5, true_if_steps_false_if_tapered=True)

calculator = LoadsCalculator('horizontal', 300)
calculator.thrust_loads()
calculator.engine_weight_loads()
shear_lift, shear_drag, moment_lift, moment_drag, normal_lift = calculator.aerodynamic_loads(
    lift=2.5 * aero.lift_gull_rh, drag=aero.drag_gull_rh)
calculator.weight_loads(2500)
shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()

tube_matrix = tube.get_tube_matrix(300)

# Extract arrays
thickness = tube_matrix[0]
Ix = tube_matrix[1]
J = tube_matrix[2]
radii_out = tube_matrix[3]
y_vals = tube_matrix[4]
area_enc = np.pi * radii_out**2


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

angle = get_angle(moment_x, Ix, 70 * 1e9, y_vals, 6.126)

def get_deflection(angle, y_values):
    deflection = cumulative_trapezoid(angle, y_values)
    return deflection

deflection = get_deflection(angle, y_vals[:-1])

plt.plot(y_vals[:-2], deflection)
plt.show()


