import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from tube_skin import TubeGeometry
from loads import LoadsCalculator
import aeroloads as aero
from stress_job_skin import StressCalculations


def correlation_factor(thickness, radius):
    theta = (1 / 16) * np.sqrt(radius / thickness)
    gamma = 1.0 - 0.731 * (1 - np.exp(-theta))
    return gamma

def axial_compression(thickness, radius, e_modulus, poisson_ratio):
    # Formulas from NASA-TM-X-73306
    # Taking the plasticity eta = 1, means that we say it will buckle before plastically deforming
    theta = (1 / 16) * np.sqrt(radius / thickness)
    gamma = 1.0 - 0.731 * (1 - np.exp(-theta))
    sigma_cr = (gamma * e_modulus * thickness / radius) / np.sqrt(3 * (1 - poisson_ratio**2))

    return sigma_cr

def curvature_parameter(length, radius, thickness, poisson_ratio):
    Z = (length ** 2 /(radius*thickness)) * np.sqrt(1 - poisson_ratio)

    return Z

def shear_torsion(t, r, l, e_modulus, gamma):
    tau_cr = (0.747 * gamma**(3/4) * e_modulus) / ((r/t)**(5/4) * (l/r)**(1/2))

    return tau_cr



# gamma = correlation_factor(0.002, 0.097)
#
# sigma_cr = axial_compression(0.002, 0.097, 73*1e9, 0.33, gamma)
# print(sigma_cr/1e6)

spanwise_change = [0, 1.5, 2.25, 3.0, 3.75, 4.5, 5.25, 6.2412]
t_change = [0, 0, 0.0002, 0.0002, 0.0004, 0.0003, 0.0008, 0]

tube = TubeGeometry(6.2412 * 2, 0.27, 0.0031, spanwise_change, t_change, taper_ratio=0.5)

matrix = tube.get_tube_matrix(300)
thicknesses = matrix[0]
inertias_Ix = matrix[1]
inertias_J = matrix[2]
radii_out = matrix[3]
y_values = matrix[4]


gamma = correlation_factor(thicknesses, radii_out)
sigma_cr = axial_compression(thicknesses, radii_out, 72.4*1e9, 0.33)

Z1 = curvature_parameter(6/25, radii_out, thicknesses, 0.33)
Z2 = curvature_parameter(1, radii_out, thicknesses, 0.33)

half_span = 12.252/2
spanwise_step = [1.5, 2, 3.5, 5]

lengths = []
for y in y_values:
    if y < spanwise_step[0]:
        lengths.append(spanwise_step[0])
    if y > spanwise_step[-1]:
        lengths.append(half_span - spanwise_step[-1])
    else:
        for i in range(len(spanwise_step) - 1):
            if y > spanwise_step[i] and y < spanwise_step[i+1]:
                lengths.append(spanwise_step[i+1] - spanwise_step[i])

# plt.plot(y_values, lengths)
# plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))
# plt.gca().yaxis.set_major_locator(MultipleLocator(0.1))
# plt.grid()
# plt.show()


tau_cr = shear_torsion(thicknesses, radii_out, lengths, 72.4*1e9, gamma)



# plt.figure(1)
# plt.plot(y_values, gamma)
# plt.figure(2)
# plt.plot(y_values, sigma_cr/1e6)
# plt.show()

# plt.figure(1)
# plt.plot(y_values, Z2*gamma)
# plt.plot(y_values, 78*(radii_out/thicknesses)**2 * (1-0.33**2))
# plt.show()
#
# ## Combined buckling
# plt.plot(y_values, tau_cr/1e6)
# plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))
# plt.gca().yaxis.set_major_locator(MultipleLocator(10))
# plt.grid()
# plt.show()


# cross_section = tube.get_tube_matrix(300)

calculator = LoadsCalculator('vertical', [0, 10000], 300)
calculator.thrust_loads()
calculator.engine_weight_loads()
calculator.weight_loads(0)
calculator.aero_moment()
shear_lift, shear_drag, moment_lift, moment_drag, normal_lift = calculator.aerodynamic_loads(lift= 2.5 * aero.lift_gull_rh, drag = np.zeros(43, dtype=int))

shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()
loads = [shear_x, shear_z, moment_x, moment_z, torque, normal]

stress = StressCalculations(matrix, loads, 360)
normal_max, normal_min, shear_max, shear_min, skin_top, skin_bottom = stress.get_combined_stresses_tube()

rb = abs(normal_min / (sigma_cr/1e6))
rt = shear_max/(tau_cr/1e6)
#
# plt.plot(y_values, rb + rt**2)
# plt.plot(y_values, np.ones(len(y_values)), color='red', linestyle='--')
# plt.grid()
# plt.show()

# matrix = tube.get_tube_matrix(624)
# radii = matrix[3]
# thickness = matrix[0]
# dy = 6.2412/624
# density = 2830
# mass = 2 * np.sum(2 * np.pi * radii * dy * thickness) * density
# print(mass)
#
combined_loads_max = np.sqrt(normal_max**2 + 3 * shear_max**2)
combined_loads_min = (normal_min**2 + 3 * shear_min**2)**0.5

# plt.plot(y_values, combined_loads_max, color='blue')
# plt.plot(y_values, combined_loads_min, color='green')
# plt.plot(y_values, 340 * np.ones(len(y_values)), color='red', linestyle='--')
# plt.grid()
# plt.show()



##### Buckling of stiffened panels

def skin_critical(kc, E, poisson, b):
    # kc = 4 for simply supported
    t = 0.001
    skin_critical = ((np.pi ** 2 * kc * E )/(12 * (1 - poisson**2))) * (t/b)**2
    return skin_critical


def effective_width(C, E, poisson, crippling_stringer):
    # C = 4 for b/t<40
    # C = 6.98 for b/t>110
    t = 0.001
    we_2 = t * np.sqrt(C*np.pi**2 / (12 * (1-poisson**2))) * np.sqrt(E/crippling_stringer)
    return we_2



E = 72.4*1e9
poisson = 0.33
thickness = 0.001 # [m] -> 1 mm

# take 25x25x3 mm L-stringer, Al 2024-T6 with sigma yield = 345 MPa
crippling_stringer = 325*1e6  # [Pa]
stringer_area = 0.000141  # [m^2]
# nr_stringers = 5

# we_2 = effective_width(6.98, E, poisson, crippling_stringer)
# width_b = 5.88 * cross_section[3]
# stringer_spacing = width_b/(nr_stringers + 1)
# new_b = stringer_spacing - we_2     # new stringer spacing taking effective width into account
#
# area1 = stringer_area * nr_stringers    # area of stringers
# area2 = we_2 * nr_stringers * thickness # area of skin, stiffened by stringers
# area3 = width_b * thickness - area2     # area of skin without effect of stringers
#
# nominator_skin_cc = area1 * crippling_stringer + area2 * crippling_stringer + area3 * skin_critical(4.0, E, poisson, new_b)
# skin_crippling_combined = nominator_skin_cc/(area1+area2+area3)



# plt.plot(cross_section[4], skin_critical(4.0, E, poisson, stringer_spacing)/1e6)
# plt.plot(cross_section[4], skin_crippling_combined/1e6)
# plt.show()

##################################################
# Critical stress per interval for given nr of stringers
# critical = []
#
# for y in y_values:
#     nr_stringers = 6  # at the root
#     stringer_steps = [0, 1, 1, 1, 1, 1, 0, 0]   # 6, 5, 4, 3, 2, 1, 1
#     widths = np.array([0.135, 0.12847952, 0.12035243, 0.11222534, 0.10409825, 0.09597116, 0.08784407]) * 5.88
#     width_b = 0.135 * 5.88
#
#     for i in range(len(spanwise_change)):
#         if y > spanwise_change[i]:
#             width_b = widths[i]
#             nr_stringers = nr_stringers - stringer_steps[i]
#
#         else:
#             nr_stringers = nr_stringers
#             width_b = width_b
#
#     we_2 = effective_width(6.98, E, poisson, crippling_stringer)
#     stringer_spacing = width_b/(nr_stringers + 1)
#     new_b = stringer_spacing - we_2     # new stringer spacing taking effective width into account
#
#     area1 = stringer_area * nr_stringers    # area of stringers
#     area2 = we_2 * nr_stringers * thickness # area of skin, stiffened by stringers
#     area3 = width_b * thickness - area2     # area of skin without effect of stringers
#
#     nominator_skin_cc = area1 * crippling_stringer + area2 * crippling_stringer + area3 * skin_critical(4.0, E, poisson, new_b)
#     skin_crippling_combined = nominator_skin_cc/(area1+area2+area3)
#     # print(skin_crippling_combined/1e6)
#     critical.append(skin_crippling_combined/1e6)
#
# plt.plot(y_values, skin_bottom, label='Stress bottom skin', color = 'red')
# # plt.plot(y_values, skin_top, label='Stress bottom skin', color = 'red')
# plt.plot(y_values, np.array(critical))
# plt.grid()
# plt.show()





##################################################
# How many stringers do we need at each section
widths = np.array([0.135, 0.12847952, 0.12035243, 0.11222534, 0.10409825, 0.09597116, 0.08784407]) * 5.88
stress_by_stringer = []

for i in range(len(widths)):
    stress_on_interval = []
    width = widths[i]
    for j in range(10):
        nr_stringers = j

        we_2 = effective_width(6.98, E, poisson, crippling_stringer)
        stringer_spacing = width/(nr_stringers + 1)
        new_b = stringer_spacing - we_2     # new stringer spacing taking effective width into account

        area1 = stringer_area * nr_stringers  # area of stringers
        area2 = we_2 * nr_stringers * thickness  # area of skin, stiffened by stringers
        area3 = width * thickness - area2  # area of skin without effect of stringers

        nominator_skin_cc = area1 * crippling_stringer + area2 * crippling_stringer + area3 * skin_critical(4.0, E, poisson, new_b)
        skin_crippling_combined = nominator_skin_cc/(area1+area2+area3)

        stress_on_interval.append(skin_crippling_combined/1e6)

    stress_by_stringer.append(stress_on_interval)


fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(16, 6), sharey=True)

# Flatten the 2D array of axes for easy iteration
axes = axes.flatten()

for i in range(len(widths)):
    axes[i].plot(np.arange(0, 10, 1), stress_by_stringer[i])
    axes[i].set_title(f"Width b = {round(widths[i], 3)} m")
    axes[i].set_xlabel("Stringers")
    axes[i].grid()
    if i % 4 == 0:  # first column of each row
        axes[i].set_ylabel("Critical Stress")

# Hide the unused 8th subplot
if len(widths) < len(axes):
    for j in range(len(widths), len(axes)):
        fig.delaxes(axes[j])

fig.suptitle("Stress by Stringer Count for Different Widths", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.savefig("stress_by_stringer.pdf", format="pdf")  # Save as PDF
plt.show()
##################################################


