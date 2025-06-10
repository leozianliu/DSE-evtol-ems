import numpy as np
import matplotlib.pyplot as plt
from tube_geometry import TubeGeometry
from loads import LoadsCalculator
import aeroloads as aero

class StressCalculations:
    def __init__(self, cross_section, loads, num_points):
        # Array with cross section at each y-location
        self.thickness = cross_section[0]
        self.inertia_Ix = cross_section[1]
        self.inertia_J = cross_section[2]
        self.radius_out = cross_section[3]

        # Number of point across y-axis and cross section
        self.y_values = cross_section[4]
        self.num_points = num_points

        # All the loads
        self.normal_load = np.array(loads[5])
        self.shear_x = np.array(loads[0])
        self.shear_z = np.array(loads[1])
        self.torque = np.array(loads[4])
        self.moment_x = np.array(loads[2])
        self.moment_z = np.array(loads[3])


    def get_torsion_stress_tube(self):
        # Calculates max stress based on max radius, tau = T*rho/J
        tau_tube = self.torque * self.radius_out / self.inertia_J

        return tau_tube, max(tau_tube)

    def get_shear_stress_tube(self):
        # Calculates shear stress for points around the cross section and returns the max value for each y location
        # Returns values in Pa

        shear_stress_max = []
        shear_stress_min = []

        for i in range(len(self.y_values)):
            radius = self.radius_out[i]
            theta = np.linspace(0, 2 * np.pi, self.num_points, endpoint=False)
            x = radius * np.cos(theta)
            z = radius * np.sin(theta)
            Vx = self.shear_x[i]
            Vz = self.shear_z[i]
            t = self.thickness[i]
            Ix = self.inertia_Ix[i]
            shear_stress_max.append(max((-Vz/Ix*(t*radius**2 - t*radius*x)-Vx/Ix*(t*radius*z))/t))
            shear_stress_min.append(min((-Vz / Ix * (t * radius ** 2 - t * radius * x) - Vx / Ix * (t * radius * z))/t))

        return shear_stress_max, shear_stress_min


    def get_normal_stress_bending_tube(self):
        # Calculates normal stress for points around the cross section and returns the max value for each y location

        normal_stress_max = []
        normal_stress_min = []

        for i in range(len(self.y_values)):
            radius = self.radius_out[i]
            theta = np.linspace(0, 2 * np.pi, self.num_points, endpoint=False)
            x = radius * np.cos(theta)
            z = radius * np.sin(theta)
            Mx = self.moment_x[i]
            Mz = self.moment_z[i]
            normal_stress_max.append(max((Mx*z + Mz*x)/self.inertia_Ix[i]))
            normal_stress_min.append(min((Mx * z + Mz * x) / self.inertia_Ix[i]))

        return normal_stress_max, normal_stress_min

    def get_normal_stress_tube(self):
        # Calculates max normal stress due to axial loads based on cross section area
        area = self.thickness * (2 * np.pi * self.radius_out)
        sigma_tube = self.normal_load / area

        #volume = np.trapz(area, x = self.y_values)

        # Stress in MPa
        return sigma_tube, max(sigma_tube, key=abs)

    def get_combined_stresses_tube(self):
        sigma_tube, max_sigma1 = self.get_normal_stress_tube()
        sigma_bending_max, sigma_bending_min = self.get_normal_stress_bending_tube()
        total_normal_max = 1.5 * (sigma_tube + sigma_bending_max) /1e6
        total_normal_min = 1.5 * (sigma_tube + sigma_bending_min) /1e6


        tau_tube, max_tau1 = self.get_torsion_stress_tube()
        tau_shear_max, tau_shear_min = self.get_shear_stress_tube()
        # tau_tube = np.zeros(300)
        total_shear_max = 1.5 * (tau_tube + tau_shear_max) /1e6
        total_shear_min = 1.5 * (tau_tube + tau_shear_min) /1e6

        # Reaturns values in MPa
        return total_normal_max, total_normal_min, total_shear_max, total_shear_min

    # def get_tube_geometry():
    #     tube = TubeGeometry(12.252, 0.26, 0.007, [2, 4, 6], [0.002, 0.003, 0.000], [1, 2, 3, 4, 5, 6], [0.026, 0.026, 0.026, 0.026, 0.026, 0.026], taper_ratio = 0.5, true_if_steps_false_if_tapered = True)
    #     return tube


# loads = [np.linspace(0, 100, 300), np.linspace(0, 100, 300), np.linspace(0, 100, 300), np.linspace(0, 100, 300), np.linspace(0, 100, 300), np.linspace(0, 100, 300)]

### TUBE STEPS
tube = TubeGeometry(12.482, 0.26, 0.007, [2, 4, 6], [0.002, 0.003, 0.000], [1, 2, 3, 4, 5, 6], [0.026, 0.026, 0.026, 0.026, 0.026, 0.026], taper_ratio = 0.5, true_if_steps_false_if_tapered = True)
#tube = TubeGeometry(12.252, 0.25, 0.007, [1, 2, 3, 4, 5, 6], [0.001, 0.001, 0.001, 0.001, 0.001, 0.002], [1, 2, 3, 4, 5, 6], [0.022, 0.022, 0.022, 0.022, 0.022, 0.022], taper_ratio = 0.5, true_if_steps_false_if_tapered = True)

### TUBE TAPER
# tube = TubeGeometry(12.252, 0.4, 0.004, taper_ratio = 0.4, true_if_steps_false_if_tapered = False)
cross_section = tube.get_tube_matrix(300)

calculator = LoadsCalculator('horizontal', [7000.0, 4000.0], 300)
calculator.thrust_loads()
calculator.engine_weight_loads()
shear_lift, shear_drag, moment_lift, moment_drag, normal_lift = calculator.aerodynamic_loads(lift= 2.5 * aero.lift_gull_rh, drag = np.zeros(43, dtype=int))
calculator.weight_loads(2500)
calculator.aero_moment()
shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()
loads = [shear_x, shear_z, moment_x, moment_z, torque, normal]

stress = StressCalculations(cross_section, loads, 360)
normal_max, normal_min, shear_max, shear_min = stress.get_combined_stresses_tube()

# print(stress.get_shear_stress_tube()[0])
# print(stress.get_shear_stress_tube()[1])


# PLOTS

plots = True

if plots:

    fig, axs = plt.subplots(2, 2, figsize=(9, 5))  # 2 rows, 3 columns

# Flatten for easier indexing
    axs = axs.flatten()

# Plot 1: Normal max
    axs[0].plot(np.linspace(0, 6.126, 300), normal_max)
    axs[0].set_title("Max normal stress")
    axs[0].set_xlabel("Spanwise location y (m)")
    axs[0].set_ylabel("Normal stress (MPa)")

# Plot 2: Normal min
    axs[1].plot(np.linspace(0, 6.126, 300), normal_min)
    axs[1].set_title("Min normal stress")
    axs[1].set_xlabel("Spanwise location y (m)")
    axs[1].set_ylabel("Normal stress (MPa)")

# Plot 3: Shear max
    axs[2].plot(np.linspace(0, 6.126, 300), shear_max)
    axs[2].set_title("Max shear stress")
    axs[2].set_xlabel("Spanwise location y (m)")
    axs[2].set_ylabel("Shear stress (MPa)")

# Plot 4: Moment Z
    axs[3].plot(np.linspace(0, 6.126, 300), shear_min)
    axs[3].set_title("Min shear stress")
    axs[3].set_xlabel("Spanwise location y (m)")
    axs[3].set_ylabel("Shear stress (MPa)")


    plt.tight_layout()
    plt.show()