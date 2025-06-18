import numpy as np
from scipy.optimize import minimize_scalar
from loads import LoadsCalculator
from tube_geometry import TubeGeometry
from stress_job import StressCalculations
import aeroloads as aero
import matplotlib.pyplot as plt

num_points = 50

def optimize_thickness_distribution(half_span, root_diameter, taper_ratio,
                                    stress_limits, flight_mode, load_factor):

    span_points = np.linspace(0, half_span, num_points)
    thicknesses = np.zeros(num_points)

    calculator = LoadsCalculator(flight_mode, thrust, num_points)
    calculator.thrust_loads()
    calculator.engine_weight_loads()
    calculator.aerodynamic_loads(lift=1.1 * load_factor * aero.lift_gull_rh, drag = np.zeros(43, dtype=int))
    calculator.aero_moment()
    calculator.weight_loads(2500)

    shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()
    load = [shear_x, shear_z, moment_x, moment_z, torque, normal]

    for i, x in enumerate(span_points):
        def optimize(t):
            tube = TubeGeometry(
                span=half_span * 2,
                root_diameter=root_diameter,
                root_thickness=t,
                taper_ratio=taper_ratio,
                true_if_steps_false_if_tapered=False
            )

            geometry = tube.get_tube_matrix(num_points)

            cross_section = [
                np.array([geometry[0][i]]),
                np.array([geometry[1][i]]),
                np.array([geometry[2][i]]),
                np.array([geometry[3][i]]),
                np.array([geometry[4][i]])
            ]

            local_load = [
                np.array([load[0][i]]),
                np.array([load[1][i]]),
                np.array([load[2][i]]),
                np.array([load[3][i]]),
                np.array([load[4][i]]),
                np.array([load[5][i]])
            ]

            stress_calc = StressCalculations(cross_section, local_load, 360)
            stress = stress_calc.get_combined_stresses_tube()

            max_stress = stress[0][0]
            min_stress = stress[1][0]
            max_shear = stress[2][0]
            min_shear = stress[3][0]

            penalty = 0.0
            if not (stress_limits['sigma_min'] <= max_stress <= stress_limits['sigma_max']):
                penalty += (abs(max_stress - np.clip(max_stress, stress_limits['sigma_min'], stress_limits['sigma_max']))) ** 2
            if not (stress_limits['sigma_min'] <= min_stress <= stress_limits['sigma_max']):
                penalty += (abs(min_stress - np.clip(min_stress, stress_limits['sigma_min'], stress_limits['sigma_max']))) ** 2
            if not (stress_limits['tau_min'] <= max_shear <= stress_limits['tau_max']):
                penalty += (abs(max_shear - np.clip(max_shear, stress_limits['tau_min'], stress_limits['tau_max']))) ** 2
            if not (stress_limits['tau_min'] <= min_shear <= stress_limits['tau_max']):
                penalty += (abs(min_shear - np.clip(min_shear, stress_limits['tau_min'], stress_limits['tau_max']))) ** 2

            return t + 1e6 * penalty

        opt = minimize_scalar(optimize, bounds=(0.001, 0.05), method='bounded')
        thicknesses[i] = opt.x

    return span_points, thicknesses

half_span = 6.2412
root_diameter = 0.27
taper_ratio = 0.5
stress_limits = {
    'sigma_max': 400,
    'sigma_min': -400,
    'tau_max': 80,
    'tau_min': -80
}

flight_cases = [
    ('horizontal', lf, t)
    for lf in [-1, 1, 2.5]
    for t in [[935, 481], [2000, 0], [0, 1000]]
] + [
    ('vertical', 1, t)
    for t in [[9000, 5000], [18000, 0], [0, 10000]]
]

final_thickness = np.zeros(num_points)

for flight_mode, load_factor, thrust in flight_cases:
    print(f"Running case: mode={flight_mode}, n={load_factor}, thrust={thrust}")
    _, thickness_per_case = optimize_thickness_distribution(
        half_span, root_diameter, taper_ratio, stress_limits, flight_mode, load_factor
    )
    final_thickness = np.maximum(final_thickness, thickness_per_case)

span_pts = np.linspace(0, half_span, num_points)

# Visualize stress for each case using final thickness
fig, axs = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

for flight_mode, load_factor, thrust in flight_cases:
    calculator = LoadsCalculator(flight_mode, thrust, num_points)
    calculator.thrust_loads()
    calculator.engine_weight_loads()
    calculator.aerodynamic_loads(lift= load_factor * aero.lift_gull_rh, drag = np.zeros(43, dtype=int))
    calculator.weight_loads(2500)
    calculator.aero_moment()
    shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()
    load = [shear_x, shear_z, moment_x, moment_z, torque, normal]

    sigma_case = np.zeros(num_points)
    tau_case = np.zeros(num_points)

    for i in range(num_points):
        tube = TubeGeometry(
            span=half_span * 2,
            root_diameter=root_diameter,
            root_thickness=final_thickness[i],
            taper_ratio=taper_ratio,
            true_if_steps_false_if_tapered=False
        )
        geometry = tube.get_tube_matrix(num_points)
        cross_section = [
            np.array([final_thickness[i]]),
            np.array([geometry[1][i]]),
            np.array([geometry[2][i]]),
            np.array([geometry[3][i]]),
            np.array([geometry[4][i]])
        ]
        local_load = [
            np.array([load[0][i]]),
            np.array([load[1][i]]),
            np.array([load[2][i]]),
            np.array([load[3][i]]),
            np.array([load[4][i]]),
            np.array([load[5][i]])
        ]
        stress_calc = StressCalculations(cross_section, local_load, 360)
        stress = stress_calc.get_combined_stresses_tube()
        sigma_case[i] = stress[0][0]
        tau_case[i] = stress[2][0]

    label = f"{flight_mode}, n={load_factor}, T={thrust}"
    axs[0].plot(span_pts, sigma_case, label=label)
    axs[1].plot(span_pts, tau_case, label=label)

axs[0].axhline(stress_limits['sigma_max'], color='gray', linestyle='--', label='σ max limit')
axs[0].axhline(stress_limits['sigma_min'], color='gray', linestyle='--', label='σ min limit')
axs[1].axhline(stress_limits['tau_max'], color='red', linestyle='--', label='τ max limit')
axs[1].axhline(stress_limits['tau_min'], color='red', linestyle='--', label='τ min limit')

axs[0].set_ylabel("Normal Stress σ [Pa]")
axs[0].set_title("Normal Stress for All Load Cases")
axs[0].legend(fontsize=8)
axs[0].grid(True)

axs[1].set_xlabel("Spanwise position [m]")
axs[1].set_ylabel("Shear Stress τ [Pa]")
axs[1].set_title("Shear Stress for All Load Cases")
axs[1].legend(fontsize=8)
axs[1].grid(True)

plt.tight_layout()
plt.show()