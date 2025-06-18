import numpy as np
from scipy.optimize import minimize_scalar
from loads import LoadsCalculator
from tube_skin import TubeGeometry
from stress_job_skin import StressCalculations
import aeroloads as aero
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

num_points = 100


def optimize_thickness_distribution(half_span, root_diameter, taper_ratio,
                                    stress_limits, flight_mode, load_factor):
    span_points = np.linspace(0, half_span, num_points)
    thicknesses = np.zeros(num_points)

    calculator = LoadsCalculator(flight_mode, thrust, num_points)
    calculator.thrust_loads()
    calculator.engine_weight_loads()
    calculator.aero_moment()
    calculator.aerodynamic_loads(lift=load_factor * aero.lift_gull_rh, drag=np.zeros(43, dtype=int))
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
                np.array([geometry[0][i]]),  # thickness
                np.array([geometry[1][i]]),  # inertia_Ix
                np.array([geometry[2][i]]),  # inertia_J
                np.array([geometry[3][i]]),  # radius_out
                np.array([geometry[4][i]]),  # y_value
                np.array([geometry[5][i]]),  # inertiax_skin
                np.array([geometry[6][i]])  # area_skin
            ]

            local_load = [
                np.array([load[0][i]]),  # shear_x
                np.array([load[1][i]]),  # shear_z
                np.array([load[2][i]]),  # moment_x
                np.array([load[3][i]]),  # moment_z
                np.array([load[4][i]]),  # torque
                np.array([load[5][i]])  # normal
            ]

            stress_calc = StressCalculations(cross_section, local_load, 360)
            stress = stress_calc.get_combined_stresses_tube()

            max_stress = stress[0][0]
            min_stress = stress[1][0]
            max_shear = stress[2][0]
            min_shear = stress[3][0]

            # Combined stress
            max_combined_stress = (max_stress ** 2 + 3 * max_shear ** 2) ** 0.5
            min_combined_stress = (min_stress ** 2 + 3 * min_shear ** 2) ** 0.5

            penalty = 0.0
            if not (stress_limits['sigma_min'] <= max_stress <= stress_limits['sigma_max']):
                penalty += (abs(max_stress - np.clip(max_stress, stress_limits['sigma_min'],
                                                     stress_limits['sigma_max']))) ** 2
            if not (stress_limits['sigma_min'] <= min_stress <= stress_limits['sigma_max']):
                penalty += (abs(min_stress - np.clip(min_stress, stress_limits['sigma_min'],
                                                     stress_limits['sigma_max']))) ** 2
            if not (stress_limits['tau_min'] <= max_shear <= stress_limits['tau_max']):
                penalty += (abs(max_shear - np.clip(max_shear, stress_limits['tau_min'],
                                                    stress_limits['tau_max']))) ** 2
            if not (stress_limits['tau_min'] <= min_shear <= stress_limits['tau_max']):
                penalty += (abs(min_shear - np.clip(min_shear, stress_limits['tau_min'],
                                                    stress_limits['tau_max']))) ** 2
            if not (stress_limits['sigma_min'] <= max_combined_stress <= stress_limits['sigma_max']):
                penalty += (abs(max_combined_stress - np.clip(max_combined_stress, stress_limits['sigma_min'],
                                                              stress_limits['sigma_max']))) ** 2
            if not (stress_limits['sigma_min'] <= min_combined_stress <= stress_limits['sigma_max']):
                penalty += (abs(min_combined_stress - np.clip(min_combined_stress, stress_limits['sigma_min'],
                                                              stress_limits['sigma_max']))) ** 2

            return t + 1e10 * penalty

        opt = minimize_scalar(optimize, bounds=(0.001, 0.1), method='bounded')
        thicknesses[i] = opt.x

    return span_points, thicknesses


half_span = 6.2412
root_diameter = 0.27
taper_ratio = 0.5
stress_limits = {
    'sigma_max': 340,
    'sigma_min': -340,
    'tau_max': 280,
    'tau_min': -280
}

flight_cases = [
                   # Horizontal flight cases
                   ('horizontal', load_factor, thickness)
                   for load_factor in [-1, 1, 2.5]
                   for thickness in [[935, 481], [0, 1000], [2000, 0]]
               ] + [
                   # Vertical flight cases
                   ('vertical', 1, thickness)
                   for thickness in [[9000, 5000], [18000, 0], [0, 10000]]
               ]

final_thickness = np.zeros(num_points)

for flight_mode, load_factor, thrust in flight_cases:
    print(f"Running case: mode={flight_mode}, n={load_factor}, thrust={thrust}")
    _, thickness_per_case = optimize_thickness_distribution(
        half_span, root_diameter, taper_ratio, stress_limits, flight_mode, load_factor
    )
    final_thickness = np.maximum(final_thickness, thickness_per_case)

span_pts = np.linspace(0, half_span, num_points)

max_sigma = np.full(num_points, -np.inf)
min_sigma = np.full(num_points, np.inf)
max_tau = np.full(num_points, -np.inf)
min_tau = np.full(num_points, np.inf)

for flight_mode, load_factor, thrust in flight_cases:
    print(f"Evaluating stress at max thickness for: mode={flight_mode}, n={load_factor}, thrust={thrust}")

    # Load calculation for the case
    calculator = LoadsCalculator(flight_mode, thrust, num_points)
    calculator.thrust_loads()
    calculator.engine_weight_loads()
    calculator.aerodynamic_loads(lift= load_factor * aero.lift_gull_rh, drag = np.zeros(43, dtype=int))
    calculator.weight_loads(2500)
    calculator.aero_moment()
    shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()
    load = [shear_x, shear_z, moment_x, moment_z, torque, normal]

    for i in range(num_points):
        # Overwrite geometry with final thickness
        # Evaluate stress using final_thickness at all span points
        tube = TubeGeometry(
            span=half_span * 2,
            root_diameter=root_diameter,
            root_thickness=final_thickness[i],
            taper_ratio=taper_ratio,
            true_if_steps_false_if_tapered=False
        )
        geometry = tube.get_tube_matrix(num_points)
        geometry[0][i] = final_thickness[i]

        cross_section = [
            np.array([geometry[0][i]]),  # thickness
            np.array([geometry[1][i]]),  # inertia_Ix
            np.array([geometry[2][i]]),  # inertia_J
            np.array([geometry[3][i]]),  # radius_out
            np.array([geometry[4][i]]),  # y_value
            np.array([geometry[5][i]]),  # inertiax_skin
            np.array([geometry[6][i]])   # area_skin
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

        sigma_max = stress[0][0]
        sigma_min = stress[1][0]
        tau_max = stress[2][0]
        tau_min = stress[3][0]

        max_sigma[i] = max(max_sigma[i], sigma_max)
        min_sigma[i] = min(min_sigma[i], sigma_min)
        max_tau[i] = max(max_tau[i], tau_max)
        min_tau[i] = min(min_tau[i], tau_min)

        max_combined = np.sqrt(max_sigma ** 2 + 3 * max_tau ** 2)
        min_combined = np.sqrt(min_sigma ** 2 + 3 * min_tau ** 2)

savethis = np.column_stack((final_thickness, span_pts))

np.savetxt('final_thickness.csv', savethis, delimiter=',', header='final thickness, span points', comments='')

fig, axs = plt.subplots(1, 2, figsize=(14, 6))

axs[0].plot(span_pts, max_sigma, label='Max Normal Stress')
axs[0].plot(span_pts, min_sigma, label='Min Normal Stress')
axs[0].plot(span_pts, max_tau, label='Max Shear Stress')
axs[0].plot(span_pts, min_tau, label='Min Shear Stress')

axs[0].plot(span_pts, max_combined, label='Max Combined Stress')
axs[0].plot(span_pts, min_combined, label='Min Combined Stress')

axs[0].axhline(stress_limits['sigma_max'], color='gray', linestyle='--', label='σ max limit')
axs[0].axhline(stress_limits['sigma_min'], color='gray', linestyle='--', label='σ min limit')
axs[0].axhline(stress_limits['tau_max'], color='red', linestyle='--', label='τ max limit')
axs[0].axhline(stress_limits['tau_min'], color='red', linestyle='--', label='τ min limit')

axs[0].set_xlabel("Spanwise position [m]")
axs[0].set_ylabel("Stress [Pa]")
axs[0].set_title("Worst-Case Stress Envelope")
axs[0].legend()
axs[0].grid(True)

axs[1].plot(span_pts, final_thickness, label='Final Min Thickness', color='tab:blue')
axs[1].set_xlabel("Spanwise position [m]")
axs[1].set_ylabel("Thickness [m]")
axs[1].set_title("Unified Thickness Distribution")
axs[1].legend()
axs[1].grid(True)

plt.tight_layout()
plt.show()

tube = TubeGeometry(
    span=half_span * 2,
    root_diameter=root_diameter,
    root_thickness=0.01,
    taper_ratio=taper_ratio,
    true_if_steps_false_if_tapered=False
)

fuselage_width = 1.8 / 2
idx = np.searchsorted(span_pts, fuselage_width)
final_thickness[:idx] = final_thickness[idx]

geometry = tube.get_tube_matrix(num_points)
r_outer = geometry[3]  # outer radius at each point
thicknesses = final_thickness
r_inner = r_outer - thicknesses

dy = span_pts[1] - span_pts[0]  # constant spacing
volume = np.sum(np.pi * (r_outer ** 2 - r_inner ** 2) * dy)  # integrate over span

print(f"Estimated wingbox volume: {volume:.4f} m³")