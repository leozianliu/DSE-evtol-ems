import numpy as np
import matplotlib.pyplot as plt
import aeroloads as aero
from scipy.interpolate import interp1d

# def calculate_thrust_loads(engine_positions_y, engine_offsets_z, engine_thrusts, half_span, num_points=100):
#     assert len(engine_positions_y) == len(engine_thrusts) == len(engine_offsets_z), "Mismatch in engine input lengths"

#     y_points = np.linspace(0, half_span, num_points)
#     shear_force = np.zeros(num_points)
#     torque = np.zeros(num_points)
#     moment = np.zeros(num_points)

#     for pos, offset_z, thrust in zip(engine_positions_y, engine_offsets_z, engine_thrusts):
#         mask = y_points <= pos
#         shear_force[mask] += thrust
#         torque[mask] += thrust * offset_z

#     dy = y_points[1] - y_points[0]
#     moment[:-1] = np.flip(np.cumsum(np.flip(shear_force[1:] * dy)))


#     return y_points, shear_force, torque, moment

# def calculate_engine_weight_loads(engine_positions_y, engine_offsets_z, engine_weights, half_span, flight_mode, num_points=100):
#     assert len(engine_positions_y) == len(engine_weights) == len(engine_offsets_z), "Mismatch in engine input lengths"

#     y_points = np.linspace(0, half_span, num_points)
#     shear_force = np.zeros(num_points)
#     torque = np.zeros(num_points)
#     moment = np.zeros(num_points)

#     for pos, offset_z, weight in zip(engine_positions_y, engine_offsets_z, engine_weights):
#         mask = y_points <= pos

#         shear_force[mask] += weight
#         if flight_mode == 'vertical':
#             torque[mask] += weight * offset_z

#     dy = y_points[1] - y_points[0]
#     moment[:-1] = np.flip(np.cumsum(np.flip(shear_force[1:] * dy)))

#     return shear_force, torque, moment

# def aeroloads_func(flight_mode, half_span, num_points=100):
#     y_points = np.linspace(0, half_span, num_points)
#     moment_lift = np.zeros(num_points)
#     lift_gull_rh = np.zeros(num_points)

#     if flight_mode == 'horizontal':
#         y_interp = y_points
#         y = np.linspace(1.05, half_span, len(aero.lift_gull_rh))

#         lift_gull_interp_func = interp1d(y, aero.lift_gull_rh, bounds_error=False, fill_value='extrapolate')
#         lift_gull_rh = lift_gull_interp_func(y_interp)

#         value_at_1_05 = lift_gull_interp_func(1.05)
#         lift_gull_rh[y_interp < 1.05] = value_at_1_05

#         dy = y_points[1] - y_points[0]
#         moment_lift[:-1] = np.flip(np.cumsum(np.flip(lift_gull_rh[1:] * dy)))

#     elif flight_mode == 'vertical':
#         y_interp = y_points
#         y = np.linspace(1.05, half_span, len(aero.drag_gull_rh))

#         drag_gull_interp_func = interp1d(y, aero.drag_gull_rh, bounds_error=False, fill_value='extrapolate')
#         lift_gull_rh = drag_gull_interp_func(y_interp)

#         value_at_1_05 = drag_gull_interp_func(1.05)
#         lift_gull_rh[y_interp < 1.05] = value_at_1_05

#         dy = y_points[1] - y_points[0]
#         moment_lift[:-1] = np.flip(np.cumsum(np.flip(lift_gull_rh[1:] * dy)))

#     return lift_gull_rh, moment_lift

# if __name__ == "__main__":
#     print(y_points)
#     # flight_mode = 'horizontal'  # Change to 'vertical' for vertical flight

#     # engine_positions_y = np.array([2.41, 4.8])  # y-locations from fuselage center (m)
#     # engine_offsets_z = np.array([-0.5, 1.0])    # z-offsets (m)
#     # engine_thrusts = np.array([8000.0, 4000.0])  # N
#     # engine_weights = np.array([50, 60]) * 9.81   # N
#     # half_span = 6.0  # m

#     # shear_force_T, torque_T, moment_T = calculate_thrust_loads(engine_positions_y, engine_offsets_z, engine_thrusts, half_span)
#     # shear_force_W, torque_W, moment_W = calculate_engine_weight_loads(engine_positions_y, engine_offsets_z, engine_weights, half_span, flight_mode)

#     # aero_lift, moment_lift = aeroloads_func(flight_mode, half_span)

#     # shear_force = shear_force_W + aero_lift
#     # torque = torque_W + torque_T
#     # moment = moment_W + moment_T + moment_lift

#     # axis_label_shear = 'Z' if flight_mode == 'horizontal' else 'X'
#     # axis_label_torque = 'Y'
#     # flight_title = 'Horizontal' if flight_mode == 'horizontal' else 'Vertical'

#     # y = np.linspace(0, half_span, len(shear_force))
#     # plt.figure(figsize=(15, 5))
#     # plt.suptitle(f'Loads in {flight_title} Flight', fontsize=16)

#     # plt.subplot(1, 3, 1)
#     # plt.plot(y, shear_force, label='Total Shear (N)')
#     # plt.xlabel('Spanwise Position y (m)')
#     # plt.ylabel(f'Shear Force (N) in {axis_label_shear}-direction')
#     # plt.title('Shear Force Distribution')
#     # plt.grid()
#     # plt.legend()

#     # plt.subplot(1, 3, 2)
#     # plt.plot(y, torque, label='Total Torque (Nm)', color='orange')
#     # plt.xlabel('Spanwise Position y (m)')
#     # plt.ylabel(f'Torque around {axis_label_torque}-axis (Nm)')
#     # plt.title('Torque Distribution')
#     # plt.grid()
#     # plt.legend()

#     # plt.subplot(1, 3, 3)
#     # plt.plot(y, moment, label='Total Bending Moment (Nm)', color='green')
#     # plt.xlabel('Spanwise Position y (m)')
#     # plt.ylabel('Bending Moment (Nm)')
#     # plt.title('Bending Moment Distribution')
#     # plt.grid()
#     # plt.legend()

#     # plt.tight_layout()
#     # plt.show()
print(np.full(29, 500))