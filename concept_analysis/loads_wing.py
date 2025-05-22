import numpy as np
import matplotlib.pyplot as plt
import aeroloads as aero
from scipy.interpolate import interp1d

def calculate_thrust_loads(engine_positions_y, engine_offsets_z, engine_thrusts, half_span, num_points=100):
    """
    Calculates shear force and torque due to thrust on a wing.

    Coordinate system:
    - y: spanwise (fuselage center to wing tip is positive)
    - z: vertical (positive downward)
    - x: out of screen

    Parameters:
    - engine_positions_y: spanwise y-positions of engines from fuselage center, in meters
    - engine_offsets_z: vertical z-offsets of engines from wing centerline in meters
    - engine_thrusts: thrust of each engine in Newtons
    - half_span: half of the total wingspan in meters
    - num_points: number of spanwise points to calculate

    Returns:
    - shear_force: Array of shear forces along the half-span (N)
    - torque: Array of torques around the spanwise axis (Nm)
    - moment: Array of cumulative bending moments along the half-span (Nm)
    """
    # Validate inputs
    assert len(engine_positions_y) == len(engine_thrusts) == len(engine_offsets_z), "Mismatch in engine input lengths"

    # Spanwise locations
    y_points = np.linspace(0, half_span, num_points)

    # Initialize output arrays
    shear_force = np.zeros(len(y_points))
    torque = np.zeros(len(y_points))
    moment = np.zeros(len(y_points))

    # Calculate shear force and torque
    for pos, offset_z, thrust in zip(engine_positions_y, engine_offsets_z, engine_thrusts):
        mask = y_points <= pos
        shear_force[mask] += thrust
        torque[mask] += thrust * offset_z

    # Calculate bending moment by integrating shear force along the span
    dy = y_points[1] - y_points[0]
    for i in range(len(moment) - 2, -1, -1):  # integrate from tip to root
        moment[i] = moment[i + 1] + shear_force[i + 1] * dy

    shear_force_T = shear_force
    torque_T = torque
    moment_T = moment 

    return shear_force_T, torque_T, moment_T

def calculate_engine_weight_loads(engine_positions_y, engine_offsets_z, engine_weights, half_span, flight_mode, num_points=100):
    """
    Calculates shear force, torque, and bending moment due to engine weight in either horizontal or vertical wing orientation.

    - engine_weights: List of weights of engines (N)
    - horizontal: Boolean, True for horizontal flight, False for vertical flight

    Returns:
    - shear_force_W: Shear force due to weight
    - torque_W: Torque around wing span axis
    - moment_W: Bending moment due to weight
    """
    assert len(engine_positions_y) == len(engine_weights) == len(engine_offsets_z), "Mismatch in engine input lengths"

    y_points = np.linspace(0, half_span, num_points)
    shear_force = np.zeros(len(y_points))
    torque = np.zeros(len(y_points))
    moment = np.zeros(len(y_points))

    for pos, offset_z, weight in zip(engine_positions_y, engine_offsets_z, engine_weights):
        mask = y_points <= pos

        if flight_mode == 'horizontal':
            shear_force[mask] += weight  # vertical weight contributes to vertical shear (z)  
            torque[mask] += 0  # No torque in z offset in horizontal flight (aligned with gravity)
        else:
            shear_force[mask] += weight  # weight acts in x direction now (shear is in x)
            torque[mask] += weight * offset_z  # torque due to offset from center of wing in z direction

    dy = y_points[1] - y_points[0]
    for i in range(len(moment) - 2, -1, -1):
        moment[i] = moment[i + 1] + shear_force[i + 1] * dy  #in horizontal flight in -x direction and in vertical flight in -z direction

    shear_force_W = shear_force 
    torque_W = torque 
    moment_W = moment

    return shear_force_W, torque_W, moment_W

def aeroloads_func(flight_mode, half_span, num_points = 100):
    y_points = np.linspace(0, half_span, num_points)
    moment_lift = np.zeros(len(y_points))

    if flight_mode == 'horizontal':
        y_interp = np.linspace(0, half_span, num_points)
        y = np.linspace(1.05, half_span, len(aero.lift_gull_rh))

        lift_gull_interp_func = interp1d(y, aero.lift_gull_rh, bounds_error=False, fill_value=0.0)

        lift_gull_rh = lift_gull_interp_func(y_interp)
        value_at_1_05 = lift_gull_interp_func(1.05)
        lift_gull_rh[y_interp < 1.05] = value_at_1_05

        dy = y_points[1] - y_points[0]
        for i in range(len(moment_lift) - 2, -1, -1):
            moment_lift[i] = moment_lift[i + 1] + lift_gull_rh[i + 1] * dy 
    else:
        lift_gull_rh = [0.0] * 100
    
    return lift_gull_rh, moment_lift


if __name__ == "__main__":
    flight_mode = 'horizontal'  # Change to 'vertical' for vertical flight
    # Inputs
    engine_positions_y = [2.41, 4.8]  # y-locations from fuselage center (m)
    engine_offsets_z = [-0.5, 1]  # z-offsets (m) from wing centerline, downward is positive
    engine_thrusts = [8000.0, 4000.0]  # N
    engine_weights = [50 * 9.81, 60 * 9.81] # N 
    half_span = 6.0  # m

    shear_force_T, torque_T, moment_T = calculate_thrust_loads(engine_positions_y, engine_offsets_z, engine_thrusts, half_span)
    shear_force_W, torque_W, moment_W = calculate_engine_weight_loads(engine_positions_y, engine_offsets_z, engine_weights, half_span, flight_mode)

    shear_force = shear_force_W + shear_force_T + aeroloads_func(flight_mode, half_span)[0]
    torque = torque_W + torque_T
    moment = moment_W + moment_T + aeroloads_func(flight_mode, half_span)[1]

    axis_label_shear = 'Z' if flight_mode == 'horizontal' else 'X'
    axis_label_torque = 'Y' if flight_mode == 'horizontal' else 'Y'
    flight_title = 'Horizontal' if flight_mode == 'horizontal' else 'Vertical'

    # Plotting
    y = np.linspace(0, half_span, len(shear_force))
    plt.figure(figsize=(15, 5))
    plt.suptitle(f'Loads in {flight_title} Flight', fontsize=16)

    plt.subplot(1, 3, 1)
    plt.plot(y, shear_force, label='Total Shear (N)')
    plt.xlabel('Spanwise Position y (m)')
    plt.ylabel(f'Shear Force (N) in {axis_label_shear}-direction')
    plt.title('Shear Force Distribution')
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(y, torque, label='Total Torque (Nm)', color='orange')
    plt.xlabel('Spanwise Position y (m)')
    plt.ylabel(f'Torque around {axis_label_torque}-axis (Nm)')
    plt.title('Torque Distribution')
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(y, moment, label='Total Bending Moment (Nm)', color='green')
    plt.xlabel('Spanwise Position y (m)')
    plt.ylabel('Bending Moment (Nm)')
    plt.title('Bending Moment Distribution')
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.show()
