import numpy as np
import matplotlib.pyplot as plt

# Constants
rho = 1.225
CD_wing = 0.035
CL = 0.8
S = 14.0
CD_fuselage = 0.144
S_fuselage = 3.43
mass = 2200.0
g = 9.81
W = mass * g
T_total = 22000

# Simulation parameters
V_target = 55.5
dt = 0.0001
max_time = 60
max_tilt_rate_deg = 45.0  # deg per second

# Initial conditions
V = 0
time = 0.0
thrust_angle_deg = 0  # Start with vertical thrust
P = 0.0
# Logging
time_history = []
V_history = []
L_history = []
D_history = []
P_history = []
T_vert_used_history = []
tilt_angle_history = []

while V < V_target and time < max_time:
    # Aerodynamic forces
    L = 0.5 * CL * rho * V**2 * S * np.cos(np.radians(thrust_angle_deg))
    D_wing = 0.5 * CD_wing * rho * V**2 * S
    D_fuselage = 0.5 * CD_fuselage * rho * V**2 * S_fuselage
    D_total = D_wing + D_fuselage

    # Solve required thrust angle to satisfy L + T_vertical = W
    T_vertical_required = W - L
    T_vertical_required = np.clip(T_vertical_required, 0.0, T_total)

    if T_vertical_required > T_total:
        raise ValueError("Thrust not sufficient to maintain vertical balance!")

    # Required angle from vertical lift
    required_angle_rad = np.arccos(T_vertical_required / T_total)
    required_angle_deg = np.degrees(required_angle_rad)

    # Smooth angle change (max 45 deg/s)
    delta_angle = required_angle_deg - thrust_angle_deg
    max_change = max_tilt_rate_deg * dt
    if abs(delta_angle) > max_change:
        thrust_angle_deg += np.sign(delta_angle) * max_change
    else:
        thrust_angle_deg = required_angle_deg

    # Final angle to use this step
    angle_rad = np.radians(thrust_angle_deg)
    T_vertical_used = T_total * np.cos(angle_rad)
    T_forward = T_total * np.sin(angle_rad)

    # Net acceleration and velocity update
    F_net = T_forward - D_total
    a = F_net / mass
    V += a * dt
    V = max(V, 0)

    P = T_forward * V

    # Logging
    time_history.append(time)
    V_history.append(V)
    L_history.append(L)
    D_history.append(D_total)
    T_vert_used_history.append(T_vertical_used)
    tilt_angle_history.append(thrust_angle_deg)
    P_history.append(P)

    time += dt

# Plotting
plt.figure(figsize=(12, 12))

plt.subplot(4, 1, 1)
plt.plot(time_history, V_history, label='Forward Speed (m/s)')
plt.axhline(V_target, linestyle='--', color='gray', label='Target Speed')
plt.ylabel("Speed (m/s)")
plt.grid()
plt.legend()

plt.subplot(4, 1, 2)
plt.plot(time_history, L_history, label='Aerodynamic Lift (N)')
plt.axhline(W, linestyle='--', color='gray', label='Weight')
plt.ylabel("Lift (N)")
plt.grid()
plt.legend()

plt.subplot(4, 1, 3)
plt.plot(time_history, T_vert_used_history, label='Vertical Thrust Used (N)')
plt.plot(time_history, D_history, '--', label='Total Drag (N)')
plt.ylabel("Forces (N)")
plt.grid()
plt.legend()

plt.subplot(4, 1, 4)
plt.plot(time_history, tilt_angle_history, label='Tilt Angle (deg)', color='orange')
plt.ylabel("Thrust Tilt (deg)")
plt.xlabel("Time (s)")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()

print(max(P_history))