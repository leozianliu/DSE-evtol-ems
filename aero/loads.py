import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

gull_wing = pd.read_csv('aero/gull_wing.csv')
straight_wing = pd.read_csv('aero/straight_wing.csv')

lift_gull = gull_wing.iloc[:,0]
drag_gull = gull_wing.iloc[:,1]
lift_vd_gull = gull_wing.iloc[:,2]
drag_vd_gull = gull_wing.iloc[:,3]

lift_straight = straight_wing.iloc[:,0]
drag_straight = straight_wing.iloc[:,1]
lift_vd_straight = straight_wing.iloc[:,2]
drag_vd_straight = straight_wing.iloc[:,3]


array_size = lift_gull.size
b = 9.66
x = np.linspace(-b/2, b/2, array_size)

n_max = 2.5 

lift_gull_max = lift_gull*n_max
lift_vd_gull_max = lift_vd_gull*n_max

lift_straight_max = lift_straight*n_max
lift_vd_straight_max = lift_vd_straight*n_max

drag_gull_max = drag_gull*n_max
drag_vd_gull_max = drag_vd_gull*n_max

drag_straight_max = drag_straight*n_max
drag_vd_straight_max = drag_vd_straight*n_max


fig, axes = plt.subplots(2, 2, figsize=(10, 8))  # 2 rows, 2 columns

# Top-left
axes[0, 0].plot(x, lift_gull_max, label='Lift Gull Wing [Vc]')
axes[0, 0].plot(x, lift_vd_gull_max, label='Lift Gull Wing [Vd]')
axes[0, 0].axvline(0, color='black', linestyle='--')
axes[0, 0].set_title("Lift Gull Wing")
axes[0, 0].set_xlabel("x [m]")
axes[0, 0].set_ylabel("Lift [N]")
axes[0, 0].legend()

# Top-right
axes[0, 1].plot(x, lift_straight_max, label='Lift Straight Wing [Vc]')
axes[0, 1].plot(x, lift_vd_straight_max, label='Lift Straight Wing [Vd]')
axes[0, 1].axvline(0, color='black', linestyle='--')
axes[0, 1].set_title("Lift Straight Wing")
axes[0, 1].set_xlabel("x [m]")
axes[0, 1].set_ylabel("Lift [N]")
axes[0, 1].legend()

# Bottom-left
axes[1, 0].plot(x, drag_gull_max, label='Drag Gull Wing [Vc]')
axes[1, 0].plot(x, drag_vd_gull_max, label='Drag Gull Wing [Vd]')
axes[1, 0].axvline(0, color='black', linestyle='--')
axes[1, 0].set_title("Drag Gull Wing")
axes[1, 0].set_xlabel("x [m]")
axes[1, 0].set_ylabel("Drag [N]")
axes[1, 0].legend()

# Bottom-right
axes[1, 1].plot(x, drag_straight_max, label='Drag Straight Wing [Vc]')
axes[1, 1].plot(x, drag_vd_straight_max, label='Drag Straight Wing [Vd]')
axes[1, 1].axvline(0, color='black', linestyle='--')
axes[1, 1].set_title("Drag Straight Wing")
axes[1, 1].set_xlabel("x [m]")
axes[1, 1].set_ylabel("Drag [N]")
axes[1, 1].legend()

plt.tight_layout()
plt.show()

total_lift_gull = np.trapz(lift_gull_max, x)
total_lift_vd_gull = np.trapz(lift_vd_gull_max, x)
total_drag_gull = np.trapz(drag_gull_max, x)
total_drag_vd_gull = np.trapz(drag_vd_gull_max, x)

total_lift_straight = np.trapz(lift_straight_max, x)
total_lift_vd_straight = np.trapz(lift_vd_straight_max, x)
total_drag_straight = np.trapz(drag_straight_max, x)
total_drag_vd_straight = np.trapz(drag_vd_straight_max, x)

total_values_gull = np.array([["total_lift_gull", total_lift_gull],
                        ["total_lift_vd_gull", total_lift_vd_gull],
                        ["total_drag_gull", total_drag_gull],
                        ["total_drag_vd_gull", total_drag_vd_gull]])

total_values_straight = np.array([["total_lift_straight", total_lift_straight],
                        ["total_lift_vd_straight", total_lift_vd_straight],
                        ["total_drag_straight", total_drag_straight],
                        ["total_drag_vd_straight", total_drag_vd_straight]])

LD_gull = total_lift_gull/total_drag_gull
LD_vd_gull = total_lift_vd_gull/total_drag_vd_gull
LD_straight = total_lift_straight/total_drag_straight
LD_vd_straight = total_lift_vd_straight/total_drag_vd_straight

LD_values = np.array([["LD_gull", LD_gull],
                        ["LD_straight", LD_straight],
                        ["LD_vd_gull", LD_vd_gull],
                        ["LD_vd_straight", LD_vd_straight]])
