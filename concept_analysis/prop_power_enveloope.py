import numpy as np
import matplotlib.pyplot as plt
import itertools

from propeller_power_calculator import Propeller
from prop_optimizer import Optimizer

energies = []
thrust_setting = 1.1
max_load_factor = 1.2

optimizer = Optimizer(2 * 2600, max_load_factor, thrust_setting)
optimizer.compute_power_energy([2444.44, 2444.44], [64.44, 64.44], [2000, 2000], [5, 5])
print(2 * optimizer.total_energy)

fig, ax = plt.subplots(figsize=(8, 5))

# Plot with better line width and color
ax.plot(np.array(optimizer.time_arr) / 60, optimizer.power_arr_total / 2,
        color='tab:red', linewidth=2, label='Total Power')

ax.plot(np.array(optimizer.time_arr) / 60, optimizer.power_arr_big,
        color='tab:blue', linewidth=2, label='Big Engines Power')

ax.plot(np.array(optimizer.time_arr) / 60, optimizer.power_arr_small,
        color='tab:green', linewidth=2, label='Small Engines Power')

# Grid and formatting
ax.grid(True, linestyle='--', alpha=0.6)
ax.set_title('Power envelope for vertical flight and cruise', fontsize=14, weight='bold')
ax.set_xlabel('Flight Time [min]', fontsize=12)
ax.set_ylabel('Total Power Usage [kW]', fontsize=12)

# Improve tick formatting
ax.tick_params(axis='both', which='major', labelsize=10)

# Optional: Add legend
ax.legend()

# Tight layout for spacing
plt.tight_layout()

# Show the plot
plt.show()

#np.savetxt("power_envelope_small_engine.csv", optimizer.power_arr_small/2, delimiter=",", fmt='%.5f')
#np.savetxt("power_envelope_big_engine.csv", optimizer.power_arr_big/2, delimiter=",", fmt='%.5f')
print(np.max(optimizer.power_arr_big/2) / (2876.6607369518106 * np.pi/30), np.max(optimizer.power_arr_small/2) / (4014.237861275711 * np.pi/30))