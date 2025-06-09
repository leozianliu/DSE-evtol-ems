import numpy as np
import matplotlib.pyplot as plt
import itertools

from propeller_power_calculator import Propeller
from prop_optimizer import Optimizer

energies = []
thrust_setting = 1.1
max_load_factor = 1.1

optimizer = Optimizer(2600, max_load_factor, thrust_setting)
optimizer.compute_power_energy([2444.44, 2444.44], [64.44, 64.44], [2000, 2000], [5, 5])
print(2 * optimizer.total_energy)

fig, ax = plt.subplots(figsize=(8, 5))

# Plot with better line width and color
ax.plot(np.array(optimizer.time_arr) / 60, optimizer.power_arr_total / 2,
        color='tab:red', linewidth=2, label='Total Power')

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
print(optimizer.power_arr_total / 2)

np.savetxt("power_envelope.csv", optimizer.power_arr_total/2, delimiter=",", fmt='%.5f')