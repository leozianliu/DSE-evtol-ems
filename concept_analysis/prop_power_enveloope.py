import numpy as np
import matplotlib.pyplot as plt
import itertools

from propeller_power_calculator import Propeller
from prop_optimizer import Optimizer

energies = []
thrust_setting = 1.2 
max_load_factor = 1.5

optimizer = Optimizer(2600, max_load_factor, thrust_setting)

optimizer.compute_power_energy([5357.142, 5357.142], [140, 140], [3555.55, 3555.55], [5, 5])
print(2 * optimizer.total_energy)

from brokenaxes import brokenaxes
import matplotlib.pyplot as plt

# Create the broken x-axis figure
fig = plt.figure(figsize=(12, 6))
bax = brokenaxes(
    xlims=((0, 0.7), (15.4, 16)),  # Show only start and end
    hspace=0.05,
    despine=False
)

# Plot the data
bax.plot(np.array(optimizer.time_arr) / 60, (optimizer.power_arr_big + optimizer.power_arr_small) / 1.05, 'r-', label='Total Power', linewidth = 4)
#bax.plot(np.array(optimizer.time_arr) / 60, optimizer.power_arr_big, 'b-', label='Big Engines Power')
#bax.plot(np.array(optimizer.time_arr) / 60, optimizer.power_arr_small, 'g-', label='Small Engines Power')

# Axis and title settings
bax.set_ylabel('Total Power Usage [kW]', fontsize=14)
bax.set_xlabel('Flight Time [min]', fontsize=14, labelpad=20)
#bax.set_title('Power envelope for vertical flight and cruise', fontsize=16, weight='bold')

# Tick label sizes
bax.tick_params(axis='both', labelsize=12)

# Legend
#bax.legend(loc='upper center', fontsize=12)
bax.grid()

#plt.tight_layout()
plt.show()
