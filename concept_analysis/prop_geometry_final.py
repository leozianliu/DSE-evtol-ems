import numpy as np
import matplotlib.pyplot as plt
import itertools

from propeller_power_calculator import Propeller
from prop_optimizer import Optimizer

T = np.linspace(1500, 15000, 8)
v = np.linspace(20, 160, 8)
rpm = np.linspace(2000, 9000, 10)

n_blades = [5]

optim = Optimizer(2600, 1.5, 2)

best_energy = np.infty

best_T = None
best_v = None
best_rpm = None
best_blades = None
all_energies = []
i = 0

for T_val, v_val, rpm_val, n_val in itertools.product(T, v, rpm, n_blades):
    i += 1
    if i % 100 == 0:
        print(i)

    try:
        optim.compute_power_energy([T_val, T_val], [v_val, v_val], [rpm_val, rpm_val], [n_val, n_val])
    except Exception as e:
        # Optional: print the error or log it for debugging
        print(f"Skipping due to error: {e}")
        continue  # Skip to the next combination
    
    if not any(optim.power_arr_big) < 0 and not any(optim.power_arr_small) < 0 and 0 < 2 * optim.total_energy < 300:
        all_energies.append(2 * optim.total_energy)
        if 0 < 2 * optim.total_energy < best_energy:
            best_energy = 2 * optim.total_energy
            best_T = T_val
            best_v = v_val
            best_rpm = rpm_val
            best_blades = n_val
        
    if i == 223:
        print(T_val, v_val, rpm_val, 2 * optim.total_energy)


print('-------------------------------------------------------------')
print("Best Energy (kWh):", best_energy)
print("Best Thrust (N):", best_T)
print("Best Velocity (m/s):", best_v)
print("Best RPM:", best_rpm)
print("Best Number of Blades:", best_blades)


plt.figure(figsize=(10, 5))
plt.plot(all_energies)       # Line plot
plt.scatter(range(len(all_energies)), all_energies, color='red', label='Test Points')  # Scatter on same axis
plt.xlabel("Model ID", fontsize=16)
plt.ylabel("Total Energy (kWh)", fontsize=16)
plt.title("Energy Consumption Across Parameter Sweep", fontsize=16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

