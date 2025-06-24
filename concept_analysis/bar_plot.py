import numpy as np
import matplotlib.pyplot as plt

init = 131.0298372552631


# Data
labels = [
    'MTOW', 'Thrust', 'Velocity', 'RPM',
    'Number of Blades', 'Load Factor', 'Drag Coefficient', 'Delta Time', 'Cruise Altitude'
]
values = np.array([
    133.08927539734205, 145.50164263652033, 142.86953190936214, 126.4031084144356,
    145.76915167665143, 131.0298372552631, 137.00574343674825, 132.28327436753062, 138.21206073854472 / 141.73001043788196 * init
]) / init*100 - 100

# Create plot
plt.figure(figsize=(10, 6))
bars = plt.bar(labels, values, color="#88BD32", edgecolor='black')

# Add labels and title
plt.ylabel('Relative Change %', fontsize=16)
plt.title('Propeller Sensitivity Study', fontsize=14)
plt.xticks(rotation=45, fontsize=16)
plt.yticks(fontsize=10)

# Add values on top of bars
# for bar in bars:
#     height = bar.get_height()
#     plt.text(bar.get_x() + bar.get_width() / 2, height + 1,
#              f'{height:.2f}', ha='center', va='bottom', fontsize=10)

# Layout adjustments
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()