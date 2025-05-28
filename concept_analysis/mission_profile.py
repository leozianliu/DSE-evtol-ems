import matplotlib.pyplot as plt
import numpy as np

# Define mission phases (time in minutes)
times = [0,1,2,20,21,22,22.16,32.16,33.16,51.16,52.16]

# Corresponding altitudes (in meters)
altitudes = [0,0,300,300,50,50,0,0,300,300,0]

# Smooth interpolation
time_smooth = np.linspace(times[0], times[-1], 500)
altitude_smooth = np.interp(time_smooth, times, altitudes)

# Plot
plt.figure(figsize=(12, 7))
plt.plot(time_smooth, altitude_smooth, label="eVTOL Ambulance Mission Profile", color="royalblue", linewidth=2)

# Font size settings
font_size = 16
plt.xlabel("Time [min]", fontsize=font_size)
plt.ylabel("Altitude [m]", fontsize=font_size)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
# plt.grid(True)  # Grid removed

# Optional: title and legend
# plt.title("Typical eVTOL Ambulance Mission Profile", fontsize=font_size + 2)
# plt.legend(fontsize=font_size)

plt.tight_layout()
plt.show()
