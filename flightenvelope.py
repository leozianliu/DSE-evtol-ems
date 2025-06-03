import numpy as np
import matplotlib.pyplot as plt

W = 22000
clmax = 1.7
S = 14

n_max = 2.5 
n_min = -1.0
Vc = 200 / 3.6 
rho = 1.225 

Vs = np.sqrt(W / (0.5 * rho * S * clmax))
Va = Vs * np.sqrt(n_max)
Vd = Vc * 1.25 
Vneg = Vs * np.sqrt(abs(n_min))

print(f"Vs: {Vs:.2f} m/s, Va: {Va:.2f} m/s, Vd: {Vd:.2f} m/s, Vneg: {Vneg:.2f} m/s")



# === Create velocity array for plotting ===
V = np.linspace(0, Vd, 500)
n_stall_pos = (V / Vs)**2
n_stall_neg = -(V / Vs)**2

# Clip to limits
n_stall_pos = np.clip(n_stall_pos, 0, n_max)
n_stall_neg = np.clip(n_stall_neg, n_min, 0)

# === Plotting ===
plt.figure(figsize=(10, 6))
plt.plot(V, n_stall_pos, label='Positive Stall Limit Curve', color='green')
plt.plot(V, n_stall_neg, label='Negative Stall Limit Curve', color='purple')

# Vertical lines for key speeds
plt.axvline(x=Vs, color='blue', linestyle='--', label='Vs')
plt.axvline(x=Va, color='cyan', linestyle='--', label='Va')
plt.axvline(x=Vneg, color='magenta', linestyle='--', label='Vneg')
plt.axvline(x=Vc, color='orange', linestyle='--', label='Vc')
plt.plot([Vd, Vd], [n_min, n_max], color='lightcoral', linestyle='--')
plt.axvline(x=Vd, color='lightcoral', label='Vd')


# Horizontal lines for load limits

plt.axhline(y=n_max, color='blue', linestyle='--')
plt.axhline(y=n_min, color='salmon', linestyle='--')
plt.axhline(y=0, color='gray', linestyle='-')
plt.plot([Va, Vd], [n_max, n_max], color='blue', linestyle='-', label='n_max')
plt.plot([Vneg, Vd], [n_min, n_min], color='salmon', linestyle='-', label='n_min')


# Plot styling
plt.title("V-n Diagram (Flight Envelope)")
plt.xlabel("Velocity (m/s)")
plt.ylabel("Load Factor (n)")
plt.legend()
plt.xlim(0, Vd + 10)
plt.ylim(n_min - 0.5, n_max + 0.5)
plt.show()
plt.close()


