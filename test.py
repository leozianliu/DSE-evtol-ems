from xfoil import XFoil
from xfoil.test import naca0012
import numpy as np
import matplotlib.pyplot as plt

# === XFOIL SETUP ===
xf = XFoil()
xf.airfoil = naca0012
xf.Re = 1e6
xf.max_iter = 40

# AOA sweep
a_deg, cl_2d, cd, cm = xf.aseq(-20, 20, 0.5)[0:4]
a_rad = np.radians(a_deg)

# === Geometry and correction factors ===
A_R = 10.28571429
beta = np.sqrt(1 - 0.15**2)  # compressibility correction

# === Compute local slope dCL/dalpha (kappa) ===
d_alpha = np.diff(a_rad)
d_cl = np.diff(cl_2d)
slope = d_cl / d_alpha
print(slope)

print('#####################')

# Pad the slope array to match a_rad
# kappa = np.concatenate(([slope[0]], slope))/(2*np.pi)
kappa = np.concatenate(([slope[0]], slope))/(2*np.pi)
print(kappa)
print('#####################')


# === Compute C_L_alpha (3D lift slope) ===
CL_alpha = (2 * np.pi * A_R) / (2 + np.sqrt((A_R * beta / kappa)**2 + 4))
print(CL_alpha)

# === Reference CL_0 at AOA = 0 ===
idx_0 = np.argmin(np.abs(a_rad))  # index of AOA ≈ 0
CL_0 = cl_2d[idx_0]

# === Final Lift Curve ===
CL_3D = CL_alpha * a_rad + CL_0

# === Plot ===
plt.figure(figsize=(8, 6))
plt.plot(a_deg, cl_2d, label='2D $C_L$ from XFOIL', color='blue')
plt.plot(a_deg, CL_3D, label='Corrected 3D $C_L$', color='red', linestyle='--')

plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Lift Coefficient $C_L$')
plt.title('Corrected 3D Lift Curve from XFOIL')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()