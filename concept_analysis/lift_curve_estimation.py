import numpy as np
import matplotlib.pyplot as plt

A_R = 10.28571429   # Aspect ratio
beta = (1-0.15**2)**0.5     # Prandtl-Glauert correction

################# This is airfoil specific ##################
kappa = (1.43737/(12.21/180*np.pi))/(2*np.pi)          # Ratio of two-dimensional lift curve slope to 2pi #############This is the input !
C_D_0 = 0.007
AOA_Critical = 17.7  # In degrees
#############################################################


AOA = (np.arange(-AOA_Critical,AOA_Critical,0.1))/180*np.pi  # In radian

C_L_alpha = 2*np.pi*A_R / (2+np.sqrt((A_R*beta/kappa)**2+4))
C_L = C_L_alpha * AOA

e = 1.78*(1-0.045*A_R**0.68)-0.64  # Oswald factor

C_D = 0.007 + C_L**2/(np.pi*A_R*e)
L_D = C_L/C_D

plt.figure(figsize=(15, 4))

plt.subplot(1, 4, 1)
plt.plot(AOA/np.pi*180, C_L, label='C_L', color='blue')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Lift Coefficient $C_L$')
plt.title('$C_L$ vs AOA')
plt.grid(True)

plt.subplot(1, 4, 2)
plt.plot(AOA/np.pi*180, C_D, label='C_D', color='red')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Drag Coefficient $C_D$')
plt.title('$C_D$ vs AOA')
plt.grid(True)

plt.subplot(1, 4, 3)
plt.plot(C_D, C_L, label='L/D', color='green')
plt.xlabel('C_D')
plt.ylabel('C_L')
plt.title('$C_L/C_D$')
plt.grid(True)

plt.subplot(1, 4, 4)
plt.plot(AOA/np.pi*180, L_D, label='L/D', color='green')
plt.xlabel('C_D')
plt.ylabel('C_L')
plt.title('$C_L/C_D$')
plt.grid(True)

plt.tight_layout()
plt.show()








