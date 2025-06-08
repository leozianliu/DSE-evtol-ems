from xfoil import XFoil
from xfoil.test import naca0012
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

xf = XFoil() 

xf.airfoil = naca0012

xf.Re = 1e6
xf.max_iter = 40
a, cl, cd, cm = xf.aseq(-20, 20, 0.5)[0:4]

print('########################################')

# print(cl)
# print(a)

list_grad_cl_a = [len(cl)-1]

for i in range(len(cl)-1):
    grad_cl_a = (cl[i+1]-cl[i])/(0.5/180*np.pi)
    list_grad_cl_a.append(grad_cl_a)

plt.plot(a, cl)
plt.show()

A_R = 10.28571429   # Aspect ratio
beta = (1-0.15**2)**0.5     # Prandtl-Glauert correction

################# This is airfoil specific ##################
kappa = [len(list_grad_cl_a)]
for i in range(len(list_grad_cl_a)):
    kappa[i] = list_grad_cl_a[i]/(2*np.pi)          # Ratio of two-dimensional lift curve slope to 2pi #############This is the input !
C_D_0 = 0.007
AOA_Critical = np.max(a)  # In degrees

#############################################################

AOA_crit_rad = AOA_Critical/180*np.pi
AOA_max_rad = np.radians(90)
AOA = np.radians(a)

C_L_alpha = 2*np.pi*A_R / (2+np.sqrt((A_R*beta/kappa)**2+4))
# print(list_grad_cl_a)
# print(kappa)
# print(C_L_alpha)

C_L_max = np.max(cl)
C_L_parabolic = C_L_max * (1 - ((AOA - AOA_crit_rad) / (AOA_max_rad - AOA_crit_rad))**2)

C_L_0 = cl[int((len(cl))/2)]
a_0 = a[int((len(a)+1)/2)]

plt.plot(AOA,list_grad_cl_a)
    
C_L = np.where(AOA <= AOA_crit_rad, C_L_linear, C_L_parabolic)
C_L = np.where(AOA >= AOA_max_rad, 0.0, C_L)  # enforce zero at 90°

e = 1.78*(1-0.045*A_R**0.68)-0.64  # Oswald factor

C_D = 0.007 + C_L_linear**2/(np.pi*A_R*e)
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








