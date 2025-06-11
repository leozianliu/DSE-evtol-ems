import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#############  Requirement definition  #############
#  Note that the drag for the main landing gear consists two wheels

V = 55.56 # [m/s] Cruise velocity
M_tow  = 2200 # [Kg] Maximum take off weight
lg_fraction = 0.05 # [-] Weight fraction of the LG
C_L_alpha = 2*np.pi # [-]  
e = 0.8 # Oswald factor [-] 
rho = 1.225
Ar = 10

d_main = 0.42418 #[m]
w_main = 0.14732 #[m]
d_nose = 0.35306 #[m]
w_nose = 0.14478 #[m]

S = 14 # [m^2] Wing reference area
C_ds_main = 1.1 # [-] 
C_ds_nose = 1.1 # [-] 

#############################################

C_d_main = (d_main*w_main)/S*(C_ds_main)
C_d_nose = (d_nose*w_nose)/S*(C_ds_nose)/2

total_dd_ut = C_d_main + C_d_nose
D_ut = 0.5 * rho * V**2 * S * total_dd_ut


################# Retractable ################


delta_cl_ret= (M_tow * lg_fraction)/ (0.5*rho*V**2*S)
total_dd_t = delta_cl_ret**2/ (np.pi*Ar*e)

if total_dd_ut > total_dd_t :
    print('Use retractable landing gear!')
    print('Drag coefficient for retractable:', total_dd_t)
    print('Drag coefficient for unretractable:', total_dd_ut)
    print(D_ut)
else:
    print('Use tractable landing gear!')
    print('Drag coefficient for retractable:', total_dd_t)
    print('Drag coefficient for unretractable:', total_dd_ut)

#################### From this part, the size and C_ds are assumed to be the same for both nose wheel and main wheels. ######


d_comb = d_main
w_comb = w_main
C_ds_comb = np.linspace(0,1,100)
# print(C_ds_comb)

C_d_combined = (d_comb*w_comb)/S*(C_ds_comb)*1.5
plt.plot(C_ds_comb,C_d_combined)
plt.show()

