import numpy as np 
import matplotlib.pyplot as plt

Wwg = 700
Wfg	= 1100

C_r = 2.0
hinge_loc = 0.4 #bound to change
xcg_fg_wrtw = C_r * hinge_loc
xcg_wg_wrtw = 0.877 #from XFLR5

MAC = 1.5
Xcglemac = MAC*0.25


Xlemac = 2.3 
Xcgwg = Xlemac + xcg_wg_wrtw
Xcgfg = Xlemac + xcg_fg_wrtw
Xlemacs = []
Xlemacs.append(Xlemac)

n = 300

for i in range(n):
    Xcgwglemac = Xcgwg - Xlemac
    Xlemac = Xcgfg - Xcglemac + Wwg/Wfg *(Xcgwglemac - Xcglemac)
    Xlemacs.append(Xlemac)


plt.plot(Xlemacs)
#plt.show()
print(Xlemacs[-1])

xcgfg_actual = Xlemacs[-1] + xcg_fg_wrtw
print('Final cg:', xcgfg_actual)
print('Xcgfg:', Xcgfg)