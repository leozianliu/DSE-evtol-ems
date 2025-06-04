import numpy as np 
import matplotlib.pyplot as plt

Wwg = 700
Wfg	= 1100

C_r = 2.0
#want the cg location of the fuselage to be right under the hinge
hinge_loc = 0.4 #bound to change
xcg_fg_wrtw = C_r * hinge_loc
#xcg_wg_wrtw = 0.877 #from XFLR5 

#cg of the wing with respect to the leading edge of the root chord 
xcg_wg_wrtw = 0.995

MAC = 1.5
Xcglemac = MAC*0.25

Xleroot = 1.5
Xcgwg = Xleroot + xcg_wg_wrtw
Xcgfg = Xleroot + xcg_fg_wrtw




'''
Xlemacs = []
Xlemacs.append(Xlemac)
n = 300

for i in range(n):
    Xcgwglemac = Xcgwg - Xlemac
    Xlemac = Xcgfg - Xcglemac + Wwg/Wfg *(Xcgwglemac - Xcglemac)
    Xlemacs.append(Xlemac)

print(Xlemacs[-1])
'''



print('Xcgfg:', Xcgfg)