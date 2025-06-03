import numpy as np 
import matplotlib.pyplot as plt
from whole_cg import Xcgfg

#weights [kg]
w_person = 80
w_stretcher = 57

#preliminary
w_fuselage = 1100 
#x_cg_struc = xcgfg_actual * 1000
x_cg_struc = Xcgfg * 1000


def cg_f(pt_x, d1_x, d2_x, pi_x, w_pat):
    x_cg_1 = (w_pat*pt_x   + w_person*d1_x + w_person*d2_x + w_person*pi_x + w_fuselage*x_cg_struc)
    x_cg_2 = (w_pat + w_person + w_person + w_person + w_fuselage)

    x_cg = x_cg_1 / x_cg_2
    return x_cg

def cg(pt_x, d1_x, d2_x, pi_x, w_pat):
    x_cg_1 = (w_pat*pt_x + w_person*d1_x + w_person*d2_x + w_person*pi_x)
    x_cg_2 = (w_pat + w_person + w_person + w_person)

    x_cg = x_cg_1 / x_cg_2
    return x_cg


#case_1 - everyone is at their spot 
x_cg_1 = cg_f(
    pt_x = 3200, 
    d1_x = 2425, 
    d2_x = 2425, 
    pi_x = 1450, 
    w_pat = w_person + w_stretcher
)

#case_2 - no patient
x_cg_2 = cg_f(
    pt_x = 3200, 
    d1_x = 2425, 
    d2_x = 2425, 
    pi_x = 1450, 
    w_pat = w_stretcher
)


#case_3 - both medics are at the end of the stretcher 
x_cg_3 = cg_f(
    pt_x = 3200, 
    d1_x = 3975, 
    d2_x = 3975, 
    pi_x = 1450, 
    w_pat = w_person + w_stretcher
)

#case_4 - one medic is at the front and the other sitting 
x_cg_4= cg_f(
    pt_x = 3200, 
    d1_x = 2425, 
    d2_x = 1450, 
    pi_x = 1450, 
    w_pat = w_person + w_stretcher
)


cabin_width = 1800
cabin_length = 3100
cabin_height = 1900
floor_height = 325
nose_length = 1200
z_loc = 747 #(above the floor)
y_graph = z_loc + floor_height 


#including fuselage skin
plt.figure(figsize=(10, 4))
plt.scatter(x_cg_struc, y_graph, label='CG total', marker='x')
plt.scatter(x_cg_1, y_graph, label='CG case 1', marker='x')
plt.scatter(x_cg_2, y_graph, label='CG case 2', marker='x')
plt.scatter(x_cg_3, y_graph, label='CG case 3', marker='x')
plt.scatter(x_cg_4, y_graph, label='CG case 4', marker='x')
#plt.axvline(x = nose_length, linestyle='--')
#plt.axvline(x = cabin_length + nose_length, linestyle='--')
#plt.axhline(y = floor_height, linestyle='--')
#plt.axhline(y = cabin_height + floor_height, linestyle='--')
plt.xlabel('Cg position along the fuselage length [mm]', fontsize = 12)
plt.yticks(color='white')
plt.legend(loc = 'lower right', fontsize = 8)
plt.show()
plt.close()

#making a list 
xcg = [x_cg_1, x_cg_2, x_cg_3, x_cg_4]
#min and max values
x_cg_min = min(xcg)
x_cg_max = max(xcg)
x_range = x_cg_max - x_cg_min

print(f"CG min: {x_cg_min} mm, CG max: {x_cg_max} mm, Range: {x_range} mm")