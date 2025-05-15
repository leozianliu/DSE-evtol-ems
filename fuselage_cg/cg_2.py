import numpy as np 
import matplotlib.pyplot as plt

#weights [kg]
w_person = 80

w_stretcher = 25 
w_oxygen = 23
w_defib = 4.5
w_suction = 2
w_imm = 7
w_medkit = 13.5

w_eq = w_oxygen + w_defib + w_suction + w_imm + w_medkit

#fuselage skin contribution
fuselage_x = 1894.27
fuselage_z = 719.34 #(this is taken from the floor)
fuselage_y = 0
w_fuselage = 625.7787 
w_total_fuselage = 2160 
w_extra = w_total_fuselage - w_fuselage

x_cg_wanted = 2180 #from the front of the cabin [mm]

x_extra = (x_cg_wanted - (fuselage_x*w_fuselage)/(w_total_fuselage))* (w_total_fuselage/w_extra)
x_cg_total = (fuselage_x*w_fuselage + x_extra*w_extra) / w_total_fuselage


def cg_f(pt_x, pt_y, eq_x, eq_y, d1_x, d1_y, d2_x, d2_y, pi_x, pi_y, w_pat):
    ##CG locations [mm]
    #patient tray 

    x_cg_1 = (w_pat*pt_x  + w_eq*eq_x + w_person*d1_x + w_person*d2_x + w_person*pi_x + w_total_fuselage*x_cg_total)
    x_cg_2 = (w_pat + w_eq + w_person + w_person + w_person + w_total_fuselage)
    y_cg_1 = (w_pat*pt_y + w_eq*eq_y + w_person*d1_y + w_person*d2_y + w_person*pi_y + w_fuselage*fuselage_y)
    y_cg_2 = (w_pat + w_eq + w_person + w_person + w_person + w_fuselage)

    x_cg = x_cg_1 / x_cg_2
    y_cg = y_cg_1 / y_cg_2
    return x_cg, y_cg

def cg(pt_x, pt_y, eq_x, eq_y, d1_x, d1_y, d2_x, d2_y, pi_x, pi_y, w_pat):
    ##CG locations [mm]
    #patient tray 

    x_cg_1 = (w_pat*pt_x + w_eq*eq_x + w_person*d1_x + w_person*d2_x + w_person*pi_x)
    x_cg_2 = (w_pat + w_eq + w_person + w_person + w_person)
    y_cg_1 = (w_pat*pt_y + w_eq*eq_y + w_person*d1_y + w_person*d2_y + w_person*pi_y )
    y_cg_2 = (w_pat + w_stretcher + w_eq + w_person + w_person + w_person)

    x_cg = x_cg_1 / x_cg_2
    y_cg = y_cg_1 / y_cg_2
    return x_cg, y_cg


#case_1 - everyone is at the right spot with fuselage skin
x_cg_1_f = cg_f(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 1244.67, 
    d1_y = 600, 
    d2_x = 1897.8 , 
    d2_y = 600, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_person + w_stretcher
)

#case_1 - everyone is at the right spot without fuselage skin
x_cg_1 = cg(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 1244.67, 
    d1_y = 600, 
    d2_x = 1897.8 , 
    d2_y = 600, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_person + w_stretcher
)

#case_2 - everyone is at the right spot with fuselage skin - patient
x_cg_2_f = cg_f(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 1244.67, 
    d1_y = 600, 
    d2_x = 1897.8 , 
    d2_y = 600, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_stretcher 
)

#case_2 - everyone is at the right spot without fuselage skin - patient
x_cg_2 = cg(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 1244.67, 
    d1_y = 600, 
    d2_x = 1897.8 , 
    d2_y = 600, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_stretcher
)

#case_3 - both medics are next to the patient - with fuselage skin
x_cg_3_f = cg_f(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 2190, 
    d1_y = -200, 
    d2_x = 2890 , 
    d2_y = -200, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_stretcher + w_person
)

#case_3 - both medics are next to the patient - without fuselage skin
x_cg_3 = cg(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 2190, 
    d1_y = -200, 
    d2_x = 2890 , 
    d2_y = -200, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_stretcher + w_person
)

#case_4 - both medics next to patient at the back - with fuselage skin
x_cg_4_f = cg_f(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 3920, 
    d1_y = -200, 
    d2_x = 3220, 
    d2_y = -200, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_stretcher + w_person
)

#case_4 - both medics next to patient at the back - without fuselage skin
x_cg_4 = cg(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 3920, 
    d1_y = -200, 
    d2_x = 3220, 
    d2_y = -200, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_stretcher + w_person
)

#case_5 - one medic in the front and one in the back - with fuselage skin
x_cg_5_f = cg_f(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 3920, 
    d1_y = -200, 
    d2_x = 2190, 
    d2_y = -200, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_stretcher + w_person
)

#case_5 - one medic in the front and one in the back - without fuselage skin
x_cg_5 = cg(
    pt_x = 3169.6, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 3920, 
    d1_y = -200, 
    d2_x = 2190, 
    d2_y = -200, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_stretcher + w_person
)

#case_6 - both medics in the back as well as the patient - with fuselage skin
x_cg_6_f = cg_f(
    pt_x = 1070 + 3200, 
    pt_y = -650, 
    eq_x = 3519, 
    eq_y = 500, 
    d1_x = 1070 + 3200, 
    d1_y = -200, 
    d2_x = 1070 + 3200, 
    d2_y = -200, 
    pi_x = 1244.6, 
    pi_y = -600,
    w_pat = w_stretcher + w_person
)


cabin_width = 2100
cabin_length = 3200
cabin_height = 1900
floor_height = 424
nose_length = 1070
z_loc = 720 #(above the floor)
y_graph = z_loc + floor_height 

#putting them all together in a list
x_cg = [x_cg_1[0], x_cg_2[0], x_cg_3[0], x_cg_4[0], x_cg_5[0]]
x_cg_f = [x_cg_total, x_cg_1_f[0], x_cg_2_f[0], x_cg_3_f[0], x_cg_4_f[0], x_cg_5_f[0], x_cg_6_f[0]]

y_cg = [x_cg_1[1], x_cg_2[1], x_cg_3[1], x_cg_4[1]]
y_cg_f = [x_cg_1_f[1], x_cg_2_f[1], x_cg_3_f[1], x_cg_4_f[1], x_cg_5_f[1]]

#min and max values
x_cg_min = min(x_cg)
x_cg_max = max(x_cg)
x_cg_f_min = min(x_cg_f)
x_cg_f_max = max(x_cg_f)
print(x_cg_f)
x_range = x_cg_f_max - x_cg_f_min
print(x_cg_f_min, x_cg_f_max, x_range)

y_cg_min = min(y_cg)
y_cg_max = max(y_cg)
y_cg_f_min = min(y_cg_f)
y_cg_f_max = max(y_cg_f)


#including fuselage skin
plt.figure(figsize=(10, 4))
#plt.scatter(fuselage_x, y_graph, label='CG skin', marker='x')
plt.scatter(x_cg_total, y_graph, label='CG total', marker='x')
plt.scatter(x_cg_1_f[0], y_graph, label='CG case 1', marker='x')
plt.scatter(x_cg_2_f[0], y_graph, label='CG case 2', marker='x')
plt.scatter(x_cg_3_f[0], y_graph, label='CG case 3', marker='x')
plt.scatter(x_cg_4_f[0], y_graph, label='CG case 4', marker='x')
plt.scatter(x_cg_5_f[0], y_graph, label='CG case 5', marker='x')
plt.scatter(x_cg_6_f[0], y_graph, label='CG case 6', marker='x')
#plt.axvline(x = nose_length, linestyle='--')
#plt.axvline(x = cabin_length + nose_length, linestyle='--')
#plt.axhline(y = floor_height, linestyle='--')
#plt.axhline(y = cabin_height + floor_height, linestyle='--')
plt.xlabel('Cg position along the fuselage length [mm]', fontsize = 12)
plt.yticks(color='white')
plt.legend(loc = 'lower right', fontsize = 8)
plt.show()
plt.close()
'''

plt.figure(figsize=(5, 6))
plt.scatter(fuselage_y, y_graph, label='CG skin', marker='x')
plt.scatter(x_cg_1_f[1], y_graph, label='CG case 1', marker='x')
plt.scatter(x_cg_2_f[1], y_graph, label='CG case 2', marker='x')
plt.scatter(x_cg_3_f[1], y_graph, label='CG case 3/4/5', marker='x')
plt.axvline(x = -cabin_width/2, linestyle='--')
plt.axvline(x = cabin_width/2, linestyle='--')
plt.axhline(y = floor_height, linestyle='--')
plt.axhline(y = cabin_height + floor_height, linestyle='--')
plt.legend(loc = 'lower right')
plt.show()
plt.close()
'''

'''
#excluding fuselage skin
plt.figure(figsize=(6, 4))
plt.scatter(fuselage_x, y_graph, label='CG skin')
plt.scatter(x_cg_1[0], y_graph, label='CG only payload case 1')
plt.scatter(x_cg_2[0], y_graph, label='CG only payload case 2')
plt.scatter(x_cg_3[0], y_graph, label='CG only payload case 3')
plt.scatter(x_cg_4[0], y_graph, label='CG only payload case 4')
plt.axvline(x = nose_length, linestyle='--')
plt.axvline(x = cabin_length + nose_length, linestyle='--')
plt.axhline(y = floor_height, linestyle='--')
plt.axhline(y = cabin_height + floor_height, linestyle='--')
plt.legend(loc = 'lower right')
plt.show()
plt.close()


plt.figure(figsize=(5, 6))
plt.scatter(fuselage_y, y_graph, label='CG skin')
plt.scatter(x_cg_1[1], y_graph, label='CG only payload case 1')
plt.scatter(x_cg_2[1], y_graph, label='CG only payload case 2')
plt.scatter(x_cg_3[1], y_graph, label='CG only payload case 3')
plt.scatter(x_cg_4[1], y_graph, label='CG only payload case 4')
plt.axvline(x = -cabin_width/2, linestyle='--')
plt.axvline(x = cabin_width/2, linestyle='--')
plt.axhline(y = floor_height, linestyle='--')
plt.axhline(y = cabin_height + floor_height, linestyle='--')
plt.legend(loc = 'lower right')
plt.show()
'''