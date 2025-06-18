import numpy as np
np.set_printoptions(suppress=True, formatter={'float_kind':'{:0.4f}'.format})

CL = 1
V = 200/3.6
rho = 1.225
taper = 0.4
area = 3.5
lift = 0.5 * rho * V**2 * area * CL
print(lift)
thickness = 0.09 # T/C
beam_width = 0.4 # Percentage of chord
AR = 5
Nr_sec = 10

SF = 3
yield_stress = 350 * 10**6 #MPa
rho_aluminium = 2780

mac = np.sqrt(area/AR)
b = AR * mac
c_r = mac * 2 / (1+taper)
t_r = c_r*thickness
w_r = c_r * beam_width
print(f"Span: {b}, root chord: {c_r}")
span_sections = np.linspace(0, b/2, Nr_sec)
load_distr = np.sqrt(1-np.linspace(0, 1, num=Nr_sec)**2) * lift / Nr_sec * 4 / np.pi / 2
moment_r = np.sum(span_sections*load_distr) * SF
print(load_distr, sum(load_distr), moment_r)
I = moment_r * (t_r / 2) / yield_stress
t = max(I * 2 / (t_r**2 * w_r), 1/1000)
area_sec = 2 * (t_r + w_r) * t
wingbox_weight = area_sec * (1 + taper) / 2 * b * rho_aluminium
skin_area = (w_r * 2 + t_r) * (1 + taper) / 2 * b
skin_weight = skin_area * 1/1000 * rho_aluminium
print(wingbox_weight, skin_weight)
total_weight = wingbox_weight + skin_weight
print('height:', t_r, 'width:', w_r)
print(t*1000, 'mm', total_weight, 'kg')


area_vtail = 2.5
lift_v = 0.5 * rho * V**2 * area_vtail * CL
h_vtail = 2
cv_t = c_r
mac_v = area_vtail / h_vtail
cv_r = (2 * mac_v) - cv_t
tv_r = cv_r * thickness
wv_r = cv_r * beam_width
tv_t = cv_t * thickness
taper_h = cv_t/cv_r
print('Horizontal tail:', cv_r, cv_t, taper_h)
print(tv_r, tv_t)
span_sections_h = np.linspace(0, h_vtail, Nr_sec)
load_distr = np.linspace(cv_r, cv_t, num=Nr_sec) * lift_v / mac_v / Nr_sec
moment_r = np.sum(span_sections_h*load_distr) * SF + moment_r
print(load_distr, sum(load_distr), moment_r, lift_v)
I = moment_r * (tv_r / 2) / yield_stress
t = max(I * 2 / (tv_r**2 * wv_r), 1/1000)
area_sec = 2 * (tv_r + wv_r) * t
wingbox_weight = area_sec * (1 + taper) / 2 * h_vtail * rho_aluminium
skin_area = (wv_r * 2 + tv_r) * (1 + taper) / 2 * h_vtail
skin_weight = skin_area * 1/1000 * rho_aluminium
print(wingbox_weight, skin_weight)
total_weight_v = wingbox_weight + skin_weight
print('height:', tv_r, 'width:', wv_r)
print(t*1000, 'mm', total_weight_v, 'kg')

print(f"Total weight: {total_weight_v + total_weight} kg")