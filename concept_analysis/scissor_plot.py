import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def C_l_alpha (A, M, eta, lambda_half):
    beta = np.sqrt(1-M**2)

    numerator = 2 * np.pi * A
    term1 = (A * beta / eta)**2
    term2 = 1 + (np.tan(lambda_half)**2) / (beta**2)
    denominator = 2 + np.sqrt(4 + term1 * term2)

    C_l_alpha = numerator / denominator
    return C_l_alpha

def C_l_alpha_aircraft (C_l_alpha_wing, b_f, b, S_net, S):
    term1 = C_l_alpha_wing * (1 + 2.15 * (b_f / b)) *  (S_net / S)
    term2 = np.pi / 2 * b_f**2 / S
    C_l_alpha_aircraft = term1 + term2 
    
    return C_l_alpha_aircraft


def x_ac_bar_wing_fuselage (x_ac_bar_wing, C_L_alpha_minus_H, b, b_f, h_f, l_fn, S, c_bar, taper, c_g, lambda_quarter):

    term1 = 1.8 / C_L_alpha_minus_H * b_f * h_f * l_fn / S / c_bar
    term2 = 0.273 / (1 + taper) * b_f * c_g * (b - b_f) / c_bar ** 2 / (b + 2.15 * b_f) * np.tan (lambda_quarter)
    x_ac_bar_wing_fuselage = x_ac_bar_wing- term1 + term2

    return x_ac_bar_wing_fuselage

def c_mac_wing_fuselage (cm0_airfoil, sweep, A, b_f, h_f, l_f, S, c_bar, c_l_0, C_l_alpha_minus_tail):
    numerator = A * np.cos(sweep)**2
    denominator = A + 2 * np.cos(sweep)
    c_mac_wing = cm0_airfoil * (numerator / denominator)

    factor = -1.8 * (1 - (2.5 * b_f / l_f))
    geometry_term = (np.pi * b_f * h_f * l_f) / (4 * S * c_bar)
    lift_ratio = c_l_0 / C_l_alpha_minus_tail
    c_mac_wing = factor * geometry_term * lift_ratio
    c_mac_wing_fuselage = c_mac_wing + c_mac_wing

    return c_mac_wing_fuselage

def c_n_beta_wing_fuselage (l_cg, l_f, h_f, S_fs, h_f1, h_f2, b_f1, b_f2, Bp, I_p_small, D_p_small, I_p_big, D_p_big, C_n_b_i):

    k_beta = 0.3 * (l_cg / l_f) + 0.75 * (h_f / l_f) - 0.105
    term1 = (S_fs * l_f) / (S * b)
    term2 = (h_f1 / h_f2) ** 0.5
    term3 = (b_f2 / b_f1) ** (1/3)
    C_n_beta_f = -k_beta * term1 * term2 * term3
    print('C_n_beta_fuselage =', C_n_beta_f)

    C_n_beta_p = -0.053 * Bp * (I_p_small * D_p_small **2 * 2/ S / b + I_p_big * D_p_big **2 * 2/ S / b )
    print('C_n_beta_propeller =', C_n_beta_p)

    print('C_n_beta_wing =', C_n_b_i)

    C_n_beta =  C_n_beta_f +  C_n_beta_p + C_n_b_i

    return C_n_beta 

######## Constants , all in SI unit and radian ##########
C_r = 1.4
taper = 0.457
b_f = 1.8
b = 12
S_net_wing = 11.272
S_wing = S_net_wing + b_f * C_r
de_da = 0.1
l_h = 4.5
S = 11.272
c_bar = 1.089
Vh_V = 0.8
s_margin = 0.05
M_wing = 0.21
A_wing = 13.565
lambda_half_wing = 0
c_g = S/b
h_f = 1.85 + 0.325
l_fn = 1.5
l_f = 6
lambda_quarter = 0
V = 55.6
rho = 1
lam_bda = 0
cm0_airfoil = -0.092
c_l_0 = 0.5
C_L_H = -0.6
C_L_minus_H = 0.8

###################### Vertical constant ##############
C_Y_b_v = np.pi/2 * 1.5
l_cg = 2.4 
S_fs = 9 
h_f1 = 1.913
h_f2 = 1.785
b_f1 = 1.730
b_f2 = 1.727
Bp = 5 # Number of blades for each propeller
I_p_small = 0.9
D_p_small = 2
I_p_big =  0.9
D_p_big = 3
C_n_b_i = -0.017 # For high wing
ds_db = 0
Vv_V = 1

A_htail = 5
lambda_half_htail = 0
eta = 0.95
l_v = 4.5

x_bar_cg = np.linspace(0,1,1000)

x_bar_ac_wing = 0.25

M_htail= M_wing * Vh_V
C_L_alpha_wing = C_l_alpha (A_wing, M_wing, eta, lambda_half_wing)
C_L_alpha_htail = C_l_alpha (A_htail, M_htail, eta, lambda_half_htail)

C_l_alpha_minus_tail = C_l_alpha_aircraft (C_L_alpha_wing, b_f, b, S_net_wing, S_wing)

x_bar_ac = x_ac_bar_wing_fuselage (x_bar_ac_wing, C_l_alpha_minus_tail, b, b_f, h_f, l_fn, S, c_bar, taper, c_g, lambda_quarter)

# Stability curve
deno_sta_1 = x_bar_cg - x_bar_ac + s_margin
deno_sta_2 = x_bar_cg - x_bar_ac
nume_sta_1 = (C_L_alpha_htail / C_l_alpha_minus_tail ) * (1 - de_da) * (l_h / c_bar) * (Vh_V ** 2)
Sh_S_sta_margined = deno_sta_1 / nume_sta_1
Sh_S_sta_un = deno_sta_2 / nume_sta_1

C_m_ac = c_mac_wing_fuselage (cm0_airfoil, lam_bda, A_wing, b_f, h_f, l_f, S, c_bar, c_l_0, c_bar)

# Controllability curve
term_cont_1 = x_bar_cg / ((C_L_H / C_L_minus_H) * (l_h / c_bar) * Vh_V ** 2)
term_cont_2_nume = (C_m_ac / C_L_minus_H) - x_bar_ac
term_cont_2_deno = (C_L_H / C_L_minus_H) * (l_h / c_bar) * Vh_V ** 2
term_cont_2 = term_cont_2_nume / term_cont_2_deno
Sh_S_cont = term_cont_1 + term_cont_2

# plt.plot(x_bar_cg, Sh_S_cont, label='Controllability')
# plt.plot(x_bar_cg, Sh_S_sta_margined, label='Stability with S.M.')
# plt.plot(x_bar_cg, Sh_S_sta_un, label='Stability')

# plt.ylim(bottom=0)  
# plt.ylabel(r"$S_h/S$")  
# plt.xlabel(r"$\bar{x_{cg}}$")  
# plt.legend()
# plt.grid(True)
# plt.show()

upper_curve = np.maximum(Sh_S_cont, Sh_S_sta_margined)

plt.plot(x_bar_cg, Sh_S_cont, label='Controllability')
plt.plot(x_bar_cg, Sh_S_sta_margined, label='Stability with S.M.')
plt.plot(x_bar_cg, Sh_S_sta_un, label='Stability')
plt.axvline(x=0.37, color='red', linestyle='--', label='x_cg = 0.37')
plt.axvline(x=0.5, color='red', linestyle='--', label='x_cg = 0.5')


y_max = np.max(upper_curve)
plt.fill_between(x_bar_cg, upper_curve, y_max, color='green', alpha=0.3, label='Safe Region')

plt.fill_between(x_bar_cg, 0, upper_curve, color='red', alpha=0.3, label='Unsafe Region')

plt.xlim(0, 1)           # Set x-axis range
plt.ylim(0, 0.55)        # Set y-axis range

plt.ylim(bottom=0)
plt.ylabel(r"$S_h/S$")
plt.xlabel(r"$\bar{x}_{cg}$")
plt.legend()
plt.grid(True)
plt.show()


################################## Vertical stabilizer ##################################

# Stability curve
Sv_S = np.linspace(0,1,1000)
C_n_beta_a_minus_v = c_n_beta_wing_fuselage (l_cg, l_f, h_f, S_fs, h_f1, h_f2, b_f1, b_f2, Bp, I_p_small, D_p_small, I_p_big, D_p_big, C_n_b_i)

C_n_beta = C_n_beta_a_minus_v + C_Y_b_v * Sv_S * l_v / b * (1 - ds_db) * Vv_V **2

# plt.plot (C_n_beta , Sv_S, label='Stability')
# plt.axvline(x=0, color='red', linestyle='--', label='x = 0')
# plt.axvspan(min(C_n_beta), 0, color='red', alpha=0.3, label='Unsafe Region' )
# plt.axvspan(0, max(C_n_beta), color='green', alpha=0.3, label='Safe Region')
# plt.ylabel(r"$S_h/S$")
# plt.xlabel(r"$C_{n_{beta}}$")
# plt.legend()
# plt.grid(True)
# plt.show()


x_min, x_max = -0.05, 0.2
y_min, y_max = 0, 1
fig, ax = plt.subplots()

ax.axvspan(x_min, 0, ymin=0, ymax=1, color='red', alpha=0.2)

ax.fill_between(C_n_beta, Sv_S, y_max, where=(C_n_beta >= 0), 
                color='green', alpha=0.3, label='Safe Region')

ax.fill_between(C_n_beta, Sv_S, y_min, where=(C_n_beta >= 0), 
                color='red', alpha=0.3, label='Unsafe Region')

ax.plot(C_n_beta, Sv_S, color='royalblue', label='Stability')

ax.axvline(x=0, color='red', linestyle='--', label='C_n_beta = 0')

ax.set_xlabel(r"$C_{n_{\beta}}$")
ax.set_ylabel(r"$S_h/S$")
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.grid(True)
ax.legend()
plt.show()












