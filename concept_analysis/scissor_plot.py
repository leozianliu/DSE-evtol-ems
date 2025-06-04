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


######## Constants ##########
C_r = 2
taper = 0.4
b_f = 1.8
b = 12
S_net_wing = 14
S_wing = 14 + b_f * C_r
de_da = 0.1
l_h = 3.5
S = 14
c_bar = 1.7
Vh_V = 1
s_margin = 0.01
M_wing = 0.16
A_wing = 10
lambda_half_wing = 0
c_g = S/b
h_f = 1.85 + 0.325
l_fn = 1.5
l_f = 6
lambda_quarter = 0
V = 55.6
rho = 1.225
lam_bda = 0
cm0_airfoil = -0.039
c_l_0 = 0.34
C_L_H = -0.6
C_L_minus_H = 0.8
M_htail= 0.16
A_htail = 10
lambda_half_htail = 0
eta = 0.95

x_bar_cg = np.linspace(0,1,100)
x_bar_ac_wing = 0.5  

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
plt.axvline(x=0.37, color='red', linestyle='--', label='x = 0.37')
plt.axvline(x=0.5, color='red', linestyle='--', label='x = 0.5')

y_max = np.max(upper_curve)
plt.fill_between(x_bar_cg, upper_curve, y_max, color='green', alpha=0.3, label='Safe Region')

plt.fill_between(x_bar_cg, 0, upper_curve, color='red', alpha=0.3, label='Unsafe Region')

plt.ylim(bottom=0)
plt.ylabel(r"$S_h/S$")
plt.xlabel(r"$\bar{x_{cg}}$")
plt.legend()
plt.grid(True)
plt.show()



















