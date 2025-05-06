import numpy as np

#Mission parameters:
m_payload = 400 #kg
range = 120 #km
landings = 3
v_cruise = 200/3.6 #m/s
v_hover = 67/3.6 #m/s
w_hover = 6 #m
l_hover = 12 #m
n_cycles = 3000
mtow_max = 3175 #kg
n_max = 2 #G
cost_miss = 80 #€/h
cost_init = 2 #M€
t_landing = 60 #sec
t_hover = landings*t_landing #sec
t_cruise = range * 1000 / v_cruise #sec

N_blades = 4 #Number of blades (between 2 and 6)
N_disks = 6 #Number of disks (between 4 and 8)
D_rotor = 3 #m (max w_hover/2)
S_rotor = D_rotor**2 * np.pi / 4 #m^2
S_disks = N_disks * S_rotor #m^2
S_wing = 20 #m^2
C_D = 0.05 #Drag coefficient

#Constants:
g = 9.81 #m/s^2
rho_air = 1.225 #kg/m^3
Density_batt = 250*0.8 #Wh/kg

#Code parameters:
battery = True
wing = False
mtow_init = 2500 #kg


#Energy calculation:

# P_parasite = 0
# if wing:
#     P_parasite += 0.5 * rho_air * (v_cruise**3) * S_wing * C_D #W

# P_climb = 0 # On average
# P_hover 
# P_hover = P_parasite + P_profile + P_climb + P_hover #W

# if wing:
#     D_cruise = 0.5 * rho_air * (v_cruise**2) * S_wing * C_D #N
#     mu = v_cruise / V_tip #m/s
#     P_cr_profile = 1/8 * rho_air * N_blades * V_tip**3 * S_disks * sigma * (1 + K * mu**2) * Cd_0 #W
#     P_cr_induced = 0 # For cruise with wings
#     P_cr_parasite = v_cruise * D_cruise #W



P_cruise = P_cr_parasite + P_cr_induced + P_cr_profile #W
E_cruise = P_cruise * t_cruise #J

T = mtow_init * g #N
P_induced = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks) #W
P_hover = P_induced + P_profile #W
E_hover = P_hover * t_hover #J

E_total = E_cruise + E_hover

#Battery calculation:
m_powersource = E_total / Density_batt #kg

#Other weigths:
m_eom = 0.2 * mtow_max #kg

if wing:
    m_wing = 0.1 * mtow_max #kg
    m_eom += m_wing #kg


m_propulsion = 0.1 * mtow_max #kg

mtow = m_payload + m_powersource + m_eom + m_propulsion #kg