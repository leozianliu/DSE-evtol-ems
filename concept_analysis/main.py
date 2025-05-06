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

LD_ratio = 10 #Lift to drag ratio


N_disks = 6 #Number of disks (between 4 and 8)
D_rotor = 3 #m (max w_hover/2)
S_rotor = D_rotor**2 * np.pi / 4 #m^2
S_disks = N_disks * S_rotor #m^2
A_front = 2 #m^2 #m^2 (area of the front of the vehicle)
C_D = 0.4 #Drag coefficient of the front of the vehicle

#Constants:
g = 9.81 #m/s^2
rho_air = 1.225 #kg/m^3
density_batt = 250*3600*0.8 #Wh/kg, 80% depth of discharge

#Code parameters:
battery = True
wing = False
mtow = 2500 * g #N
mtow_prev = 0 #N

while abs(mtow_prev - mtow) > 0.1:
    print('Hey, I am in the loop')
    #Energy calculation:
    if wing:
        T = mtow / LD_ratio #N
        P_cruise = T * v_cruise #W
    else:
        F_sideways = 0.5 * rho_air * (v_cruise**2) * A_front * C_D #N
        T = np.sqrt(mtow + F_sideways)
        P_cruise = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks) #W
    E_cruise = P_cruise * t_cruise #J

    T = mtow #N
    P_induced = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks) #W
    P_profile = 0 # Neglect for now, but can be added later
    P_hover = P_induced + P_profile #W
    E_hover = P_hover * t_hover #J

    E_total = E_cruise + E_hover #J

    # Energy source mass as a function of type, energy capacity and output power
    m_powersource = E_total / density_batt #kg

    if wing:
        # Wing mass, area, drag as a function of mtow, v_cruise
        m_eom = 2200 / 3900 * mtow #kg (ratio for winged vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
    else:
        m_eom = 1700/3450 * mtow #kg (ratio for non-winged vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)

    mtow_prev = mtow #N
    mtow = m_payload + m_powersource + m_eom #N

print("Final MTOW: ", mtow / g, "kg")
print("Final Power source ratio: ", m_powersource/mtow*100, "%")