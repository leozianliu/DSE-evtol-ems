import numpy as np
from engine import *

#Mission parameters:
m_payload = 400 #kg
range = 120 #km
landings = 2
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

LD_ratio = 12 #Lift to drag ratio
eff_motor = 0.95 #Efficiency of the motor
eff_propeller = 0.85 #Efficiency of the propeller

# Note: If integrated propulsion is used, N_disks_cruise is not used
N_disks_cruise = 2 #Number of disks (between 2 and 4)
N_disks_takeoff = 8 #Number of disks (between 4 and 8)
D_rotor = 3 #m (max w_hover/2)
S_rotor = D_rotor**2 * np.pi / 4 #m^2
S_disks = N_disks_takeoff * S_rotor #m^2
A_front = 2 #m^2 #m^2 (area of the front of the vehicle)
C_D = 0.4 #Drag coefficient of the front of the vehicle

#Constants:
g = 9.81 #m/s^2
rho_air = 1.225 #kg/m^3
density_batt = 300*3600*0.8*0.85 #300Wh/kg (Chinese), 80% depth of discharge

#Code parameters:
battery = True
wing = True
integrated_prop = False #Use the same motors for hover and cruise
mtow = 2500 * g #N
mtow_prev = 0 #N
n = 1

print("Winged: ", wing)

if not wing:
    S_disks = 120 #m^2 (for now, but can be changed later)
diameter_rotor = 2 * np.sqrt(S_disks/N_disks_takeoff / np.pi) #m (diameter of the rotor)
print("Diameter of the rotor: ", diameter_rotor, "m")

while abs(mtow_prev - mtow) > 0.1 and n<1000:
    n+=1
    #Energy calculation:
    if wing:
        if integrated_prop:
            # Cruise power calculation
            T = mtow / LD_ratio #N
        if not integrated_prop:
            # Cruise power calculation
            T = mtow / (LD_ratio * 0.9)
        P_cruise = T * v_cruise
    else:
        F_sideways = 0.5 * rho_air * (v_cruise**2) * A_front * C_D #N
        T = np.sqrt(mtow**2 + F_sideways**2)
        P_cruise = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks) #W
    E_cruise = P_cruise * t_cruise #J

    T = mtow #N
    P_induced = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks) #W
    P_profile = 0 # Neglect for now, but can be added later
    P_hover = (P_induced + P_profile) / eff_motor / eff_propeller #W
    E_hover = P_hover * t_hover #J

    E_total = E_cruise + E_hover #J

    # Energy source mass as a function of type, energy capacity and output power
    m_powersource = E_total / density_batt #kg

    if wing and integrated_prop:
        # Wing mass, area, drag as a function of mtow, v_cruise
        # oem = battery + equipment + propulsion + structure
        m_equipment = 450/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
        m_structure = 1000/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
        P_cruise_single = P_cruise / N_disks_takeoff / 1000 #KW, again for integrated propulsion takeoff=cruise for N
        P_hover_single = P_hover / N_disks_takeoff / 1000 #KW (energy per disk)
        m_propulsion = calculatePropulsionMass(P_cruise_single, P_hover_single, N_disks_takeoff) #kg (mass of the propulsion system)
        m_propulsion = m_propulsion + N_disks_takeoff * 10 # mass of motor tilter system (20kg per motor)
        print("Propulsion mass: ", m_propulsion, "kg")
        m_eom = m_powersource + m_equipment + m_structure + m_propulsion
        print("Power mass:", m_powersource, "kg")
    elif wing and (not integrated_prop):
        # Wing mass, area, drag as a function of mtow, v_cruise
        # oem = battery + equipment + propulsion + structure
        m_equipment = 450/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
        m_structure = 1000/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
        P_cruise_single = P_cruise / N_disks_cruise / 1000 #KW, again for integrated propulsion takeoff=cruise for N
        P_hover_single = P_hover / N_disks_takeoff / 1000 #KW (energy per disk)
        m_propulsion = calculatePropulsionMass(P_cruise_single, 0, N_disks_cruise) + calculatePropulsionMass(0, P_hover_single, N_disks_takeoff) #kg (mass of the propulsion system)
        print("Propulsion mass cruise: ", calculatePropulsionMass(P_cruise_single, 0, N_disks_cruise), "kg")
        print("Propulsion mass hover: ", calculatePropulsionMass(0, P_hover_single, N_disks_takeoff), "kg")
        m_eom = m_powersource + m_equipment + m_structure + m_propulsion
    elif (not wing):
        m_eom = 1500/3450 * mtow/g #kg (ratio for non-winged vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
    else:
        print("Error: No mass ratio for this configuration")
        m_eom = 0

    mtow_prev = mtow #N
    mtow = (m_payload + m_powersource + m_eom) * g #N

print("Number of iterations: ", n)
print("Final MTOW: ", mtow / g, "kg")
print("Final Power source ratio: ", m_powersource*g/mtow*100, "%")
print("Final eom mass ratio: ", m_eom*g/mtow*100, "%")
print(P_hover/1000)