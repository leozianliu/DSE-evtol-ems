import numpy as np
from engine import *

#Mission parameters:
m_payload = 400 #kg
range = 120 #km
landings = 4
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

LD_ratio = 15 #Lift to drag ratio
eff_motor = 0.95 #Efficiency of the motor
eff_propeller_takeoff = 0.85 #Efficiency of the propeller

# Note: If integrated propulsion is used, N_disks_cruise is not used
S_disks = 25.13 #m^2 (total disk area based on size requirements, can be changed later)
N_disks_cruise = 1 #Number of disks (between 2 and 4)
N_disks_takeoff = 8 #Number of disks (between 4 and 8)
#D_rotor = 4 #m (max w_hover/2)
#S_rotor = D_rotor**2 * np.pi / 4 #m^2
S_rotor = S_disks / N_disks_takeoff #m^2 (disk area per rotor)
D_rotor = 2 * np.sqrt(S_rotor / np.pi) #m (diameter of the rotor)
#S_disks = N_disks_takeoff * S_rotor #m^2
A_front = 5 #m^2 #m^2 (area of the front of the vehicle)
C_D = 0.4 #Drag coefficient of the front of the vehicle

#Configuration parameters:
battery = True
wing = True
integrated_prop = False #Use the same motors for hover and cruise
tilt_wing = False

#MTOW
mtom = 2500 #kg #Maximum takeoff mass (MTOM) of the vehicle (GUESS)

#Constants:
g = 9.81 #m/s^2
rho_air = 1.225 #kg/m^3
density_batt_whkg = 300 #300Wh/kg (Chinese),
density_batt = density_batt_whkg*3600*0.80 #degrade to 80% of the battery capacity
separate_prop_extra_drag_factor = 1.15  0 # Skin friction of motor booms increases skin drag by 30% and total drag by 15% (skin drag is 50% of total drag)
LD_reduction_factor = 1/separate_prop_extra_drag_factor
blockage_factor_tiltwing = 0.90 # Free area over total area for propellers in tilt-wing configuration
blockage_factor_tiltrotor = 0.78 # Free area over total area for propellers in tilt-rotor configuration
blockage_factor_sepprop = 0.86 # Free area over total area for propellers in separate propulsion configuration
m_tilt_mech = 60 #kg for the tilt-wing mechanism

#Code parameters
mtow_prev = 0 #N
n = 1

mtow = mtom * g #N (initial guess for MTOW)

if tilt_wing and not integrated_prop:
    print("Warning: Tilt-wing configuration must use integrated propulsion in this program.")
    raise ValueError("Warning: Tilt-wing configuration must use integrated propulsion in this program.")

print("Winged: ", wing)
print("Integrated propulsion: ", integrated_prop)
print("Tilt-wing: ", tilt_wing)

if integrated_prop:
    if tilt_wing:
        blockage_factor = blockage_factor_tiltwing
    else:
        blockage_factor = blockage_factor_tiltrotor
elif not integrated_prop:
    blockage_factor = blockage_factor_sepprop
else:
    raise ValueError("Warning: No propeller blockage ratio for this configuration. Cry about it.")

if not wing:
    raise ValueError("Warning: Non-winged configuration not implemented yet. Cry about it.")
#     S_disks = 120 #m^2 (for now, but can be changed later)
# diameter_rotor = 2 * np.sqrt(S_disks/N_disks_takeoff / np.pi) #m (diameter of the rotor)
# print("Diameter of the rotor: ", diameter_rotor, "m")

print("Disk area original: ", S_disks, "m^2")
S_disks = blockage_factor * S_disks #m2 (blockage factor for tilt-wing configuration)
print('Disk area adjusted for blockage factor: ', S_disks, "m^2")

while abs(mtow_prev - mtow) > 0.1 and n<1000:
    n+=1
    #Energy calculation:
    if wing:
        # Cruise power calculation
        LD_ratio = 13
        if tilt_wing:
            LD_ratio = 15
        T = mtow / LD_ratio #N
        if not integrated_prop:
            LD_ratio = 15
            # Cruise power calculation
            T = mtow / (LD_ratio * LD_reduction_factor) # Skin friction of motor booms increases skin drag by 30% and total drag by 15% (skin drag is 50% of total drag)
        P_cruise = T * v_cruise# / eff_motor
    else:
        F_sideways = 0.5 * rho_air * (v_cruise**2) * A_front * C_D #N
        T = np.sqrt(mtow**2 + F_sideways**2)
        P_cruise = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks) #W
    E_cruise = P_cruise * t_cruise #J

    T = mtow #N
    P_induced = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks) #W
    P_profile = 0 # Neglect for now, but can be added later
    P_hover = (P_induced + P_profile) / eff_motor / eff_propeller_takeoff #W
    E_hover = P_hover * t_hover #J

    E_total = E_cruise + E_hover #J

    # Energy source mass as a function of type, energy capacity and output power
    m_powersource = E_total / density_batt #kg

    if wing and integrated_prop:
        if not tilt_wing:
            # Wing mass, area, drag as a function of mtow, v_cruise
            # oem = battery + equipment + propulsion + structure
            m_equipment = 450/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
            m_structure = 1000/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
            P_cruise_single = P_cruise / N_disks_takeoff / 1000 #KW, again for integrated propulsion takeoff=cruise for N
            P_hover_single = P_hover / N_disks_takeoff / 1000 #KW (energy per disk)
            # print('Single motor power takeoff (kW): ', P_hover_single)
            m_propulsion = calculatePropulsionMass(P_cruise_single, P_hover_single, N_disks_takeoff, D_rotor) #kg (mass of the propulsion system)
            m_propulsion = m_propulsion + N_disks_takeoff * 15 # mass of motor tilter actuator (10kg per motor)
            # print("Propulsion mass: ", m_propulsion, "kg")
            m_eom = m_powersource + m_equipment + m_structure + m_propulsion
            # print("Power mass:", m_powersource, "kg")
        if tilt_wing:
            # Wing mass, area, drag as a function of mtow, v_cruise
            # oem = battery + equipment + propulsion + structure
            m_equipment = 450/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
            m_structure = 1000/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
            m_structure = m_structure + m_tilt_mech #kg for the tilt-wing mechanism 
            P_cruise_single = P_cruise / N_disks_takeoff / 1000 #KW, again for integrated propulsion takeoff=cruise for N
            P_hover_single = P_hover / N_disks_takeoff / 1000 #KW (energy per disk)
            # print('Single motor power takeoff (kW): ', P_hover_single)
            m_propulsion = calculatePropulsionMass(P_cruise_single, P_hover_single, N_disks_takeoff, D_rotor) #kg (mass of the propulsion system)
            # print("Propulsion mass: ", m_propulsion, "kg")
            m_eom = m_powersource + m_equipment + m_structure + m_propulsion
            # print("Power mass:", m_powersource, "kg")
    elif wing and (not integrated_prop):
        # Wing mass, area, drag as a function of mtow, v_cruise
        # oem = battery + equipment + propulsion + structure
        m_equipment = 450/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
        m_structure = 1000/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
        P_cruise_single = P_cruise / N_disks_cruise / 1000 #KW, again for integrated propulsion takeoff=cruise for N
        P_hover_single = P_hover / N_disks_takeoff / 1000 #KW (energy per disk)
        # print('Single motor power takeoff (kW): ', P_hover_single)
        m_propulsion = calculatePropulsionMass(P_cruise_single, 0, N_disks_cruise, D_rotor) + calculatePropulsionMass(0, P_hover_single, N_disks_takeoff, D_rotor) #kg (mass of the propulsion system)
        # print("Propulsion mass cruise: ", calculatePropulsionMass(P_cruise_single, 0, N_disks_cruise, D_rotor), "kg")
        # print("Propulsion mass hover: ", calculatePropulsionMass(0, P_hover_single, N_disks_takeoff, D_rotor), "kg")
        m_eom = m_powersource + m_equipment + m_structure + m_propulsion
    elif (not wing):
        m_eom = 1500/3450 * mtow/g #kg (ratio for non-winged vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
    else:
        raise ValueError("Warning: No mass ratio for this configuration. Cry about it.")

    mtow_prev = mtow #N
    mtow = (m_payload + m_powersource + m_eom) * g #N

print("---------------------------------")
print("Rotor diameter: ", D_rotor, "m")
print('Single motor power takeoff: ', P_hover_single, "kW")
# print("Power mass:", m_powersource, "kg")
print("P_hover: ", P_hover / 1000, "kW")
print("P_cruise: ", P_cruise / 1000, "kW")
print("Number of iterations: ", n)
print("Final MTOW: ", mtow / g, "kg")
print("Final Propulsion mass: ", m_propulsion, "kg")
print("Final Power Source mass: ", m_powersource, "kg")
print("Final EOM mass ratio: ", m_eom*g/mtow*100, "%")
#print(P_hover/1000)