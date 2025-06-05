import numpy as np
from engine import *

#Mission parameters:
m_payload = 400 #kg
flight_radius = 60 #km (radius * 2)
v_cruise = 200/3.6 #m/s
v_takeoff = 5 #m/s

flights = 2 #Number of flights per mission
# One way: takeoff, (transition1), climb, cruise, descent, (transition2), hover, landing
t_takeoff = 30 #sec
t_transition1 = 20 
t_climb = 40
t_cruise = flight_radius * 1000 / v_cruise #sec
t_descent = 60 
t_transition2 = 20
t_hover = 10
t_landing = 30 #sec

t_takeoff = flights * t_takeoff #sec
t_transition1 = flights * t_transition1 #sec
t_climb = flights * t_climb #sec
t_cruise = flights * t_cruise #sec
t_descent = flights * t_descent #sec
t_transition2 = flights * t_transition2 #sec
t_hover = flights * t_hover #sec
t_landing = flights * t_landing #sec


LD_ratio = 14 #Lift to drag ratio
eff_motor = 0.95 #Efficiency of the motor
eff_propeller_takeoff = 0.85 #Efficiency of the propeller

D_rotor_big = 3.1
D_rotor_small = 2.1
S_disks = ((D_rotor_big / 2)**2 * np.pi + (D_rotor_small / 2)**2 * np.pi) * 2 #m^2 (total disk area based on size requirements, can be changed later)
N_disks_takeoff = 4 #Number of disks

A_front = 4 #m^2 #m^2 (area of the front of the vehicle)
C_D = 0.4 #Drag coefficient of the front of the vehicle

#Configuration parameters:
tilt_wing = True

#MTOW guess
mtom = 2500 #kg #Maximum takeoff mass (MTOM) of the vehicle (GUESS)

#Constants:
g = 9.81 #m/s^2
rho_air = 1.225 #kg/m^3
density_batt_whkg = 350 #300Wh/kg (Chinese),
density_batt = density_batt_whkg*3600 #J/kg (density of the battery in J/kg)
DoD = 0.85 #Depth of discharge for the battery (100% to 15%)
blockage_factor_tiltwing = 1.0  #0.90 # Free area over total area for propellers in tilt-wing configuration
m_tilt_mech = 60 #kg for the tilt-wing mechanism
figure_of_merit = 0.8 # converting induced to total power for hovering

#Code parameters
mtow_prev = 0 #N
n = 1

mtow = mtom * g #N (initial guess for MTOW)

blockage_factor = blockage_factor_tiltwing

print("Disk area original: ", S_disks, "m^2")
S_disks = blockage_factor * S_disks #m2 (blockage factor for tilt-wing configuration)
print('Disk area adjusted for blockage factor: ', S_disks, "m^2")

def singleMotorHoverPower(T, S_disk):
    P_indueced = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disk)
    P_hover = P_induced / figure_of_merit
    P_hover = P_hover / eff_motor / eff_propeller_takeoff #W, with efficiency applied

while abs(mtow_prev - mtow) > 0.1 and n<1000:
    n+=1
    
    # Cruise power calculation
    T = mtow / LD_ratio #N
    P_cruise = T * v_cruise # / eff_motor

    T = mtow #N

    # Now consider each rotor separately for power calculation
    S_rotor_big = (D_rotor_big / 2)**2 * np.pi #m^2 (disk area per rotor)
    S_rotor_small = (D_rotor_small / 2)**2 * np.pi #m^2 (disk area per rotor)


    P_induced = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks) #W
    P_hover = P_induced / figure_of_merit
    P_hover = P_hover / eff_motor / eff_propeller_takeoff #W, with efficiency applied

    vihv = P_hover  / mtow #m/s (velocity induced by hover power)
    VcbVihv = v_takeoff / vihv
    PcbPhv = 0.5 * VcbVihv + np.sqrt((0.5 * VcbVihv)**2 + 1)
    P_takeoff = P_hover * PcbPhv #W (power required for takeoff)
    P_landing = P_hover #W (power required for landing)
    P_transition1 = P_takeoff #W (power required for transition 1), valid assumption
    P_transition2 = P_takeoff #W (power required for transition 2), valid assumption
    P_climb = P_cruise # NEED TO BE CHANGED WAITING FOR NUMBER FROM VIKKI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    P_descent = 0 # unpowered descent

    # Mission energy calculation
    # takeoff, (transition), climb, cruise, (transition), descent, hover, landing
    # All energies in Joules
    E_takeoff = P_hover * t_takeoff #J
    E_transition1 = P_hover * t_transition1 #J
    E_climb = P_climb * t_climb #J
    E_cruise = P_cruise * t_cruise #J
    E_descent = P_descent * t_descent #J
    E_transition2 = P_hover * t_transition2 #J
    E_hover = P_hover * t_landing #J
    E_landing = P_hover * t_landing #J

    E_total = E_takeoff + E_transition1 + E_climb + E_cruise + E_descent + E_transition2 + E_hover + E_landing #J
    E_copter_mode = (E_takeoff + E_landing + E_hover + E_transition1 + E_transition2) / 3600 / 1000#kWh (energy for the copter mode, takeoff, landing, hover and transitions)
    E_plane_mode = (E_climb + E_cruise + E_descent) / 3600 / 1000 #kWh (energy for the plane mode, climb, cruise and descent)

    # Energy source mass as a function of type, energy capacity and output power
    m_powersource = E_total / density_batt / DoD #kg

    m_equipment = 450/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
    m_structure = 1000/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
    m_structure = m_structure + m_tilt_mech #kg for the tilt-wing mechanism
    P_cruise_single = P_cruise / N_disks_takeoff / 1000 #KW, again for integrated propulsion takeoff=cruise for N
    P_hover_single = P_hover / N_disks_takeoff / 1000 #KW (energy per disk)
    m_propulsion = calculatePropulsionMass(P_cruise_single, P_hover_single, N_disks_takeoff, D_rotor_big) #kg (mass of the propulsion system)
    m_eom = m_powersource + m_equipment + m_structure + m_propulsion

    mtow_prev = mtow #N
    mtow = (m_payload + m_powersource + m_eom) * g #N

print("---------------------------------")
print("Number of iterations: ", n)
print('Single motor power takeoff: ', P_hover_single, "kW")
print("P_hover: ", P_hover / 1000, "kW")
print("P_cruise: ", P_cruise / 1000, "kW")
print("Energy in multicopter mode:", E_copter_mode, "Wh")
print("Energy in plane mode:", E_plane_mode, "Wh")
print("Total energy: ", E_total / 3600 / 1000, "kWh")

print("Final MTOW: ", mtow / g, "kg")
print("Final Propulsion mass: ", m_propulsion, "kg")
print("Final Power Source mass: ", m_powersource, "kg")
print("Final EOM mass ratio: ", m_eom*g/mtow*100, "%")
#print(P_hover/1000)