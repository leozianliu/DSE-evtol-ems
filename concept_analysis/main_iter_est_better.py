import numpy as np

#Mission parameters:
m_payload = 400 #kg
flight_radius = 50 #km (radius)
v_cruise = 200/3.6 #m/s
v_takeoff = 5 #m/s

flights = 2 #Number of flights per mission
# One way: takeoff, (transition1), climb, cruise, descent, (transition2), hover, landing
t_takeoff = 30 #sec
t_transition1 = 20 
t_climb = 40
t_descent = 40
t_cruise = (flight_radius * 1000 / v_cruise) - t_climb - t_descent #sec
t_transition2 = 20
t_hover = 10
t_landing = 30 #sec

dh_climb = 150 #m (climb height in fixed wing mode)
roc = dh_climb / t_climb #m/s (rate of climb)

t_takeoff = flights * t_takeoff #sec
t_transition1 = flights * t_transition1 #sec
t_climb = flights * t_climb #sec
t_cruise = flights * t_cruise #sec
t_descent = flights * t_descent #sec
t_transition2 = flights * t_transition2 #sec
t_hover = flights * t_hover #sec
t_landing = flights * t_landing #sec

#Aircraft parameters:
eff_motor = 0.95 #Efficiency of the motor
eff_small_propeller_takeoff = 0.75 #Efficiency of the propeller during takeoff
eff_small_propeller_cruise = 0.75 #Efficiency of the propeller
eff_big_propeller_takeoff = 0.80 #Efficiency of the propeller during takeoff
eff_big_propeller_cruise = 0.80 #Efficiency of the propeller during cruise

D_rotor_big = 3.20
D_rotor_small = 2.25
S_disks = ((D_rotor_big / 2)**2 * np.pi + (D_rotor_small / 2)**2 * np.pi) * 2 #m^2 (total disk area based on size requirements, can be changed later)

# Aerodynamic parameters:
Cl_wing_cruise = 0.8 #Lift coefficient of the wing in cruise
S_wing = 14
Cd_wing = 0.035
S_fuselage = 3.43
Cd_fuselage = 0.144
S_tail = 0
Cd_tail = 0
S_strut = 0
Cd_strut = 0
S_gear = 14
Cd_gear = 0.007 #Drag coefficient of the landing gear

LD_ratio_aircraft = (S_wing * Cl_wing_cruise) / (S_wing * Cd_wing + S_fuselage * Cd_fuselage + S_tail * Cd_tail + S_strut * Cd_strut + S_gear * Cd_gear) #Lift to drag ratio of the aircraft
print("Lift to drag ratio of the aircraft: ", LD_ratio_aircraft)

#Configuration parameters:
density_batt_whkg = 300 #300Wh/kg (Chinese),
density_batt = density_batt_whkg*3600 #J/kg (density of the battery in J/kg)
DoD = 0.8 #Depth of discharge (DoD), same effect as battery degradation, 80% of the battery capacity is used
# blockage_factor_tiltwing = 1.0  #0.90 # Free area over total area for propellers in tilt-wing configuration
m_tilt_mech = 100 #kg for the tilt-wing mechanism
figure_of_merit = 0.8 # converting induced to total power for hovering

#Weight safety factor:
safety_factor_battery_weight = 1.1 #Safety factor for MTOW estimation
safety_factor_propulsion_weight = 1 #Safety factor for MTOW estimation

#Code parameters
mtow_prev = 0 #N
n = 1

#MTOW guess
mtom = 2500 #kg #Maximum takeoff mass (MTOM) of the vehicle (GUESS)

#Constants:
g = 9.81 #m/s^2
rho_air = 1.225 #kg/m^3

mtow = mtom * g #N (initial guess for MTOW)

# blockage_factor = blockage_factor_tiltwing

print("Disk area original: ", S_disks, "m^2")
# S_disks = blockage_factor * S_disks #m2 (blockage factor for tilt-wing configuration)
# print('Disk area adjusted for blockage factor: ', S_disks, "m^2")

def singleMotorCruisePower(weight, v_cruise, eff_motor, eff_propeller):
    T_cruise = weight / LD_ratio_aircraft
    P_cruise = T_cruise * v_cruise
    P_cruise = P_cruise / eff_motor / eff_propeller #W, with efficiency applied
    return P_cruise

def singleMotorClimbPower(roc, weight, v_air):
    P_climb = roc * weight + 0.5 * rho_air * v_air**3 * (S_wing * Cd_wing + S_fuselage * Cd_fuselage + S_tail * Cd_tail + S_strut * Cd_strut + S_gear * Cd_gear)
    return P_climb

def singleMotorHoverPower(T, S_disk, eff_motor, eff_propeller):
    P_induced = T**(3/2) / np.sqrt(2 * rho_air * S_disk)
    P_hover = P_induced / figure_of_merit
    P_hover = P_hover / eff_motor / eff_propeller #W, with efficiency applied
    return P_hover, P_induced

def singleMotorTakeoffPower(T, v_takeoff, S_disk, eff_motor, eff_propeller):
    P_hover, P_induced = singleMotorHoverPower(T, S_disk, eff_motor, eff_propeller)
    vihv = P_induced / T #m/s (velocity induced by hover power)
    VcbVihv = v_takeoff / vihv
    PcbPhv = 0.5 * VcbVihv + np.sqrt((0.5 * VcbVihv)**2 + 1)
    P_takeoff = P_hover * PcbPhv #W (power required for takeoff)
    return P_takeoff

def powerSingleRotorAllPhases(T, v_takeoff, v_cruise, S_rotor_single, eff_motor, eff_propeller_takeoff, eff_propeller_cruise):
    # Power calculations for all phases
    P_takeoff = singleMotorTakeoffPower(T, v_takeoff, S_rotor_single, eff_motor, eff_propeller_takeoff) #W (power required for takeoff for single rotor)
    P_transition1 = None #W (power required for transition 1 for single rotor), valid assumption
    P_climb = None # NEED TO BE CHANGED WAITING FOR NUMBER FROM VIKKI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    P_cruise = singleMotorCruisePower(T, v_cruise, eff_motor, eff_propeller_cruise)
    P_descent = 0 # unpowered descent
    P_transition2 = None #W (power required for transition 2 for single rotor), valid assumption
    P_hover, _ = singleMotorHoverPower(T, S_rotor_single, eff_motor, eff_propeller_takeoff) #W (power required for hover for single rotor)
    P_landing = P_hover #W (power required for landing for single rotor)

    # Overwrite
    P_climb = singleMotorClimbPower(roc, T, v_cruise)
    P_transition1 = P_transition2 = P_hover

    return [P_takeoff, P_transition1, P_climb, P_cruise, P_descent, P_transition2, P_hover, P_landing]

def calculateMissionEnergy(P_takeoff, P_transition1, P_climb, P_cruise, P_descent, P_transition2, P_hover, P_landing, t_takeoff, t_transition1, t_climb, t_cruise, t_descent, t_transition2, t_hover, t_landing):
    # Mission energy calculation
    # takeoff, (transition), climb, cruise, (transition), descent, hover, landing
    E_takeoff = P_takeoff * t_takeoff #J
    E_transition1 = P_transition1 * t_transition1 #J
    E_climb = P_climb * t_climb #J
    E_cruise = P_cruise * t_cruise #J
    E_descent = P_descent * t_descent #J
    E_transition2 = P_transition2 * t_transition2 #J
    E_hover = P_hover * t_hover #J
    E_landing = P_landing * t_landing #J

    E_total_J = E_takeoff + E_transition1 + E_climb + E_cruise + E_descent + E_transition2 + E_hover + E_landing #J
    E_copter_mode_Wh = (E_takeoff + E_landing + E_hover + E_transition1 + E_transition2) / 3600 #Wh (energy for the copter mode, takeoff, landing, hover and transitions)
    E_plane_mode_Wh = (E_climb + E_cruise + E_descent) / 3600 #Wh (energy for the plane mode, climb, cruise and descent)
    return E_total_J, E_copter_mode_Wh, E_plane_mode_Wh

def calculateSingleRotorMass(cont_power, D_rotor):#def calculatePropulsionMass(cont_power, peak_power, cont_torque, peak_torque, n_motor):
    # Motor & driver mass as a function of power and torque
    #power and torque inputs PER MOTOR
    #power [W]
    #torque [Nm]
    #masse [kg]
    #diameter [m]

    cont_power = cont_power / 1000 # W to kW

    m_motor_cont_power = 0.1156 * cont_power + 4.52
    m_motor_cont_torque = 0
    m_motor = max(m_motor_cont_power, m_motor_cont_torque)#, m_motor_cont_torque, m_motor_peak_torque)
    m_inverter = 0 #0.0187*cont_power + 0.433 # motor mass includes inverter mass
    m_propeller = 10 * (D_rotor/2.3) # Volocopter's rotor mass is 7.5 kg for 2.3 m diameter rotor, assumes a linear scaling with diameter

    return (m_motor + m_inverter + m_propeller)


while abs(mtow_prev - mtow) > 0.1 and n<1000:
    # Iterations
    n+=1
    T = mtow #N
    T_small_single = T / 6
    T_big_single = T_small_single * 2

    # Now consider each rotor separately for power calculation
    S_rotor_big = (D_rotor_big / 2)**2 * np.pi #m^2 (disk area per rotor)
    S_rotor_small = (D_rotor_small / 2)**2 * np.pi #m^2 (disk area per rotor)

    # P_arr = [P_takeoff, P_transition1, P_climb, P_cruise, P_descent, P_transition2, P_hover, P_landing]
    # Single small rotor power calculations
    P_arr_small = powerSingleRotorAllPhases(
        T_small_single, v_takeoff, v_cruise, S_rotor_small, eff_motor, eff_small_propeller_takeoff, eff_small_propeller_cruise)
    P_takeoff_small, P_transition1_small, P_climb_small, P_cruise_small, P_descent_small, P_transition2_small, P_hover_small, P_landing_small = P_arr_small

    # Single big rotor power calculations
    P_arr_big = powerSingleRotorAllPhases(
        T_big_single, v_takeoff, v_cruise, S_rotor_big, eff_motor, eff_big_propeller_takeoff, eff_big_propeller_cruise)
    P_takeoff_big, P_transition1_big, P_climb_big, P_cruise_big, P_descent_big, P_transition2_big, P_hover_big, P_landing_big = P_arr_big

    P_takeoff_tot = P_arr_small[0] * 2 + P_arr_big[0] * 2 #W (total power required for takeoff)
    P_transition1_tot = P_arr_small[1] * 2 + P_arr_big[1] * 2 #W (total power required for transition)
    P_climb_tot = P_arr_small[2] * 2 + P_arr_big[2] * 2 #W (total power required for climb)
    P_cruise_tot = P_arr_small[3] * 2 + P_arr_big[3] * 2 #W (total power required for cruise)
    P_descent_tot = P_arr_small[4] * 2 + P_arr_big[4] * 2 #W (total power required for descent)
    P_transition2_tot = P_arr_small[5] * 2 + P_arr_big[5] * 2 #W (total power required for transition)
    P_hover_tot = P_arr_small[6] * 2 + P_arr_big[6] * 2 #W (total power required for hover)
    P_landing_tot = P_arr_small[7] * 2 + P_arr_big[7] * 2 #W (total power required for landing)

    # E_total, E_copter_mode_Wh, E_plane_mode_Wh, FOR SINGLE ROTOR SMALL
    E_total, E_copter_mode_Wh, E_plane_mode_Wh = calculateMissionEnergy(
        P_takeoff_tot, P_transition1_tot, P_climb_tot,
        P_cruise_tot, P_descent_tot, P_transition2_tot, P_hover_tot, P_landing_tot, t_takeoff, t_transition1, t_climb,
        t_cruise, t_descent, t_transition2, t_hover, t_landing
    ) #J, Wh, Wh
    
    if n <= 5 or n % 10 == 0:  # Print first 5 iterations then every 10th
        print(f"Iteration {n}: MTOW = {mtow / g:.1f} kg, Total Energy = {E_total/ 3600 / 1000:.2f} kWh")

    P_takeoff_small2xThrust = singleMotorTakeoffPower(T_small_single * 2, v_takeoff, S_rotor_small, eff_motor, eff_small_propeller_takeoff) # Size for 2x hover thrust
    P_takeoff_big2xThrust = singleMotorTakeoffPower(T_big_single * 2, v_takeoff, S_rotor_big, eff_motor, eff_big_propeller_takeoff) # Size for 2x hover thrust

    # Energy source mass as a function of type, energy capacity and output power
    m_powersource = E_total / density_batt / DoD #kg
    m_powersource = m_powersource * safety_factor_battery_weight # for cooling and casing masses
    m_equipment = 450/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
    m_structure = 1000/3700 * mtow/g #kg (ratio for battery powered vehicles from https://www.researchgate.net/publication/318235979_A_Study_in_Reducing_the_Cost_of_Vertical_Flight_with_Electric_Propulsion/figures)
    m_structure = m_structure + m_tilt_mech #kg for the tilt-wing mechanism
    
    m_propulsion = calculateSingleRotorMass(P_takeoff_small2xThrust, D_rotor_small) * 2\
        + calculateSingleRotorMass(P_takeoff_big2xThrust, D_rotor_big) * 2 #kg (mass of the propulsion system)
    m_propulsion = m_propulsion * safety_factor_propulsion_weight # Safety factor for propulsion mass
    m_eom = m_powersource + m_equipment + m_structure + m_propulsion

    mtow_prev = mtow #N
    mtow = (m_payload + m_eom) * g #N


print("---------------------------------")
print("Number of iterations: ", n)
print('Single small motor power takeoff: ', P_takeoff_small / 1000, "kW")
print('Single big motor power takeoff: ', P_takeoff_big / 1000, "kW")
print('Single small motor thrust 2x hover: ', T_small_single / 1000, "kN")
print('Single big motor thrust 2x hover: ', T_big_single / 1000, "kN")
print("\n")
print("P_takeoff: ", P_takeoff_tot / 1000, "kW")
print("P_transition1: ", P_transition1_tot / 1000, "kW")
print("P_climb: ", P_climb_tot / 1000, "kW")
print("P_cruise: ", P_cruise_tot / 1000, "kW")
print("P_descent: ", P_descent_tot / 1000, "kW")
print("P_transition2: ", P_transition2_tot / 1000, "kW")
print("P_hover: ", P_hover_tot / 1000, "kW")
print("P_landing: ", P_landing_tot / 1000, "kW")
print("\n")
print("Energy in multicopter mode:", E_copter_mode_Wh, "Wh")
print("Energy in plane mode:", E_plane_mode_Wh, "Wh")
print("Total energy: ", E_total / 3600 / 1000, "kWh")
print("\n")
print("Final MTOW: ", mtow / g, "kg")
print("Final Propulsion mass: ", m_propulsion, "kg")
print("Final Power Source mass: ", m_powersource, "kg")
print("Final EOM mass ratio: ", m_eom*g/mtow*100, "%")

# Convergence check
if n >= 1000:
    print("WARNING: Maximum iterations reached without convergence!")
else:
    print(f"Converged after {n} iterations")