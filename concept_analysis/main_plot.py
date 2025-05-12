import numpy as np
import matplotlib.pyplot as plt
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

LD_ratio = 12 #Lift to drag ratio
eff_motor = 0.95 #Efficiency of the motor
eff_propeller = 0.85 #Efficiency of the propeller

# Note: If integrated propulsion is used, N_disks_cruise is not used
S_disks = 35.33 #m^2 (total disk area based on size requirements, can be changed later)
N_disks_cruise = 0 #Number of disks (between 2 and 4)
N_disks_takeoff = 6 #Number of disks (between 4 and 8)
S_rotor = S_disks / N_disks_takeoff #m^2 (disk area per rotor)
D_rotor = 2 * np.sqrt(S_rotor / np.pi) #m (diameter of the rotor)
A_front = 5 #m^2 #m^2 (area of the front of the vehicle)
C_D = 0.4 #Drag coefficient of the front of the vehicle

#Configuration parameters:
battery = True
wing = True
integrated_prop = True #Use the same motors for hover and cruise
tilt_wing = True

# Constants
g = 9.81  # m/s^2
rho_air = 1.225  # kg/m^3
separate_prop_extra_drag_factor = 1.10
LD_reduction_factor = 1 / separate_prop_extra_drag_factor
blockage_factor_tiltwing = 0.90
blockage_factor_tiltrotor = 0.78
blockage_factor_sepprop = 0.86
m_tilt_mech = 60  # kg

# Setup blockage factor
if integrated_prop:
    if tilt_wing:
        blockage_factor = blockage_factor_tiltwing
    else:
        blockage_factor = blockage_factor_tiltrotor
else:
    blockage_factor = blockage_factor_sepprop

S_disks = blockage_factor * S_disks  # Adjusted disk area

# Battery energy density range (Wh/kg)
battery_densities = np.linspace(200, 1000, 50)  # From 200 to 1000 Wh/kg
mtow_results = []

for density_batt_whkg in battery_densities:
    # Calculate effective energy density (J/kg) with DoD and efficiency
    density_batt = density_batt_whkg*3600*0.80 #degrade to 80% of the battery capacity

    
    # Iterative MTOW calculation
    mtow_prev = 0
    mtow = 2500 * g  # Initial guess in N
    n = 0
    
    while abs(mtow_prev - mtow) > 0.1 and n < 1000:
        n += 1
        mtow_prev = mtow
        
        # Energy calculations
        if wing:
            if integrated_prop:
                LD_ratio = 13
                if tilt_wing:
                    LD_ratio = 15
                T = mtow / LD_ratio
            else:
                T = mtow / (LD_ratio * LD_reduction_factor)
            P_cruise = T * v_cruise
        else:
            F_sideways = 0.5 * rho_air * (v_cruise**2) * A_front * C_D
            T = np.sqrt(mtow**2 + F_sideways**2)
            P_cruise = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks)
        E_cruise = P_cruise * t_cruise

        T = mtow
        P_induced = 1.15 * T**(3/2) / np.sqrt(2 * rho_air * S_disks)
        P_hover = (P_induced) / eff_motor / eff_propeller
        E_hover = P_hover * t_hover

        E_total = E_cruise + E_hover
        m_powersource = E_total / density_batt

        # Mass calculations
        if wing and integrated_prop:
            m_equipment = 450/3700 * mtow/g
            m_structure = 1000/3700 * mtow/g
            if tilt_wing:
                m_structure += m_tilt_mech
            P_cruise_single = P_cruise / N_disks_takeoff / 1000
            P_hover_single = P_hover / N_disks_takeoff / 1000
            m_propulsion = calculatePropulsionMass(P_cruise_single, P_hover_single, N_disks_takeoff, D_rotor)
            if not tilt_wing:
                m_propulsion += N_disks_takeoff * 10
            m_eom = m_powersource + m_equipment + m_structure + m_propulsion
        elif wing and (not integrated_prop):
            m_equipment = 450/3700 * mtow/g
            m_structure = 1000/3700 * mtow/g
            P_cruise_single = P_cruise / N_disks_cruise / 1000
            P_hover_single = P_hover / N_disks_takeoff / 1000
            m_propulsion = (calculatePropulsionMass(P_cruise_single, 0, N_disks_cruise, D_rotor) + 
                            calculatePropulsionMass(0, P_hover_single, N_disks_takeoff, D_rotor))
            m_eom = m_powersource + m_equipment + m_structure + m_propulsion
        else:
            m_eom = 1500/3450 * mtow/g

        mtow = (m_payload + m_powersource + m_eom) * g

    mtow_results.append(mtow / g)  # Convert to kg for results

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(battery_densities, mtow_results, 'b-', linewidth=2)
plt.axhline(y=mtow_max, color='r', linestyle='--', label='Max MTOW Limit')
plt.xlabel('Battery Energy Density (Wh/kg)')
plt.ylabel('MTOW (kg)')
plt.title('MTOW vs Battery Energy Density')
plt.grid(True)
plt.legend()
plt.show()