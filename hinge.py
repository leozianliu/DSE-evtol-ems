import numpy as np

### CALCULATING GEAR 
#inputs
#rotation around the wingbox
# Constants
angular = 90  # degrees
wingbox_diameter = 25 #cm
radius = 25 / 2 * 10  # mm
gear_thickness = 1.5 * 10  # mm
module = 6  #for heavy duty gears, it is around 6 to 10 mm 
helix_angle = 20 * np.pi / 180  # radians (not used)
R = radius + gear_thickness  # mm

def gear_sizing():
    d = R * 2  # pitch diameter in mm
    teeth = round(d / module)

    teeth_for_part = round(teeth * (90 / 360))

    pd = teeth / d  # diametral pitch (1/mm)
    pc = np.pi / pd  # circular pitch (mm)

    a = 1 / pd       # addendum - top part of the tooth 
    b = 1.25 / pd    # dedendum - bottom part of the tooth
    t = 0.5 * pc     # tooth thickness

    clearance = 0.25 * pd  # clearance

    #helical considerations
    axial_pitch = pc / np.cos(helix_angle)  # axial pitch (mm)

    return teeth, teeth_for_part, d, pd, pc, a, b, t, axial_pitch, clearance

# Run
gear = gear_sizing()
teeth, teeth_for_part, d, pd, pc, a, b, t, axial_pitch, clearance = gear

print(f'HELICAL GEAR -- Teeth: {teeth}, Teeth for part: {teeth_for_part:.2f}, Pitch Diameter: {d:.2f} mm, '
      f'Diametral Pitch: {pd:.4f}, Circular Pitch: {pc:.2f} mm, '
      f'Addendum: {a:.2f} mm, Dedendum: {b:.2f} mm, Tooth Thickness: {t:.2f} mm, '
      f'Clearance: {clearance:.2f}, Axial Pitch: {axial_pitch:.2f} mm, Helix Angle: 20 degrees')
print()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#WORM GEAR
circular_pitch = pc

circ =np.pi * R * 2
angle = 135
segment = circ *(angle/360)
#print(f'Segment: {segment:.2f} mm')

width = 100 #mm 
off = width* np.tan(helix_angle)
#print(f'Offset: {off:.2f} mm')

total_height = clearance + a + b 

worm_gear_thickness = gear_thickness/2
diameter_worm = width - total_height
#print(f'Diameter: {diameter:.2f} mm')

circ_worm = np.pi * diameter_worm
#print(f'Circumference of the worm: {circ_worm:.2f} mm')

off_worm = circ_worm* np.tan(helix_angle)
#print(f'Offset of the worm: {off_worm:.2f} mm')

pitch = circ / teeth
#print(f'Pitch: {pitch:.2f} mm')

print(f'WORM GEAR -- Segment: {segment:.2f} mm, Offset: {off:.2f} mm, Diameter: {diameter_worm:.2f} mm, \n'
      f'Circumference of the worm: {circ_worm:.2f} mm,Offset of the worm: {off_worm:.2f} mm, \n'
      f'Pitch: {pitch:.2f} mm')
print()


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##TORQUE CALCULATIONS
##Calculating the forcessssss
g = 9.81
m_wing = 250
W_wing = m_wing * g  # N
ang_acc = 20*np.pi/180  # rad/s^2

#Inertial torque
l_wing = 12
Inertial_wing = 1/3 * m_wing * l_wing**2  # kg*m^2
tau_I = Inertial_wing * ang_acc  # N*m

#Aero torque
V = 200/3.6 # m/s
rho = 1.225  # kg/m^3
q = 0.5 * rho * V**2  # dynamic pressure
cm = 0.044  # moment coefficient
c_r = 2  # root chord in m
tau_aero = q * c_r* cm  # N*m

#Grav torque
w_cg = 0.88 #from the leading edge
hinge = 0.4*c_r
tau_grav = W_wing * abs((w_cg - hinge))  # N*m
control_torque = 18000 #Nm #from Leonhard's calculations

total_torque = (tau_grav + control_torque + tau_I + tau_aero)/1000
print(f'TORQUE IT NEEDS TO SUSTAIN -- Gravitational Torque: {tau_grav:.2f} Nm, Control Torque: {control_torque:.2f} Nm, \n'
      f'Inertial Torque: {tau_I:.2f} Nm, Aerodynamic Torque: {tau_aero:.2f} Nm, Total Torque: {total_torque:.2f} kNm')
print()

stress = 650 # MPa
width = 50 # mm
b = width/1000 #face widht in m
m = module / 1000 # module in m
Y = (0.399+0.408)/2 #Lewis form factor for total number of teeth 
#Y = 0.21 #actual number of teeth for actual number of teeth 

W_t = stress * 10**6 * b * m * Y  # Tangential load in N
T = W_t * R /1000 / 1000# Torque in kNm first for m second fir k

print(f'TORQUE OF GEAR -- Tangential Load: {W_t:.2f} N, Torque Capacity of the Gear: {T:.2f} kNm')
print()


'''
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### MOTOR REQUIREMENTS
ang = 45
ang_r = 45/180 * np.pi
segment_v = (circ * ang / 360)/1000 #segement length for the 45 deg turn in m 
print(segment_v)

#power of the motor
time = 1 #in one second
ang_vel = ang_r / time #rad/s
gear_efficiency = 0.9
power_required = total_torque* ang_vel #in kW, here total torque is in kNm and 
power_from_motor = power_required / gear_efficiency


#RPM for the motor
N_helical = ang_vel * 30 / np.pi #RPM of the helical gear 
helical_gear_d = R * 2 + total_height
worm_gear_d = diameter_worm + total_height
circ_helical = helical_gear_d * np.pi 
circ_worm = worm_gear_d
circ_ratio = circ_helical / circ_worm #Ratio between the circumfernce of the helical and worm gear
gear_ratio = teeth / 1 
# N_motor = N_helical * circ_ratio #RPM requirement for the motor from helical rotation and ratio
N_motor = N_helical * gear_ratio

torque_motor = power_from_motor / (2*np.pi * N_motor/60) * 1000

print(f'MOTOR REQUIREMENTS-- Total Power Required: {power_from_motor:.2f} kW, for two gears: {power_from_motor/2:.2f} kW, \n'
      f'Helical spped: {N_helical:.2f} RPM, Speed Requirement for the Motor: {N_motor:.2f} RPM, \n '
      f'Torque Requirement: {torque_motor:.2f} Nm, for half: {torque_motor/2:.2f} Nm')
print()

print(f'FOR ELECTRO HYDRAULIC ACTUATOR -- Total Power Required: {power_from_motor:.2f} kW, for two gears: {power_from_motor/2:.2f} kW, \n'
      f'Sped requirement for the actuator: {segment_v:.2f} m/s')
'''

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#MOTOR REQUIREMENTS #2 
#Assuming there will be multiple gears 
n_gears = 2 
gear_ratio = teeth/1 
torque = total_torque / n_gears  #in kNm 
angle_to_rotate = 45 #degrees
w_out = angle_to_rotate * np.pi / 180 # angle of the helical gear in radians 
time = 1 #second
vel_ang = w_out/time 

power_output = torque * vel_ang # in kW
gear_efficiency = 0.9
power_required = power_output/gear_efficiency

w_in = w_out * gear_ratio # ang vel for worm
RPM_in = w_in * 60 / (2*np.pi) # worm rpm
torque_input = power_required / w_in * 1000 #Nm 

print(f'Helical gear velocity: {w_out:.2f} rad/s, Worm velocity: {w_in:.2f} rad/s, \n'
      f'RPM of the worm: {RPM_in:.2f} RPM, \n'
      f'Motor Power Required: {power_required:.2f} kW, \n'
      f'Motor Torque Required: {torque_input:.2f} kNm')
