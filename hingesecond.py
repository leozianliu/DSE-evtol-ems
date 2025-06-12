import numpy as np

def RPM(v): 
    RPM = v * 60 / (2* np.pi)
    return RPM

def deg_to_rad(deg):
    rad = deg * np.pi /180
    return rad

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#HELICAL WORM GEAR 
angular = 90  # degrees
wingbox_diameter = 27 #cm
radius = wingbox_diameter / 2 * 10  # mm
gear_thickness = 1 * 10  # mm
module = 6  #for heavy duty gears, it is around 6 to 10 mm 
helix_angle = deg_to_rad(20)
R = radius + gear_thickness  # mm (Pitch radius basically)
pressure_angle = 20 #degrees

#calculating gear dimensions
d = R * 2  # pitch diameter in mm (where the teeth of the gears starti)
teeth = round(d / module)

teeth_for_part = round(teeth * (90 / 360)) #don't really use this

pd = teeth / d  # diametral pitch (1/mm)
pc = np.pi / pd  # circular pitch (mm)

a = 1 / pd       # addendum - top part of the tooth 
b = 1.25 / pd    # dedendum - bottom part of the tooth
t = 0.5 * pc     # tooth thickness

clearance = 0.25 * pd  # clearance

#helical considerations
axial_pitch = pc / np.cos(helix_angle)  # axial pitch (mm)

width = 80 #mm of the helical gear
off = width* np.tan(helix_angle)

print(f'HELICAL GEAR -- \n'
      f'Teeth: {teeth}, \n'
      f'Teeth for part: {teeth_for_part:.2f}, \n'
      f'Pitch Diameter: {d:.2f} mm, \n'
      f'Diametral Pitch: {pd:.4f}, \n' 
      f'Circular Pitch: {pc:.2f} mm, \n'
      f'Offset: {off:.2f} mm, \n'
      f'Axial Pitch: {axial_pitch:.2f} mm, \n'
      f'Helix Angle: {helix_angle:.2f} degrees')
print()
print(f'TOOTH DIMENSIONS: --'
      f'Addendum: {a:.2f} mm, \n'
      f'Dedendum: {b:.2f} mm, \n'
      f'Tooth Thickness: {t:.2f} mm, \n'
      f'Clearance: {clearance:.2f}, \n'
      f'Totah Height: {a + b + clearance:.2f} mm')
print()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#WORM
circular_pitch = pc
circ =np.pi * R * 2
angle = 135 #degrees
segment = circ *(angle/360)

total_height = clearance + a + b 

worm_gear_thickness = gear_thickness/2
diameter_worm = width - total_height*2

circ_worm = np.pi * diameter_worm

off_worm = circ_worm* np.tan(helix_angle)

pitch = circ / teeth

print(f'WORM GEAR -- \n' 
      f'Segment: {segment:.2f} mm, \n' 
      f'Diameter: {diameter_worm:.2f} mm, \n'
      f'Circumference of the worm: {circ_worm:.2f} mm, \n'
      f'Offset of the worm: {off_worm:.2f} mm, \n'
      f'Pitch: {pitch:.2f} mm')
print()


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##TORQUE REQUIREMENTS
g = 9.81
m_wing = 190

#there is no angular acceleration -- no inertial torque
#Inertial torque
#r_inertia = R/1000 #turn to m
#ang_acc = # rad/s^2 #assuming this
#inertial_wing = 1/3 * m_wing * r_inertia**2  # kg*m^2
#torque_inertial = inertial_wing * ang_acc  # Nm -- this turns out to be negligible in the end

#aerodynamic torque
c_r = 1.4 #root chord [m]
c_hinge = 0.309 * c_r #center of the hinge from the leading edge [m]
c_cg = 0.584 # from the leading edge [m]
arm = abs(c_hinge - c_cg)

Cm = 0.491 
rho = 1.225 #kg/m^3
air_velocity = 35 #m/s
rotational_velocity = 45 * np.pi /180 * arm #m/s
aero_velocity = np.sqrt(air_velocity**2 + rotational_velocity**2) #m/s
S = 11.3
torque_aero = Cm * 1/2 * rho * aero_velocity**2 * S * c_r

#gravitational torque
w_wing = m_wing * g  # N
torque_gravitational = w_wing * arm  # N*m

control_torque = 21000 #Nm #from Leonhard's calculations

total_torque = (torque_gravitational + control_torque + torque_aero)/1000
print(f'TORQUE IT NEEDS TO SUSTAIN -- \n'
      f'Gravitational Torque: {torque_gravitational:.2f} Nm, \n'
      f'Control Torque: {control_torque:.2f} Nm, \n'
      f'Aerodynamic Torque: {torque_aero:.2f} Nm, \n'
      f'Total Torque: {total_torque:.2f} kNm')
print()


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#TORQUE ON GEAR
stress = 650 # MPa
face = width/1000 #face widht in m #width defined in WORM section
m = module / 1000 # module in m
Y = (0.399+0.408)/2 #Lewis form factor for total number of teeth -- 47 teeth

tangential_load = stress * 10**6 * face * m * Y  # Tangential load in N
torque_on_gear = tangential_load * R /1000 / 1000# Torque in kNm first for m second fir k

print(f'TORQUE OF ONE GEAR -- \n'
      f'Tangential Load: {tangential_load:.2f} N, \n'
      f'Torque Capacity of the Gear: {torque_on_gear:.2f} kNm, \n'
      f'Torque for Two Gears: {torque_on_gear*2:.2f} kNm')
print()


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#MOTOR REQUIREMENTS 
#Assuming there will be multiple gears 
n_gears = 2 #total number of worm gears
threads_worm = 3
gear_ratio = teeth / threads_worm
torque = total_torque / n_gears  #in kNm 
angle_to_rotate = 45 #degrees
helical_angle_rad = deg_to_rad(angle_to_rotate) # angle of the helical gear in radians 
time = 1 #second
helical_angular_velocity = helical_angle_rad/time #rad/s for the helical gear 
worm_angle_rad = helical_angle_rad * gear_ratio #how many radians does it need to turn for 45 deg in helical
worm_angular_velocity = worm_angle_rad/time #rad/s for the worm 
worm_velocity = worm_angular_velocity * diameter_worm #deg/s for the worm for the sim

power_requirement = torque * helical_angular_velocity #torque in kNm, so power in kW
gear_efficiency = 0.9  
motor_power_requirement = power_requirement / gear_efficiency
worm_RPM = RPM(worm_angle_rad) #RPM of the worm 
motor_torque_requirement = motor_power_requirement / worm_angular_velocity  #the power is in kW so torque will be in kNmw


print(f'MOTOR REQUIREMENTS: -- \n'
      f'Helical velocity: {helical_angular_velocity:.2f} rad/s, \n'
      f'Worm velocity: {worm_velocity:.2f} deg/s, \n'
      f'Worm velocity: {worm_angular_velocity:.2f} rad/s, \n'   
      f'RPM of the motor: {worm_RPM:.2f} RPM \n'
      f'Motor Power Requirement: {motor_power_requirement:.2f} kW, \n'
      f'Motor Torque Requirement: {motor_torque_requirement:.2f} kNm, \n')


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#GEARBOX sizing
diameter_motor = 159 #mm
RPM_motor = 840 #RPM 
torque_motor = 174 #Nm 
diameter_worm = diameter_worm #mm 
RPM_worm = worm_RPM

gearbox_ratio = RPM_motor / worm_RPM #7:1 RPM motor to RPM worm
print(f'GEARBOX SIZING: -- \n'
      f'Gearbox ratio: {gearbox_ratio:.2f} \n')


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Total hinge sizing
height = 20 + wingbox_diameter*10 + total_height*2 + gear_thickness*2 + diameter_worm + total_height*2
print(f'Total hinge height: {height:.2f} mm')