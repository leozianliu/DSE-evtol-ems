import numpy as np


g = 9.81
m_wing = 100
d_wingbox = 25 / 100 #m 

# weight force of the wing
Fg = m_wing * g
r_g = d_wingbox / 2

#torque the wing creates about the hinge 
M_wing = Fg * r_g

#torque the hinge needs to counteract
M_hinge = M_wing

'''will need to find 
M_hydraulic = F_hydraulic * length of piston arm
F_hydraulic = hydraulic pressure * area of piston cross section'''


### CALCULATING GEAR 
#inputs
#rotation around the wingbox
# Constants
angular = 90  # degrees
radius = 25 / 2 * 10  # mm
gear_thickness = 1.5 * 10  # mm
module = 8  #for heavy duty gears, it is around 6 to 10 mm 
helix_angle = 20 * np.pi / 180  # radians (not used)
R = radius + gear_thickness  # mm

def gear_sizing():
    d = R * 2  # pitch diameter in mm
    teeth = round(d / module)

    teeth_for_part = round(teeth * (90 / 360))

    pd = teeth / d  # diametral pitch (1/mm)
    pc = np.pi / pd  # circular pitch (mm)

    a = 1 / pd       # addendum
    b = 1.25 / pd    # dedendum
    t = 0.5 * pc     # tooth thickness

    return teeth, teeth_for_part, d, pd, pc, a, b, t

# Run
gear = gear_sizing()
teeth, teeth_for_part, d, pd, pc, a, b, t = gear

print(f'Teeth: {teeth}, Teeth for part: {teeth_for_part:.2f}, Pitch Diameter: {d:.2f} mm, '
      f'Diametral Pitch: {pd:.4f}, Circular Pitch: {pc:.2f} mm, '
      f'Addendum: {a:.2f} mm, Dedendum: {b:.2f} mm, Tooth Thickness: {t:.2f} mm')