import numpy as np
from wing_design import design_wing

exp = True
lift = 2600 * 9.81  # Maximum Takeoff Weight in Newtons
V = 200 / 3.6  # Velocity in m/s
b = 12 # Span in meters
tratio = 1
t_r = 0.28 * tratio # m
t_m = 0.256 * tratio # Thickness at first engine in meters
t_t = 0.14 * tratio # m
anhedral = 20
stallmargin = 0
# These values should result in an L/D of 89.0976735576122

area, lift, drag, wingdata = design_wing(exp, lift, V, b, t_r, t_m, t_t, anhedral, stallmargin)
# Contains for each section: Airfoil name, cl_(clcdmax), chord, t(verification), alpha_clcd, cd_(clcd_max), stall_margin (deg), alpha_induced, real_alpha, real_CD
print(wingdata)
print(f"Area: {area} m2, lift: {lift} N, drag: {drag} N, L/D: {lift/drag}")



# Idea; try a couple of thickness and velocity values
# For each combination, find the L/D and the weight of the wing
# Using this values, find MTOW
# Select the combination that gave the lowest MTOW
