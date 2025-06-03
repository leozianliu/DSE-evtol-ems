import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#############  Requirement definition  #############


p_required = np.pi/6   # [Rad/s]
V = 55.56    # [m/s]
b2 = 6 #  [m], aileron tip position 
b = 12 # [m],span
taper = 0.41 # Taper ratio [-]
C_r = 2 # Root chord [m]
S = 14 # Wing area [m^2]
c_la = np.pi*2 # Constant
c_d0 = 0.01 # Constant
d_a = np.pi/12 # Aileron deflection [Rad]
cl_delta_a = np.pi/2 # Constant 

######################################################

pb2v_ng = - p_required * b / (2*V)
C_lp = -(c_la+c_d0)*C_r*b*(1+3*taper)/(24*S)

def aileron_solver (b1, pb2v_ng,C_lp,d_a,cl_delta_a):
    Cl_da = (cl_delta_a*C_r)/(S*b)*((b2**2-b1**2)+ 4*(taper-1)*(b2**3-b1**3)/(3*b))
    f = pb2v_ng*C_lp/d_a - Cl_da

    return f

initial_value = 1
b1 = fsolve(aileron_solver, initial_value, args = (pb2v_ng,C_lp,d_a,cl_delta_a))

print('The position of the starting point of the ailelron is:', b1[0])
    

