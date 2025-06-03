import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#############  Requirement definition  #############


p_required = 1   # [Rad/s]
V = 1    # [m/s]
b2 = 1 #  [m], aileron tip position 
b = 1 # [m],span
taper = 0.41 # Taper ratio [-]
C_r = 1 # Root chord [m]
S = 14 # Wing area [m^2]
c_la = 1 # Constant
c_d0 = 1 # Constant
d_a = 1 # Aileron deflection [Rad/s]
cl_delta_a = 1 # Constant 

######################################################

pb2v_ng = - p_required * b / (2*V)
C_lp = -(c_la+c_d0)*C_r*b*(1+3*taper)/(24*S)

def aileron_solver (b1, pb2v_ng,C_lp,d_a,cl_delta_a):
    Cl_da = (cl_delta_a*C_r)/(S*b)*((b2**2-b1**2)+ 4*(taper-1)*(b2**3-b1**3)/(3*b))
    f = pb2v_ng*C_lp/d_a - Cl_da

    return f

initial_value = 1
b1 = fsolve(aileron_solver, initial_value, args = (pb2v_ng,C_lp,d_a,cl_delta_a))

print(b1[0])
    

