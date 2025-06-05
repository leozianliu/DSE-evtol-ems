import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#############  Requirement definition  #############


p_required = np.pi/6   # [Rad/s]
q_required = np.pi/6   # [Rad/s]
r_required = np.pi/6   # [Rad/s]

V = 55.56    # [m/s]
Vh_V = 0.8 # [-] 
Vv_V = 0.8    #  [-] 
b2_wing = 6 #  [m], aileron tip position 
b_wing = 12 # [m],span
taper_wing = 0.41 # Taper ratio [-]
C_r_wing = 2 # Root chord [m]
S_wing = 14 # Wing area [m^2]
c_la_wing = np.pi*2 # Constant
c_d0_wing = 0.01 # Constant
d_a_wing = np.pi/12 # Aileron deflection [Rad]
cl_delta_a_wing = np.pi/2 # Constant


b2_rudder = 2.5
b_rudder = 3
taper_rudder = 1
C_r_rudder = 2
S_rudder = 5
c_n_beta_rudder = np.pi*2 # Constant
c_d0_rudder = 0.01 # Constant
d_r_rudder = np.pi/12 # Aileron deflection [Rad]
cn_delta_r_rudder = np.pi/2 # Constant

b2_elevator = 2.5
b_elevator = 3
taper_elevator = 1
C_r_elevator = 2
S_elevator = 5
c_m_a_elevator = np.pi*2 # Constant
c_d0_elevator = 0.01 # Constant
d_e_elevator = np.pi/12 # Aileron deflection [Rad]
cm_delta_e_elevator = np.pi/2 # Constant
c_elevator = 1


######################################################

pb2v_ng_wing = - p_required * b_wing / (2*V)
C_lp_wing = -(c_la_wing+c_d0_wing)*C_r_wing*b_wing*(1+3*taper_wing)/(24*S_wing)

rb2v_ng_rudder = - r_required * b_rudder / (2*V*Vv_V)
C_nr_rudder = -(c_n_beta_rudder+c_d0_rudder)*C_r_rudder*b_rudder*(1+3*taper_rudder)/(24*S_rudder)

v_c_ng_rudder = - c_elevator / (V*Vh_V)
C_mq_elevator = -(c_m_a_elevator+c_d0_elevator)*C_r_elevator*b_elevator*(1+3*taper_elevator)/(24*S_elevator)

def Solver (b1, pb2v_ng,C_lp,d_a,cl_delta_a):
    Cl_da = (cl_delta_a*C_r_wing)/(S_wing*b_wing)*((b2_wing**2-b1**2)+ 4*(taper_wing-1)*(b2_wing**3-b1**3)/(3*b_wing))
    f = pb2v_ng*C_lp/d_a - Cl_da

    return f

initial_value = 1
b1_wing = fsolve(Solver, initial_value, args = (pb2v_ng_wing,C_lp_wing,d_a_wing,cl_delta_a_wing))
b1_rudder = fsolve(Solver, initial_value, args = (rb2v_ng_rudder ,C_nr_rudder,d_r_rudder,cn_delta_r_rudder))
b1_elevator = fsolve(Solver, initial_value, args = (v_c_ng_rudder ,C_mq_elevator ,d_e_elevator,cm_delta_e_elevator))

print('The position of the starting point of the ailelron is:', b1_wing[0])
print('The position of the starting point of the rudder is:', b1_rudder[0])
print('The position of the starting point of the elevator is:', b1_elevator[0])

#####################################################
    

