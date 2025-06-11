import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#############  Requirement definition  #############


p_required = np.pi/6   # [Rad/s]
q_required = np.pi/6   # [Rad/s]
r_required = np.pi/6   # [Rad/s]

V = 55.56/2    # [m/s]
Vh_V = 1 # [-] 
Vv_V = 1    #  [-] 
b2_wing = 6 #  [m], aileron tip position 
b_wing = 12 # [m],span
taper_wing = 0.41 # Taper ratio [-]
C_r_wing = 2 # Root chord [m]
S_wing = 14 # Wing area [m^2]
c_la_wing =  2 * np.pi# Constant
c_d0_wing = 0.01 # Constant
d_a_wing = np.radians(30) # Aileron deflection [Rad]
cl_delta_a_wing = 3.81 * np.cos(np.radians(15)) # Constant


b2_rudder = 1.67
b_rudder = 1.67
taper_rudder = 0.67
C_r_rudder = 1.794
S_rudder =  2.5
c_n_beta_rudder = 2 * np.pi # Constant
c_d0_rudder = 0.01 # Constant
d_r_rudder = np.radians(30) # Aileron deflection [Rad]
cn_delta_r_rudder = 3.81 # Constant

b2_elevator = 4.14/2
b_elevator = 4.14
taper_elevator = 0.4
C_r_elevator = 1.196
S_elevator = 3.5
c_m_a_elevator = np.pi*2 # Constant
c_d0_elevator = 0.01 # Constant
d_e_elevator = np.radians(30) # Aileron deflection [Rad]
cm_delta_e_elevator = 3.81 # Constant
c_elevator = 1


######################################################

pb2v_ng_wing = - p_required * b_wing / (2*V)
C_lp_wing = -((c_la_wing+c_d0_wing)*C_r_wing*b_wing*(1+3*taper_wing)/(24*S_wing))

rb2v_ng_rudder = - r_required * b_rudder / (2*V*Vv_V)
C_nr_rudder = -(c_n_beta_rudder+c_d0_rudder)*C_r_rudder*b_rudder*(1+3*taper_rudder)/(24*S_rudder)*10*(6/3.5)**2  # To be replaced by simulation result!!!!!!!!!

v_c_ng_elevator = - c_elevator / (V*Vh_V)
C_mq_elevator = -(c_m_a_elevator+c_d0_elevator)*C_r_elevator*b_elevator*(1+3*taper_elevator)/(24*S_elevator)*10*(6/3.5)**2 # To be replaced by simulation result!!!!!!!

def Solver_wing (b1, pb2v_ng_wing,C_lp_wing,d_a_wing,cl_delta_a_wing):
    Cl_da = (cl_delta_a_wing*C_r_wing)/(S_wing*b_wing)*((b2_wing**2-b1**2)+ 4*(taper_wing-1)*(b2_wing**3-b1**3)/(3*b_wing))
    f = pb2v_ng_wing*C_lp_wing /d_a_wing - Cl_da

    return f

def Solver_rudder (b1, rb2v_ng_rudder ,C_nr_rudder,d_r_rudder,cn_delta_r_rudder):
    Cl_da = (cn_delta_r_rudder*C_r_rudder)/(S_rudder*b_rudder)*((b2_rudder**2-b1**2)+ 4*(taper_rudder-1)*(b2_rudder**3-b1**3)/(3*b_rudder))
    f = rb2v_ng_rudder*C_nr_rudder/d_r_rudder - Cl_da

    return f

def Solver_elevator (b1, v_c_ng_rudder, C_mq_elevator, d_e_elevator,cm_delta_e_elevator):
    Cl_da = (cm_delta_e_elevator*C_r_elevator)/(S_elevator*b_elevator)*((b2_elevator**2-b1**2)+ 4*(taper_elevator-1)*(b2_elevator**3-b1**3)/(3*b_elevator))
    f = v_c_ng_rudder *C_mq_elevator/d_e_elevator - Cl_da

    return f

def area_cal_double (C_r, taper, starting_point, length, span):
    end_point = starting_point + length
    half_span = span / 2
    C_t = C_r * taper
    upper = C_r - (C_r-C_t)/half_span*end_point
    lower = C_r - (C_r-C_t)/half_span*starting_point
    area = (upper+lower)*length/2 # One side!

    return area

def area_cal_single (C_r, taper, starting_point, length, span):
    end_point = starting_point + length
    half_span = span
    C_t = C_r * taper
    upper = C_r - (C_r-C_t)/half_span*end_point
    lower = C_r - (C_r-C_t)/half_span*starting_point
    area = (upper+lower)*length/2 # Whole side!

    return area


initial_value = 1
b1_wing = fsolve(Solver_wing, initial_value, args = (pb2v_ng_wing,C_lp_wing,d_a_wing,cl_delta_a_wing))
b1_rudder = fsolve(Solver_rudder, initial_value, args = (rb2v_ng_rudder ,C_nr_rudder,d_r_rudder,cn_delta_r_rudder))
b1_elevator = fsolve(Solver_elevator, initial_value, args = (v_c_ng_elevator ,C_mq_elevator ,d_e_elevator,cm_delta_e_elevator))

print('The position of the starting point of the aileron is:', b1_wing[0], '\\\\\\\\\\\  The length of the aileron is:',b2_wing-b1_wing[0])
print('The position of the starting point of the rudder is:', b1_rudder[0],'\\\\\\\\\\  The length of the rudder is:',b2_rudder-b1_rudder[0])
print('The position of the starting point of the elevator is:',b1_elevator[0], '\\\\\\\\\\  The length of the elevator is:',b2_elevator-b1_elevator[0])


#Calculating area needed for control
area_aileron = area_cal_double (C_r_wing, taper_wing, b1_wing[0], b2_wing-b1_wing[0], b_wing)/4
area_elevator = area_cal_double (C_r_elevator, taper_elevator, b1_elevator[0], b2_elevator-b1_elevator[0],b_elevator)/4
area_rudder = area_cal_single (C_r_rudder, taper_rudder, b1_rudder[0], b2_rudder-b1_rudder[0],b_rudder)/4

print('The area of the aileron is:',area_aileron)
print('The area of the rudder is:',area_rudder)
print('The area of the elevator is:',area_elevator)





