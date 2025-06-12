import numpy as np 
import matplotlib.pyplot as plt 

#weights [kg]
w_person = 90
w_stretcher = 35
w_skin = 210 
w_frame = 150
w_side_door = 30
w_back_door = 20
w_windshield = 30
w_floor = 60

#center of gravity from the nose in [m]
x_cabin = 1.2
x_skin = 1.5 + x_cabin
x_frame = x_skin
x_side_door = 3.3
x_back_door = 4.8
x_windshield = 0.8
x_floor = 2.2

top = w_skin * x_skin + w_frame * x_frame + w_side_door * x_side_door + w_back_door * x_back_door + w_windshield * x_windshield + w_floor * x_floor
w_fuselage = w_skin + w_frame + w_side_door + w_back_door + w_windshield + w_floor
x_cg_fuselage = top/w_fuselage
print(f'Center of gravity of the fuselage structure: {x_cg_fuselage:.2f} m from the nose \n')

w_HV_battery = 630 + 50 #including ac 
w_LV_battery = 20
w_battery_MS = 20

x_LV_battery = 0.8
x_HV_battery = 1.87 #orignal 1.95
top_with_battery = top + x_LV_battery*w_LV_battery + w_HV_battery * x_HV_battery 
w_with_battery = w_fuselage + w_HV_battery + w_LV_battery
x_cg_fuselage_battery = top_with_battery / w_with_battery
print(f'Including the battery weights: {x_cg_fuselage_battery:.2f} m from nose \n')

w_pilot_seat = 12
w_doctor_seat = 12

x_pilot_seat = 1.45
x_doctor_seat = 2.25 
top_with_seats = top_with_battery + w_pilot_seat * x_pilot_seat + w_doctor_seat * x_pilot_seat * 2
w_with_seats = w_with_battery + w_pilot_seat + 2 * w_doctor_seat
x_cg_fuselage_battery_seats = top_with_seats / w_with_seats
print(f'Including the seats: {x_cg_fuselage_battery_seats:.2f} m from nose \n')

w_cockpit_instruments = 20 
w_equipment = 61
w_stretcher = 35

x_cockpit_instruments = 0.45
x_eq = 1.6 
x_stretcher = 3.2

top_with_more = top_with_seats + w_cockpit_instruments * x_cockpit_instruments + w_equipment * x_eq + w_stretcher * x_stretcher
w_with_more = w_with_seats + w_cockpit_instruments + w_equipment + w_stretcher
x_cg_fuselage_more = top_with_more / w_with_more
print(f'No passengers fuselage center of gravity : {x_cg_fuselage_more:.2f} m from nose') 


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### CASES for cg range 
w_pilot = 90
x_pilot = 1.45
w_doctor = 90


#including the pilot
top_with_pilot = top_with_more + w_pilot * x_pilot
w_with_pilot = w_with_more + w_pilot
x_cg_pilot = top_with_pilot/w_with_pilot

def cg(x_d1, x_d2, x_patient, w_pat):
    top = top_with_pilot + x_d1 * w_doctor + x_d2 * w_doctor + x_patient * w_pat
    weight = w_with_pilot + w_doctor * 2 + w_pat
    cg = top/weight
    return cg

#case 1 - everyone is on the plane at thier spots 
x_cg_1 = cg(
    x_d1 = 2.25,
    x_d2 = 2.25,
    x_patient = 3.2,
    w_pat = 90
)

#case 2 - no patient 
x_cg_2 = cg(
    x_d1 = 2.25,
    x_d2 = 2.25,
    x_patient = 3.2,
    w_pat = 0
)

#case 3 - medics are at the end of the stretcher 
x_cg_3 = cg(
    x_d1 = 3.7,
    x_d2 = 3.7,
    x_patient = 3.2,
    w_pat = 90
)

#case 4 - one medic is next to the pilot
x_cg_4 = cg(
    x_d1 = 1.45,
    x_d2 = 2.25,
    x_patient = 3.2,
    w_pat = 90
)

root_chord = 1.4 
x_hinge_loc = 0.309 * root_chord #from the LE of the root chord
x_leroot = 1.7
x_hinge = x_hinge_loc + x_leroot #this is where we want the fuselage center of gravity to be 
print(f'Location fo the hinge: {x_hinge:.2f} m from the nose \n')

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#PLOTTING
y_graph = 0
plt.figure(figsize=(6,2))
#plt.plot([1700, 3100], [y_graph, y_graph], color='lightgrey', label='root chord')
#plt.scatter(x_leroot, y_graph, marker='|')
#plt.axvline(x = x_hinge*1000)
plt.scatter(x_cg_fuselage_more, y_graph, label = 'Structural CG')
plt.scatter(x_cg_1, y_graph, label='CG case 1', marker='x')
plt.scatter(x_cg_2, y_graph, label='CG case 2', marker='x')
plt.scatter(x_cg_3, y_graph, label='CG case 3', marker='x')
plt.scatter(x_cg_4, y_graph, label='CG case 4', marker='x')
plt.yticks(color='white')
plt.xlabel('Location along the fuselage')
plt.legend(loc = 'lower right', fontsize = 10)
plt.show()


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#LIST

cg = [x_cg_1, x_cg_2, x_cg_3, x_cg_4]
print(f'CENTER of gravity variations: -- \n'
      f'Case 1: {x_cg_1:.2f} m, \n'
      f'Case 2: {x_cg_2:.2f} m, \n'
      f'Case 3: {x_cg_3:.2f} m, \n'
      f'Case 4: {x_cg_4:.2f} m, \n')

cg_min = min(cg)
cg_max = max(cg)
cg_range = cg_max - cg_min
print(f'Range of the cg: {cg_range:.2f} m')