import numpy as np
import matplotlib.pyplot as plt
from propeller_power_calculator import Propeller

cd_fuselage = 0.1444908713
cd_wing = 0.035
cd_total = (cd_fuselage + cd_wing) * 1.5
wing_surface = 14
frontal_area = 3.436449
g = 9.80665
mtow = 2200 * g
big_diameter = 3.33 
small_diameter = 2.353
engine_rotation_speed = 5 # deg/s
max_load_factor = 1.1 # gs
lift_over_drag = 10


def drag(velocity):
    return 0.5 * 1.225 * velocity ** 2 * (cd_fuselage * frontal_area + cd_wing * wing_surface) * 1.5


prop_area_big = np.pi * (big_diameter/2) ** 2
prop_area_small = np.pi * (small_diameter/2) ** 2

hover_thrust_big = mtow / 2  * prop_area_big / (prop_area_big + prop_area_small)
hover_thrust_small = mtow / 2  * prop_area_small / (prop_area_big + prop_area_small)

big_prop = Propeller(hover_thrust_big / lift_over_drag, 3, big_diameter, 4500, 5)
small_prop = Propeller(hover_thrust_small / lift_over_drag, 30, small_diameter, 4500, 5)

big_prop.compute_maximum_RPM(use_maximum_RPM=True)
small_prop.compute_maximum_RPM(use_maximum_RPM=True)

big_prop.design_propeller_geometry()
small_prop.design_propeller_geometry()

total_energy_prop1 = 0
total_energy_prop2 = 0

take_off_time = 120 # seconds per round trip
landing_time = 120 # seconds per round trip
cruise_time = 2160 # seconds per round trip
velocity = 0
dt = 2 
mass = mtow / g

hcruise = 300
htransition = 300
height_current = 0

flight_altitude_arr = np.array([0])
power_arr_big = np.array([0])
power_arr_small = np.array([0])
time_arr = np.array([0])
energy = np.array([0])
velocity_arr = []

## Take-Off
while height_current < htransition:
    
    flight_altitude_arr = np.append(flight_altitude_arr, height_current)
    if velocity < 5:
        thrust_big = max_load_factor * (mtow / 2 + drag(velocity) / 2) * prop_area_big / (prop_area_big + prop_area_small)
        thrust_small = max_load_factor * (mtow / 2 + drag(velocity) / 2) * prop_area_small / (prop_area_big + prop_area_small)
        acceleration = ((thrust_big + thrust_small) * 2 - (mtow + drag(velocity))) / mass
        velocity += acceleration * dt
    else:
        thrust_big = (mtow / 2 + drag(5) / 2) * prop_area_big / (prop_area_big + prop_area_small)
        thrust_small = (mtow / 2 + drag(5) / 2) * prop_area_small / (prop_area_big + prop_area_small)
        velocity = 5

    velocity_arr.append(velocity)

    big_prop.change_flight_regime(velocity, thrust_big)
    small_prop.change_flight_regime(velocity, thrust_small)

    power_arr_big = np.append(power_arr_big, big_prop.shaft_power/1000)
    power_arr_small = np.append(power_arr_small, small_prop.shaft_power/1000)

    height_current += velocity * dt
    time_arr = np.append(time_arr, time_arr[-1] + dt)
    energy = np.append(energy, 2 * (big_prop.shaft_power + small_prop.shaft_power)/1000 * dt / 3600)
    

## Cruise

v_cruise = 55.56 # m/s
t_cruise = 18 * 60 # seconds

t_start_cruise = time_arr[-1]

thrust_big = (drag(v_cruise + big_prop.vi) / 2) * prop_area_big / (prop_area_big + prop_area_small)
thrust_small = (drag(v_cruise + big_prop.vi) / 2) * prop_area_small / (prop_area_big + prop_area_small)

big_prop.change_flight_regime(v_cruise, thrust_big)
small_prop.change_flight_regime(v_cruise, thrust_small)

while time_arr[-1] <= t_start_cruise + t_cruise:
    
    time_arr = np.append(time_arr, time_arr[-1] + dt)

    power_arr_big = np.append(power_arr_big, big_prop.shaft_power/1000)
    power_arr_small = np.append(power_arr_small, small_prop.shaft_power/1000)

    energy = np.append(energy, 2 * (big_prop.shaft_power + small_prop.shaft_power)/1000 * dt / 3600)


print(small_prop.thrust, big_prop.thrust)
print(sum(energy) * 2)
plt.plot(time_arr/60, 2 * power_arr_big, color='red')
plt.plot(time_arr/60, 2 * power_arr_small, color='green')
plt.plot(time_arr/60, 2 * power_arr_big + 2 * power_arr_small, color='black')
plt.ylabel('Power Required [kW]')
plt.xlabel('Time [min]')
plt.grid()
plt.show()