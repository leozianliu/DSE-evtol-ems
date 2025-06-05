import numpy as np
import matplotlib.pyplot as plt
from propeller_power_calculator import Propeller
class Optimizer:

    def __init__(self, mtom, lift_over_drag, max_load_factor):

        g = 9.80665
        self.mtow = mtom * g
        self.mass = mtom
        self.max_load_factor = max_load_factor

        self.big_diameter = 3.33 
        self.small_diameter = 2.353

        self.lift_over_drag = lift_over_drag

        self.engine_rotation_speed = 5 # deg/s
        self.max_load_factor = max_load_factor # gs
        lift_over_drag = 10

    def drag(velocity):

        cd_fuselage = 0.1444908713
        cd_wing = 0.035
        wing_surface = 14
        frontal_area = 3.436449

        return 0.5 * 1.225 * velocity ** 2 * (cd_fuselage * frontal_area + cd_wing * wing_surface) * 1.5

    def compute_power_energy(self, design_thrust : list, design_velocity : list, design_rpm : list, design_number_blades : list, climb_speed=5):

        prop_area_big = np.pi * (self.big_diameter/2) ** 2
        prop_area_small = np.pi * (self.small_diameter/2) ** 2

        self.big_loading_factor = prop_area_big / (prop_area_big + prop_area_small)
        self.small_loading_factor = prop_area_small / (prop_area_big + prop_area_small)

        big_prop = Propeller(design_thrust[0], design_velocity[0], self.big_diameter, design_rpm[0], design_number_blades[0])
        small_prop = Propeller(design_thrust[1], design_velocity[1], self.small_diameter, design_rpm[1], design_number_blades[1])

        big_prop.compute_maximum_RPM(use_maximum_RPM=True)
        small_prop.compute_maximum_RPM(use_maximum_RPM=True)

        big_prop.design_propeller_geometry()
        small_prop.design_propeller_geometry()

        velocity = 0
        dt = 0.2 

        htransition = 300
        height_current = 0

        flight_altitude_arr = np.array([0])
        power_arr_big = np.array([0])
        power_arr_small = np.array([0])
        time_arr = np.array([0])
        energy = np.array([0])
        velocity_arr = []

        ## Take-Off
        #-------------------------------------------------------------------------------------------------------------------------------------------
        while height_current < htransition:
            
            flight_altitude_arr = np.append(flight_altitude_arr, height_current)
            if velocity < 5:
                thrust_big = self.max_load_factor * (self.mtow / 2 + self.drag(velocity) / 2) * self.big_loading_factor
                thrust_small = self.max_load_factor * (self.mtow / 2 + self.drag(velocity) / 2) * self.small_loading_factor
                acceleration = ((thrust_big + thrust_small) * 2 - (self.mtow + self.drag(velocity))) / self.mass
                velocity += acceleration * dt
            else:
                thrust_big = (self.mtow / 2 + self.drag(climb_speed) / 2) * self.big_loading_factor
                thrust_small = (self.mtow / 2 + self.drag(climb_speed) / 2) * self.small_loading_factor
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
        #----------------------------------------------------------------------------------------------------------------------------------------

        v_cruise = 55.56 # m/s
        t_cruise = 18 * 60 # seconds

        t_start_cruise = time_arr[-1]

        thrust_big = (self.drag(v_cruise + big_prop.vi) / 2) * self.big_loading_factor
        thrust_small = (self.drag(v_cruise + big_prop.vi) / 2) * self.small_loading_factor

        big_prop.change_flight_regime(v_cruise, thrust_big)
        small_prop.change_flight_regime(v_cruise, thrust_small)

        while time_arr[-1] <= t_start_cruise + t_cruise:
            
            time_arr = np.append(time_arr, time_arr[-1] + dt)

            power_arr_big = np.append(power_arr_big, big_prop.shaft_power/1000)
            power_arr_small = np.append(power_arr_small, small_prop.shaft_power/1000)

            energy = np.append(energy, 2 * (big_prop.shaft_power + small_prop.shaft_power)/1000 * dt / 3600)

