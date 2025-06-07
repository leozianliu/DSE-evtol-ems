import numpy as np
import matplotlib.pyplot as plt
from propeller_power_calculator import Propeller
class Optimizer:

    def __init__(self, mtom, max_load_factor):

        g = 9.80665
        self.mtow = mtom * g
        self.mass = mtom
        self.max_load_factor = max_load_factor

        self.big_diameter = 3.19 
        self.small_diameter = 2.286

        self.engine_rotation_speed = 5 # deg/s
        self.max_load_factor = max_load_factor # gs

    def drag(self, velocity_drag):

        cd_fuselage = 0.1444908713
        cd_wing = 0.035
        wing_surface = 14
        frontal_area = 3.436449

        return 0.5 * 1.225 * velocity_drag ** 2 * (cd_fuselage * frontal_area + cd_wing * wing_surface) * 1.5

    def compute_power_energy(self, design_thrust : list, design_velocity : list, design_rpm : list, design_number_blades : list, climb_speed=5):

        prop_area_big = np.pi * (self.big_diameter/2) ** 2
        prop_area_small = np.pi * (self.small_diameter/2) ** 2

        self.big_loading_factor = prop_area_big / (prop_area_big + prop_area_small)
        self.small_loading_factor = prop_area_small / (prop_area_big + prop_area_small)

        big_prop = Propeller(design_thrust[0], design_velocity[0], self.big_diameter, design_rpm[0], design_number_blades[0])
        small_prop = Propeller(design_thrust[1], design_velocity[1], self.small_diameter, design_rpm[1], design_number_blades[1])

        big_prop.design_propeller_geometry()
        small_prop.design_propeller_geometry()

        big_prop.compute_maximum_RPM(use_maximum_RPM=True)
        small_prop.compute_maximum_RPM(use_maximum_RPM=True)

        velocity_drag = 0
        dt = 1 # s

        htransition = 300
        height_current = 0

        flight_altitude_arr = [0]
        self.power_arr_big = [0]
        self.power_arr_small = [0]
        self.time_arr = [0]
        energy = [0]
        velocity_arr = [0]

        ## Take-Off
        #-------------------------------------------------------------------------------------------------------------------------------------------
        while height_current < htransition:
            
            flight_altitude_arr = np.append(flight_altitude_arr, height_current)
            if velocity_drag < 5:
                thrust_big = self.max_load_factor * (self.mtow / 2 + self.drag(velocity_drag) / 2) * self.big_loading_factor
                thrust_small = self.max_load_factor * (self.mtow / 2 + self.drag(velocity_drag) / 2) * self.small_loading_factor
                acceleration = ((thrust_big + thrust_small) * 2 - (self.mtow + self.drag(velocity_drag))) / self.mass
                velocity_drag += acceleration * dt
            else:
                thrust_big = (self.mtow / 2 + self.drag(climb_speed) / 2) * self.big_loading_factor
                thrust_small = (self.mtow / 2 + self.drag(climb_speed) / 2) * self.small_loading_factor
                velocity_drag = 5

            velocity_arr.append(velocity_drag)

            big_prop.change_flight_regime(velocity_drag, thrust_big)
            small_prop.change_flight_regime(velocity_drag, thrust_small)

            self.power_arr_big.append(big_prop.shaft_power/1000)
            self.power_arr_small.append(small_prop.shaft_power/1000) 

            height_current += velocity_drag * dt
            self.time_arr.append(self.time_arr[-1] + dt)
            energy.append(1.7 *2 * (big_prop.shaft_power + small_prop.shaft_power) * dt / 3600 / 1000)

            

        ## Cruise
        #----------------------------------------------------------------------------------------------------------------------------------------

        v_cruise = 55.56 # m/s
        t_cruise = 18 * 60 # seconds

        t_start_cruise = self.time_arr[-1]

        thrust_big = (self.drag(v_cruise + big_prop.vi) / 2) * self.big_loading_factor
        thrust_small = (self.drag(v_cruise + big_prop.vi) / 2) * self.small_loading_factor

        big_prop.change_flight_regime(v_cruise, thrust_big)
        small_prop.change_flight_regime(v_cruise, thrust_small)

        while self.time_arr[-1] <= t_start_cruise + t_cruise:
            
            self.time_arr.append(self.time_arr[-1] + dt)
            self.power_arr_big.append(big_prop.shaft_power/1000)
            self.power_arr_small.append(small_prop.shaft_power/1000)

            energy.append(2 * (big_prop.shaft_power + small_prop.shaft_power) * dt / 3600 / 1000)


        ## Transition
        #----------------------------------------------------------------------------------------------------------------------------------------



        self.total_energy = np.sum(energy)