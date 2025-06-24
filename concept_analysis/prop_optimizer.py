import numpy as np
import matplotlib.pyplot as plt
from propeller_power_calculator import Propeller
class Optimizer:

    def __init__(self, mtom, max_load_factor, thrust_setting):

        g = 9.80665
        self.mtow = mtom * g
        self.mass = mtom
        self.max_load_factor = max_load_factor
        self.thrust_setting = thrust_setting

        self.big_diameter = 3.19 
        self.small_diameter = 2.286

        self.engine_rotation_speed = 5 # deg/s
        self.max_load_factor = max_load_factor # gs

    def drag(self, velocity_wing, velocity_fuselage):

        cd_fuselage = 0.09
        cd_wing = 0.035
        wing_surface = 11.27
        frontal_area = 3.436449

        return 0.5 * 1.225 * (velocity_fuselage ** 2 * cd_fuselage * frontal_area + velocity_wing ** 2 * cd_wing * wing_surface) * 1.5

    def compute_power_energy(self, design_thrust : list, design_velocity : list, design_rpm : list, design_number_blades : list, climb_speed=10, plot=False):

        g = 9.80665

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
        #print(2 * big_prop.maximum_RPM, 2 * small_prop.maximum_RPM)

        dt = 1 # s

        htransition = 300
        height_current = 0

        flight_altitude_arr = [0]
        self.power_arr_big = [0]
        self.power_arr_small = [0]
        self.power_arr_total = [0]
        self.time_arr = [0]
        energy = [0]
        velocity_arr = [0]

        big_prop.change_flight_regime(0, self.mtow / 2 * self.big_loading_factor)
        small_prop.change_flight_regime(0, self.mtow / 2 * self.small_loading_factor)
        velocity_fuselage = 0
        velocity_wing = big_prop.vi

        ## Take-Off
        #-------------------------------------------------------------------------------------------------------------------------------------------
        while height_current < htransition:
            
            flight_altitude_arr.append(height_current)

            if velocity_fuselage < climb_speed:
                thrust_big = self.thrust_setting * (self.mtow + self.drag(velocity_wing, velocity_fuselage)) / 2 * self.big_loading_factor
                thrust_small = self.thrust_setting * (self.mtow + self.drag(velocity_wing, velocity_fuselage)) / 2 * self.small_loading_factor
                acceleration = ((thrust_big + thrust_small) * 2 - self.mtow - self.drag(velocity_wing, velocity_fuselage)) / self.mass
                

                if acceleration > self.max_load_factor * g:
                    acceleration = self.max_load_factor * g
                    thrust_big = (self.mtow +  self.drag(velocity_wing, velocity_fuselage) + acceleration * self.mass) / 2 * self.big_loading_factor
                    thrust_small = (self.mtow +  self.drag(velocity_wing, velocity_fuselage) + acceleration * self.mass) / 2 * self.small_loading_factor
                    

                velocity_fuselage += acceleration * dt
                velocity_wing = velocity_fuselage 
            else:
                velocity_fuselage = climb_speed
                print(self.time_arr[-1] / 60)
                velocity_wing = velocity_fuselage
                thrust_big = (self.mtow + self.drag(velocity_wing, velocity_fuselage)) / 2 * self.big_loading_factor
                thrust_small = (self.mtow + self.drag(velocity_wing, velocity_fuselage)) / 2 * self.small_loading_factor


            velocity_arr.append(velocity_fuselage)

            big_prop.change_flight_regime(velocity_fuselage, thrust_big)
            small_prop.change_flight_regime(velocity_fuselage, thrust_small)

            self.power_arr_big.append(big_prop.shaft_power/1000)
            self.power_arr_small.append(small_prop.shaft_power/1000) 

            height_current += velocity_fuselage * dt
            self.time_arr.append(self.time_arr[-1] + dt)
            energy.append(2 * (big_prop.shaft_power + small_prop.shaft_power) * dt / 3600 / 1000)
        
        #print(big_prop.vi)
        ## Cruise
        #----------------------------------------------------------------------------------------------------------------------------------------

        v_cruise = 55.56 * 1 # m/s
        t_cruise = 18 * 60 # seconds

        t_start_cruise = self.time_arr[-1]

        big_prop.change_flight_regime(v_cruise , self.mtow / 11.3 * self.big_loading_factor / 2)
        small_prop.change_flight_regime(v_cruise , self.mtow / 11.3 * self.small_loading_factor / 2)

        preq = 0.5 * 1.225 * v_cruise**3 * 11.27 * 0.014 + (2600*9.81)**2 / (0.5 * 1.225 * v_cruise * 11.27 * np.pi * 0.9 * 13.1) + v_cruise * (0.09 * 0.5 * 1.225 * v_cruise**2 * 3.436449)
        
        t = preq / v_cruise
        
        #big_prop.change_flight_regime(v_cruise, t * 0.9 / 2)
        #small_prop.change_flight_regime(v_cruise, t * 0.1 / 2)


        velocity_fuselage = v_cruise
        velocity_wing = v_cruise + big_prop.vi
        #velocity_wing = v_cruise
        height_current = htransition

        while self.time_arr[-1] <= t_start_cruise + 50*1000 / v_cruise:

            flight_altitude_arr.append(height_current)
            self.time_arr.append(self.time_arr[-1] + dt)
            self.power_arr_big.append(big_prop.shaft_power/1000)
            self.power_arr_small.append(small_prop.shaft_power/1000)

            energy.append(2 * (big_prop.shaft_power + small_prop.shaft_power) * dt / 3600 / 1000)


        ## Landing
        #----------------------------------------------------------------------------------------------------------------------------------------

        velocity_fuselage = 0.
        velocity_wing = velocity_fuselage - big_prop.vi

        thrust_big = (self.mtow + self.drag(velocity_wing, 0)) / 2 * self.big_loading_factor
        thrust_small = (self.mtow + self.drag(velocity_wing, 0)) / 2 * self.small_loading_factor
        total_thrust = 2 * (thrust_big + thrust_small) 
        hover_thrust = total_thrust

        big_prop.change_flight_regime(0, thrust_big)
        small_prop.change_flight_regime(0, thrust_small)

        brake_distance = 0.5 * self.max_load_factor * g * (climb_speed / (self.max_load_factor * g)) ** 2
        
        while height_current > 0:
            
            flight_altitude_arr.append(height_current)

            if height_current > brake_distance:

                if velocity_fuselage > - climb_speed:

                    acceleration = (total_thrust - total_thrust * self.thrust_setting) / self.mass 
                    total_thrust /= self.thrust_setting
                    thrust_big = total_thrust / 2 * self.big_loading_factor
                    thrust_small = total_thrust / 2 * self.small_loading_factor

                    if acceleration < - self.max_load_factor * g:
                        acceleration = - self.max_load_factor * g
                        total_thrust = hover_thrust
                        thrust_big = total_thrust / 2 * self.big_loading_factor
                        thrust_small = total_thrust / 2 * self.small_loading_factor

                else:
                    #velocity_fuselage = - climb_speed
                    acceleration = 0
                    thrust_big = hover_thrust / 2 * self.big_loading_factor
                    thrust_small = hover_thrust / 2 * self.small_loading_factor

            else:
                #height_current = brake_distance
                #print(self.time_arr[-1] / 60)
                acceleration = self.max_load_factor * g
                thrust_big = hover_thrust / 2 * self.big_loading_factor
                thrust_small = hover_thrust / 2 * self.small_loading_factor
            
            velocity_fuselage += acceleration * dt
            velocity_wing = velocity_fuselage - big_prop.vi
            velocity_arr.append(velocity_fuselage)
            total_thrust = 2 * (thrust_big + thrust_small)

            big_prop.change_flight_regime(np.abs(velocity_fuselage), thrust_big)
            small_prop.change_flight_regime(np.abs(velocity_fuselage), thrust_small)

            self.power_arr_big.append(big_prop.shaft_power/1000)
            self.power_arr_small.append(small_prop.shaft_power/1000) 

            height_current += velocity_fuselage * dt
            self.time_arr.append(self.time_arr[-1] + dt)
            energy.append(2 * (big_prop.shaft_power + small_prop.shaft_power) * dt / 3600 / 1000)
            

        self.total_energy = np.sum(energy)
        self.flight_altitude_arr = flight_altitude_arr
        self.power_arr_small = np.array(self.power_arr_small)
        self.power_arr_big = np.array(self.power_arr_big)

        self.power_arr_total = 2 * (np.array(self.power_arr_big) + np.array(self.power_arr_small))
        self.energy = energy

        if plot:
            plt.plot(self.time_arr/60, flight_altitude_arr)
            plt.title('Flight envelope without transitioning phase')
            plt.grid()
            plt.tight_layout()
            plt.xlabel('Flight Time [min]')
            plt.ylabel('Altitude')
            plt.show()