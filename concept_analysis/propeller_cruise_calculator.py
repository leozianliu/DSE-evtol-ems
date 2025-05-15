import numpy as np

class Propeller:
    def __init__(self, shaft_power, freestream_velocity, freestream_density, freestream_temperature, propeller_diameter, RPM):
        
        # Inputs, all are floats
        self.shaft_power = shaft_power # shaft power after motor efficiency in Watts
        self.freestream_velocity = freestream_velocity # velocity in m/s
        self.freestream_density = freestream_density # ambient density in kg/m^3
        self.freestream_temperature = freestream_temperature # ambient temperature in Kelvin
        self.propeller_diameter = propeller_diameter # Propeller diameter in m
        self.propeller_area = np.pi * propeller_diameter ** 2 / 4 # Propeller area in m^2
        self.gamma_specific_speeds_ratio = 1.4
        self.R_air_constant = 287 # Air constant in J/(kg Kelvin)
        self.maximum_blade_mach = 0.8 # Maximum Mach number to avoid shockwaves on the blade tips
        self.speed_of_sound = np.sqrt(self.gamma_specific_speeds_ratio * self.freestream_temperature * self.R_air_constant)
        self.RPM = RPM # Rotations per minute of propeller blades
        self.advance_ratio = self.freestream_velocity / self.propeller_diameter / (self.RPM / 60)

        # These will be computed
        self.downstream_velocity = None
        self.propeller_thrust = None
        self.ideal_efficiency = None
        self.real_efficiency = None
        self.maximum_RPM = None
        self.maximum_tip_speed = None
        self.induced_velocity = None
        self.power_required = None
        

    def compute_thrust(self):

        if self.freestream_velocity < 1e-2: # This edge case needs to be tested
            # Handle hover separately:
            self.induced_velocity = (self.shaft_power / (2 * self.freestream_density * self.propeller_area)) ** (1 / 3)
            self.downstream_velocity = 2 * self.induced_velocity
        else:
            a = 0.5 / self.freestream_velocity ** 3
            b = 1.5 / self.freestream_velocity ** 2
            c = -0.5 / self.freestream_velocity
            d = -1.5 - self.shaft_power / (0.5 * self.freestream_density * self.freestream_velocity ** 3 * self.propeller_area)

            cubic_coefficients = [a, b, c, d]
            cubic_roots = np.roots(cubic_coefficients)

            valid_roots = [
                v.real for v in cubic_roots
                if np.isreal(v) and v.real > 0 and self.freestream_velocity < v.real < self.freestream_velocity + 50
            ]

            if len(valid_roots) == 1:
                self.downstream_velocity = float(valid_roots[0])
                propeller_airspeed = (self.freestream_velocity + self.downstream_velocity) / 2
                self.propeller_thrust = (
                    self.freestream_density * propeller_airspeed * self.propeller_area *
                    (self.downstream_velocity - self.freestream_velocity))
            else:
                print("Multiple or no real and positive downstream velocity solutions found in Propeller.compute_thrust, returning zero thrust")
                self.propeller_thrust = 0
                return self.propeller_thrust
        
        return self.propeller_thrust
    

    def compute_power_required(self):

        # Estimate power required from thust, formulas need to be checked

        if self.freestream_velocity < 1e-2:
            self.induced_velocity = np.sqrt(self.propeller_thrust / (2 * self.freestream_density * self.propeller_area))
            self.power_required = self.propeller_thrust * self.induced_velocity
        else:
            self.induced_velocity = self.propeller_thrust / (2 * self.freestream_density * self.propeller_area * self.freestream_velocity)
            self.power_required = self.propeller_thrust * self.freestream_velocity + 0.5 * self.propeller_thrust * self.induced_velocity  # induced + profile

        return self.power_required
    

    def compute_efficiency(self):
        if self.propeller_thrust is None:
            raise RuntimeError("Must compute thrust before efficiency.")

        self.ideal_efficiency = 2 * self.freestream_velocity / (self.freestream_velocity + self.downstream_velocity)
        self.real_efficiency = self.propeller_thrust * self.freestream_velocity / self.shaft_power

        # Note: This approximation method estimates maximum real efficiency as 50% for some reason

        return self.ideal_efficiency, self.real_efficiency
    
    
    def compute_maximum_tip_speed(self):

        self.maximum_tip_speed = self.maximum_blade_mach * self.speed_of_sound
        return self.maximum_tip_speed
    

    def compute_maximum_RPM(self):
        
        self.maximum_RPM = 60 * self.maximum_tip_speed / np.pi / self.propeller_diameter
        return self.maximum_RPM