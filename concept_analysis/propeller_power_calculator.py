import numpy as np

class Propeller:
    def __init__(self, shaft_power, freestream_velocity, freestream_density, freestream_temperature, propeller_diameter, RPM, number_blades):
        
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
        self.number_blades = number_blades # Number of blades of the propeller

        # These will be computed
        self.downstream_velocity = None
        self.propeller_thrust = None 
        self.propeller_efficiency = None
        self.maximum_RPM = None
        self.maximum_tip_speed = None
        self.induced_velocity = None
        self.power_required = None
        
    def compute_power_required(self, true_if_compute_power_false_if_compute_thrust=True):

        # === Step 1: Initialize ===

        N = 100 # Number of panels
        epsilon = 0 # initial guess
        Cl = 0.7 # Assumed performance metrics of the airfoil
        Cd = 0.014
        s = Cd/Cl
    
        r_hub = 0.15 * self.propeller_diameter/2

        omega = 0.10472 * self.RPM

        converged = False
        max_iterations = 1e6
        epsilon_tol = 1e-6

        r = np.linspace(r_hub, self.propeller_diameter/2, N)
        dr = (self.propeller_diameter/2 - r_hub) / N
        lamda = self.freestream_velocity / (omega * (self.propeller_diameter/2))  # speed ratio


        xi = r / (self.propeller_diameter/2)
        x = xi / lamda

        if true_if_compute_power_false_if_compute_thrust:
            power_coefficient = 2 * self.power_required / (self.freestream_density * self.freestream_velocity ** 3 * np.pi * (self.propeller_diameter/2) ** 2)
        else: 
            thrust_coefficient = 2 * self.propeller_thrust / (self.freestream_density * self.freestream_velocity ** 2 * np.pi * (self.propeller_diameter/2) ** 2)

        # Estimate power required from thust
        while not converged and max_iterations:

            tan_phi_tip = lamda * (1 + epsilon / 2)
            phi_tip = np.arctan(tan_phi_tip)

            phi = np.arctan(tan_phi_tip / xi)

            F = (2/np.pi) * np.rad2deg(np.arccos(np.exp(-(self.number_blades/2) * (1 - xi) / (np.sin(phi_tip)))))
            G = F * np.cos(phi) * np.sin(phi) * x

            Wc = (4 * np.pi * lamda * G * self.freestream_velocity * (self.propeller_diameter/2) * epsilon) / (Cl * self.number_blades)
            
            a = (epsilon/2) * np.cos(phi)**2 * (1 - s * np.tan(phi))
            W = self.freestream_velocity * (1 + a) / np.sin(phi)

            chord = Wc / W
            alpha = np.radians(5)  # Assume 5 deg fixed angle of attack of blades
            beta = alpha + phi

            # === Eq. 11 integrands ===
            I1_integrand = 4 * G * xi * (1 - s * np.tan(phi)) 
            I2_integrand = lamda * (I1/ xi / 2) * (1 + s / np.tan(phi)) * np.sin(phi) * np.cos(phi)
            J1_integrand = 4 * xi * G * (1 + s / np.tan(phi))
            J2_integrand = (J1 / 2) * (1 - s * np.tan(phi)) * np.cos(phi)**2
            

            I1 = np.trapz(I1_integrand, xi)
            I2 = np.trapz(I2_integrand, xi)
            J1 = np.trapz(J1_integrand, xi)
            J2 = np.trapz(J2_integrand, xi)
            
            # Integrate (Eq. 12 & 13)
            if true_if_compute_power_false_if_compute_thrust:
                epsilon_new = (I1 / I2 / 2) - np.sqrt((I1 / I2 / 2)**2 - thrust_coefficient/I2)
                power_coefficient = J1 * epsilon + J2 * epsilon**2
            else:
                epsilon_new = - (J1 / 2 / J2) + np.sqrt((J1 / 2 / J2) ** 2 + power_coefficient/J2)
                thrust_coefficient = I1 * epsilon_new - I2 * epsilon_new ** 2

            # Update epsilon
            max_iter -= 1

            if abs(epsilon_new - epsilon) < epsilon_tol:
                converged = True
            
            epsilon = epsilon_new

        self.propeller_thrust = 0.5 * thrust_coefficient * self.freestream_velocity * self.freestream_velocity ** 2 * np.pi * (self.propeller_diameter/2) ** 2
        self.power_required = power_coefficient * 0.5 * self.freestream_density * self.freestream_velocity ** 3 * np.pi * (self.propeller_diameter/2) ** 2

        CT = self.propeller_thrust / (self.freestream_density * (self.RPM / 60)**2 * (2 * (self.propeller_diameter/2)) ** 4)
        CP = self.power_required / (self.freestream_density * (self.RPM / 60)**3 * (2 * (self.propeller_diameter/2)) ** 5)

        self.propeller_efficiency = self.advance_ratio * CT / CP  

        return self.power_required
        
    
    def compute_maximum_tip_speed(self):

        self.maximum_tip_speed = self.maximum_blade_mach * self.speed_of_sound
        return self.maximum_tip_speed
    

    def compute_maximum_RPM(self):
        
        self.maximum_RPM = 60 * self.maximum_tip_speed / np.pi / self.propeller_diameter
        return self.maximum_RPM
    
