import numpy as np
import matplotlib.pyplot as plt

class Propeller:
    def __init__(self, shaft_power, propeller_thrust, freestream_velocity, freestream_density, freestream_temperature, propeller_diameter, RPM, number_blades, ducted_fan=False):
        
        # Inputs, all are floats
        self.shaft_power = shaft_power # shaft power after motor efficiency in Watts
        self.propeller_thrust = propeller_thrust # motor 
        self.freestream_velocity = freestream_velocity # velocity in m/s
        self.freestream_density = freestream_density # ambient density in kg/m^3
        self.freestream_temperature = freestream_temperature # ambient temperature in Kelvin
        self.propeller_diameter = propeller_diameter # Propeller diameter in m
        self.propeller_area = np.pi * propeller_diameter ** 2 / 4 # Propeller area in m^2
        self.gamma_specific_speeds_ratio = 1.4
        self.R_air_constant = 287 # Air constant in J/(kg Kelvin)
        self.maximum_blade_mach = 0.7 # Maximum Mach number to avoid shockwaves on the blade tips
        self.speed_of_sound = np.sqrt(self.gamma_specific_speeds_ratio * self.freestream_temperature * self.R_air_constant)
        self.RPM = RPM # Rotations per minute of propeller blades
        self.number_blades = number_blades # Number of blades of the propeller
        self.hub_ratio = 0.15
        self.ducted_fan = ducted_fan
        self.assumed_ducted_efficiency_increase = 1

        if self.ducted_fan:
            self.assumed_ducted_efficiency_increase += 0.15

        # These will be computed
        self.downstream_velocity = None

        self.propeller_efficiency = None
        self.maximum_RPM = None
        self.maximum_tip_speed = None
        self.advance_ratio = None
        self.Reynolds_number = None
        
    def compute_power_required(self, true_if_compute_power_false_if_compute_thrust=True):

        # === Step 1: Initialize ===

        N = 100 # Number of panels
        epsilon = 0 # initial guess
        Cl = 0.7 # Assumed performance metrics of the airfoil
        Cd = 0.07
        s = Cd/Cl
    
        r_hub = self.hub_ratio * self.propeller_diameter/2

        omega = 0.10472 * self.RPM

        converged = False
        max_iterations = 10000
        epsilon_tol = 1e-6

        r = np.linspace(r_hub, self.propeller_diameter/2, N)
        
        lamda = self.freestream_velocity / (omega * (self.propeller_diameter/2))  # speed ratio


        xi = r / (self.propeller_diameter/2)
        x = xi / lamda
        self.radial_station_positions = r

        if true_if_compute_power_false_if_compute_thrust:

            thrust_coefficient = 2 * self.propeller_thrust / (self.freestream_density * self.freestream_velocity ** 2 * np.pi * (self.propeller_diameter/2) ** 2)
        else: 
            
            power_coefficient = 2 * self.shaft_power / (self.freestream_density * self.freestream_velocity ** 3 * np.pi * (self.propeller_diameter/2) ** 2)
        
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
            I2_integrand = lamda * (I1_integrand/ xi / 2) * (1 + s / np.tan(phi)) * np.sin(phi) * np.cos(phi)
            J1_integrand = 4 * xi * G * (1 + s / np.tan(phi))
            J2_integrand = (J1_integrand / 2) * (1 - s * np.tan(phi)) * np.cos(phi)**2
            

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
            max_iterations -= 1

            if abs(epsilon_new - epsilon) < epsilon_tol:
                converged = True
            
            epsilon = epsilon_new

        if true_if_compute_power_false_if_compute_thrust:

            power_coefficient = J1 * epsilon + J2 * epsilon ** 2
            self.shaft_power = (1/self.assumed_ducted_efficiency_increase) * 0.5 * power_coefficient * self.freestream_density * self.freestream_velocity ** 3 * np.pi * (self.propeller_diameter/2) ** 2
        else:
            thrust_coefficient = I1 * epsilon - I2 * epsilon ** 2
            self.propeller_thrust = self.assumed_ducted_efficiency_increase * 0.5 * thrust_coefficient * self.freestream_density * self.freestream_velocity ** 2 * np.pi * (self.propeller_diameter/2) ** 2


        self.chord_distribution = chord
        self.twist_distribution = np.rad2deg(beta)

        CT = self.propeller_thrust / (self.freestream_density * (self.RPM / 60)**2 * self.propeller_diameter ** 4)
        CP = self.shaft_power / (self.freestream_density * (self.RPM / 60)**3 * self.propeller_diameter ** 5)

        self.advance_ratio = lamda * np.pi
        self.propeller_efficiency = self.advance_ratio * CT / CP  
        self.power_loading = self.propeller_thrust / self.shaft_power
        self.disk_loading = self.propeller_thrust / self.propeller_area

        return              

    def compute_maximum_RPM(self, use_maximum_RPM=True):
        
        self.maximum_tip_speed = self.maximum_blade_mach * np.sqrt(self.speed_of_sound ** 2 - self.freestream_velocity ** 2)
        self.maximum_RPM = 60 * self.maximum_tip_speed / np.pi / self.propeller_diameter

        if use_maximum_RPM:
            self.RPM = self.maximum_RPM

        return 
    
    def adjust_chord_for_desired_Reynolds(self, desired_Reynolds_number=0):
        if desired_Reynolds_number:

            self.Reynolds_number = (self.freestream_density * self.radial_station_positions * self.RPM * 0.10472) * self.chord_distribution / 1.81e-5

            low_reynolds_element = round(((0.4 - self.hub_ratio) / (1 - (self.hub_ratio / self.number_blades))) * self.number_blades)
            high_reynolds_element = round(((0.95 - self.hub_ratio) / (1 - (self.hub_ratio / self.number_blades))) * self.number_blades)
            low_reynolds_number = self.Reynolds_number[high_reynolds_element]
            high_reynolds_number = self.Reynolds_number[low_reynolds_element]
            min_reynolds_number = min(low_reynolds_number,high_reynolds_number)


            self.chord_distribution = (desired_Reynolds_number/min_reynolds_number) * self.chord_distribution
            self.Reynolds_number = (self.freestream_density * (self.radial_station_positions * self.RPM * 0.10472) * self.chord_distribution) / 1.81e-5

        else:

            chord_factor = 0.2 / np.max(self.chord_distribution)
            self.chord_distribution = self.chord_distribution * chord_factor
            self.Reynolds_number = (self.freestream_density * self.radial_station_positions * self.RPM * 0.10472) * self.chord_distribution / 1.81e-5

        return
    

    def summary_parameters(self, plot_chord_twist_distributions=False):

        print("\nSummary of Propeller Parameters:")
        print("----------------------")
        
        print(f"Freestream velocity: {self.freestream_velocity} [m/s]")
        print(f"Propeller diameter: {self.propeller_diameter} [m]")
        print(f"Propeller area: {np.round(self.propeller_area, 2)} [m^2]")
        print(f"Rotational speed: {np.round(self.RPM, 2)} [RPM], or {np.round(self.RPM * 0.10472, 2)} [rad/s]")
        print(f"Propeller thrust: {np.round(self.propeller_thrust, 2)} [N]")
        print(f"Propeller power input: {np.round(self.shaft_power/1000, 2)} [kW]")
        print(f"Number of blades: {self.number_blades}")
        print(f"Propeller advance ratio: {np.round(self.advance_ratio, 5)}")
        print(f"Propeller efficiency: {np.round(self.propeller_efficiency*100, 2)}%")
        
        print()  # Extra space for readability

        if plot_chord_twist_distributions:
                fig, axs = plt.subplots(1, 2, figsize=(12, 5))

                # Plot 1
                axs[0].plot(self.radial_station_positions / (self.propeller_diameter/2), self.chord_distribution, color='red', marker='x')
                axs[0].set_title('Chord Distribution')
                axs[0].set_xlabel('Position along blade radius')
                axs[0].set_ylabel('Chord lenght [m]')
                axs[0].grid(True)



                # Plot 2
                axs[1].plot(self.radial_station_positions / (self.propeller_diameter/2), self.twist_distribution, color='red', marker='x')
                axs[1].set_title('Twist Distribution')
                axs[1].set_xlabel('Position along blade radius')
                axs[1].set_ylabel('Twist angle [deg]')
                axs[1].grid(True)



                plt.show()

        return


