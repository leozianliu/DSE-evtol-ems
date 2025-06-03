import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar

_polar_df = pd.read_csv("concept_analysis/polar.txt", sep=r'\s+', skiprows=2,
                         names=["AOA","cl","cd","cm"])
cl_interp = interp1d(_polar_df.AOA, _polar_df.cl, kind='cubic', fill_value="extrapolate")
cd_interp = interp1d(_polar_df.AOA, _polar_df.cd, kind='cubic', fill_value="extrapolate")

def get_aero_coefficients(alpha):
    alpha = np.array(alpha)
    CL = np.array([])
    CD = np.array([])
    for angle in alpha:
        CL = np.append(CL, cl_interp(angle))
        CD = np.append(CD, cd_interp(angle))
    return CL , CD


class Propeller:
    def __init__(self, thrust, velocity, diameter, RPM, number_blades, density=1.225, temperature=293.15, ducted_fan=False):
        
        ### Inputs, all are floats

        ### Geometry and performance:

        self.diameter = diameter # Propeller diameter in m
        self.hub_ratio = 0.15
        self.rtip = diameter / 2
        self.rhub = self.hub_ratio * self.rtip
        self.propeller_area = np.pi * diameter ** 2 / 4 # Propeller area in m^2

        self.r = np.linspace(self.rhub, self.rtip, 1000) # Split blades into sections
        self.dr = self.r[1] - self.r[0]

        self.number_blades = number_blades # Number of blades of the propeller
        self.thrust = thrust # Thrust required
        self.rpm = RPM # Rotations per minute of propeller blades
        self.omega = 2 * np.pi * self.rpm / 60 # Same but converted to radians per second

        ### Ambient:

        self.velocity = velocity # velocity in m/s
        self.density = density # ambient density in kg/m^3
        self.temperature = temperature # ambient temperature in Kelvin
        self.advance_ratio = self.velocity / (self.rpm / 60) / self.diameter

        ### Constants:
        self.gamma = 1.4 # ratio of specific heats of air (Cp/Cv)
        self.R_air_constant = 287 # Air constant in J/(kg Kelvin)
        self.maximum_blade_mach = 0.7 # Maximum Mach number to avoid shockwaves on the blade tips
        self.speed_of_sound = np.sqrt(self.gamma * self.temperature * self.R_air_constant) # please tell me you know what the speed of sound is
        
        ### Ducted fan assumption:
        self.ducted_fan = ducted_fan
        self.assumed_ducted_efficiency_increase = 1

        if self.ducted_fan:
            self.assumed_ducted_efficiency_increase += 0.15 # Source: miki's hopes and dreams

        ### These will be computed, they are initialized as None for easier debugging (yes it did help ok???)

        self.efficiency = None
        self.maximum_RPM = None
        self.maximum_tip_speed = None
        self.Reynolds_number = None
        self.power_coefficient = None
        self.thrust_coefficient = None
        

    def compute_maximum_RPM(self, use_maximum_RPM=False): # Use this to prevent blades from going supersonic at the tips
        
        self.maximum_tip_speed = self.maximum_blade_mach * self.speed_of_sound 
        self.maximum_RPM = 60 * self.maximum_tip_speed / np.pi / self.diameter

        if use_maximum_RPM:
            self.rpm = self.maximum_RPM
            self.advance_ratio = self.velocity / (self.rpm / 60) / self.diameter 

        return 
            
    def blade_element__momentum_method(self):

        ### this is very innacurate and probably doesn't work, it assumes that induced velocity is perpendicular to 
        ### total velocity, which only holds somewhere under 200 N/mm^2 disk loading, which is not really doable for eVTOLs

        a = np.ones_like(self.r) * 0.3
        aline = np.ones_like(self.r) * 0.01
        Niterations = 1000
        Erroriterations = 1e-7
        sigma = self.number_blades * (self.chord + 0.0000001) / (2 * np.pi * self.r)
        anew = np.zeros_like(a)

        for _ in range(Niterations):

            Uaxial = self.velocity * (1 - a)
            Utan = (1 + aline) * self.omega * self.r
            inflowangle = np.arctan(Uaxial / Utan)
            vrel = np.sqrt(Uaxial ** 2 + Utan ** 2)

            alpha = self.twist - np.degrees(inflowangle)

            cl, cd = get_aero_coefficients(alpha) 
            cnormal = cl * np.cos(inflowangle) - cd * np.sin(inflowangle)
            ctang = cl * np.sin(inflowangle) + cd * np.cos(inflowangle)

            f = (self.number_blades / 2) * (self.diameter/2 - self.r) / (self.r * np.sin(inflowangle) + 1e-6)  # add epsilon to avoid div by zero
            F = (2 / np.pi) * np.arccos(np.exp(-f))

            anew = 1 / (1 + (4 * F * np.sin(inflowangle)**2) / (sigma * cnormal + 1e-6))
            alinenew = 1 / (1 + (4 * F * np.sin(inflowangle) * np.cos(inflowangle)) / (sigma * ctang + 1e-6))

            if np.mean(np.abs(a - anew)) < Erroriterations \
            and np.mean(np.abs(aline - alinenew)) < Erroriterations:
                break

            
            a = anew
            aline = alinenew
            
        self.dL = 0.5 * self.density * vrel ** 2 * self.chord * cl
        self.dD = 0.5 * self.density * vrel ** 2 * self.chord * cd

        self.dT = self.number_blades * (self.dL * np.cos(inflowangle) - self.dD * np.sin(inflowangle))
        self.dQ = self.r * self.number_blades * (self.dL * np.sin(inflowangle) + self.dD * np.cos(inflowangle))
        self.thrust = np.trapz(self.dT, self.r)
        self.shaft_power = np.trapz(self.dQ, self.r) * self.omega

        self.thrust_coefficient =  self.thrust / (self.density * (self.rpm / 60)**2 * self.diameter ** 4)
        self.power_coefficient = self.shaft_power / (self.density * (self.rpm / 60)**3 * self.diameter ** 5)

        self.advance_ratio = self.velocity / (self.rpm / 60) / self.diameter
        self.efficiency = self.assumed_ducted_efficiency_increase * self.advance_ratio * self.thrust_coefficient / self.power_coefficient

        self.power_loading = self.thrust / self.shaft_power
        self.disk_loading = self.thrust / self.propeller_area
        
        self.figure_of_merit = self.thrust ** 1.5 / (self.shaft_power * np.sqrt(2 * self.density * self.propeller_area))

        return
    
    def design_propeller_geometry(self):

        ### assumes optimal circulation along the blade

        def thrust_func(vline): ## Ignore the verbose, the root_scalar functions works in interesting ways
            x1 = 2 * np.pi * (self.rpm/60) * self.r / self.velocity
            x2 = 2 * np.pi * (self.rpm/60) * (self.rhub**2 / self.r) / self.velocity
            x3 = 2 * np.pi * (self.rpm/60) * self.rhub / self.velocity
            gamma1 = self.velocity * vline * (x1**2 / (1 + x1**2)) / (self.rpm/60) / self.number_blades
            gamma2 = self.velocity * vline * (x2**2 / (1 + x2**2)) / (self.rpm/60) / self.number_blades
            gamma3 = self.velocity * vline * (x3**2 / (1 + x3**2)) / (self.rpm/60) / self.number_blades
            j = self.velocity / 2 / np.pi / (self.rpm / 60)
            f = (self.number_blades / 2) * ((self.rtip - self.r) / self.rtip) * (np.sqrt(1 + j**2) / j)
            F = (2/np.pi) * np.arccos(np.exp(-f))
            total_circulation = (gamma1 + gamma2 + gamma3) * F
            t_integrand = self.density * total_circulation * (self.omega * self.r) - self.number_blades * total_circulation/4/np.pi/self.r * self.dr
            return np.trapz(t_integrand, self.r) - self.thrust

        solution = root_scalar(thrust_func, bracket=[0.01, 150], method='brentq') 

        vline = solution.root * 1.06

        x1 = 2 * np.pi * (self.rpm/60) * self.r / self.velocity
        x2 = 2 * np.pi * (self.rpm/60) * (self.rhub**2 / self.r) / self.velocity
        x3 = 2 * np.pi * (self.rpm/60) * self.rhub / self.velocity

        gamma1 = self.velocity * vline * (x1**2 / (1 + x1**2)) / (self.rpm/60) / self.number_blades
        gamma2 = self.velocity * vline * (x2**2 / (1 + x2**2)) / (self.rpm/60) / self.number_blades
        gamma3 = self.velocity * vline * (x3**2 / (1 + x3**2)) / (self.rpm/60) / self.number_blades

        j = self.velocity / 2 / np.pi / (self.rpm / 60)
        f = (self.number_blades / 2) * ((self.rtip - self.r) / self.rtip) * (np.sqrt(1 + j**2) / j)
        F = (2/np.pi) * np.arccos(np.exp(-f))

        total_circulation = (gamma1 + gamma2 + gamma3) * F / self.number_blades

        circulation_function = interp1d(self.r[:-100], total_circulation[:-100], kind='cubic', fill_value="extrapolate")
        total_circulation = circulation_function(self.r)

        phi = np.arctan((self.velocity + vline) / (self.omega * self.r))
        alpha = np.array([5]) # degrees

        v_tang = total_circulation * self.number_blades / (4 * np.pi * self.r)
        v_rel = np.sqrt((self.velocity + vline) ** 2 + (self.omega * self.r) ** 2) - v_tang / np.cos(phi)

        cl, cd = get_aero_coefficients(alpha)
        
                    # these values I pulled out of my imagination because the paper doesn't provide airfoils, 
                    # don't change them and don't try to improve them, there are better things to do in life

        chord = 2 * total_circulation / v_rel / cl # in milimeters for visibility
        twist = alpha + np.degrees(phi)

        v_axial = vline - v_tang * np.tan(phi)

        self.dT = self.number_blades * self.density * total_circulation * (self.omega * self.r - v_tang)
        self.dQ = self.number_blades * self.density * total_circulation * self.r * v_axial

        self.error_percentage = 100 - self.thrust / (np.trapz(self.dT, self.r)) * 100
        self.shaft_power = np.trapz(self.dQ, self.r) * self.omega / self.number_blades
        
        return chord, twist
    
    def change_flight_regime(self, new_velocity, new_thrust):

        self.velocity = new_velocity
        self.thrust = new_thrust
        self.design_propeller_geometry()
        
        return
    
    def compute_coefficients(self, print_coeff=False):

        self.thrust_coefficient = self.thrust / (self.density * (self.rpm/60) ** 2 * self.diameter ** 4)
        self.power_coefficient = self.shaft_power / (self.density * (self.rpm/60) ** 3 * self.diameter ** 5)
        self.efficiency = self.advance_ratio * self.thrust_coefficient / self.power_coefficient
        if print_coeff:
            print("\nSummary of coefficients:")
            print("------------------------------------")
            
            print(f"Thrust Coefficient: {np.round(self.thrust_coefficient, 7)}")
            print(f"Power Coefficient: {np.round(self.power_coefficient, 7)}")
            print(f"New Efficiency: {np.round(self.efficiency*100, 2)}")
            print()


        return
    
    def summary_parameters(self, plot_chord_twist_distributions=False):

        print("\nSummary of Propeller Parameters:")
        print("------------------------------------")
        
        print(f"Freestream velocity: {self.velocity} [m/s]")
        print(f"Propeller diameter: {self.diameter} [m]")
        print(f"Propeller area: {np.round(self.propeller_area, 2)} [m^2]")
        print(f"Rotational speed: {np.round(self.rpm, 2)} [RPM], or {np.round(self.rpm * 2 * np.pi/60, 2)} [rad/s]")
        print(f"Propeller thrust: {np.round(self.thrust, 2)} [N]")
        print(f"Propeller power input: {np.round(self.shaft_power/1000, 2)} [kW]")
        print(f"Number of blades: {self.number_blades}")
        print(f"Propeller advance ratio: {np.round(self.advance_ratio, 5)}")
        print(f"Error Percentage: {np.round(100 - self.thrust / (np.trapz(self.dT, self.r)) * 100, 2)}%")

        if self.efficiency:
            print(f"Propeller efficiency: {np.round(self.efficiency*100, 2)}%")

        
        print()  # Extra space for readability

        if plot_chord_twist_distributions:
                fig, axs = plt.subplots(1, 2, figsize=(12, 5))

                # Plot Chord
                axs[0].plot(self.r / self.rtip, self.chord * 1000, color='red')
                axs[0].set_title('Chord Distribution')
                axs[0].set_xlabel('Position along blade radius')
                axs[0].set_ylabel('Chord lenght [mm]')
                axs[0].grid(True)



                # Plot Twist
                axs[1].plot(self.r / self.rtip, self.twist, color='red')
                axs[1].set_title('Twist Distribution')
                axs[1].set_xlabel('Position along blade radius')
                axs[1].set_ylabel('Twist angle [deg]')
                axs[1].grid(True)



                plt.show()

        return


