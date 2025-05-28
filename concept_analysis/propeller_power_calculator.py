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

def circulation_distribution(v0, vline, rpm, number_blades, r):
    
    x = 2 * np.pi * (rpm/60) * r / v0
    return v0 * vline / (rpm/60) / number_blades * (x**2 / (1 + x**2))


class Propeller:
    def __init__(self, shaft_power, thrust, velocity, diameter, RPM, number_blades, density=1.225, temperature=293.15, ducted_fan=False):
        
        # Inputs, all are floats
        self.shaft_power = shaft_power # shaft power after motor efficiency in Watts
        self.thrust = thrust # motor 
        self.velocity = velocity # velocity in m/s
        self.density = density # ambient density in kg/m^3
        self.temperature = temperature # ambient temperature in Kelvin
        self.diameter = diameter # Propeller diameter in m
        self.propeller_area = np.pi * diameter ** 2 / 4 # Propeller area in m^2
        self.gamma_specific_speeds_ratio = 1.4
        self.R_air_constant = 287 # Air constant in J/(kg Kelvin)
        self.maximum_blade_mach = 0.8 # Maximum Mach number to avoid shockwaves on the blade tips
        self.speed_of_sound = np.sqrt(self.gamma_specific_speeds_ratio * self.temperature * self.R_air_constant)
        self.RPM = RPM # Rotations per minute of propeller blades
        self.number_blades = number_blades # Number of blades of the propeller
        self.hub_ratio = 0.15
        self.ducted_fan = ducted_fan
        self.assumed_ducted_efficiency_increase = 1
        self.r = np.linspace(self.hub_ratio * self.diameter/2, self.diameter/2, 1000)
        self.dr = self.r[1] - self.r[0]
        self.omega = 0.10472 * self.RPM

        if self.ducted_fan:
            self.assumed_ducted_efficiency_increase += 0.15

        # These will be computed

        self.efficiency = None
        self.maximum_RPM = None
        self.maximum_tip_speed = None
        self.advance_ratio = None
        self.Reynolds_number = None
        self.power_coefficient = None
        self.thrust_coefficient = None
        
        
    def design_propeller_geometry(self, true_if_compute_power_false_if_compute_thrust=True):

        # === Assumptions: ===
        # RPMs higher than 1000
        # No tip losses (might add later)
        # Addition of duct increases efficiency by exactly 15%
        # Error is around +- 5-10 % 

        # === Step 1: Initialize ===

        N = 1000 # Number of panels
        epsilon = 0 # initial guess

        alpha = 0  # Assume 5 deg fixed angle of attack of blades
        Cl, Cd = get_aero_coefficients([alpha])
        Cl = 0.4
        Cd = 0.02
        s = Cd/Cl
        
        r_hub = self.hub_ratio * self.diameter/2
       
        converged = False
        max_iterations = 1000
        epsilon_tol = 1e-6


        self.dr = self.r[1] - self.r[0]
        
        
        lamda = self.velocity / (self.omega * (self.diameter/2))  # speed ratio

        self.advance_ratio = self.velocity / (self.RPM / 60) / self.diameter   


        xi = self.r / (self.diameter/2)
        x = xi / lamda

        if true_if_compute_power_false_if_compute_thrust:

            self.thrust_coefficient =  2 * self.thrust / (self.density * self.velocity ** 2 * np.pi * (self.diameter/2) ** 2)
        else: 
            
            self.power_coefficient = 2 * self.shaft_power / (self.density * self.velocity ** 3 * np.pi * (self.diameter/2) ** 2)
        
        # Estimate power required from thust
        while not converged and max_iterations:

            tan_phi_tip = lamda * (1 + epsilon / 2)
            phi_tip = np.arctan(tan_phi_tip)

            phi = np.arctan(tan_phi_tip / xi)
            F = (2 / np.pi) * np.degrees(np.arccos(np.exp(-(self.number_blades/2) * (1 - xi) / (np.sin(phi_tip)))))
            G = F * np.cos(phi) * np.sin(phi) * x

            Wc = (4 * np.pi * lamda * G * self.velocity * (self.diameter/2) * epsilon) / (Cl * self.number_blades)
            
            a = (epsilon/2) * np.cos(phi)**2 * (1 - s * np.tan(phi))

            W = self.velocity * (1 + a) / np.sin(phi)
            a_prime = 0.5 * epsilon * np.cos(phi) * np.sin(phi) * (1 + s / np.tan(phi)) / x

            chord = Wc / W
            
            beta = np.radians(alpha) + phi

            # === Eq. 11 integrands ===
            I1_integrand = 4 * G * xi * (1 - s * np.tan(phi)) 
            I2_integrand = lamda * (I1_integrand/ xi / 2) * (1 + s / np.tan(phi)) * np.sin(phi) * np.cos(phi)
            J1_integrand = 4 * xi * G * (1 + s / np.tan(phi))
            J2_integrand = (J1_integrand / 2) * (1 - s * np.tan(phi)) * np.cos(phi)**2
            
            self.dL = self.density * W * (2 * np.pi * self.velocity**2 * epsilon * G / self.omega) * self.dr
            self.dT = self.dL * np.cos(phi) * (1 - s * np.tan(phi))
            self.dQ = 2 * np.pi * xi * (self.diameter/2) * self.density * self.velocity * (1 + a) * (2 * self.omega * xi * (self.diameter/2) * a_prime * F) * self.dr

            # Integrate (Eq. 12 & 13)
            I1 = np.trapz(I1_integrand, xi)
            I2 = np.trapz(I2_integrand, xi)
            J1 = np.trapz(J1_integrand, xi)
            J2 = np.trapz(J2_integrand, xi)
            
            
            if true_if_compute_power_false_if_compute_thrust:
                epsilon_new = (I1 / I2 / 2) - np.sqrt((I1 / I2 / 2)**2 - self.thrust_coefficient/I2)
            else:
                epsilon_new = - (J1 / 2 / J2) + np.sqrt((J1 / 2 / J2) ** 2 + self.power_coefficient/J2)

            max_iterations -= 1

            if abs(epsilon_new - epsilon) < epsilon_tol:
                converged = True
            
            epsilon = epsilon_new
        
        if true_if_compute_power_false_if_compute_thrust:

            self.power_coefficient = J1 * epsilon + J2 * epsilon ** 2
            self.shaft_power = (1/self.assumed_ducted_efficiency_increase) * 0.5 * self.power_coefficient * self.density * self.velocity ** 3 * np.pi * (self.diameter/2) ** 2
            
        else:
            self.thrust_coefficient = I1 * epsilon - I2 * epsilon ** 2
            self.thrust = self.assumed_ducted_efficiency_increase * 0.5 * self.thrust_coefficient * self.density * self.velocity ** 2 * np.pi * (self.diameter/2) ** 2

        self.chord_distribution = chord
        self.twist_distribution = np.rad2deg(beta)

        return              

    def compute_maximum_RPM(self, use_maximum_RPM=False):
        
        self.maximum_tip_speed = self.maximum_blade_mach * self.speed_of_sound 
        self.maximum_RPM = 60 * self.maximum_tip_speed / np.pi / self.diameter

        if use_maximum_RPM:
            self.RPM = self.maximum_RPM
            self.advance_ratio = self.velocity / (self.RPM / 60) / self.diameter 

        return 
    
    def adjust_chord_for_desired_Reynolds(self, desired_Reynolds_number=1e5):
        if desired_Reynolds_number:

            self.Reynolds_number = (self.density * self.r * self.omega) * self.chord_distribution / 1.81e-5

            low_reynolds_element = round(((0.4 - self.hub_ratio) / (1 - (self.hub_ratio / self.number_blades))) * self.number_blades)
            high_reynolds_element = round(((0.95 - self.hub_ratio) / (1 - (self.hub_ratio / self.number_blades))) * self.number_blades)
            low_reynolds_number = self.Reynolds_number[high_reynolds_element]
            high_reynolds_number = self.Reynolds_number[low_reynolds_element]
            min_reynolds_number = min(low_reynolds_number,high_reynolds_number)


            self.chord_distribution = (desired_Reynolds_number/min_reynolds_number) * self.chord_distribution
            self.Reynolds_number = (self.density * (self.r * self.RPM * 0.10472) * self.chord_distribution) / 1.81e-5

        else:

            chord_factor = 0.2 / np.max(self.chord_distribution)
            self.chord_distribution = self.chord_distribution * chord_factor
            self.Reynolds_number = (self.density * self.r * self.RPM * 0.10472) * self.chord_distribution / 1.81e-5

        return
        
    def compute_performance_parameters(self):

        a = np.ones_like(self.r) * 0.3
        aline = np.ones_like(self.r) * 0.01
        Niterations = 100
        Erroriterations = 1e-7
        sigma = self.number_blades * (self.chord_distribution + 0.0000001) / (2 * np.pi * self.r)
        anew = np.zeros_like(a)

        for _ in range(Niterations):

            Uaxial = self.velocity * (1 - a)
            Utan = (1 + aline) * self.omega * self.r
            inflowangle = np.arctan(Uaxial / Utan)
            vrel = np.sqrt(Uaxial ** 2 + Utan ** 2)

            alpha = self.twist_distribution - np.degrees(inflowangle)

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
            
        self.dL = 0.5 * self.density * vrel ** 2 * self.chord_distribution * cl
        self.dD = 0.5 * self.density * vrel ** 2 * self.chord_distribution * cd

        self.dT = self.number_blades * (self.dL * np.cos(inflowangle) - self.dD * np.sin(inflowangle))
        self.dQ = self.r * self.number_blades * (self.dL * np.sin(inflowangle) + self.dD * np.cos(inflowangle))
        self.thrust = np.sum(self.dT) * self.dr
        self.shaft_power = np.sum(self.dQ) * self.dr * self.omega

        self.thrust_coefficient =  self.thrust / (self.density * (self.RPM / 60)**2 * self.diameter ** 4)
        self.power_coefficient = self.shaft_power / (self.density * (self.RPM / 60)**3 * self.diameter ** 5)

        self.advance_ratio = self.velocity / (self.RPM / 60) / self.diameter
        self.efficiency = self.assumed_ducted_efficiency_increase * self.advance_ratio * self.thrust_coefficient / self.power_coefficient

        self.power_loading = self.thrust / self.shaft_power
        self.disk_loading = self.thrust / self.propeller_area
        
        self.figure_of_merit = self.thrust ** 1.5 / (self.shaft_power * np.sqrt(2 * self.density * self.propeller_area))

        return
    
    def high_loading_design(self):

        # Constants
        J = self.velocity / 2 / np.pi / (self.RPM / 60) / (self.diameter/2)
        f = self.number_blades / 2 * ((self.diameter/2 - self.r) / (self.diameter/2)) * (np.sqrt(1 + J ** 2) / J)
        F = 2 / np.pi * (np.arccos(np.exp(-f)))
        assumed_cl = 1.5

        # Pre-compute constant multiplier
        omega = 2 * np.pi * self.RPM / 60  # Angular velocity in rad/s

        def compute_thrust(vline):
            circulation = (
                circulation_distribution(self.velocity, vline, self.RPM, self.number_blades, self.r) +
                circulation_distribution(self.velocity, vline, self.RPM, self.number_blades, (self.hub_ratio * self.diameter / 2) ** 2 / self.r) -
                circulation_distribution(self.velocity, vline, self.RPM, self.number_blades, self.hub_ratio * self.diameter / 2 * np.ones_like(self.r)) 
            ) * F

            thrust_distribution = (
                self.density * circulation * (omega * self.r - self.number_blades * circulation / (4 * np.pi * self.r)) * self.dr
            )
            return np.sum(thrust_distribution)

        def residual(vline):
            return compute_thrust(vline) - self.thrust

        # Solve for vline that makes thrust match desired value
        solution = root_scalar(residual, bracket=[0.01, 200], method='brentq')  # adjust bracket as needed

        if not solution.converged:
            print("Failed to converge to a vline value.")
        print(solution.root)
        vline = solution.root

        # Recalculate final circulation with solved vline
        total_circulation = (
            circulation_distribution(self.velocity, vline, self.RPM, self.number_blades, self.r) +
            circulation_distribution(self.velocity, vline, self.RPM, self.number_blades, (self.hub_ratio * self.diameter / 2) ** 2 / self.r) - 
            circulation_distribution(self.velocity, vline, self.RPM, self.number_blades, self.hub_ratio * self.diameter / 2) 
        ) * F

        # Store results
        self.vline = vline
        self.total_circulation = total_circulation
        
        # interpolate away from blade tip and extrapolate up to blade tip
        circulation_function = interp1d(self.r[:800], total_circulation[:800], kind='cubic', fill_value='extrapolate')
        total_circulation = np.array([circulation_function(R) for R in self.r])

        # Find W
        phi = np.arctan((self.velocity + vline) / self.omega / self.r)
        v_tang = total_circulation * self.number_blades / 4 / np.pi / self.r
        W = np.sqrt((self.velocity + vline)**2 + (self.omega * self.r)**2) - v_tang / np.cos(phi)

        # Find initial chord
        #self.chord_distribution = 2 * total_circulation / W / assumed_cl

        # Find optimum angle of attack and lift coefficient from reynolds, mach and chord, and update chord
        alpha_optim = 5
        cl, cd = get_aero_coefficients([alpha_optim])
        self.chord_distribution = 2 * total_circulation / W / cl

        # Find twist distribution
        self.twist_distribution = np.degrees(phi) + alpha_optim

        # Find lift and drag
        self.dT = 0.5 * self.density * W ** 2 * self.chord_distribution * (cl * np.cos(phi) -  cd * np.sin(phi)) * self.dr
        self.dQ = 0.5 * self.density * W ** 2 * self.chord_distribution * self.omega * (cl * np.sin(phi) + cd * np.cos(phi)) * self.r * self.dr

        #self.thrust = np.sum(self.dT)
        correction = self.thrust / np.sum(self.dT)
        self.shaft_power = correction * self.omega * np.sum(self.dQ) / 1000 

        self.advance_ratio = self.velocity / (self.RPM / 60) / self.diameter
        self.efficiency = self.thrust * self.velocity / self.shaft_power
        return
    
    def compute_noise_signature(self, n_harmonics=np.array([1, 2, 3, 4, 5])):
        self.sound_pressure_level = np.array([])
        
        for n in n_harmonics:
            
            observer_distance = np.linspace(0, 200, 200)
            Ri = 100
            r = np.linspace(self.hub_ratio * self.diameter/2, self.diameter/2, 1000)
            dr = r[1] - r[0]

            omega = n * self.number_blades * 0.10472 * self.RPM
            decay = 1 / (1 + (omega / (self.number_blades * omega/n))**2)
            self.circulation = self.dL / self.density / self.velocity
            self.source_strenghts = self.circulation * self.chord_distribution

            k     = omega / self.speed_of_sound
            p_L   = np.sum( decay * 1j * omega / (4 * np.pi * self.speed_of_sound) * self.circulation * dr * np.exp(-1j * k * Ri) / Ri)
            p_T   = np.sum( decay * -1j * omega * self.density/(4*np.pi*self.speed_of_sound) * self.source_strenghts * dr * np.exp(-1j * k * Ri) / Ri)
            p_hat = p_L + p_T
            self.sound_pressure_level = np.append(self.sound_pressure_level , 20 * np.log10(abs(p_hat)/2e-6))

        plt.scatter(n_harmonics, self.sound_pressure_level)
        plt.grid()
        plt.show()
        return self.sound_pressure_level

    def summary_parameters(self, plot_chord_twist_distributions=False):

        print("\nSummary of Propeller Parameters:")
        print("------------------------------------")
        
        print(f"Freestream velocity: {self.velocity} [m/s]")
        print(f"Propeller diameter: {self.diameter} [m]")
        print(f"Propeller area: {np.round(self.propeller_area, 2)} [m^2]")
        print(f"Rotational speed: {np.round(self.RPM, 2)} [RPM], or {np.round(self.RPM * 0.10472, 2)} [rad/s]")
        print(f"Propeller thrust: {np.round(self.thrust, 2)} [N]")
        print(f"Propeller power input: {np.round(self.shaft_power/1000, 2)} [kW]")
        print(f"Number of blades: {self.number_blades}")
        print(f"Propeller advance ratio: {np.round(self.advance_ratio, 5)}")
        print(f"Propeller efficiency: {np.round(self.efficiency*100, 2)}%")
        #print(f"Propeller tonal noise: {np.round(self.sound_pressure_level, 2)} [dB]")

        
        print()  # Extra space for readability

        if plot_chord_twist_distributions:
                fig, axs = plt.subplots(1, 2, figsize=(12, 5))

                # Plot Chord
                axs[0].plot(self.r / (self.diameter/2), self.chord_distribution/ (self.diameter/2), color='red')
                axs[0].set_title('Chord Distribution')
                axs[0].set_xlabel('Position along blade radius')
                axs[0].set_ylabel('Chord lenght [m]')
                axs[0].grid(True)



                # Plot Twist
                axs[1].plot(self.r / (self.diameter/2), self.twist_distribution, color='red')
                axs[1].set_title('Twist Distribution')
                axs[1].set_xlabel('Position along blade radius')
                axs[1].set_ylabel('Twist angle [deg]')
                axs[1].grid(True)



                plt.show()

        return


