from const_and_convers import *
import numpy as np

class Propulsion_powertrain:
    """
    Class to represent a propulsion powertrain.
    All units are SI units unless otherwise specified.
    """

    def __init__(self, mto_mass, l_over_d, cruise_speed, equipment_power=0):
        self.mto_weight = mto_mass * g
        self.l_over_d = l_over_d
        self.cruise_speed = cruise_speed
        self.equipment_power = equipment_power
        self.pitch_fan_power = 0 # for now, later add ducted fan power for pitch control
        self.thrust_ratio_inner_outer = 2 # thrust ratio for inner (higher) and outer propellers (lower)

        self.bat_energy_dens_wh = 280 # Wh/kg
        self.bat_volme_dens = 450 # Wh/m^3
        self.motor_power_dens = 9000 # W/kg
        self.motor_number = 4
        self.propeller_eff_vtol = 1#0.85
        self.propeller_eff_cruise = 1#0.85
        self.motor_eff = 0.95

        self.mtow_margin = 1
        self.factor_k = 1.15 # induced power from thrust for hovering
        self.figure_of_merit = 0.8 # converting induced to total power for hovering
        self.cruise_power_SF = 1 # safety factor for cruise power
        self.vtol_power_SF = 1 # safety factor for VTOL power
        self.take_off_time = 120 # seconds per round trip
        self.landing_time = 120 # seconds per round trip
        self.cruise_time = 2160 # seconds per round trip
        self.climb_speed = 5 # m/s

        self.propulsion_eff_vtol = self.propeller_eff_vtol * self.motor_eff
        self.propulsion_eff_cruise = self.propeller_eff_cruise * self.motor_eff


    def calculate_hover_power(self):
        hover_power_induced = self.factor_k * (self.mto_weight * self.mtow_margin)**1.5 / np.sqrt(2 * rho_air * self.propulsion_eff_vtol)
        hover_power_total = hover_power_induced / self.figure_of_merit
        induced_velo = hover_power_induced / (self.mto_weight * self.mtow_margin)
        clb_hov_velo_ratio = self.climb_speed / induced_velo
        clb_hov_power_ratio = 0.5 * clb_hov_velo_ratio + np.sqrt((clb_hov_velo_ratio / 2)**2 + 1)
        climb_power = hover_power_total * clb_hov_power_ratio
        self.climb_power = climb_power * self.vtol_power_SF / self.propulsion_eff_vtol
        self.hover_power = hover_power_total * self.vtol_power_SF / self.propulsion_eff_vtol
        self.descend_power = self.hover_power

    def calculate_cruise_power(self):
        self.cruise_power = (self.mto_weight * self.mtow_margin) * self.cruise_power_SF  / (self.l_over_d * self.propulsion_eff_cruise)
    
    def calculate_energy(self):
        self.calculate_hover_power()
        self.calculate_cruise_power()

        self.mission_energy = self.climb_power * self.take_off_time + self.descend_power * self.landing_time + self.cruise_power * self.cruise_time

    def calculate_single_motor_inverter_mass(single_motor_power):
        return 0.0913 * single_motor_power * 1000 + 4.63

    def calculate_proptrain_mass(self):
        self.calculate_energy()
        self.battery_mass = self.mission_energy / (self.bat_energy_dens_wh * 3600)
        self.battery_volume = self.mission_energy / (self.bat_volme_dens * 3600)
        self.inner_motor_power_nominal = self.climb_power * self.thrust_ratio_inner_outer / (1 + self.thrust_ratio_inner_outer)
        self.outer_motor_power_nominal = self.climb_power / (1 + self.thrust_ratio_inner_outer)
        self.inner_motor_power_max = self.hover_power * self.thrust_ratio_inner_outer / (1 + self.thrust_ratio_inner_outer)
        self.outer_motor_power_max = self.hover_power / (1 + self.thrust_ratio_inner_outer)
        self.inner_motor_esc_mass = self.calculate_single_motor_inverter_mass(self.inner_motor_power_nominal) * (self.motor_number / 2) # 2/4 big motors
        self.outer_motor_esc_mass = self.calculate_single_motor_inverter_mass(self.outer_motor_power_nominal) * (self.motor_number / 2) # 2/4 small motors
        