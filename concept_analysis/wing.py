# import sys
# from pathlib import Path
# ROOT_DIR = Path(__file__).parent.parent
# sys.path.append(ROOT_DIR) # Replace with your own path so you can use .py files from other directories

import numpy as np
from helper import *

# Wing mass, area, drag as a function of mtow, v_cruise
class Wing:
    def __init__(self, mtow, v_cruise, b_span):
        self.mtow = mtow
        self.mtom = N_2_kg(mtow)
        self.v_cruise = v_cruise
        self.b_span = b_span

        self.wing_area = None
        self.calculate_wing_area()
        self.wing_loading = None
        self.calculate_wing_loading()
        self.wing_mass = None
        self.calculate_wing_mass()

    def calculate_wing_area(self):
        CL_cruise = 0.8
        self.wing_area = (self.mtow / (0.5 * CL_cruise * 1.225 * self.v_cruise ** 2))  # m^2 (using lift equation)
        return self.wing_area
    
    def calculate_wing_loading(self):
        # Wing loading is the ratio of weight to wing area
        self.wing_loading = self.mtom / self.wing_area
        return self.wing_loading

    def calculate_wing_mass(self):
        # Constants for wing mass estimation (these values are placeholders and should be replaced with actual data)
        b_ref = 1.905 # m
        n_ult = 3.8 * 1.5 # ultimate load factor (example value)
        cantilever_ratio = 40 # example value for cantilever ratio
        
        self.wing_mass = self.mtom * 4.9e-3 * self.b_span**0.75 * (1 + np.sqrt(b_ref / self.b_span)) * n_ult**0.55 * (cantilever_ratio / self.wing_loading)**0.3 # kg (example formula for wing mass estimation)
        return self.wing_mass

    # def calculate_drag_coefficient(self):
    #     # Constants for drag coefficient estimation (these values are placeholders and should be replaced with actual data)
    #     aspect_ratio = 10  # example value for aspect ratio
    #     C_D0 = 0.02  # zero-lift drag coefficient
    #     C_D = C_D0 + (1 / aspect_ratio) * (self.mtow / self.v_cruise) ** 2  # drag coefficient estimation
    #     return C_D