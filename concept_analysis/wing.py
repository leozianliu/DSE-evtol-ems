# Wing mass, area, drag as a function of mtow, v_cruise
class Wing:
    def __init__(self, mtow, v_cruise, S_disks, rho_air, g):
        self.mtow = mtow
        self.v_cruise = v_cruise
        self.S_disks = S_disks
        self.rho_air = rho_air
        self.g = g

    def calculate_wing_mass(self):
        # Constants for wing mass estimation (these values are placeholders and should be replaced with actual data)
        wing_area = 20  # m^2 (example value for wing area)
        wing_loading = self.mtow / wing_area  # N/m^2 (wing loading)
        wing_mass = 0.1 * wing_area * (wing_loading ** 0.5)  # kg (example formula for wing mass estimation)
        return wing_mass

    def calculate_drag_coefficient(self):
        # Constants for drag coefficient estimation (these values are placeholders and should be replaced with actual data)
        aspect_ratio = 10  # example value for aspect ratio
        C_D0 = 0.02  # zero-lift drag coefficient
        C_D = C_D0 + (1 / aspect_ratio) * (self.mtow / self.v_cruise) ** 2  # drag coefficient estimation
        return C_D