import numpy as np
import matplotlib.pyplot as plt

class TubeGeometry:

    def __init__(self, span, root_diameter, root_thickness, taper_ratio = 1, true_if_steps_false_if_tapered = True):
        # This generates the cross_section at the root
        self.half_span = span/2     # [m]
        self.taper_ratio = taper_ratio # 1 unless specified otherwise
        self.when_steps = true_if_steps_false_if_tapered # true unless specified otherwise
        self.thickness_ratio = (1-taper_ratio)*root_thickness/self.half_span
        self.diameter_ratio = (1-taper_ratio)*root_diameter/self.half_span

        # Values at the root
        self.root_diameter = root_diameter          # [mm]
        self.root_thickness = root_thickness    # [mm]



    def get_diameter(self, y, diameter_spanwise_steps, diameter_steps):
        # Changes the diameter based on y-location, steps in diameter are from one to the other, not from the root
        diameter = self.root_diameter

        if self.when_steps:
            for i in range(len(diameter_spanwise_steps)):
                if y > diameter_spanwise_steps[i]:
                    diameter = diameter - diameter_steps[i]
                else:
                    diameter = diameter

        else:
            diameter = diameter - self.diameter_ratio * y/self.half_span

        return diameter

    def get_thickness(self, y, thickness_spanwise_steps, thickness_steps):
        # Changes the thickness based on y-location, steps in thickness are from one to the other, not from the root
        thickness = self.root_thickness

        if self.when_steps:
            for i in range(len(thickness_spanwise_steps)):
                if y > thickness_spanwise_steps[i]:
                    thickness = thickness - thickness_steps[i]
                else:
                    thickness = thickness

        else:
            thickness = thickness - self.thickness_ratio * y/self.half_span

        return thickness


    def inertia_tube(self, diameter, thickness):
        # For a thin-walled tube Ixx = Izz = pi*t*d^3/8, Ixz = 0, J = pi*t*d^3/4
        # Otherwise, solid circular section Ixx = Izz = pi*d^4/64, Ixz = 0, J = pi*d^4/32

        # Check if thin-walled assumption is valid at this point
        if thickness/diameter <= 0.1:
            true_if_thinwalled = True
        else:
            true_if_thinwalled = False

        # Returns values in [mm^4]
        if true_if_thinwalled:
            moment_Ix = np.pi * thickness * diameter ** 3 /8
            moment_J = np.pi * thickness * diameter ** 3 /4

        else:
            out_diameter = diameter + thickness
            in_diameter = diameter - thickness

            moment_Ix = np.pi * (out_diameter ** 4 - in_diameter ** 4) /64
            moment_J = np.pi * (out_diameter ** 4 - in_diameter ** 4) /32

        return moment_Ix, moment_J





