import numpy as np
import matplotlib.pyplot as plt

class TubeGeometry:

    def __init__(self, span, root_diameter, root_thickness, thickness_spanwise_steps = [0], thickness_steps = [0], diameter_spanwise_steps = [0], diameter_steps = [0], taper_ratio = 1, true_if_steps_false_if_tapered = True):
        # This generates the cross_section at the root
        self.half_span = span/2     # [m]
        self.taper_ratio = taper_ratio # 1 unless specified otherwise
        self.when_steps = true_if_steps_false_if_tapered # true unless specified otherwise
        self.thickness_ratio = (1-taper_ratio)*root_thickness/self.half_span
        self.diameter_ratio = (1-taper_ratio)*root_diameter/self.half_span

        # Values at the root
        self.root_diameter = root_diameter      # [m]
        self.root_thickness = root_thickness    # [m]

        # Steps in thickness are from one to the other, not from the root, all in mm
        self.thickness_spanwise_steps = thickness_spanwise_steps
        self.thickness_steps = thickness_steps

        # Steps in diameter are from one to the other, not from the root, all in mm
        self.diameter_spanwise_steps = diameter_spanwise_steps
        self.diameter_steps = diameter_steps



    def get_diameter(self, y):
        # Changes the diameter based on y-location
        diameter = self.root_diameter

        # if self.when_steps:
        #     for i in range(len(self.diameter_spanwise_steps)):
        #         if y > self.diameter_spanwise_steps[i]:
        #             diameter = diameter - self.diameter_steps[i]
        #         else:
        #             diameter = diameter

        #else:
        diameter = diameter - self.diameter_ratio * (y - 0.9)

        return diameter

    def get_thickness(self, y):
        # Changes the thickness based on y-location
        thickness = self.root_thickness

        if self.when_steps:
            for i in range(len(self.thickness_spanwise_steps)):
                if y > self.thickness_spanwise_steps[i]:
                    thickness = thickness - self.thickness_steps[i]
                else:
                    thickness = thickness

        else:
            thickness = thickness - self.thickness_ratio * y

        return thickness


    def get_inertia_tube(self, diameter, thickness):
        # For a thin-walled tube Ixx = Izz = pi*t*d^3/8, Ixz = 0, J = pi*t*d^3/4
        # Otherwise, solid circular section Ixx = Izz = pi*d^4/64, Ixz = 0, J = pi*d^4/32

        # Check if thin-walled assumption is valid at this point
        if thickness/diameter <= 0.1:
            true_if_thinwalled = True
        else:
            true_if_thinwalled = False

        # Returns values in [m^4]
        if true_if_thinwalled:
            moment_Ix = (np.pi * thickness * diameter ** 3 /8)
            moment_J = (np.pi * thickness * diameter ** 3 /4)
            radius = diameter /2

        else:
            out_diameter = diameter + thickness
            in_diameter = diameter - thickness

            moment_Ix = (np.pi * (out_diameter ** 4 - in_diameter ** 4) /64)
            moment_J = (np.pi * (out_diameter ** 4 - in_diameter ** 4) /32)

            radius = out_diameter /2

        return moment_Ix, moment_J, radius

    def get_tube_matrix(self, num_points_span):

        y_values = np.linspace(0, self.half_span, num_points_span)

        thicknesses = []
        inertias_Ix = []
        inertias_J = []
        radii_out = []

        for y in y_values:
            t = self.get_thickness(y)
            thicknesses.append(t)
            d = self.get_diameter(y)
            Ix, J, r = self.get_inertia_tube(d, t)
            inertias_Ix.append(Ix)
            inertias_J.append(J)
            radii_out.append(r)

        matrix = np.array([thicknesses, inertias_Ix, inertias_J, radii_out, y_values])

        return matrix
    

# tube = TubeGeometry(12, 50, 5, [1, 2, 3, 4, 5], [1, 1, 1, 1, 0], [1, 2, 3, 4, 5], [1, 1, 1, 1, 0])
# tube.get_diameter()

# matrix = tube.get_tube_matrix(100)
# radii_out = matrix[3]
# print(radii_out)





