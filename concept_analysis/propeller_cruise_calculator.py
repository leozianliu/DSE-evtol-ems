import numpy as np

class Propeller:
    def __init__(self, shaft_power, freestream_velocity, freestream_density, propeller_diameter):
        
        # Inputs, all are floats
        self.shaft_power = shaft_power
        self.freestream_velocity = freestream_velocity
        self.freestream_density = freestream_density
        self.propeller_diameter = propeller_diameter
        self.propeller_area = np.pi * propeller_diameter ** 2 / 4

        # These will be computed
        self.downstream_velocity = None
        self.propeller_thrust = None
        self.ideal_efficiency = None
        self.real_efficiency = None

    def compute_thrust(self):
        if not self.propeller_thrust:
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
            else:
                #raise ValueError("Multiple or no real and positive downstream velocity solutions found in Propeller.compute_thrust")
                self.downstream_velocity = self.freestream_velocity
                return 0

            propeller_airspeed = (self.freestream_velocity + self.downstream_velocity) / 2
            self.propeller_thrust = (
                self.freestream_density * propeller_airspeed * self.propeller_area *
                (self.downstream_velocity - self.freestream_velocity)
            )

        return self.propeller_thrust

    def compute_efficiency(self):
        if self.propeller_thrust is None:
            raise RuntimeError("Must compute thrust before efficiency.")

        self.ideal_efficiency = 2 * self.freestream_velocity / (self.freestream_velocity + self.downstream_velocity)
        self.real_efficiency = self.propeller_thrust * self.freestream_velocity / self.shaft_power

        return self.ideal_efficiency, self.real_efficiency
