class Propulsion_powertrain:
    """
    Class to represent a propulsion powertrain.
    """

    def __init__(self, name, power, efficiency):
        """
        Initialize the propulsion powertrain with a name, power, and efficiency.

        :param name: Name of the propulsion powertrain
        :param power: Power of the propulsion powertrain in kW
        :param efficiency: Efficiency of the propulsion powertrain (0-1)
        """
        self.name = name
        self.power = power
        self.efficiency = efficiency

    def 