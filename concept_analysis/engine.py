# Motor & driver mass as a function of power and torque


def calculatePropulsionMass(cont_power, peak_power, cont_torque, peak_torque, n_motor):
    #Returns mass of all motors and inverters combined
    #power and torque inputs PER MOTOR
    
    #power [kW]
    #torque [Nm]
    #masse [kg]

    m_motor_peak_power = 0.12*peak_power+3.65
    m_motor_cont_power = 0.249*cont_power+1.56
    m_motor_peak_torque = 0.0411*peak_torque+3.94
    m_motor_cont_torque = 0.0995*cont_torque-0.117


    m_motor = max(m_motor_cont_power, m_motor_peak_power, m_motor_cont_torque, m_motor_peak_torque)

    m_inverter = 0.0187*cont_power + 0.433

    return (m_motor + m_inverter)*n_motor

