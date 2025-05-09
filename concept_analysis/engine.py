# Motor & driver mass as a function of power and torque


def calculatePropulsionMass(cont_power, peak_power, n_motor, D_rotor):#def calculatePropulsionMass(cont_power, peak_power, cont_torque, peak_torque, n_motor):
    #Returns mass of all motors and inverters combined
    #power and torque inputs PER MOTOR
    
    #power [kW]
    #torque [Nm]
    #masse [kg]
    #diameter [m]

    m_motor_peak_power = 0.12*peak_power+3.65
    m_motor_cont_power = 0.249*cont_power+1.56
    # m_motor_peak_torque = 0.0411*peak_torque+3.94
    # m_motor_cont_torque = 0.0995*cont_torque-0.117


    m_motor = max(m_motor_cont_power, m_motor_peak_power)#, m_motor_cont_torque, m_motor_peak_torque)

    m_inverter = 0.0187*cont_power + 0.433

    m_propeller = 7.5 * (D_rotor/2.3)**2 # Volocopter's rotor mass is 7.5 kg for 2.3 m diameter rotor, assumes a linear scaling with diameter

    return (m_motor + m_inverter + m_propeller)*n_motor

