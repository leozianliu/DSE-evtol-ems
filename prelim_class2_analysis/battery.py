import numpy as np

# Constants (based on INR18650-30Q)
Cp = 830  # J/kg.C
Wcell = 0.048  # kg
C_cell = 3.0  # Ah
VOC_const = 3.6  # Simplified average OCV
CTh = 2000  # Thevenin capacitance in F
eta_loss = 0.78

# Constraints
MAX_CRATE = 5.0
MAX_DOD = 0.8
MAX_TEMP = 50.0  # Celsius

# Air properties
k_air = 0.026  # W/m.C
rho_air = 1.2  # kg/m3
Cp_air = 1005  # J/kg.C
mu_air = 1.8e-5  # Pa.s

# Cooling parameters
A_cell = 0.01  # m2 (approximate surface area)
D_cell = 0.018  # m
hn = 16.25  # W/m2.C (natural convection)

# Calibration parameters (example values)
alpha = 5e-3
beta = -0.5


def voc_soc(soc):
    return 3.43 * np.exp(0.1938 * soc) - 0.5517 * np.exp(-5.0321 * soc)

def r0_soc(soc, T):
    return 0.0697 * np.exp(-12.4547 * soc) + 0.0219 * np.exp(0.1674 * soc)

def rth_soc(soc, T):
    return -0.0356 * soc + 0.0592

def calibrate_r0(r0, c_rate):
    return r0 + alpha / c_rate

def calibrate_rth(rth, c_rate):
    return rth * np.exp(beta * c_rate)

def battery_voltage(soc, i_out, vth):
    return voc_soc(soc) - vth - i_out * calibrate_r0(r0_soc(soc, 25), c_rate=i_out / C_cell)

def update_vth(vth, i_out, soc, dt):
    Rth = calibrate_rth(rth_soc(soc, 25), c_rate=i_out / C_cell)
    return vth + dt * ((i_out / CTh) - (vth / (Rth * CTh)))

def heat_generation(voc, vout, i_out, soc):
    entropy_term = i_out * (25 * dVoc_dT(soc))
    return i_out * (voc - vout) - entropy_term

def dVoc_dT(soc):
    return 1e-3 * (-36.827 * soc**6 + 109.94 * soc**5 - 119.77 * soc**4 +
                  58.206 * soc**3 - 13.498 * soc**2 + 2.4143 * soc - 0.3405)

def update_temperature(Tb, Qbat, dt):
    return Tb + (Qbat * dt) / (Wcell * Cp)

def air_cooling(Tb, Ta):
    return hn * A_cell * (Tb - Ta)

def check_constraints(crate, dod, temp):
    return crate <= MAX_CRATE and dod <= MAX_DOD and temp <= MAX_TEMP

# Simulate one mission step
def simulate_mission_step(power_demand_W, duration_s, Ncells, soc_init, Tb_init, Ta=25):
    dt = 1  # 1 second step
    steps = int(duration_s / dt)
    soc = soc_init
    Tb = Tb_init
    vth = 0.0

    for t in range(steps):
        p_cell = power_demand_W / Ncells
        i_out = p_cell / battery_voltage(soc, i_out=1, vth=vth)  # Initial guess
        for _ in range(3):  # Iterate for convergence
            vout = battery_voltage(soc, i_out, vth)
            i_out = p_cell / vout
            vth = update_vth(vth, i_out, soc, dt)

        dod = 1 - soc
        Qheat = heat_generation(voc_soc(soc), vout, i_out, soc)
        Qcool = air_cooling(Tb, Ta)
        Tb = update_temperature(Tb, Qheat - Qcool, dt)
        soc -= (i_out * dt) / (C_cell * 3600)

        crate = i_out / C_cell
        if not check_constraints(crate, dod, Tb):
            return False, crate, dod, Tb, soc

    return True, crate, dod, Tb, soc

# Simulate full mission profile
def simulate_mission_profile(profile, Ncells, soc_init=1.0, Tb_init=25):
    soc = soc_init
    Tb = Tb_init

    for idx, segment in enumerate(profile):
        name, power, duration = segment
        print(f"Segment {idx+1} ({name}): P = {power/1000:.1f} kW, duration = {duration}s")
        success, crate, dod, Tb, soc = simulate_mission_step(power, duration, Ncells, soc, Tb)
        if not success:
            print(f"FAILED at {name}: C-rate = {crate:.2f}, DoD = {dod:.2f}, Temp = {Tb:.2f}C")
            return
        else:
            print(f"  PASS: C-rate = {crate:.2f}, DoD = {dod:.2f}, Temp = {Tb:.2f}C, SoC = {soc:.2f}")

    print("Mission completed successfully.")

# Example usage:
if __name__ == "__main__":
    mission = [
        ("Takeoff", 150000, 60),
        ("Transition", 200000, 30),
        ("Cruise", 100000, 300),
        ("Descent", 80000, 60),
        ("Landing", 120000, 45)
    ]
    Ncells = 16000
    simulate_mission_profile(mission, Ncells)



