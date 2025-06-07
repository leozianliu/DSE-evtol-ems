import json
import numpy as np
from xfoil.xfoil import XFoil
from naca import naca
np.set_printoptions(suppress=True, formatter={'float_kind':'{:0.4f}'.format})

def t(b_sec, b_m, b): # b is the spanwise fraction of the wing, from 0 at the root to b/2 at the tip
    if b_sec <= b_m:
        return t_r - (t_r - t_m) * (b_sec / b_m) # Linear interpolation between root and first engine
    else:
        return t_m - (t_m - t_t) * ((b_sec - b_m) / (b/2 - b_m)) # Linear interpolation between first engine and tip

def analyse(cl, code, Re):
    xf = XFoil()
    xf.airfoil = naca(code, 200, finite_TE = False, half_cosine_spacing = True, rotation = 0)
    xf.Re = Re
    xf.max_iter = 40
    a, cd = xf.cl(cl)[0:2]
    return a, cd

def average_array(array):
    return (array[:1]+array[1:])/2

def design_wing(lift, V, b, t_r, t_m, t_t, stalldelay_tip):

    Nr_sections = 9 # Number of sections per half-span, including the tip section.
    min_analyses = 3 # Number of xfoil analyses per section
    b_m = 0.2 * b # First engine span in meters

    rootsec_stall = 0.2 # Spanwise fraction that stalls first
    midsec_stall = 0.5 # Spanwise fraction that stalls next

    # Model parameters deciding what thicknesses must be able to end up in every section:
    min_t_max = 0.21
    max_t_min = 0.09

    # Section calculations
    rho = 1.225
    L = 1.0  # Characteristic length in m
    Mu = 1.81e-5  # Dynamic viscosity of air in kg/(m·s)
    Re = rho * V * L / Mu
    Total_circ = lift / (0.5 * rho * V**2 * b) # Circulation per meter of span in c * Cl
    span_sections = np.linspace(0, b/2, num=Nr_sections)  # Spanwise sections from root to tip
    circ_distr = np.sqrt(1-np.linspace(0, 1, num=Nr_sections)**2) * Total_circ * 4 / np.pi # Quarter ellipse with area equal to total_circ
    t_sections = np.array([t(b_sec, b_m, b) for b_sec in span_sections])  # Thickness at each section

    with open("airfoildata.json", "r") as fp:
        # Load the dictionary from the file
        airfoildata = json.load(fp)

    airfoildata_array = np.array(list(airfoildata.values()))
    # Add a column for airfoil names
    airfoildata_array = np.hstack((airfoildata_array, (np.array([int('1'+k+'1') for k in airfoildata.keys()])).reshape(-1, 1)))
    # Add a column for t/c (thickness ratio)
    airfoildata_array = np.hstack((airfoildata_array, (np.array([int(k[2:4]) for k in airfoildata.keys()])/100).reshape(-1, 1)))
    # Add a column for Cl_cruise / (t/c)
    airfoildata_array = np.hstack((airfoildata_array, (airfoildata_array[:, 1]/airfoildata_array[:, -1]).reshape(-1, 1)))
    # Add a column for Alpha_max - Alpha_cruise
    airfoildata_array = np.hstack((airfoildata_array, (airfoildata_array[:, 3] - airfoildata_array[:, 0]).reshape(-1, 1)))
    # Add a column for Cl/Cd_max
    airfoildata_array = np.hstack((airfoildata_array, (airfoildata_array[:, 1]/airfoildata_array[:, 2]).reshape(-1, 1)))

    t_max = min(np.max(airfoildata_array[:, -4]), min_t_max)  # Maximum thickness ratio
    t_min = max(np.min(airfoildata_array[:, -4]), max_t_min)  # Minimum thickness ratio

    # Sort the airfoildata_array by stall margin
    sorted_indices = np.argsort(airfoildata_array[:, -2])
    airfoildata_array_sorted = airfoildata_array[sorted_indices]

    t_max = airfoildata_array_sorted[np.where(airfoildata_array_sorted[:, -4] == t_max)[0], -2]
    t_min = airfoildata_array_sorted[np.where(airfoildata_array_sorted[:, -4] == t_min)[0], -2]
    stallangle = max(t_max[0], t_min[0])
    stallangle_stop = min(t_max[-1], t_min[-1]) - stalldelay_tip

    assert stallangle_stop > stallangle, "Decrease stalldelay or provide more airfoils!"
        
    LD_min = np.min(airfoildata_array_sorted[:, -1])  # Minimum L/D ratio

    previous_score = -1
    final_stallangle = stallangle_stop
    while stallangle <= stallangle_stop:
        # Score the airfoils in each category based L/D
        LD_r = airfoildata_array_sorted[np.where(airfoildata_array_sorted[:, -2] <= stallangle)[0], -1]
        LD_m = airfoildata_array_sorted[np.where((stallangle <= airfoildata_array_sorted[:, -2]) & (airfoildata_array_sorted[:, -4] <= stallangle+stalldelay_tip))[0], -1]
        LD_t = airfoildata_array_sorted[np.where(airfoildata_array_sorted[:, -2] >= stallangle+stalldelay_tip)[0], -1]
        score = np.sum(LD_r)**2 + np.sum(LD_m)**2 + np.sum(LD_t)**2
        if score < previous_score:
            final_stallangle = stallangle
        stallangle += 0.2

    print("Final stall angle:", final_stallangle)

    aoa_margin = 0.5  # Margin in degrees for AoA tilt

    airfoils_root = airfoildata_array_sorted[np.where(airfoildata_array_sorted[:, -2] <= final_stallangle + aoa_margin)[0], :]
    airfoils_mid = airfoildata_array_sorted[np.where((final_stallangle  - aoa_margin <= airfoildata_array_sorted[:, -2]) & (airfoildata_array_sorted[:, -2] <= final_stallangle + stalldelay_tip + aoa_margin))[0], :]
    airfoils_tip = airfoildata_array_sorted[np.where(airfoildata_array_sorted[:, -2] >= final_stallangle + stalldelay_tip - aoa_margin)[0], :]

    circ_over_t_req = circ_distr / t_sections  # Circulation per meter of span divided by thickness ratio

    # Contains for each section: Airfoil name, cl_(clcdmax), chord, t(verification), alpha_clcd, clcd_max, stall_margin
    wingdata = np.zeros((Nr_sections, 7))

    for section, spanpos in enumerate(span_sections):
        if spanpos <= rootsec_stall * b / 2:
            airfoils = airfoils_root
            min_aoa = -20
            max_aoa = final_stallangle
        elif spanpos <= midsec_stall * b / 2:
            airfoils = airfoils_mid
            min_aoa = final_stallangle
            max_aoa = final_stallangle + stalldelay_tip
        else:
            airfoils = airfoils_tip
            min_aoa = final_stallangle + stalldelay_tip
            max_aoa = 30

        print(airfoils)

        angle_error = (airfoils[:, -3] - circ_over_t_req[section]) * airfoils[:, -4] / (2 * np.pi)**2 * 180 # Required change in AoA to match required circulation
        
        stallmargin_new = airfoils[:, -2] + angle_error
        airfoils = np.hstack((airfoils, angle_error.reshape(-1, 1)))  # Add angle error to airfoils array

        airfoils = np.squeeze(airfoils[np.argwhere((stallmargin_new >= min_aoa) & (stallmargin_new <= max_aoa))])  # Find airfoils that are within the stall margin range

        assert airfoils.size != 0, "No available airfoils for the requested lift, thickness and stallmargin. Expand the range of thicknesses of your airfoils."

        # Add a column for expected Cl/Cd after reorientation (ClCd is modeled as a cone from max to 0 in the origin)
        airfoils = np.hstack((airfoils, (airfoils[:, -2] * (1 - np.abs(airfoils[:, -1])/airfoils[:, 0])).reshape(-1, 1))) # Cl/Cd_new = Cl/Cd_old * (1- error_angle/max_ld_angle)
        final_order = np.argsort(airfoils[:, -1])[::-1] # From highest to lowest L/D

        info = []
        ld_ratio = 0
        for i, argument in enumerate(final_order):
            code = str(airfoils[argument, -7])[1:5]
            cl = circ_over_t_req[section] * airfoils[argument, -6]
            a, cd = analyse(cl, code, Re)
            # info.append(f"Analysis completed, Cl/Cd = {cl/cd}, alpha_margin = {airfoils[argument, 3] - a}, expected value was {airfoils[argument, -4] + airfoils[argument, -2]}")
            alpha_margin = airfoils[argument, 3] - a
            info.append(f"Iteration: {i} - {code}, Cl/Cd = {cl/cd}, alpha_margin = {alpha_margin}, {cl}, {argument}")
            if min_aoa <= alpha_margin <= max_aoa:
                if cl/cd > ld_ratio:
                    ld_ratio = cl/cd
                    airfoil = argument
                    airfoils[argument, 0] = a
                    airfoils[argument, 1] = cl
                    airfoils[argument, 2] = cd
                    airfoils[argument, -4] = alpha_margin
                if i+1 >= min_analyses:
                    break
        # print(info)
        # print(f"Selected airfoil: {airfoils[argument]}, L/D = {ld_ratio}, alpha_margin = {alpha_margin}")

        # Put airfoil data in an array
        cl_cruise = circ_over_t_req[section] * airfoils[airfoil, 6]
        wingdata[section] = [airfoils[airfoil, 5], cl_cruise, t_sections[section]/airfoils[airfoil, 6],
                                t_sections[section], airfoils[airfoil, 0],
                                airfoils[airfoil, 1]/airfoils[airfoil, 2], airfoils[airfoil, 8]]

    dccl = np.diff(wingdata[:, 1]*wingdata[:, 2])
    dccl = np.concatenate((np.zeros((1)), average_array(dccl), np.array([dccl[-1]])))

    alpha_induced = []
    for i, span_pos in enumerate(span_sections):
        denominator = span_pos - np.append(-span_sections[::-1], span_sections[1:])
        integral = np.sum(dccl[i] / denominator[:np.argwhere(-denominator >= denominator[-1])[0, 0]])
        alpha_induced.append(1/ (8 * np.pi) * integral)

    # Add a column for induced angle of attack
    wingdata = np.hstack((wingdata, np.rad2deg(alpha_induced).reshape(-1, 1)))
    # Calculate CD and add it to the wingdata
    CD = wingdata[:, 1] * (1 / wingdata[:, 5] - np.sin(np.deg2rad(wingdata[:, -1])))
    wingdata = np.hstack((wingdata, CD.reshape(-1, 1)))

    area = np.sum(average_array(wingdata[:, 2])) * b/(Nr_sections-1)
    lift = 0.5 * rho * V**2 * np.sum(average_array(wingdata[:, 1] * wingdata[:, 2])) * b/(Nr_sections-1)
    drag = 0.5 * rho * V**2 * np.sum(average_array(wingdata[:, -1] * wingdata[:, 2])) * b/(Nr_sections-1)

    print(wingdata)
    print(f"Area: {area} m2, lift: {lift} N, drag: {drag} N, L/D: {lift/drag}")

    with open("wingdata.json", "w") as fp:
        json.dump({'design velocity':V, 'span':b, 'root thickness':t_r, 'mid_thickness':t_m, 'tip_thickness':t_t, 'tip_stalldelay':stalldelay_tip, 'wingdata':wingdata.tolist()}, fp)  # encode dict into JSON
    print("Saved wing data to wingdata.json!")

    return area, lift, drag, wingdata

if __name__ == "__main__":
    lift = 2200 * 9.81  # Maximum Takeoff Weight in Newtons
    V = 200 / 3.6  # Velocity in m/s
    b = 12 # Span in meters
    t_r = 0.24 # m
    t_m = 0.18 # Thickness at first engine in meters
    t_t = 0.08 # m
    stallmargin = 2
    # These values should result in an L/D of 89.0976735576122

    area, lift, drag, wingdata = design_wing(lift, V, b, t_r, t_m, t_t, stallmargin)
    print(wingdata)
    print(f"Area: {area} m2, lift: {lift} N, drag: {drag} N, L/D: {lift/drag}")

# Simplify code (assess all airfoil by off-design clcd)
# Run xfoil multiple times (multiple nr of points 160, 170 etc.) and take average
# Include more detailed ClCd approximation (third order polynomial?)