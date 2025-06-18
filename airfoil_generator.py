from xfoil.xfoil import XFoil
from naca import naca
import numpy as np
from scipy.interpolate import PchipInterpolator
import json
import matplotlib.pyplot as plt

unittest = False
plot = False

def generate_data(nacacode, Re):
    """Generate data for a given airfoil and Reynolds number."""
    xf = XFoil()
    airfoil = naca(nacacode, 200, finite_TE = False, half_cosine_spacing = True, rotation = 0)
    xf.airfoil = airfoil
    xf.Re = Re
    xf.max_iter = 40
    a_raw, cl_raw, cd_raw = xf.aseq(-5, 20, 0.5)[0:3]
    a_cd = a_raw[~np.isnan(cd_raw)]
    cd_cleaned = cd_raw[~np.isnan(cd_raw)]
    a_cl = a_raw[~np.isnan(cl_raw)]
    cl_cleaned = cl_raw[~np.isnan(cl_raw)]
    cl_interp = PchipInterpolator(a_cl, cl_cleaned, extrapolate=False)
    cd_interp = PchipInterpolator(a_cd, cd_cleaned, extrapolate=False)
    a = np.arange(max(-5, min(a_cl), min(a_cd)), min(20.5, max(a_cl), max(a_cd)), 0.5)
    cl = cl_interp(a)
    cd = cd_interp(a)
    # Take the average of different nr of datapoint airfoils
    # Add a plot for the actual generated data
    # Do a polynomial fit on ld/alpha
    return a, cl, cd



camber_values = ['00', '14', '24', '34', '44', '54', '64']  # Example camber values
thickness_ratios = ['09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30']  # Example thickness ratios
V = 200/3.6  # Velocity in m/s
L = 1.0  # Characteristic length in m
Rho = 1.225  # Density of air in kg/m^3
Mu = 1.81e-5  # Dynamic viscosity of air in kg/(m·s)
Re = Rho * V * L / Mu
codeset = [f'{camber}{thickness}' for camber in camber_values for thickness in thickness_ratios]

with open("airfoildata.json", "r") as fp:
    # Load the dictionary from the file
    airfoildata = json.load(fp)

airfoildata_new = {code: [] for code in codeset if code not in airfoildata.keys()}

for nacacode in airfoildata_new.keys():
    a, cl, cd = generate_data(nacacode, 1e6)
    clcd_max = np.argmax(cl/cd)
    cl_max = np.argmax(cl)
    if 0 > a[clcd_max] > 15:
        print(f"Warning: Optimal angle of attack {a[clcd_max]} for NACA{nacacode} is outside expected range (0 to 15 degrees).")
    if a[cl_max] < 5:
        print(f"Warning: Alpha CL/CD max {a[clcd_max]} for NACA{nacacode} is outside expected range.")
    airfoildata_new[nacacode] = [a[clcd_max], cl[clcd_max], cd[clcd_max], a[cl_max], cl[cl_max]]
    if plot:
        fig, ax1 = plt.subplots()
        ax1.plot(a, cl, color='tab:red')
        ax1.plot(a[cl_max], cl[cl_max], 'o')
        ax1.tick_params(axis='y', labelcolor='tab:red')
        ax1.grid(True)
        ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis
        ax2.plot(a, cl/cd, color='tab:blue')
        ax2.plot(a[clcd_max], cl[clcd_max]/cd[clcd_max], 'o')
        ax2.tick_params(axis='y', labelcolor='tab:blue')
        plt.title('NACA ' + nacacode + ' at Re = ' + str(int(Re)))
        plt.show()

airfoildata.update(airfoildata_new)  # update the dictionary with new data
with open("airfoildata.json", "w") as fp:
    json.dump(airfoildata, fp)  # encode dict into JSON
print("Saved airfoil data to airfoildata.json!")

if __name__ == "__main__" and unittest:
    # Test if naca module is sensitive to second digit
    naca0412 = naca('0412', 160, finite_TE = False, half_cosine_spacing = True, rotation = 0)
    naca0012 = naca('0012', 160, finite_TE = False, half_cosine_spacing = True, rotation = 0)
    assert naca0412 == naca0012, "naca0412 and naca0012 should be equal"

    # Compare naca code generated airfoil to testdata airfoil
    xf = XFoil()
    xf.Re = 1e6
    xf.max_iter = 40
    
    naca0012_generated = naca('0012', 160, finite_TE = False, half_cosine_spacing = True, rotation = 0)
    xf.airfoil = naca0012_generated
    a_gen, cl_gen, cd_gen, cm_gen = xf.aseq(-20, 20, 0.5)[0:4]

    naca0012_test = ([+1.0000e+00, +9.9168e-01, +9.8037e-01, +9.6727e-01, +9.5272e-01, +9.3720e-01, +9.2112e-01, +9.0474e-01,
                    +8.8821e-01, +8.7160e-01, +8.5494e-01, +8.3827e-01, +8.2158e-01, +8.0488e-01, +7.8817e-01, +7.7146e-01,
                    +7.5475e-01, +7.3803e-01, +7.2132e-01, +7.0460e-01, +6.8789e-01, +6.7118e-01, +6.5447e-01, +6.3777e-01,
                    +6.2108e-01, +6.0440e-01, +5.8772e-01, +5.7106e-01, +5.5441e-01, +5.3778e-01, +5.2116e-01, +5.0456e-01,
                    +4.8798e-01, +4.7143e-01, +4.5489e-01, +4.3839e-01, +4.2191e-01, +4.0546e-01, +3.8905e-01, +3.7268e-01,
                    +3.5635e-01, +3.4007e-01, +3.2383e-01, +3.0766e-01, +2.9154e-01, +2.7550e-01, +2.5953e-01, +2.4366e-01,
                    +2.2788e-01, +2.1222e-01, +1.9670e-01, +1.8135e-01, +1.6619e-01, +1.5127e-01, +1.3666e-01, +1.2246e-01,
                    +1.0877e-01, +9.5752e-02, +8.3582e-02, +7.2423e-02, +6.2395e-02, +5.3537e-02, +4.5806e-02, +3.9101e-02,
                    +3.3294e-02, +2.8256e-02, +2.3867e-02, +2.0028e-02, +1.6658e-02, +1.3692e-02, +1.1080e-02, +8.7858e-03,
                    +6.7811e-03, +5.0484e-03, +3.5772e-03, +2.3627e-03, +1.4037e-03, +6.9909e-04, +2.4355e-04, +2.6000e-05,
                    +2.6000e-05, +2.4355e-04, +6.9909e-04, +1.4037e-03, +2.3627e-03, +3.5772e-03, +5.0484e-03, +6.7811e-03,
                    +8.7858e-03, +1.1080e-02, +1.3692e-02, +1.6658e-02, +2.0028e-02, +2.3867e-02, +2.8256e-02, +3.3294e-02,
                    +3.9101e-02, +4.5806e-02, +5.3537e-02, +6.2395e-02, +7.2423e-02, +8.3582e-02, +9.5752e-02, +1.0877e-01,
                    +1.2246e-01, +1.3666e-01, +1.5127e-01, +1.6619e-01, +1.8135e-01, +1.9670e-01, +2.1222e-01, +2.2788e-01,
                    +2.4366e-01, +2.5953e-01, +2.7550e-01, +2.9154e-01, +3.0766e-01, +3.2383e-01, +3.4007e-01, +3.5635e-01,
                    +3.7268e-01, +3.8905e-01, +4.0546e-01, +4.2191e-01, +4.3839e-01, +4.5489e-01, +4.7143e-01, +4.8798e-01,
                    +5.0456e-01, +5.2116e-01, +5.3778e-01, +5.5441e-01, +5.7106e-01, +5.8772e-01, +6.0440e-01, +6.2108e-01,
                    +6.3777e-01, +6.5447e-01, +6.7118e-01, +6.8789e-01, +7.0460e-01, +7.2132e-01, +7.3803e-01, +7.5475e-01,
                    +7.7146e-01, +7.8817e-01, +8.0488e-01, +8.2158e-01, +8.3827e-01, +8.5494e-01, +8.7160e-01, +8.8821e-01,
                    +9.0474e-01, +9.2112e-01, +9.3720e-01, +9.5272e-01, +9.6727e-01, +9.8037e-01, +9.9168e-01, +1.0000e+00]
                , [+1.2600e-03, +2.4215e-03, +3.9814e-03, +5.7619e-03, +7.7062e-03, +9.7433e-03, +1.1815e-02, +1.3886e-02,
                    +1.5936e-02, +1.7956e-02, +1.9943e-02, +2.1894e-02, +2.3810e-02, +2.5689e-02, +2.7532e-02, +2.9338e-02,
                    +3.1107e-02, +3.2839e-02, +3.4534e-02, +3.6190e-02, +3.7807e-02, +3.9384e-02, +4.0921e-02, +4.2415e-02,
                    +4.3865e-02, +4.5271e-02, +4.6630e-02, +4.7940e-02, +4.9199e-02, +5.0406e-02, +5.1557e-02, +5.2650e-02,
                    +5.3683e-02, +5.4652e-02, +5.5554e-02, +5.6385e-02, +5.7143e-02, +5.7822e-02, +5.8419e-02, +5.8930e-02,
                    +5.9349e-02, +5.9671e-02, +5.9891e-02, +6.0004e-02, +6.0002e-02, +5.9879e-02, +5.9628e-02, +5.9241e-02,
                    +5.8710e-02, +5.8027e-02, +5.7181e-02, +5.6164e-02, +5.4967e-02, +5.3580e-02, +5.1997e-02, +5.0216e-02,
                    +4.8243e-02, +4.6095e-02, +4.3805e-02, +4.1422e-02, +3.9000e-02, +3.6592e-02, +3.4237e-02, +3.1957e-02,
                    +2.9760e-02, +2.7644e-02, +2.5598e-02, +2.3613e-02, +2.1675e-02, +1.9770e-02, +1.7888e-02, +1.6017e-02,
                    +1.4147e-02, +1.2270e-02, +1.0381e-02, +8.4792e-03, +6.5676e-03, +4.6570e-03, +2.7615e-03, +9.0564e-04,
                    -9.0564e-04, -2.7615e-03, -4.6570e-03, -6.5676e-03, -8.4792e-03, -1.0381e-02, -1.2270e-02, -1.4147e-02,
                    -1.6017e-02, -1.7888e-02, -1.9770e-02, -2.1675e-02, -2.3613e-02, -2.5598e-02, -2.7644e-02, -2.9760e-02,
                    -3.1957e-02, -3.4237e-02, -3.6592e-02, -3.9000e-02, -4.1422e-02, -4.3805e-02, -4.6095e-02, -4.8243e-02,
                    -5.0216e-02, -5.1997e-02, -5.3580e-02, -5.4967e-02, -5.6164e-02, -5.7181e-02, -5.8027e-02, -5.8710e-02,
                    -5.9241e-02, -5.9628e-02, -5.9879e-02, -6.0002e-02, -6.0004e-02, -5.9891e-02, -5.9671e-02, -5.9349e-02,
                    -5.8930e-02, -5.8419e-02, -5.7822e-02, -5.7143e-02, -5.6385e-02, -5.5554e-02, -5.4652e-02, -5.3683e-02,
                    -5.2650e-02, -5.1557e-02, -5.0406e-02, -4.9199e-02, -4.7940e-02, -4.6630e-02, -4.5271e-02, -4.3865e-02,
                    -4.2415e-02, -4.0921e-02, -3.9384e-02, -3.7807e-02, -3.6190e-02, -3.4534e-02, -3.2839e-02, -3.1107e-02,
                    -2.9338e-02, -2.7532e-02, -2.5689e-02, -2.3810e-02, -2.1894e-02, -1.9943e-02, -1.7956e-02, -1.5936e-02,
                    -1.3886e-02, -1.1815e-02, -9.7433e-03, -7.7062e-03, -5.7619e-03, -3.9814e-03, -2.4215e-03, -1.2600e-03]
                )
    xf.airfoil = naca0012_test
    a_test, cl_test, cd_test, cm_test = xf.aseq(-20, 20, 0.5)[0:4]

    plt.plot(a_gen, cl_gen, color='blue', label='NACA 0012 from naca generator')
    plt.plot(a_test, cl_test, color='red', label='NACA 0012 from model test')
    plt.legend()
    plt.show()