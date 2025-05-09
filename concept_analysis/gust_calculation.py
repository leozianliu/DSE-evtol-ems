import numpy as np
import matplotlib.pyplot as plt

#Euler integration to calculate required angular acceleration of vehicle to comply with gust requirement

#toggle plotting
plotting = True


#Vehicle parameters
m = 3200 #vehicle mass [kg]
S = 5 #reference area for drag [m^2]
CD0 = 0.5 #fuselage zero lift drag coefficient (because of assumed hover, induced lift is 0)

#time
dt = 0.001 #integration step size [s]
t = 0 #time [s]

#state
x = 0 #horizontal vehicle position [m]
x_d= 0 #horizontal vehicle velocity [m/s]
x_dd = 0 #horizontal vehicle acceleration [m/s/s]
theta = 0 #vehicle angle[rad]
theta_d = 0 #angular velocity [rad/s]
theta_dd = 0 #angular acceleration [rad/s/s]


#environment
g = 9.81 #earth acceleration [m/s/s]
rho = 1.225 # density [kg/m^3]
VG = 9.144 #gust velocity


#vehicle behavior
tspinup = 0.132 #time for props to change their velocity to provide desired thrust differential
maxAngularAccel = 0.235 #maximum angular acceleration [rad/s/s]
thold = 0  #time holding max angular accel. if too high: overshoot linear velocity, if too low: undershoot linear velocity, if extremely low: simulation becomes unstable

theta_steady = np.arctan((rho*np.power(VG, 2)*S*CD0)/(2*m*g))#steady

#set up arrays
time_arr = []
theta_arr = []
theta_d_arr = []
theta_dd_arr = []
x_arr = []
x_d_arr = []
x_dd_arr = []

time_arr.append(t)
theta_arr.append(theta)
theta_d_arr.append(theta_d)
theta_dd_arr.append(theta_dd)
x_d_arr.append(x_d)
x_dd_arr.append(x_dd)

x_arr.append(x)




while((theta < theta_steady) or (not np.isclose(np.abs(theta_d),0))):
    #determine angular acceleration from angular acceleration profile
    theta_ddd = maxAngularAccel/tspinup
    
    if(t<tspinup):
        theta_dd = theta_ddd*t
    elif((t>=tspinup) and (t<(tspinup + thold))):
        theta_dd = maxAngularAccel
    elif((t >= (tspinup + thold)) and (t < (3*tspinup + thold))):
        theta_dd = maxAngularAccel - theta_ddd * (t - tspinup - thold)
    elif((t >= (3*tspinup + thold)) and t < (3*tspinup + 2*thold)):
        theta_dd = -maxAngularAccel
    else:
        theta_dd = -maxAngularAccel + theta_ddd * (t - 3*tspinup - 2*thold)

    #perform forward euler on theta
    theta_d += theta_dd * dt
    theta += theta_d * dt

    #calculte linear acceleration
    x_dd = (rho*np.power((VG + x_d),2)*S*CD0)/(2*m) - g*np.tan(theta)
    
    #perform forward euler on position
    x_d += x_dd*dt
    x += x_d+dt

    #perform time step
    t += dt


    #append to arrays
    time_arr.append(t)
    theta_arr.append(theta)
    theta_d_arr.append(theta_d)
    theta_dd_arr.append(theta_dd)
    x_arr.append(x)
    x_d_arr.append(x_d)
    x_dd_arr.append(x_dd)

    print("theta: ", theta, " ; x: " , x)


    #break if sim is broken
    if(theta >= np.pi):
        break

print("Max x: ", np.max(x_arr))

if(plotting):
    plt.figure(figsize=(8, 8))

    # Plot 1: theta vs time
    plt.subplot(2, 3, 1)
    plt.plot(time_arr, theta_arr, label='θ (theta)', color='blue')
    plt.xlabel('Time')
    plt.ylabel('Theta (θ)')
    plt.title('Theta vs Time')
    plt.grid(True)
    plt.legend()

    # Plot 2: theta_d vs time
    plt.subplot(2, 3, 2)
    plt.plot(time_arr, theta_d_arr, label='theta_d', color='orange')
    plt.xlabel('Time')
    plt.ylabel('Theta_d')
    plt.title('Theta_d vs Time')
    plt.grid(True)
    plt.legend()

    # Plot 3: theta_dd vs time
    plt.subplot(2, 3, 3)
    plt.plot(time_arr, theta_dd_arr, label='theta_dd', color='green')
    plt.xlabel('Time')
    plt.ylabel('Theta_dd')
    plt.title('Theta_dd vs Time')
    plt.grid(True)
    plt.legend()

    # Plot 4: x vs time

    plt.subplot(2, 3, 4)
    plt.plot(time_arr, x_arr, label='x', color='red')
    plt.xlabel('Time')
    plt.ylabel('x')
    plt.title('x vs Time')
    plt.grid(True)
    plt.legend()

    # Plot 5: x_d vs time
    plt.subplot(2, 3, 5)
    plt.plot(time_arr, x_d_arr, label='x_d', color='red')
    plt.xlabel('Time')
    plt.ylabel('x_d')
    plt.title('x_d vs Time')
    plt.grid(True)
    plt.legend()

    # Plot 6: x_dd vs time

    plt.subplot(2, 3, 6)
    plt.plot(time_arr, x_dd_arr, label='x_dd', color='red')
    plt.xlabel('Time')
    plt.ylabel('x_d')
    plt.title('x_dd vs Time')
    plt.grid(True)
    plt.legend()

        

    plt.tight_layout()
    plt.show()
