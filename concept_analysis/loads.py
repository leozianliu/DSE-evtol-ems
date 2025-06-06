import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import aeroloads as aero

class LoadsCalculator:
    def __init__(self, flight_mode, thrust, num_points):
        self.engine_positions_y = np.array([2.415, 5.029])  # [2.41, 4.884] y-locations from along wingbox (m) 4.884
        self.engine_offsets_x = np.array([1, 1])    # x-offsets from wingbox(m)
        self.engine_thrusts = np.array(thrust)  #  N cruise: [2800, 1400]   vertical: [7000.0, 4000.0] 
        self.engine_weights = np.array([100, 100]) * 9.81   # [50, 60] N
        self.half_span = 6.241                             # m  6.126    
        self.num_points = num_points
        self.flight_mode = flight_mode
        self.y_points = np.linspace(0, self.half_span, self.num_points)
        self.gull_location = 2.241                        # m 
        self.gull_angle = np.radians(-20)                            # degrees
        self.engine_offsets_z = np.array([-0.696, 0.447])    # z-offsets (m)
        self.fuselage_width = 1.8    /2 

    def thrust_loads(self):
        shear_t = np.zeros(self.num_points)
        torque_t = np.zeros(self.num_points)
        moment_t = np.zeros(self.num_points)

        for pos, thrust in zip(self.engine_positions_y, self.engine_thrusts):
            mask = self.y_points <= pos
            shear_t[mask] += -thrust
            
        idx_kink = np.argmin(np.abs(self.y_points - self.gull_location))
        torque_t[self.y_points <= self.engine_positions_y[1]] += self.engine_thrusts[1] * -self.engine_offsets_z[1]
        torque_t[self.y_points <= self.gull_location] *= math.cos(self.gull_angle)
        moment_due_torque_t = torque_t[idx_kink] * math.sin(self.gull_angle) 
        torque_t[self.y_points <= self.gull_location] += self.engine_thrusts[1] * (self.engine_positions_y[1] - self.gull_location) * math.sin(self.gull_angle)
        torque_t[self.y_points <= self.engine_positions_y[0]] += self.engine_thrusts[0] * -self.engine_offsets_z[0]

        # Assume the kink starts at index `kink_index`
        dy_before_kink = self.y_points[1] - self.y_points[0]            # Before kink
        dy_after_kink = dy_before_kink * np.cos(self.gull_angle)        # Adjusted dy after kink

        dy_array = np.where(
            np.arange(len(self.y_points) - 1) < idx_kink,
            dy_before_kink,
            dy_after_kink
        )
        moment_t[:-1] = np.flip(np.cumsum(np.flip(shear_t[:-1] * dy_array)))

        moment_t[self.y_points <= self.gull_location] += moment_due_torque_t

        self.shear_t = shear_t      
        self.torque_t = torque_t   
        self.moment_t = moment_t    
        
        return self.shear_t, self.torque_t, self.moment_t
    
    def engine_weight_loads(self):
        shear_we = np.zeros(self.num_points)
        torque_we = np.zeros(self.num_points)
        moment_we = np.zeros(self.num_points)
        normal_we = np.zeros(self.num_points)

        idx_kink = np.argmin(np.abs(self.y_points - self.gull_location))
        dy_before_kink = self.y_points[1] - self.y_points[0]            # Before kink
        dy_after_kink = dy_before_kink * np.cos(self.gull_angle)        # Adjusted dy after kink

        dy_array = np.where(
            np.arange(len(self.y_points) - 1) < idx_kink,
            dy_before_kink,
            dy_after_kink
        )

        for pos, offset_z, offset_x, weight in zip(self.engine_positions_y, self.engine_offsets_z, self.engine_offsets_x, self.engine_weights):
            mask = self.y_points <= pos        
            if self.flight_mode == 'horizontal':
                torque_we[mask] += -weight * offset_x
                shear_we[mask] += -weight

                # Rotating shear to reference frame after kink
                shear_we[self.y_points > self.gull_location] *= np.cos(self.gull_angle)
                normal_we[self.y_points > self.gull_location] += shear_we[self.y_points > self.gull_location] * math.tan(self.gull_angle)  
                moment_we[:-1] = np.flip(np.cumsum(np.flip(shear_we[1:] * dy_array))) * -1

            else:
                torque_we[mask] += -weight * offset_z
                shear_we[mask] += weight
                moment_we[:-1] = np.flip(np.cumsum(np.flip(shear_we[:-1] * dy_array)))    
        
        if self.flight_mode == 'vertical':
            torque_we[self.y_points < self.gull_location] += self.engine_weights[1] * (self.engine_positions_y[1] - self.gull_location) * math.sin(self.gull_angle) 

        self.shear_we = shear_we      
        self.torque_we = torque_we    
        self.moment_we = moment_we    
        self.normal_we = normal_we

        return self.shear_we, self.torque_we, self.moment_we, self.normal_we

    def aerodynamic_loads(self, lift, drag):
        shear_lift = np.zeros(self.num_points)
        moment_lift = np.zeros(self.num_points)
        shear_drag = np.zeros(self.num_points)
        moment_drag = np.zeros(self.num_points) 
        normal_lift = np.zeros(self.num_points)

        y = np.linspace(1.05, self.half_span, len(lift))

        lift_interp_func = interp1d(y, lift, bounds_error=False, fill_value='extrapolate')
        lift = lift_interp_func(self.y_points)
        drag_interp_func = interp1d(y, drag, bounds_error=False, fill_value='extrapolate')
        drag = drag_interp_func(self.y_points)

        #lift_value_at_105 = lift_interp_func(1.05)
        #drag_value_at_105 = drag_interp_func(1.05)
        lift[self.y_points <= self.fuselage_width] = 0 #lift_value_at_105
        drag[self.y_points <= self.fuselage_width] = 0 #drag_value_at_105

        dy = self.y_points[1] - self.y_points[0]
        shear_lift[:-1] = np.flip(np.cumsum(np.flip(lift[1:] * dy))) 
        shear_drag[:-1] = np.flip(np.cumsum(np.flip(drag[1:] * dy)))

        shear_lift[self.y_points > self.gull_location] = shear_lift[self.y_points > self.gull_location] / math.cos(self.gull_angle)
        shear_drag[self.y_points > self.gull_location] = shear_drag[self.y_points > self.gull_location]
        
        moment_lift[:-1] = np.flip(np.cumsum(np.flip(shear_lift[1:] * dy))) * -1
        moment_drag[:-1] = np.flip(np.cumsum(np.flip(shear_drag[1:] * dy))) 

        # Normal before kink from rotated lift after kink
        idx_kink = np.argmin(np.abs(self.y_points - self.gull_location))
        normal_lift[self.y_points  <= self.gull_location] = shear_lift[idx_kink] * math.sin(self.gull_angle) 

        self.shear_lift = shear_lift    
        self.shear_drag = shear_drag     
        self.moment_lift = moment_lift
        self.moment_drag = moment_drag  
        self.normal_lift = normal_lift      

        return self.shear_lift, self.shear_drag, self.moment_lift, self.moment_drag, self.normal_lift
    
    def weight_loads(self, MTOW):
        shear_weight = np.zeros(self.num_points)
        shear_weight[self.y_points < self.fuselage_width] = MTOW * 9.81 / 2 

        self.shear_weight = shear_weight

        return self.shear_weight
    
    def aero_moment(self):
        mask = self.y_points <= 5.1
        root_torque = 150 * math.ceil(5.1)
        aerodynamic_moment = np.zeros_like(self.y_points)
        aerodynamic_moment[mask] = np.linspace(root_torque, 150, mask.sum())
        self.aerodynamic_moment = aerodynamic_moment
        return self.aerodynamic_moment
    
    def combined_loads(self):
        dy = self.y_points[1] - self.y_points[0]
        if self.flight_mode == 'horizontal':
            shear_x = self.shear_t * -1 # Neglected: + self.shear_drag  
            shear_z = self.shear_weight + ((self.shear_lift + self.shear_we) * -1)
            moment_z = self.moment_t  # Neglected: + self.moment_drag
            moment_x = (self.moment_we + self.moment_lift) * -1
            torque = -self.torque_t + self.torque_we - self.aerodynamic_moment
            normal = - self.normal_lift + self.normal_we
        
        else:
            shear_x = (self.shear_t + self.shear_we) * -1 - self.shear_weight
            shear_z = np.zeros(self.num_points)
            moment_z = self.moment_we + self.moment_t 
            moment_x = np.zeros(self.num_points)
            torque =  self.torque_we - self.torque_t
            normal = np.zeros(self.num_points) 

        # fixing all degrees of freedom at fuselage interface
        shear_x[self.y_points < self.fuselage_width] = 0 
        shear_z[self.y_points < self.fuselage_width] = 0 
        moment_z[self.y_points < self.fuselage_width] = 0 
        moment_x[self.y_points < self.fuselage_width] = 0 
        torque[self.y_points < self.fuselage_width] = 0 
        normal[self.y_points < self.fuselage_width] = 0

        return shear_x, shear_z, moment_x, moment_z, torque, normal
    
calculator = LoadsCalculator('horizontal', [2800, 1400], 300)
calculator.thrust_loads()
calculator.engine_weight_loads()
calculator.weight_loads(2500)
calculator.aero_moment()
shear_lift, shear_drag, moment_lift, moment_drag, normal_lift = calculator.aerodynamic_loads(lift= 2.5 * aero.lift_gull_rh, drag= 2.5 * aero.drag_gull_rh)

shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()

if __name__ == "__main__":
    import matplotlib.pyplot as plt 

    # calculator = LoadsCalculator('horizontal', 300)
    # calculator.thrust_loads()
    # calculator.engine_weight_loads()
    # shear_lift, shear_drag, moment_lift, moment_drag, normal_lift = calculator.aerodynamic_loads(lift= aero.lift_gull_rh, drag= aero.drag_gull_rh)

    # shear_x, shear_z, moment_x, moment_z, torque, normal = calculator.combined_loads()

    # print("Shear X:\n", shear_x)
    # print("Shear Z:\n", shear_z)
    # print("Moment X:\n", moment_x)
    # print("Moment Z:\n", moment_z)
    # print("Torque:\n", torque)

    fig, axs = plt.subplots(2, 3, figsize=(18, 10))  # 2 rows, 3 columns

# Flatten for easier indexing
    axs = axs.flatten()

# Plot 1: Shear X
    axs[0].plot(calculator.y_points, shear_x)
    axs[0].set_title("Shear Force X-direction")
    axs[0].set_xlabel("Spanwise location y (m)")
    axs[0].set_ylabel("Shear X (N)")

# Plot 2: Shear Z
    axs[1].plot(calculator.y_points, shear_z)
    axs[1].set_title("Shear Force Z-direction")
    axs[1].set_xlabel("Spanwise location y (m)")
    axs[1].set_ylabel("Shear Z (N)")

# Plot 3: Moment X
    axs[2].plot(calculator.y_points, moment_x)
    axs[2].set_title("Bending Moment around X-axis")
    axs[2].set_xlabel("Spanwise location y (m)")
    axs[2].set_ylabel("Moment X (Nm)")

# Plot 4: Moment Z
    axs[3].plot(calculator.y_points, moment_z)
    axs[3].set_title("Bending Moment around Z-axis")
    axs[3].set_xlabel("Spanwise location y (m)")
    axs[3].set_ylabel("Moment Z (Nm)")

# Plot 5: Torque
    axs[4].plot(calculator.y_points, torque)
    axs[4].set_title("Torsional Moment (Torque)")
    axs[4].set_xlabel("Spanwise location y (m)")
    axs[4].set_ylabel("Torque (Nm)")

# Plot 6: Normal
    axs[5].plot(calculator.y_points, normal)
    axs[5].set_title("Normal")
    axs[5].set_xlabel("Spanwise location y (m)")
    axs[5].set_ylabel("Normal load (N)")

    plt.tight_layout()
    plt.show()
