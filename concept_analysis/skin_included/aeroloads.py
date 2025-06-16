import numpy as np
import pandas as pd

gull_wing = pd.read_csv('gull_wing_new.csv')

lift_gull = gull_wing.iloc[:,0]
# drag_gull = gull_wing.iloc[:,1]
# lift_vd_gull = gull_wing.iloc[:,2]
# drag_vd_gull = gull_wing.iloc[:,3]

array_size = lift_gull.size
b = 12
x = np.linspace(-b/2, b/2, array_size)

n_max = 2.5

lift_gull_max = lift_gull*n_max
# lift_vd_gull_max = lift_vd_gull*n_max

# drag_gull_max = drag_gull*n_max
# drag_vd_gull_max = drag_vd_gull*n_max

#right portion for loading diagram
mask = x >= 1.05
x_rh = x[mask]
lift_gull_rh = lift_gull[mask]
# drag_gull_rh = drag_gull[mask]

lift_gull_max_rh = lift_gull_max[mask]
# drag_gull_max_rh = drag_gull_max[mask]

# print(lift_gull_rh)

