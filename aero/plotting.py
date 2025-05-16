import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df_NACA0015 = pd.read_csv('aero/NACA0015.csv')
df_NACA3412 = pd.read_csv('aero/NACA3412.csv')
df_S4233 = pd.read_csv('aero/S4233.csv')

x_NACA0015 = df_NACA0015.iloc[:,0]
y_NACA0015 = df_NACA0015.iloc[:,1]

x_NACA3412 = df_NACA3412.iloc[:,0]
y_NACA3412 = df_NACA3412.iloc[:,1]

x_S4233 = df_S4233.iloc[:,0]
y_S4233 = df_S4233.iloc[:,1]


plt.figure(figsize=(10,2))
plt.plot(x_NACA0015, y_NACA0015, label='NACA 0015')
plt.plot(x_NACA3412, y_NACA3412, label='NACA 3412')
plt.plot(x_S4233, y_S4233, label='S4233')
plt.axhline(y=0, color='black', linestyle='--')
plt.legend()
plt.show()