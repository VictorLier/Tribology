
# Problem 3

import numpy as np
import matplotlib.pyplot as plt

# Example data
T = np.array([27, 37.6, 49.0, 60.7, 76.4]) + 273  # Temperature in Kelvin
v = np.array([53.7, 31.4, 21.4, 13.9, 9.05])        # Kinematic viscosity in cSt

# Plot using log(T) and log(log(v))
plt.figure(1)
plt.clf()

plt.plot(np.log(T), np.log(np.log(v)), '-o')
plt.grid(True)

# Make y-axis
labels = np.array([5, 10, 15, 20, 30, 100, 200])
Ytick = np.log(np.log(labels))
plt.gca().set_yticks(Ytick)
plt.gca().set_yticklabels(labels)

# Make x-axis
xlabels = np.arange(-20, 101, 10) + 273  # Temperature labels from -20 to 100 with a step of 10
Xtick = np.log(xlabels)
plt.gca().set_xticks(Xtick)
plt.gca().set_xticklabels(xlabels)

# Set range of axis
plt.axis([np.log(293), np.log(353), np.log(np.log(8)), np.log(np.log(100))])

# Labels
plt.xlabel('Temperature [K]')
plt.ylabel('Kinematic viscosity [cSt = mm^2/s]')

# Show the plot
plt.show()

# 2 - fit the data

modify_nu = np.log(np.log(v+0.8))
modify_T = np.log(T)
a=np.polyfit(modify_T, modify_nu, 1)
m = a[0]
C = a[1]

print(f"Following the ASTM walther relation, the constants are: m = {m:.3f} and C = {C:.3f}")
print("-----")
m_simple = (np.log(np.log(v[4]+0.8)) - np.log(np.log(v[0]+0.8)))/(np.log(T[4])-np.log(T[0]))
C_simple = np.log(np.log(v[0]+0.8)) - m_simple*np.log(T[0])
print(f"Following the simple relation, the constants are: m = {m_simple:.3f} and C = {C_simple:.3f}")

# 3 - how much should the temperature increase to halft the kinematic viscosity

nu_30 = 30
nu_15 = nu_30/2
T_30 = np.exp((np.log(np.log(nu_30 + 0.8)) - C_simple)/m_simple)
T_15 = np.exp((np.log(np.log(nu_15 + 0.8)) - C_simple)/m_simple)

print(f"The temperature should increase by {T_15-T_30:.2f} K to halft the kinematic viscosity")