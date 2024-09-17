
# Problem 3

import math
import matplotlib.pyplot as plt
import numpy as np
from math import log

# variables

T = np.array([27, 37.6, 49.0, 60.7, 76.4])+np.ones(5)*273.25 # [K] - temperature
nu = np.array([53.7, 31.4, 21.4, 13.9, 9.05]) # [cSt] - kinematic viscosity

# 1 - plot the data

# Make double logarithmic on y axis
yticks = [np.exp(np.exp(i)) for i in range(0, 2)]
xticks = [np.exp(i) for i in range(6, 7)]

plt.plot(T, nu, 'o-')
plt.yticks(yticks)
plt.xticks(xticks)
plt.xlabel('Temperature [degree celcius]')
plt.ylabel('Kinematic viscosity [cSt]')
plt.title('Temperature vs Kinematic viscosity')
plt.grid()
plt.show()

# 2 - fit the data

modify_nu = np.log(np.log(nu+0.8))
modify_T = np.log(T)
a=np.polyfit(modify_T, modify_nu, 1)
m = a[0]
C = a[1]

print(f"Following the ASTM walther relation, the constants are: m = {m:.3f} and C = {C:.3f}")
print("-----")
m_simple = (np.log(np.log(nu[4]+0.8)) - np.log(np.log(nu[0]+0.8)))/(np.log(T[4])-np.log(T[0]))
C_simple = np.log(np.log(nu[0]+0.8)) - m_simple*np.log(T[0])
print(f"Following the simple relation, the constants are: m = {m_simple:.3f} and C = {C_simple:.3f}")

# 3 - how much should the temperature increase to halft the kinematic viscosity

nu_30 = 30
nu_15 = nu_30/2
T_30 = np.exp((np.log(np.log(nu_30 + 0.8)) - C_simple)/m_simple)
T_15 = np.exp((np.log(np.log(nu_15 + 0.8)) - C_simple)/m_simple)

print(f"The temperature should increase by {T_15-T_30:.2f} K to halft the kinematic viscosity")