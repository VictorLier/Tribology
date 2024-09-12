
# Problem 3

import math
import matplotlib.pyplot as plt
import numpy as np
from math import log

# variables

T = np.array([27, 37.6, 49.0, 60.7, 76]) # [degree celcius] - temperature
nu = np.array([53.7, 31.4, 21.4, 13.9, 9.05]) # [cSt] - kinematic viscosity

# 1 - plot the data

# Make double logarithmic on y axis
yticks = [np.exp(np.exp(i)) for i in range(0, 2)]
xticks = [np.exp(i) for i in range(3, 5)]

plt.plot(T, nu, 'o-')
plt.yticks(yticks)
plt.xticks(xticks)
plt.xlabel('Temperature [degree celcius]')
plt.ylabel('Kinematic viscosity [cSt]')
plt.title('Temperature vs Kinematic viscosity')
plt.grid()
plt.show()

# 2 - fit the data

modify = np.log(np.log(nu+0.8))
a=np.polyfit(T, modify, 1)
m = a[0]
C = a[1]

print(f"Following the ASTM walther relation, the constants are: m = {m:.3f} and C = {C:.3f}")

# 3 - how much should the temperature increase to halft the kinematic viscosity
