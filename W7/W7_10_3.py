
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# lambda = 2 rb / b 
# b / rb = 2 / lambda

lamb = np.array([1, 2, 4])
eps = np.linspace(0, 1, 100)

plt.figure()
for i in range(len(lamb)):
    upper = (2 + eps**2)*(16*eps**2 + np.pi**2*(1-eps**2))**(1/2)
    lower = 24*(1-eps**2)*(np.pi**2 - eps**2*(np.pi**2 - 4))**(1/2)
    term = (2/lamb[i])**2 * upper / lower
    plt.plot(eps, term, label="lambda = " + str(lamb[i]))

plt.legend()
plt.grid(True)
plt.xlabel('Eccentricity Ratio')
plt.ylabel('Term')
plt.show()


