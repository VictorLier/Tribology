import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Someya import Table1
import scipy.sparse as sps

b = 50.8e-3
c = 102e-6 
eta_0 = 0.08
omega_b = 30 * 2 * np.pi
epsilon = 0.1
e = c * epsilon

# plotting 3d plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
phi = np.linspace(0, np.pi, 100)
y = np.linspace(-b/2, b/2, 100)
phi, y = np.meshgrid(phi, y)
p = 3 * eta_0 * omega_b * epsilon / c**2 * (b**2/4 - y**2) * np.sin(phi)/(1+epsilon*np.cos(phi))**3
ax.plot_surface(phi, y, p)
ax.set_xlabel('phi')
ax.set_ylabel('y')
ax.set_zlabel('p')
plt.show()

# making finite difference method

nx = 25
dphi = abs(phi[0] - phi[1])
dy = abs(y[0] - y[1])

h = np.zeros((nx, nx))
for i in range(nx):
    h[i,:] = c*(1+epsilon*np.cos(phi[i]))

M = sps.eye(nx**2)
M = M.tocsr()
rhs = np.zeros(nx**2)

for i in range(1, nx-1):
    for j in range(1, nx-1):
        c = j + nx * (i)  # Current index
        n = c + 1  # North neighbor
        s = c - 1  # South neighbor
        e = j + nx * (i+1)  # East neighbor
        w = j + nx * (i - 1)  # West neighbor

        # Filling the M matrix
        M[c, c] = -2 * h[j, i]**3 * (1/(dx**2) + 1/(dy**2))
        M[c, n] = 3 * h.flatten('F')[c]**2 * (h.flatten('F')[n] - h.flatten('F')[s]) / (2 * dy**2) + h.flatten('F')[c]**3 / (dy**2)
        M[c, s] = -3 * h.flatten('F')[c]**2 * (h.flatten('F')[n] - h.flatten('F')[s]) / (2 * dy**2) + h.flatten('F')[c]**3 / (dy**2)
        M[c, e] = 3 * h.flatten('F')[c]**2 * (h.flatten('F')[e] - h.flatten('F')[w]) / (2 * dx**2) + h.flatten('F')[c]**3 / (dx**2)
        M[c, w] = -3 * h.flatten('F')[c]**2 * (h.flatten('F')[e] - h.flatten('F')[w]) / (2 * dx**2) + h.flatten('F')[c]**3 / (dx**2)

        

        # Filling the rhs vector
        rhs[c] = 6 * omega_b * eta_0 * (- e * np.sin(phi[c])/(c**3*(1-epsilon*np.cos(phi[c]))) )

p2 = sps.linalg.spsolve(M, rhs)

self.p2 = p2.reshape(nx, nx, order = 'F')

if plot:
    X, Y = np.meshgrid(self.X, self.Y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, self.p2)
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.plot_surface(X, Y, h)
    plt.show()
    plt.show()

    np.savetxt('W5/data/pressureDistributionShroudHeight.txt', np.array([X.flatten('F'), Y.flatten('F'), h.flatten('F')]).T)


if printbol:
    self.F2 = np.trapezoid(np.trapezoid(self.p2, self.X), self.Y)
    print(f"The load capacity is: {self.F2:.3g} N")











