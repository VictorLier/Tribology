import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Someya import Table1
import scipy.sparse as sps
from Someya import Table1

b = 50.8e-3
C = 102e-6 
eta_0 = 0.08
omega_b = 30 * 2 * np.pi
R_b = 50.902e-3

def analytical(nx, epsilon):
    phi = np.linspace(0, np.pi, nx)
    y = np.linspace(-b/2, b/2, nx)
    phi, y = np.meshgrid(phi, y)
    p = 3 * eta_0 * omega_b * epsilon / C**2 * (b**2/4 - y**2) * np.sin(phi)/(1+epsilon*np.cos(phi))**3
    return phi, y, p

def finite_difference(nx, epsilon):
    phi = np.linspace(0, np.pi, nx)
    y = np.linspace(-b/2, b/2, nx)
    dphi = abs(phi[0] - phi[1])
    dy = abs(y[0] - y[1])

    h = np.zeros((nx, nx))
    for i in range(nx):
        h[i,:] = C*(1+epsilon*np.cos(phi[i]))

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
            M[c, c] =   2 / dphi**2 + 2 / dy**2
            M[c, e] = - 1 / dphi**2 + ( 3 * epsilon * np.sin(phi[i]) / (1 + epsilon * np.cos(phi[i]))) * 1 / (2*dphi)
            M[c, w] = - 1 / dphi**2 - ( 3 * epsilon * np.sin(phi[i]) / (1 + epsilon * np.cos(phi[i]))) * 1 / (2*dphi)
            M[c, n] = - 1 / dy**2
            M[c, s] = - 1 / dy**2

            # Filling the rhs vector
            rhs[c] = 6 * omega_b * eta_0 * epsilon * np.sin(phi[i])/(C**2*(1 + epsilon*np.cos(phi[i]))**3)

    p2 = sps.linalg.spsolve(M, rhs)
    p2 = p2.reshape(nx, nx, order = 'F')
    phi, y = np.meshgrid(phi, y)
    return phi, y, p2

nx = 25
epsilon = 0.75

phi, y, p2 = finite_difference(nx,epsilon)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(phi*180/np.pi, y*1e3, p2*1e-6)
ax.set_xlabel('phi [deg]')
ax.set_ylabel('y [mm]')
ax.set_zlabel('p [MPa]')

phi, y, p = analytical(nx,epsilon)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(phi*180/np.pi, y*1e3, p*1e-6)
ax.set_xlabel('phi [deg]')
ax.set_ylabel('y [mm]')
ax.set_zlabel('p [MPa]')

plt.show()
plt.show()

# np.savetxt('W5/data/pressureDistributionShroudHeight.txt', np.array([X.flatten('F'), Y.flatten('F'), h.flatten('F')]).T)
# if printbol:
#     self.F2 = np.trapezoid(np.trapezoid(self.p2, self.X), self.Y)
#     print(f"The load capacity is: {self.F2:.3g} N")

# Q2

epsilon_table = Table1[:,1]
Phi_table = Table1[:,2]

Phi_a = []
w_r_a = []
S_a = []

Phi_n = []
w_r_n = []
S_n = []

for epsilon in Table1[:,1]:
    # analytical results
    w_x_a = eta_0 * omega_b * R_b * b**3 / (4 * C**2) * np.pi*epsilon / ((1-epsilon**2)**(3/2))
    w_z_a = eta_0 * omega_b * R_b * b**3 / (C**2)     * epsilon**2    / ((1-epsilon**2)**2) 
    w_r_a.append(np.sqrt(w_x_a**2 + w_z_a**2))
    Phi_a.append(np.arctan(w_x_a / w_z_a))
    S_a.append(eta_0 * omega_b/2*np.pi * b * 2*R_b / w_r_a[-1] * (R_b / C) ** 2)

    # numerical results
    nx = 25
    phi, y, p2 = finite_difference(nx,epsilon)
    dphi = abs(phi[0,0] - phi[0,1])
    dy = abs(y[0,0] - y[1,0])
    w_x_n = 0
    w_z_n = 0
    for i in range(nx): # loop over phi
        for j in range(nx): # loop over y
                w_x_n = w_x_n + p2[i,j] * R_b * np.sin(phi[i,j]) * dphi * dy
                w_z_n = w_z_n + p2[i,j] * R_b * np.cos(phi[i,j]) * dphi * dy
    w_r_n.append(np.sqrt(w_x_n**2 + w_z_n**2))
    Phi_n.append(np.arctan(w_x_n / w_z_n))
    S_n.append(eta_0 * omega_b/2*np.pi * b * 2*R_b / w_r_n[-1] * (R_b / C) ** 2)

# "Table of Sommerfeld numbers; epsilon, S, S_a, %error_a, S_n, %error_n"
table = np.zeros((len(S_a),6))
table[:,0] = epsilon_table
table[:,1] = Table1[:,0]
table[:,2] = S_a
table[:,3] = (S_a - Table1[:,0])/Table1[:,0]*100
table[:,4] = S_n
table[:,5] = (S_n - Table1[:,0])/Table1[:,0]*100

# "Table of epsilon and Phi; epsilon, Phi, Phi_a, %error_a, Phi_n, %error_n"
table2 = np.zeros((len(S_a),6))
table2[:,0] = epsilon_table
table2[:,1] = Phi_table
table2[:,2] = np.array(Phi_a)*180/np.pi
table2[:,3] = (np.array(Phi_a)*180/np.pi - Phi_table)/Phi_table*100
table2[:,4] = -180/np.pi*np.array(Phi_n)
table2[:,5] = (-180/np.pi*np.array(Phi_n) - Phi_table)/Phi_table*100

# print tables with 3 decimals and not in scientific notation
np.set_printoptions(precision=3, suppress=True)
print("Table of Sommerfeld numbers; epsilon, S, S_a, %error_a, S_n, %error_n")
print(table)
print(f"Mean error for S_a: {np.mean(table[:,3]):.3f} and for S_n: {np.mean(table[:,5]):.3f}")
print("Table of epsilon and Phi; epsilon, Phi, Phi_a, %error_a, Phi_n, %error_n")
print(table2)
print(f"mean error for Phi_a: {np.mean(table2[:,3]):.3f} and for Phi_n: {np.mean(table2[:,5]):.3f}")

