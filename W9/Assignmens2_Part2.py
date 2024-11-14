import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Someya import Table1
import scipy.sparse as sps
from Someya import Table1

# Constants
mass = 460          # Mass           [kg]
R_b = 50.902e-3     # Bearing radius [m]
D = 2*50.8e-3       # Shaft diameter [m]
Cp = 102e-6         # Clearance      [m]
W = mass*9.82/2     # Load           [N]
psi = Cp/R_b        # Eccentricity ratio
pf = 1e5            # Feed pressure  [N]

m_p = 0
LD_fac = 0.5
L = D * LD_fac      # Bearing length [m]
Cb = (1 - m_p) * Cp # Bearing clearance [m]
Qf = 2 * 0.15

# ISO VG 32
w_air = 1 # [m/s] air velocity of surroundings
alpha = 9.807 * (0.7 + 1.2 * w_air**(1/2))
A = 9 * D * D**(1/2)
rho = 876 # [kg/m^3]
cp = 1964 # [J/kg*K]
lamb = 1/3

nu_40 = 0.08/rho*1e6
t_40 = 40
nu_100 = 0.007/rho*1e6
t_100 = 100
m_lub = (np.log(np.log(nu_100 + 0.8)) - np.log(np.log(nu_40 + 0.8))) / (np.log((t_40+273.15)/(t_100+273.15)))
k_lub = np.log(np.log(nu_40 + 0.8)) + m_lub*np.log(t_40+273.15)
t_1 = 30
t_0 = 20
p_f = 2e5 # [Pa]
damping = 0.4
chi = 1 # !

def eta_i(temp):
    return rho*1e-6*(np.exp(np.exp(-m_lub*np.log(temp+273.15) + k_lub)) - 0.8)

def calculate_lub_temp(N):
    t = [50]  # Initial temperature
    omega = 2 * np.pi * N
    table = np.flip(Table1, axis=0)
    
    i = 0
    while True:
        eta = eta_i(t[i])
        S_current = eta * N * L * D / W * (R_b / Cp) ** 2
        Qs = np.interp(S_current, table[:, 0], table[:, 3])
        f_J = psi * np.interp(S_current, table[:, 0], table[:, 5])
        epsi_current = np.interp(S_current, table[:, 0], table[:, 1])

        h1 = Cp - epsi_current * Cb
        qf = 8 * h1 ** 3 / eta * pf * Qf
        q = R_b * omega * Cp * L * Qs + qf

        t_new = ((1 - lamb) * (f_J * R_b * W * omega + alpha * A * t_0) + cp * rho * q * t_1) / (cp * rho * q + alpha * A * (1 - lamb))
        t.append(t[i] + damping * (t_new - t[i]))

        if i > 50:
            t[-1] = np.mean(t[-2:])
            break

        if np.abs((t[i + 1] - t[i])/t[i]) < 1e-3:
            break
        
        i += 1
    
    lub_temp = t[-1]

    return lub_temp



b = 50.8e-3
C = 102e-6 
N = 100
eta_0 = eta_i(calculate_lub_temp(N))
print(eta_0)
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

nx_test = 40
epsilon_test = 0.3

phi, y, p2 = finite_difference(nx_test,epsilon_test)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(phi*180/np.pi, y*1e3, p2*1e-6, cmap='viridis', edgecolor='none')
ax.set_xlabel('$\phi$ [deg]')
ax.set_ylabel('y [mm]')
ax.set_zlabel('Pressure $p$ [MPa]')
ax.view_init(elev=30, azim=-120)  # Change the elevation and azimuth to adjust the orientation

phi, y, p = analytical(nx_test,epsilon_test)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(phi*180/np.pi, y*1e3, p*1e-6, cmap='viridis', edgecolor='none')
ax.set_xlabel('$\phi$ [deg]')
ax.set_ylabel('y [mm]')
ax.set_zlabel('Pressure $p$ [MPa]')
ax.view_init(elev=30, azim=-120)  # Change the elevation and azimuth to adjust the orientation

plt.show()
plt.show()

# calculating a mean sum squared error between the analytical and numerical solution
def squared_error(nx, epsilon):
    phi, y, p2 = finite_difference(nx, epsilon)
    phi, y, p = analytical(nx, epsilon)
    return np.mean((p2-p)**2)

print(f"Mean squared error between analytical and numerical solution: {squared_error(nx_test, epsilon_test):.3g}")


# np.savetxt('W5/data/pressureDistributionShroudHeight.txt', np.array([X.flatten('F'), Y.flatten('F'), h.flatten('F')]).T)
# if printbol:
#     self.F2 = np.trapezoid(np.trapezoid(self.p2, self.X), self.Y)
#     print(f"The load capacity is: {self.F2:.3g} N")

# Q2

Phi_a, S_a = [], []
Phi_n, S_n = [], []

Table1 = np.flip(Table1, axis=0)
epsilon_table = Table1[:,1]

for epsi in Table1[:,1]:
    # analytical results
    w_x_a = eta_0 * omega_b * R_b * b**3 / (4 * C**2) * np.pi*epsi / ((1-epsi**2)**(3/2))
    w_z_a = eta_0 * omega_b * R_b * b**3 / (C**2)     * epsi**2    / ((1-epsi**2)**2) 
    Phi_a.append(np.abs(np.arctan(w_x_a / w_z_a)*180/np.pi))
    S_a.append(np.interp(Phi_a[-1], Table1[:,2], Table1[:,0]))

    # numerical results
    nx_number = 25
    phi, y, p2 = finite_difference(nx_number,epsi)
    dphi = abs(phi[0,0] - phi[0,1])
    dy = abs(y[0,0] - y[1,0])
    w_x_n = 0
    w_z_n = 0
    for i in range(nx_number): # loop over phi
        for j in range(nx_number): # loop over y
                w_x_n = w_x_n + p2[i,j] * R_b * np.sin(phi[i,j]) * dphi * dy
                w_z_n = w_z_n + p2[i,j] * R_b * np.cos(phi[i,j]) * dphi * dy
    Phi_n.append(np.abs(np.arctan(w_x_n / w_z_n)*180/np.pi))
    S_n.append((np.interp(Phi_n[-1], Table1[:,2], Table1[:,0])))

# "Table of Sommerfeld numbers; epsilon, S, S_a, %error_a, S_n, %error_n"
table = np.zeros((len(S_a),5))
table[:,0] = Table1[:,0]
table[:,1] = S_a
table[:,2] = (S_a - Table1[:,0])/Table1[:,0]*100
table[:,3] = S_n
table[:,4] = (S_n - Table1[:,0])/Table1[:,0]*100

# print tables with 3 decimals and not in scientific notation
np.set_printoptions(precision=3, suppress=True)
print("Table of Sommerfeld numbers; epsilon, S, S_a, %error_a, S_n, %error_n")
print(table)
print(f"Mean error for S_a: {np.mean(table[:,2]):.3f} and for S_n: {np.mean(table[:,4]):.3f}")

plt.figure()
plt.plot(epsilon_table, S_a, label='Analytical')
plt.plot(epsilon_table, S_n, label='Numerical')
plt.plot(epsilon_table, Table1[:,0], label='Table1')
plt.xlabel('epsilon')
plt.ylabel('Sommerfeld number')
plt.legend()

plt.figure()
plt.plot(epsilon_table, np.array(Phi_a), label='Analytical')
plt.plot(epsilon_table, np.array(Phi_n), label='Numerical')
plt.plot(epsilon_table, Table1[:,2], label='Table1')
plt.xlabel('epsilon')
plt.ylabel('Phi')
plt.legend()
plt.show()

plt.show()

np.savetxt('ASSIGN2_AUST/Sa.txt', np.array([epsilon_table, S_a]).T)
np.savetxt('ASSIGN2_AUST/Sn.txt', np.array([epsilon_table, S_n]).T)
np.savetxt('ASSIGN2_AUST/Sbench.txt', np.array([epsilon_table, Table1[:,0]]).T)
np.savetxt('ASSIGN2_AUST/Phia.txt', np.array([epsilon_table, Phi_a]).T)
np.savetxt('ASSIGN2_AUST/Phin.txt', np.array([epsilon_table, Phi_n]).T)
np.savetxt('ASSIGN2_AUST/Phibench.txt', np.array([epsilon_table, Table1[:,2]]).T)