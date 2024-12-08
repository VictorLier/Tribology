import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import cmath
from scipy.optimize import fsolve
import scipy.sparse as sps

# Define the variables
d_i = 69.5e-3    # [m] - Inner race diameter
d_o = 85.5e-3    # [m] - Outer race diameter
d   = 8e-3       # [m] - Diameter of cylindrical rollers
l   = 8e-3       # [m] - Axial length of cylindrical rollers
n   = 20         # [-] - Numer of rollers
w_z = 8e3        # [N] - Maximum radial load
omega_i = 837    # [rad/s] - Inner race angular speed
omega_o = 0      # [rad/s] - Outer race angular speed

E   = 2.1e11     # [Pa] - Young's modulus
nu  = 0.3        # [-] - Poisson's ratio

eta_0 = 0.01     # [Pa*s] - Viscosity of lubricant
xi  = 2e-8       # [1/Pa/s] - Pressure-viscosity coefficient

print("-------------------------------------------")
# --------------------- Question 1 ------------------------------------

# 1) ------------------------------------------------------------------
# w_z = n * w_z_max / 4   # Eq. 21.47 =>
w_z_max = 4 * w_z / n

print(f"Q1.1 - Load of the roller with the highest load: {w_z_max*1e-3:.2f} kN")
print("-------------------------------------------")

# 2) ------------------------------------------------------------------

Lambda = 3
Rqa = 0.2e-6  * np.pi/(2*np.sqrt(2))
Rqb = 0.05e-6 * np.pi/(2*np.sqrt(2))

h_min = Lambda * (np.sqrt(Rqa**2 + Rqb**2))
print(f"Q1.2 - Minimum film thickness: {h_min*1e6:.2f} um")
print("-------------------------------------------")


# --------------------- Question 2 ------------------------------------

# 1) ------------------------------------------------------------------

r_ax = np.array([d/2, d/2])
r_bx = np.array([d_i/2, -d_o/2])
R_x = ( 1 / r_ax + 1 / r_bx )**(-1) # Eq. 17.4
print(f"Q2.1 - Equivalent radius for inner: {R_x[0]*1e3:.2f} mm")
print(f"Q2.1 - Equivalent radius for outer: {R_x[1]*1e3:.2f} mm")

w_z_mark = w_z_max / l
E_mark = E / (1 - nu**2) # Eq. 17.17
W_mark = w_z_mark / (E_mark * R_x) # Eq. 17.38

p_m = E_mark * np.sqrt(W_mark/(2*np.pi)) # Eq. 17.40
print(f"Q2.1 - Maximum contact pressure for inner: {p_m[0]*1e-9:.2f} GPa")
print(f"Q2.1 - Maximum contact pressure for outer: {p_m[1]*1e-9:.2f} GPa")

delta_m = 2*W_mark*R_x / np.pi * (np.log(2*np.pi/W_mark) - 1) # Eq. 17.39
print(f"Q2.1 - Maximum deflection for inner: {delta_m[0]*1e6:.2f} um")
print(f"Q2.1 - Maximum deflection for outer: {delta_m[1]*1e6:.2f} um")
print("-------------------------------------------")

# 2) ------------------------------------------------------------------

# The deflection is of the same order of magnitude as the minimum film thickness

# 3) ------------------------------------------------------------------

Nx = 200

b = R_x * np.sqrt(8*W_mark/np.pi) # Eq. 17.37

# -------- For the inner ------------

D_x = 2*b[0] # ?
x = np.linspace(-2*D_x, 2*D_x, Nx)
p = p_m[0] * np.array([cmath.sqrt(1 - (2 * xi / D_x) ** 2) for xi in x]).real # Eq. 17.6, for Dy -> infinity

P = np.real(p / p_m[0]) # Eq. 18.23
delta_dash = 0
delta_dash = np.zeros(Nx)
X = np.linspace(-4, 4, Nx)
Delta = X[1] - X[0]
for i in range(1, Nx-1):
    for j in range(1, Nx-1):
        term1 = np.abs( (X[i+1] + X[i])/2 - X[j] )
        term2 = np.abs( (X[i-1] + X[i])/2 - X[j] )
        delta_dash[i] = delta_dash[i] - Delta/(2*np.pi)*P[j]*np.log( term1 * term2 )   # Eq. 18.31

S = (2*x)**2 / (2 *R_x[0]) # Eq. 12.52
delta = (D_x)**2 * delta_dash / R_x[0] # Eq. 18.23

# plotting
plt.figure()
plt.plot(x*1e6, p*1e-9)
plt.title("Pressure distribution")
plt.xlabel("x [um]")
plt.ylabel("p [GPa]")
plt.grid()

plt.figure()
plt.plot(x*1e6, delta*1e6)
plt.title("Deformation field")
plt.xlabel("x [um]")
plt.ylabel("delta' [um]")
plt.grid()

plt.figure()
plt.plot(x*1e6, S*1e6, label="Undeformed roller", linestyle='--')
plt.plot(x*1e6, (S + delta)*1e6, label = "Deformed roller", linestyle='-')
plt.xlabel("x [um]")
plt.ylabel("y [um]")
plt.legend()
plt.grid()

np.savetxt("Report_3/data/Q2_Si.txt", np.array([x*1e6, S*1e6]).T)
np.savetxt("Report_3/data/Q2_defi.txt", np.array([x*1e6, (S + delta)*1e6]).T)
np.savetxt("Report_3/data/Q2_pi.txt", np.array([x*1e6, p*1e-9]).T)

# --------- for the outer ------------

D_x = 2*b[1] # ?
x = np.linspace(-2*D_x, 2*D_x, Nx)
p = p_m[1] * np.array([cmath.sqrt(1 - (2 * xi / D_x) ** 2) for xi in x]).real # Eq. 17.6, for Dy -> infinity
P = np.real(p / p_m[1]) # Eq. 18.23
delta_dash = 0
delta_dash = np.zeros(Nx)
for i in range(1, Nx-1):
    for j in range(1, Nx-1):
        term1 = np.abs( (X[i+1] + X[i])/2 - X[j] )
        term2 = np.abs( (X[i-1] + X[i])/2 - X[j] )
        delta_dash[i] = delta_dash[i] - Delta/(2*np.pi)*P[j]*np.log( term1 * term2 )   # Eq. 18.31

S = (2*x)**2 / (2 *R_x[1]) # Eq. 12.52
delta = (D_x)**2 * delta_dash / R_x[1] # Eq. 18.23

# plotting
plt.figure()
plt.plot(x*1e6, p*1e-9)
plt.title("Pressure distribution")
plt.xlabel("x [um]")
plt.ylabel("p [GPa]")
plt.grid()

plt.figure()
plt.plot(x*1e6, delta*1e6)
plt.title("Deformation field")
plt.xlabel("x [um]")
plt.ylabel("delta' [um]")
plt.grid()

plt.figure()
plt.plot(x*1e6, S*1e6, label="Undeformed roller", linestyle='--')
plt.plot(x*1e6, (S + delta)*1e6, label = "Deformed roller", linestyle='-')
plt.xlabel("x [um]")
plt.ylabel("y [um]")
plt.legend()
plt.grid()

plt.show()

np.savetxt("Report_3/data/Q2_So.txt", np.array([x*1e6, S*1e6]).T)
np.savetxt("Report_3/data/Q2_defo.txt", np.array([x*1e6, (S + delta)*1e6]).T)
np.savetxt("Report_3/data/Q2_po.txt", np.array([x*1e6, p*1e-9]).T)

# 4) ------------------------------------------------------------------

# ...

# 5) ------------------------------------------------------------------

d_o_c = 85.51e-3    # [m] - Outer race diameter
c_d_c = d_o_c - d_i - 2*d # Eq. 21.2
print(f"Q2.5 - Clearance: {c_d_c*1e6:.2f} um")

r_bx_c = np.array([d_i/2, -d_o_c/2])
R_x_c = ( 1 / r_ax + 1 / r_bx_c )**(-1) # Eq. 17.4
# W_mark = w_z_mark / (E_mark * R_x) # Eq. 17.38
# delta_m = 2*W_mark*R_x / np.pi * (np.log(2*np.pi/W_mark) - 1) # Eq. 17.39

# Z_w_top = np.pi * ( 1 - c_d / (2*delta_m) )**(3/2)
# Z_w_bottom = 2.491 * ( np.sqrt( 1 + (( 1 - c_d / (2*delta_m) ) / 1.23 )**2) - 1) 
# Z_w = Z_w_top / Z_w_bottom # Eq. 21.46

# # w_z = n * w_z_max_clearance / Z_w # Eq. 21.44
# w_z_max_clearance = w_z_max * Z_w / n

delta_m_c = [12e-6] # initial guess
w_z_max_c = []

while True:
    psi = np.arccos(c_d_c / (2 * delta_m_c[-1])) # Eq. 21.41
    w_z_max_c.append(2 * np.pi * w_z * (1 - (c_d_c / (2 * delta_m_c[-1]))) / ((psi - (c_d_c / (2 * delta_m_c[-1])) * np.sin(psi)) * n)) # Eq. 21.44 

    W_mark_i = (w_z_max_c[-1] / l) / (E_mark * R_x_c[0])  # Eq. 17.38
    W_mark_o = (w_z_max_c[-1] / l) / (E_mark * R_x_c[1])  # Eq. 17.38
    
    delta_m_c_i = 2 * W_mark_i * R_x_c[0] / np.pi * (np.log((2 * np.pi) / W_mark_i) - 1) # Eq. 17.39
    delta_m_c_o = 2 * W_mark_o * R_x_c[1] / np.pi * (np.log((2 * np.pi) / W_mark_o) - 1) # Eq. 17.39
    delta_m_c.append(delta_m_c_o + delta_m_c_i) # Eq. 21.36

    if np.abs((delta_m_c[-2] - delta_m_c[-1] ) / delta_m_c[-2]) < 1e-6: # convergence criteria
        psi = np.arccos(c_d_c / (2 * delta_m_c[-1])) # Eq. 21.41
        w_z_max_c.append(2 * np.pi * w_z * (1 - (c_d_c / (2 * delta_m_c[-1]))) / ((psi - (c_d_c / (2 * delta_m_c[-1])) * np.sin(psi)) * n)) # Eq. 21.44
        break

plt.figure()
plt.plot(w_z_max_c)
plt.title("Load with clearance")
plt.xlabel("Iterations")
plt.ylabel("Load [N]")
plt.grid()

plt.figure()
plt.plot(delta_m_c)
plt.title("Deflection with clearance")
plt.xlabel("Iterations")
plt.ylabel("Deflection [m]")
plt.grid()

plt.show()

print(f"Q2.5 - Load with clearance: {w_z_max_c[-1]*1e-3:.2f} kN")
print(f"Q2.5 - total deflection with clearance: {delta_m_c[-1]*1e6:.2f} um")
print("-------------------------------------------")


# --------------------- Question 3 ------------------------------------

# At first consider a case where the contact is assumed to be between rigid elements (cf.
# Chapter 16) and consider a hydrodynamic solution of the problem.

# 1) ------------------------------------------------------------------
# Using both the short and infinite solutions find the corresponding h0 that balances
# the external load. Is the short/long bearing solution representative?

omega_o = 0

# w_za_mark = 2.44*eta_0*(u_b+u_a)*R/h_0 # eq.xxx
r = d/2
d_e = (d_o+d_i)/2
u_a=abs(omega_o-omega_i)*d_e/4*(1-d**2/d_e**2)
u_b=abs(omega_i-omega_o)*d_e/4*(1-d**2/d_e**2)

print(f"Q3.1 - u_a: {u_a:.2f} m/s") 
print(f"Q3.1 - u_b: {u_b:.2f} m/s")

# long bearing solution
h_0_long = 2.44 * eta_0 * (u_a+u_b) * r / (w_z_max/l) # Eq. 16.36

# short bearing solution
h_0_short = np.sqrt((eta_0*(u_a+u_b) * l**3) / (4 * w_z_max)) # Eq. 16.35

print(f"Q3.1 - Long bearing solution: {h_0_long*1e6:.4f} um")
print(f"Q3.1 - Short bearing solution: {h_0_short*1e6:.4f} um")
print("-------------------------------------------")

# 2) ------------------------------------------------------------------
# Adapt your FD model (e.g. using (16.39) or (16.45)) that represents the finite-width
# HL solution and discuss pressure level and film thickness.

lambda_j = l / (2*r)

def finite_difference(nx, h_0):
    X = np.linspace(0, 1, nx)
    Y = np.linspace(-1, 1, nx)
    dX = abs(X[0] - X[1])
    dY = abs(Y[0] - Y[1])

    H = np.zeros((nx, nx))
    dH = np.zeros((nx, nx))
    ddH = np.zeros((nx, nx))
    P = np.zeros((nx, nx))

    for i in range(nx):
        H[i,:] = 1 + r/(2*h_0) * (X[i] - 1)**2
        dH[i,:] = r/h_0 * (X[i] - 1)
        ddH[i,:] = (np.sqrt(2)*(r*(X[i]-1)**2+h_0)*r)/(h_0**2*np.sqrt(r/h_0*(X[i]-1)**2+2))

    H, dH, ddH = H.flatten('F'), dH.flatten('F'), ddH.flatten('F')

    M = sps.eye(nx**2)
    M = M.tocsr()
    rhs = np.zeros(nx**2)
    for i in range(1, nx-1):
        for j in range(1, nx-1):
            c = j + nx * (i)     # Current index
            n = c + 1            # North neighbor
            s = c - 1            # South neighbor
            e = j + nx * (i+1)   # East neighbor
            w = j + nx * (i - 1) # West neighbor

            # Filling the M matrix
            M[c, c] = - 2 / dX**2 - 1 / lambda_j**2 * 2 / dY**2 - 3 / 2 * 1 / H[c]**(3/2) * ddH[c]
            M[c, e] =   1 / dX**2
            M[c, w] =   1 / dX**2
            M[c, n] =               1 / lambda_j**2 * 1 / dY**2
            M[c, s] =               1 / lambda_j**2 * 1 / dY**2

            # Filling the rhs vector
            rhs[c] = 1 / H[c]**(3/2) * dH[c]

    Gamma = sps.linalg.spsolve(M, rhs)
    P = Gamma / H**(3/2)
    # Gamma = Gamma.reshape(nx, nx, order = 'F')
    P = P.reshape(nx, nx, order = 'F')
    H = H.reshape(nx, nx, order = 'F')
    X, Y = np.meshgrid(X, Y)
    return X, Y, P, H

nx_test = 10

h_0_test = np.array([1])*h_0_long

for h_0_t in h_0_test:

    X, Y, P, H = finite_difference(nx_test, h_0_t)
    x = X*r - r
    y = l*Y/2
    p = 6*eta_0*(u_a+u_b)*r*P/h_0_t**2 # Eq. 16.38
    W_tot = np.trapezoid(np.trapezoid(p, x[0, :], axis=0), y[:, 0], axis=0)
    print(f"Q3.2 - Total load for h_0 = {h_0_t*1e6:.4f} um: {W_tot:.2f} N")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x*1e3, y*1e3, p*1e-6, cmap='viridis', edgecolor='none')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('p [MPa]')
    ax.view_init(elev=30, azim=-120)  # Change the elevation and azimuth to adjust the orientation

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x*1e3, y*1e3, p*1e-6, cmap='viridis', edgecolor='none')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('p [MPa]')
    ax.set_ylim(3, 4)
    ax.view_init(elev=30, azim=-120)  # Change the elevation and azimuth to adjust the orientation

    plt.show()

    # Exercise 3.3

    p_p = - 1 / xi * np.log( 1 - xi * p) # Eq. 18.6

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x*1e3, y*1e3, p_p*1e-6, cmap='viridis', edgecolor='none')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('p_p [MPa]')
    ax.view_init(elev=30, azim=-120)  # Change the elevation and azimuth to adjust the orientation

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x*1e3, y*1e3, p_p*1e-6, cmap='viridis', edgecolor='none')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('p_p [MPa]')
    ax.set_ylim(3, 4)
    ax.view_init(elev=30, azim=-120)  # Change the elevation and azimuth to adjust the orientation

    eta_eta0 = np.exp(xi * p) # Eq. 18.3

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x*1e3, y*1e3, eta_eta0, cmap='viridis', edgecolor='none')
    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Y [mm]')
    ax.set_zlabel('$\eta_p / \eta_0$ [-]')
    ax.view_init(elev=30, azim=-120)  # Change the elevation and azimuth to adjust the orientation

    plt.show()


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# surf = ax.plot_surface(x*1e3, y*1e3, H*h_0_short*1e3, cmap='viridis', edgecolor='none')
# ax.set_xlabel('X [mm]')
# ax.set_ylabel('Y [mm]')
# ax.set_zlabel('h [mm]')
# ax.view_init(elev=30, azim=-120)  # Change the elevation and azimuth to adjust the orientation

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# surf = ax.plot_surface(X, Y, P, cmap='viridis', edgecolor='none')
# ax.set_xlabel('X [-]')
# ax.set_ylabel('Y [-]')
# ax.set_zlabel('P [-]')
# ax.view_init(elev=30, azim=-120)  # Change the elevation and azimuth to adjust the orientation


# 3) ------------------------------------------------------------------
# With the pressure levels experienced, would it have an influence on the viscosity?
# Could a transformation in the style used in (18.2)-(18.6) be applied and present
# useful results?




# --------------------- Question 4 ------------------------------------
# We now consider the EHL solution to the roller-race contact problem as presented in
# Chapter 18

# 1) ------------------------------------------------------------------
# For the EHL solution, find the minimum and central film thicknesses and determine
# the peak pressure using the interpolation relations given by Hamrock.

U = eta_0 * (u_a + u_b)/2 / ( E_mark * R_x )
G = xi * E_mark

# Minimum film thickness 
h_e_min_dash_inner = 1.806*(w_z/l)**(-0.128) * (eta_0*omega_i*d_i/2)**0.694 * xi**0.568 * R_x[0]**0.434  # Eq. 18.73
h_e_min_dash_outer = 1.806*(w_z/l)**(-0.128) * (eta_0*omega_i*d_i/2)**0.694 * xi**0.568 * R_x[1]**0.434  # Eq. 18.73

# Center film thickness
H_ec_dash_inner = 2.922*W_mark[0]**(-0.166) * U[0]**0.692 * G**0.470 # Eq. 18.74
h_c_dash_inner = H_ec_dash_inner * R_x[0]
H_ec_dash_outer = 2.922*W_mark[1]**(-0.166) * U[1]**0.692 * G**0.470 # Eq. 18.74
h_c_dash_outer = H_ec_dash_outer * R_x[1]

# Peak pressure
P_es_dash_inner = 0.648*W_mark[0]**0.185 * U[0]**0.275 * G**0.391 # Eq. 18.70
p_sk_dash_inner = P_es_dash_inner * E_mark 
P_es_dash_outer = 0.648*W_mark[1]**0.185 * U[1]**0.275 * G**0.391 # Eq. 18.70
p_sk_dash_outer = P_es_dash_outer * E_mark
# Pressure spike location
X_es_dash_inner = 1.111*W_mark[0]**0.606 * U[0]**(-0.021) * G**0.077 # Eq. 18.71
X_es_dash_outer = 1.111*W_mark[1]**0.606 * U[1]**(-0.021) * G**0.077 # Eq. 18.71
x_dash_inner = X_es_dash_inner * R_x[0]
x_dash_outer = X_es_dash_outer * R_x[1]

# Center of pressure
X_ecp_dash_inner = -3.595*W_mark[0]**(-1.019) * U[0]**0.638 * G**(-0.358) # Eq. 18.76
X_ecp_dash_outer = -3.595*W_mark[1]**(-1.019) * U[1]**0.638 * G**(-0.358) # Eq. 18.76
x_cp_dash_inner = X_ecp_dash_inner * R_x[0]
x_cp_dash_outer = X_ecp_dash_outer * R_x[1]

# Minimum film thickness
H_e_min_dash_inner = 1.714*W_mark[0]**(-0.128) * U[0]**0.694 * G**0.568 # Eq. 18.72
h_e_min_dash_inner = H_e_min_dash_inner * R_x[0]
H_e_min_dash_outer = 1.714*W_mark[1]**(-0.128) * U[1]**0.694 * G**0.568 # Eq. 18.72
h_e_min_dash_outer = H_e_min_dash_outer * R_x[1]
# Minimum film thickness location
X_e_min_dash_inner = 1.439*W_mark[0]**0.548 * U[0]**(-0.011) * G**0.026 # Eq. 18.75
X_e_min_dash_outer = 1.439*W_mark[1]**0.548 * U[1]**(-0.011) * G**0.026 # Eq. 18.75
x_e_min_dash_inner = X_e_min_dash_inner * R_x[0]
x_e_min_dash_outer = X_e_min_dash_outer * R_x[1]

print(f"Q4.1 - Dimensionless U: {U}")
print(f"Q4.1 - Dimensionless G: {G}")
print(f"Q4.1 - Dimensionless W_mark: {W_mark}")

print(f"Q4.1 - Minimum film thickness inner: {h_e_min_dash_inner*1e6:.2f} um")
print(f"Q4.1 - Minimum film thickness outer: {h_e_min_dash_outer*1e6:.2f} um")
print(f"Q4.1 - Central film thickness inner: {h_c_dash_inner*1e6:.2f} um")
print(f"Q4.1 - Central film thickness outer: {h_c_dash_outer*1e6:.2f} um")
print(f"Q4.1 - Peak pressure inner: {p_sk_dash_inner*1e-9:.2f} GPa")
print(f"Q4.1 - Peak pressure outer: {p_sk_dash_outer*1e-9:.2f} GPa")
print(f"Q4.1 - Pressure spike location inner: {x_dash_inner*1e6:.2f} um")
print(f"Q4.1 - Pressure spike location outer: {x_dash_outer*1e6:.2f} um")
print(f"Q4.1 - Center of pressure inner: {x_cp_dash_inner*1e6:.2f} um")
print(f"Q4.1 - Center of pressure outer: {x_cp_dash_outer*1e6:.2f} um")
print(f"Q4.1 - Minimum film thickness inner: {h_e_min_dash_inner*1e6:.2f} um")
print(f"Q4.1 - Minimum film thickness outer: {h_e_min_dash_outer*1e6:.2f} um")
print(f"Q4.1 - Minimum film thickness location inner: {x_e_min_dash_inner*1e6:.2f} um")
print(f"Q4.1 - Minimum film thickness location outer: {x_e_min_dash_outer*1e6:.2f} um")
print("-------------------------------------------")

# 2) ------------------------------------------------------------------
# Sketch the pressure distribution and film thickness in the contact zone and mark
# the important values. See Fig. 18.8 for inspiration. Compare results for quater,
# half and full load. Does it act as expected from Fig. 18.11?

D_x = 2*b
x_inner = np.linspace(-2*D_x[0], 2*D_x[0], Nx)
x_outer = np.linspace(-2*D_x[1], 2*D_x[1], Nx)
p_inner = p_m[0] * np.array([cmath.sqrt(1 - (2 * xi / D_x[0]) ** 2) for xi in x_inner]).real # Eq. 17.6, for Dy -> infinity
P_inner = np.real(p_inner / p_m[0]) # Eq. 18.23
p_outer = p_m[1] * np.array([cmath.sqrt(1 - (2 * xi / D_x[1]) ** 2) for xi in x_outer]).real # Eq. 17.6, for Dy -> infinity
P_outer = np.real(p_outer / p_m[1]) # Eq. 18.23

# Defining point of interest - for pressure and film thickness

# mimimum film thickness points
P_h0_min_inner = [2*x_e_min_dash_inner/D_x[0], h_e_min_dash_inner/h_min]
P_h0_min_outer = [2*x_e_min_dash_outer/D_x[1], h_e_min_dash_outer/h_min]

# central film thickness points
P_h0_c_inner = [2*x_cp_dash_inner/D_x[0], h_c_dash_inner/h_min]
P_h0_c_outer = [2*x_cp_dash_outer/D_x[1], h_c_dash_outer/h_min]

# spike pressure points
P_p_es_inner = [2*x_dash_inner/D_x[0], p_sk_dash_inner/p_m[0]]
P_p_es_outer = [2*x_dash_outer/D_x[1], p_sk_dash_outer/p_m[1]]

# peak pressure (center of pressure) points
P_p_cp_inner = [2*x_cp_dash_inner/D_x[0], max(P_inner)]
P_p_cp_outer = [2*x_cp_dash_outer/D_x[1], max(P_outer)]

plt.figure()
plt.plot(P_h0_min_inner[0], P_h0_min_inner[1], 'ro', label='Minimum film thickness inner')
plt.plot(P_h0_min_outer[0], P_h0_min_outer[1], 'ro', label='Minimum film thickness outer')
plt.plot(P_h0_c_inner[0], P_h0_c_inner[1], 'go', label='Central film thickness inner')
plt.plot(P_h0_c_outer[0], P_h0_c_outer[1], 'go', label='Central film thickness outer')
plt.plot(P_p_es_inner[0], P_p_es_inner[1], 'bo', label='Peak pressure inner')
plt.plot(P_p_es_outer[0], P_p_es_outer[1], 'bo', label='Peak pressure outer')
plt.plot(P_p_cp_inner[0], P_p_cp_inner[1], 'yo', label='Center of pressure inner')
plt.plot(P_p_cp_outer[0], P_p_cp_outer[1], 'yo', label='Center of pressure outer')

plt.plot(2*x_inner/D_x[0], P_inner, label='Pressure inner')
plt.plot(2*x_outer/D_x[1], P_outer, label='Pressure outer')
plt.xlim([-2, 2])
plt.legend()

# make plots of the pressure for W_mark[0]/4, W_mark[0]/2, W_mark[0]
W_mark_values = [W_mark[0], W_mark[0]/2, W_mark[0]/4]
colors = ['r', 'g', 'b']
plt.figure()
i = 0
for W_val in W_mark_values:
    print(f"Calculating for W_mark = {W_val}")
    b = R_x[0] * np.sqrt(8*W_val/np.pi) # Eq. 17.37
    D_x = 2*b
    x_inner = np.linspace(-D_x, D_x, Nx)
    p_m = E_mark * np.sqrt(W_val/(2*np.pi)) # Eq. 17.40
    p_inner = p_m * np.array([cmath.sqrt(1 - (2 * xi / D_x) ** 2) for xi in x_inner]).real # Eq. 17.6, for Dy -> infinity
    P_inner = np.real(p_inner / p_m) # Eq. 18.23

    P_es_dash_inner = 0.648*W_val**0.185 * U[0]**0.275 * G**0.391 # Eq. 18.70
    X_ecp_dash_inner = -3.595*W_val**(-1.019) * U[0]**0.638 * G**(-0.358) # Eq. 18.76
    x_cp_dash_inner = X_ecp_dash_inner * R_x[0]
    H_ec_dash_inner = 2.922*W_val**(-0.166) * U[0]**0.692 * G**0.470 # Eq. 18.74
    h_c_dash_inner = H_ec_dash_inner * R_x[0]
    X_es_dash_inner = 1.111*W_val**0.606 * U[0]**(-0.021) * G**0.077 # Eq. 18.71

    P_p_es_inner = [X_es_dash_inner, P_es_dash_inner]
    P_p_cp_inner = [X_ecp_dash_inner, p_m/E_mark]

    plt.plot(P_p_es_inner[0], P_p_es_inner[1], colors[i] + 'o', label='Peak pressure inner')
    plt.plot(P_p_cp_inner[0], P_p_cp_inner[1], colors[i] + 'o', label='Center of pressure inner')
    plt.plot(x_inner/R_x[0] + X_ecp_dash_inner , p_inner/E_mark, colors[i] , label='Pressure inner')
    i = i + 1
plt.xlim([-0.03, 0.03])
plt.ylim([0, 0.008])
plt.legend()

plt.show()