import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Define the variables
d_i = 69.5e-3    # [m] - Inner race diameter
d_o = 85.5e-3    # [m] - Outer race diameter
d   = 8e-3       # [m] - Diameter of cylindrical rollers
l   = 8e-3       # [m] - Axial length of cylindrical rollers
n   = 20         # [-] - Numer of rollers
w_z_max = 8e3    # [N] - Maximum radial load
omega_i = 837    # [rad/s] - Inner race angular speed
omega_o = 0      # [rad/s] - Outer race angular speed

E   = 2.1e11     # [Pa] - Young's modulus
nu  = 0.3        # [-] - Poisson's ratio

eta_0 = 0.01     # [Pa*s] - Viscosity of lubricant
xi  = 2e-8       # [1/Pa/s] - Pressure-viscosity coefficient

# --------------------- Question 1 ------------------------------------

# 1) ------------------------------------------------------------------
w_z = n * w_z_max / 4   # Eq. 21.47
print(f"Load of the roller with the highest load: {w_z*1e-3:.2f} kN")

# 2) ------------------------------------------------------------------

Lambda = 3
Rqa = 0.2e-6  * np.pi/(2*np.sqrt(2))
Rqb = 0.05e-6 * np.pi/(2*np.sqrt(2))

h_min = Lambda * (np.sqrt(Rqa**2 + Rqb**2))
print(f"Minimum film thickness: {h_min*1e6:.2f} um")

# --------------------- Question 2 ------------------------------------

# 1) ------------------------------------------------------------------

r_ax = np.array([d/2, d/2])
r_bx = np.array([d_i/2, -d_o/2])
R_x = ( 1 / r_ax + 1 / r_bx )**(-1) # Eq. 17.4
print(f"Equivalent radius for inner: {R_x[0]*1e3:.2f} mm")
print(f"Equivalent radius for outer: {R_x[1]*1e3:.2f} mm")

w_zmark = w_z_max / l
E_mark = E / (1 - nu**2) # Eq. 17.17
W_mark = w_zmark / (E_mark * R_x) # Eq. 17.38

p_m = E_mark * np.sqrt(W_mark/(2*np.pi)) # Eq. 17.40
print(f"Maximum contact pressure for inner: {p_m[0]*1e-9:.2f} GPa")
print(f"Maximum contact pressure for outer: {p_m[1]*1e-9:.2f} GPa")

delta_m = 2*W_mark*R_x / np.pi * (np.log(2*np.pi/W_mark) - 1) # Eq. 17.39
print(f"Maximum deflection for inner: {delta_m[0]*1e6:.2f} um")
print(f"Maximum deflection for outer: {delta_m[1]*1e6:.2f} um")


# 2) ------------------------------------------------------------------

# The deflection is of the same order of magnitude as the minimum film thickness

# 3) ------------------------------------------------------------------

Nx = 200

b = R_x * np.sqrt(8*W_mark/np.pi) # Eq. 17.37

# -------- For the inner ------------

D_x = 2*b[0] # ?
x = np.linspace(-D_x/2, D_x/2, Nx)
x_extended = np.linspace(-2*D_x, 2*D_x, 2*Nx) 

p = p_m[0] * np.sqrt(1-(2*x/D_x)**2) # Eq. 17.6, for Dy -> infinity
P = p / p_m[0] # Eq. 18.23
delta_dash = 0
delta_dash = np.zeros(Nx)
X = np.linspace(-1, 1, Nx)
Delta = X[1] - X[0]
for i in range(1, Nx-1):
    for j in range(1, Nx-1):
        term1 = np.abs( (X[i+1] + X[i])/2 - X[j] )
        term2 = np.abs( (X[i-1] + X[i])/2 - X[j] )
        delta_dash[i] = delta_dash[i] - Delta/(2*np.pi)*P[j]*np.log( term1 * term2 )   # Eq. 18.31

S = x**2 / (2 *R_x[0]) # Eq. 12.52
S_extended = x_extended**2 / (2 *R_x[0])

delta = D_x**2 * delta_dash / R_x[0] # Eq. 18.23
delta0_inner = np.zeros(len(x_extended))
delta0_inner[Nx//2 : Nx + Nx//2] = delta

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
plt.axis('equal')
plt.grid()

plt.figure()
plt.plot(x_extended*1e6, S_extended*1e6, label="Undeformed roller", linestyle='--')
plt.plot(x_extended*1e6, (S_extended + delta0_inner)*1e6, label = "Deformed roller", linestyle='-')
plt.xlabel("x [um]")
plt.ylabel("y [um]")
plt.legend()
plt.axis('equal')
plt.grid()

p_extended = np.zeros(2*Nx) - 1e9
p_extended[Nx//2 : Nx + Nx//2] = p

np.savetxt("Report_3/data/Q2_Si.txt", np.array([x_extended*1e6, S_extended*1e6]).T)
np.savetxt("Report_3/data/Q2_defi.txt", np.array([x_extended*1e6, (S_extended + delta0_inner)*1e6]).T)
np.savetxt("Report_3/data/Q2_pi.txt", np.array([x_extended*1e6, p_extended*1e-9]).T)

# --------- for the outer ------------

D_x = 2*b[1] # ?
x = np.linspace(-D_x/2, D_x/2, Nx)
x_extended = np.linspace(-2*D_x, 2*D_x, 2*Nx) 
p = p_m[1] * np.sqrt(1-(2*x/D_x)**2) # Eq. 17.6, for Dy -> infinity
P = p / p_m[1] # Eq. 18.23
delta_dash = 0
delta_dash = np.zeros(Nx)
for i in range(1, Nx-1):
    for j in range(1, Nx-1):
        term1 = np.abs( (X[i+1] + X[i])/2 - X[j] )
        term2 = np.abs( (X[i-1] + X[i])/2 - X[j] )
        delta_dash[i] = delta_dash[i] - Delta/(2*np.pi)*P[j]*np.log( term1 * term2 )   # Eq. 18.31

S = x**2 / (2 *R_x[1]) # Eq. 12.52
S_extended = x_extended**2 / (2 *R_x[1])

delta = D_x**2 * delta_dash / R_x[1] # Eq. 18.23
delta0_outer = np.zeros(len(x_extended))
delta0_outer[Nx//2 : Nx + Nx//2] = delta

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
plt.axis('equal')
plt.grid()

plt.figure()
plt.plot(x_extended*1e6, S_extended*1e6, label="Undeformed roller", linestyle='--')
plt.plot(x_extended*1e6, (S_extended + delta0_outer)*1e6, label = "Deformed roller", linestyle='-')
plt.xlabel("x [um]")
plt.ylabel("y [um]")
plt.legend()
plt.axis('equal')
plt.grid()

plt.show()

p_extended = np.zeros(2*Nx) - 1e9
p_extended[Nx//2 : Nx + Nx//2] = p

np.savetxt("Report_3/data/Q2_So.txt", np.array([x_extended*1e6, S_extended*1e6]).T)
np.savetxt("Report_3/data/Q2_defo.txt", np.array([x_extended*1e6, (S_extended + delta0_outer)*1e6]).T)
np.savetxt("Report_3/data/Q2_po.txt", np.array([x_extended*1e6, p_extended*1e-9]).T)

# 4) ------------------------------------------------------------------

# ...

# 5) ------------------------------------------------------------------

d_o = 85.51e-3    # [m] - Outer race diameter
c_d = d_o - d_i - 2*d # Eq. 21.2
print(f"Clearance: {c_d*1e6:.2f} um")

r_bx = np.array([d_i/2, -d_o/2])
R_x = ( 1 / r_ax + 1 / r_bx )**(-1) # Eq. 17.4
W_mark = w_zmark / (E_mark * R_x) # Eq. 17.38
delta_m = 2*W_mark*R_x / np.pi * (np.log(2*np.pi/W_mark) - 1) # Eq. 17.39

Z_w_top = np.pi * ( 1 - c_d / (2*delta_m) )**(3/2)
Z_w_bottom = 2.491 * ( np.sqrt( 1 + (( 1 - c_d / (2*delta_m) ) / 1.23 )**2) - 1) 
Z_w = Z_w_top / Z_w_bottom # Eq. 21.46

w_z_clearance = n * w_z_max / Z_w # Eq. 21.44
print(f"Load with clearance for inner: {w_z_clearance[0]*1e-3:.2f} kN")
print(f"Load with clearance for outer: {w_z_clearance[1]*1e-3:.2f} kN")
print(f"Using Z_w for inner: {Z_w[0]:.2f}")
print(f"Using Z_w for outer: {Z_w[1]:.2f}")
