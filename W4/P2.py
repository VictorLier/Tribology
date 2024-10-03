# Importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Define the symbols
ub = 10        # velocity of the slider
h0 = 60e-6     # minimum film thickness
l = 0.1        # length of the slider
eta = 0.1      # viscosity of the lubricant at operation temperature
cp = 2000      # specific heat capacity of the lubricant
rho = 850      # density of the lubricant

# 1 - Shoulder height
x, sh = sp.symbols('x sh')
h_expr = (h0 + sh) * sp.exp(-x / l)

# Solve for shoulder height 'sh' using boundary condition: h(l) = h0
shoulder_height = sp.solve(h_expr.subs(x, l) - h0, sh)[0]

print(f"Shoulder height: {shoulder_height}")
sh = shoulder_height

# Redefine h as a function of x using the solved shoulder height
h = (h0 + sh) * sp.exp(-x / l)

# 2 - Pressure within the bearing
p = sp.Function('p')(x)  # pressure as a function of x
A = sp.symbols('A')

# Differential equation for pressure
diff_eq = sp.Eq(p.diff(x), 6 * ub * eta / h**2 + eta * A / h**3)

# Boundary conditions
# p(0) = 0 (pressure at x = 0)
# p(l) = 0 (pressure at x = l)
boundary_conditions = {p.subs(x, l): 0}

# Solve the differential equation with the boundary conditions
pressure_solution = sp.dsolve(diff_eq, p, ics=boundary_conditions)

# Display the solution for pressure
print(f"Pressure solution: {pressure_solution}")

# Use solution and boundary condition 1 to solve for A

A = sp.solve(pressure_solution.rhs.subs(x, 0), A)[0]

# Substitute A into the pressure solution
pressure_solution = pressure_solution.subs('A', A)

x_vals = np.linspace(0, l, 100)
p_vals = [pressure_solution.rhs.subs(x, val) for val in x_vals]

# Plot the pressure distribution
plt.plot(x_vals, p_vals)
plt.xlabel("x [m]")
plt.ylabel("P [Pa]")
plt.title("Pressure distribution")
# plt.show()

# 3 - Normal load component and shear force components

# Normal load component
w_xam = sp.integrate(pressure_solution.rhs*h, (x, 0, l))
print(w_xam)
w_zam = sp.integrate(pressure_solution.rhs, (x, 0, l))
print(w_zam)
f_am = sp.integrate(-h/2*sp.diff(pressure_solution.rhs, x)+eta*ub/h, (x, 0, l))
print(f_am)

# 4 - friction factor and powerloss

w_xam = w_xam*10 # !!!!!!!!!!!!!!!

mu = (f_am + w_xam)/(w_zam*ub)
print(f"Friction factor: {mu}")

h_pm = (f_am + w_xam)*ub
print(f"Power loss: {h_pm}")

q_m = -h.subs(x, 0)**3/12*eta*sp.diff(pressure_solution.rhs, x).subs(x, 0)+ub*h.subs(x,0)/2
print(f"Flow rate: {q_m}")
# !!!!!!!!!!!!!!!!!

Dtm = h_pm/(rho*cp*q_m)
print(f"Thermal diffusivity: {Dtm}")