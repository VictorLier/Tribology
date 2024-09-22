
import sympy as sp


# 1)

# sp.diff(rho*h**3/12*eta*sp.diff(p,x),x) = h*(u_a+u_b)/2*sp.diff(rho,x)
# sp.diff((rho_0 * (1-alpha*((DeltaT*x/l + T_0)-T_0)))*h**3/12*eta*sp.diff(p,x),x) = h*(u_a)/2*sp.diff(rho_0 * (1-alpha*((DeltaT*x/l + T_0)-T_0)),x)
# (1-alpha*((DeltaT*x/l + T_0)-T_0)))*h**3/12*eta*sp.diff(p,x) = h*(u_a)/2*(- alpha*DeltaT/l)*x + C

import sympy as sp

# Define symbols
x, l, h, eta, u_a, u_b, T, T_0, DeltaT, rho, rho_0, alpha = sp.symbols('x l h eta u_a u_b T T_0 DeltaT rho rho_0 alpha')

# Define function p(x)
p = sp.Function('p')(x)

# Boundary conditions
u_b = 0

# Define temperature and density
T = DeltaT*x/l + T_0
rho = rho_0 * (1 - alpha*(T - T_0))

# Define the differential equation
eq1 = sp.Eq(sp.diff(rho*h**3/(12*eta)*sp.diff(p, x), x) - h*(u_a + u_b)/2 * sp.diff(rho, x), 0)

# Solve the differential equation with boundary conditions
SOL1 = sp.dsolve(eq1, p, ics={p.subs(x, 0): 0, p.subs(x, l): 0})

# Extract the solution
p_funcy = sp.simplify(SOL1.rhs)

# Output the solution
print(f"p(x) = {p_funcy}")

# 2) Compute maximum pressure

# values
u_a_val = 1
rho_0_val = 860
alpha_val = 6.4e-4
eta_val = 0.1
h_val = 10e-6
l_val = 40e-3
DeltaT_val = 5

# substitute values
p_funcy_val = p_funcy.subs({u_a: u_a_val, rho_0: rho_0_val, alpha: alpha_val, eta: eta_val, h: h_val, l: l_val, DeltaT: DeltaT_val})

# Find the maximum pressure
p_max = sp.solve(sp.diff(p_funcy_val, x), x)
print(p_max)
# Output the maximum pressure
print(f"Maximum pressure is {p_funcy_val.subs(x, p_max[0])}")

print(3*eta_val*u_a_val*alpha_val*DeltaT_val*l_val/(4*h_val**2))

# 3) Compute the load capacity per unit width of the bearing

# integrating p(x) from 0 to l to get the load capacity per unit width
load_capacity = sp.integrate(p_funcy, (x, 0, l))
print(sp.simplify(load_capacity))

# 4) compute the flow rate at inlet, outlet and the point of maximum pressure

# flow rate at inlet
q_mark_x_inlet =  - h_val**3/(12*eta_val)*sp.diff(p_funcy_val, x).subs(x, 0) + u_a_val*h_val/2
print(f"Flow rate at the inlet is {q_mark_x_inlet}")

# flow rate at outlet
q_mark_x_outlet =  - h_val**3/(12*eta_val)*sp.diff(p_funcy_val, x).subs(x, l_val) + u_a_val*h_val/2
print(f"Flow rate at the outlet is {q_mark_x_outlet}")

# flow rate at the point of maximum pressure
q_mark_x_max =  - h_val**3/(12*eta_val)*sp.diff(p_funcy_val, x).subs(x, p_max[0]) + u_a_val*h_val/2
print(f"Flow rate at the point of maximum pressure is {q_mark_x_max}")

# 5) sketch the flow profiles at the computed instances

import matplotlib.pyplot as plt
import numpy as np

# Define the flow profile
x_vals = np.linspace(0, l_val, 100)
p_vals = [p_funcy_val.subs(x, x_val) for x_val in x_vals]

# Plot the flow profile
plt.plot(x_vals, p_vals)
plt.xlabel('x')
plt.ylabel('p(x)')
plt.title('Flow profile')
plt.show()





