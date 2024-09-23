
import sympy as sp

# defining simplified differential equation Eq. 7.46

# defining the symbols
x, eta, rho, w_a, r_i, r_o, h = sp.symbols('x eta rho w_a r_i r_o h')

# defining the function p(x)
p = sp.Function('p')(x)

# defining the differential equation
eq1 = sp.diff(rho*h**3/(12*eta)*sp.diff(p, x), x) - rho*w_a

# defining the boundary conditions
p_i = 0
p_o = 0

# solving the differential equation

SOL1 = sp.dsolve(eq1, p, ics={p.subs(x, r_i): p_i, p.subs(x, r_o): p_o})

print(sp.simplify(SOL1))

# # defining the variables
r_i_val = 95e-3
r_o_val = (250e-3)/2
h_val = 50e-6
eta_val = 0.1
rho_0_val = 860
# w_z_val = 96000.0

