
import sympy as sp
import numpy as np

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

print(sp.simplify(SOL1.rhs))

# defining the variables
r_i_val = 95e-3
r_o_val = (250e-3)/2
h0_val = 50e-6
eta_val = 0.1
rho_0_val = 860
w_a_val = -0.1
# w_z_val = 96000.0

# make a plot of the pressure distribution, using the given values of the variables
# but for different film thicknesses; h = {h0, h0/2, h0/4, h0/8}
import matplotlib.pyplot as plt
# defining the values of the film thickness
h_vals = [h0_val, h0_val/2, h0_val/4, h0_val/8]

# defining the values of w_a
w_a_vals = [-0.1, -0.05, -0.01, -0.005]

# defining the range for x
x_vals = np.linspace(r_i_val, r_o_val, 100)

# plotting the pressure distribution for different film thicknesses and w_a values
plt.figure(figsize=(12, 8))

for w_a_val in w_a_vals:
    for h_val in h_vals:
        p_expr = SOL1.rhs.subs({eta: eta_val, w_a: w_a_val, r_i: r_i_val, r_o: r_o_val, h: h_val})
        p_vals = [p_expr.subs(x, x_val).evalf() for x_val in x_vals]
        plt.plot(x_vals, p_vals, label=f'h = {h_val:.2e} m, w_a = {w_a_val:.2e} m/s^2')

plt.xlabel('x (m)')
plt.ylabel('Pressure distribution (Pa)')
plt.title('Pressure distribution vs x for different film thicknesses and w_a values')
plt.legend()
plt.show()

# Derive an expression for the flow rate and sketch the  flow profile at ri, ro and the point of maximum pressure.

# ask if it should be made by hand
# EQ 7.28
# EQ 7.36

r_mean_val = (r_i_val + r_o_val)/2
z = sp.symbols('z')
u = sp.Function('u')(x,z)


p = sp.simplify(SOL1.rhs)

u = sp.simplify(-z*((h-z)/(2*eta))*sp.diff(p, x))

q_mark = sp.integrate(u, (z, 0, h))

# point of maximum pressure by differentiating with respect to x and setting to zero
p_diff = sp.diff(p, x)
p_max = sp.solve(p_diff, x)

# Compute the velocity wa that balance the external load and the time it may take
# for the bearing surfaces to get into contact at this speed (assuming it is constant).

wz = sp.symbols('wz')

eq2 = wz - sp.integrate(p, (x, r_i_val, r_o_val))*(2*sp.pi*r_mean_val)

SOL2 = sp.solve(eq2, w_a)

w_a = SOL2[0]

print(w_a.subs({eta: eta_val, r_i: r_i_val, r_o: r_o_val, h: h0_val}))

# the time it takes for the bearing surfaces to get into contact at this speed
t = h0_val/w_a

# plotting the time as a function of the load w_z

time = t.subs({eta: eta_val, r_i: r_i_val, r_o: r_o_val, h: h0_val})

print("Time expression:", time)

wz_vals = np.linspace(0, 100000, 100)

Coeff = -float(sp.simplify(time*wz))
print("Coeff:", Coeff)

plt.figure(figsize=(12, 8))
plt.plot(wz_vals, Coeff/wz_vals)
plt.xlabel('Load (N)')
plt.ylabel('Time (s)')
plt.title('Time vs Load')
plt.show()



# USE WZ FROM PROBLEM 1