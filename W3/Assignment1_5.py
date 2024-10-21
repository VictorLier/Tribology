
import sympy as sp
import os
import numpy as np
import matplotlib.pyplot as plt

# defining simplified differential equation Eq. 7.46

x, eta, rho, w_a, r_i, r_o, h = sp.symbols('x eta rho w_a r_i r_o h')
p = sp.Function('p')(x)
eq1 = sp.diff(rho*h**3/(12*eta)*sp.diff(p, x), x) - rho*w_a # defining the differential equation

# defining the boundary conditions
p_i = 0
p_o = 0

# solving the differential equation
SOL1 = sp.dsolve(eq1, p, ics={p.subs(x, r_i): p_i, p.subs(x, r_o): p_o})

print(sp.simplify(SOL1.rhs))

# defining the variables
r_i_val = 95e-3         # inner radius
r_o_val = (250e-3)/2    # outer radius
h0_val = 50e-6          # film thickness
rho_0_val = 850          # density
eta_val = 32e-6 * rho_0_val # dynamic viscosity
w_a_val = -0.1          # velocity
W_z_val = 4e5           # load capacity

h_vals = [h0_val, h0_val/2, h0_val/4, h0_val/8]
w_a_vals = [-0.1, -0.05, -0.01, -0.005]
x_vals = np.linspace(r_i_val, r_o_val, 100)

# Create directory if it does not exist
output_dir = 'W3'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# plotting the pressure distribution for different film thicknesses and w_a values
plt.figure()
for w_a_val in w_a_vals:
    for h_val in h_vals:
        p_expr = SOL1.rhs.subs({eta: eta_val, w_a: w_a_val, r_i: r_i_val, r_o: r_o_val, h: h_val})
        p_vals = [p_expr.subs(x, x_val).evalf() for x_val in x_vals]
        filename = os.path.join(output_dir, f'1_5_wa_{w_a_val:.2e}_h_{h_val:.2e}.txt')
        np.savetxt(filename, np.array([x_vals, p_vals]).T)
        plt.plot(x_vals, p_vals, label=f'h = {h_val:.2e} m')

    plt.xlabel('x (m)')
    plt.ylabel('Pressure distribution (Pa)')
    plt.title(f'Pressure distribution vs x for w_a = {w_a_val:.2e} m/s^2')
    plt.legend()
    plt.show()


# Derive an expression for the flow rate and sketch the  flow profile at ri, ro and the point of maximum pressure.

# ask if it should be made by hand
r_mean_val = (r_i_val + r_o_val)/2
z = sp.symbols('z')
u = sp.Function('u')(x,z)
p = sp.simplify(SOL1.rhs)
u = sp.simplify(-z*((h-z)/(2*eta))*sp.diff(p, x))
print(f"The velocity profile is given by: {u}")
q_mark = sp.integrate(u, (z, 0, h))
print("Flow rate")
print(q_mark)


# point of maximum pressure by differentiating with respect to x and setting to zero
p_diff = sp.diff(p, x)
p_max = sp.solve(p_diff, x)
print("Point of maximum pressure")
print(p_max)
print("maximum pressure")
print((p.subs(x, p_max[0])).subs({eta: eta_val, r_i: r_i_val, r_o: r_o_val, h: h0_val, w_a: w_a_val}).evalf())
print(sp.simplify(p.subs(x, r_i/2 + r_o/2)))

print("load capacity per unit width")
Load_capacity = W_z_val/(2*sp.pi*r_mean_val)
# printing the load capacity per unit width as a float
print(float(Load_capacity.evalf()))

# Compute the velocity wa that balance the external load and the time it may take
# for the bearing surfaces to get into contact at this speed (assuming it is constant).

r_m = (r_i + r_o)/2
W_z = sp.symbols('W_z')

eq2 = W_z/(2*sp.pi*r_m) - sp.integrate(p, (x, r_i, r_o))
SOL2 = sp.solve(eq2, w_a)
w_a = SOL2[0]
print("Velocity wa that balance the external load")
print(sp.simplify(w_a))
w_a = w_a.subs({eta: eta_val, r_i: r_i_val, r_o: r_o_val, W_z: W_z_val})
print(float(w_a.subs(h, h0_val).evalf()))

# the time it takes for the bearing surfaces to get into contact at this speed
t = h/w_a
print("Time it takes for the bearing surfaces to get into contact at this speed")
print(float(t.subs({h: h0_val}).evalf()))

print(-t.subs({h: h0_val, w_a: w_a.subs({eta: eta_val, r_i: r_i_val, r_o: r_o_val})}))
print(-t.subs({h: h0_val/2, w_a: w_a.subs({eta: eta_val, r_i: r_i_val, r_o: r_o_val})}))
print(-t.subs({h: h0_val/4, w_a: w_a.subs({eta: eta_val, r_i: r_i_val, r_o: r_o_val})}))
print(-t.subs({h: h0_val/8, w_a: w_a.subs({eta: eta_val, r_i: r_i_val, r_o: r_o_val})}))

# Plot the time it takes for the bearing surfaces to get into contact as a function of the film thickness.
h_vals = np.linspace(h0_val, h0_val/8, 100)
t_vals = [-t.subs({h: h_val, w_a: w_a.subs({eta: eta_val, r_i: r_i_val, r_o: r_o_val})}) for h_val in h_vals]
h_vals_micro = [h_val*1e6 for h_val in h_vals]
t_vals_mili = [t_val*1e3 for t_val in t_vals]   
plt.figure()
plt.plot(h_vals_micro, t_vals_mili)
plt.xlabel('Film thickness (µm)')
plt.ylabel('Time (ms)')
plt.show()

# expoerting data to latex
np.savetxt('W3/Assignment1_5.txt', np.array([h_vals_micro, t_vals_mili]).T)