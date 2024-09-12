from math import pi

# Variables
dia = 0.004 # [m]
width = 0.004 # [m]
visc = 0.17 # [Pa*s]
rpm = 1800 # [rpm]
load = 2200 # [N]
h_min = 45.6e-6 # [m]
R_q = 1.3e-6 # [m]
friction = 0.021 # [-]

# 1

# Film paramter
Lambda = h_min/ ((R_q**2 + R_q**2)**0.5) # [-]

if Lambda < 1:
    print(f"With a film parameter of {Lambda:.2f} the bearing is hydrodynamic")
if Lambda > 1 and Lambda < 5:
    print(f"With a film parameter of {Lambda:.2f} the bearing is partial lubrication")
if Lambda > 3 and Lambda < 10:
    print(f"With a film parameter of {Lambda:.2f} the bearing is elastohydrodynamic lubrication")
if Lambda > 5 and Lambda < 100:
    print(f"With a film parameter of {Lambda:.2f} the bearing is hydrodynamic lubrication")

# Hersey number

# Pressure
Pres = load / (width * dia * 0.7) # [Pa]  -  0.7 er fundet på
v = pi * dia * rpm / 60 # [m/s] - Velocity
Hersey = (visc * v) / (Pres) # [-] - Hersey number

print(f"Hersey number is {Hersey:.3G}")


# 2
# When running down, from fig. 3.13 the friction and wear will first fall, until a film parameter of 10 is reached, then the friction and wear will rise again.

# 3 - Determine the film thickness at Lambda = 4
import sympy as sp
lam = 5
h = sp.symbols('h')
eq1 = -lam + h / ((R_q**2 + R_q**2)**0.5)
sol = sp.solve(eq1, h)
print(sol)

# in the given hydrodynamic regime, the friction scales proportionally with the film height:

# making the sol variable a float
sol = float(sol[0])

mu_3 = sol/h_min*friction
print(f"The new scaled friction is {mu_3:.4f}")

# 4 - At this transition point, the balance between load and speed have changed consider:
#     • The load that can be carried if the rotational speed is constant

# a relation between the minimum film thickness , speed and load is given by the following proportion:
# h proportional to (u/w)**0.5
# where u is the speed and w is the load
# We use this relation to determine the load that can be carried if the rotational speed is constant
# recall that the viscocity and film thickness are proportional
# mu1/mu2 = (u1/u2)**0.5 / (w1/w2)**0.5 => mu1/mu2 = 1 / (w1/w2)**0.5 => w2 = (mu1/mu2)**2 * w1

w1 = load
mu1 = friction
mu2 = mu_3

w2 = (mu1/mu2)**2 * w1

print(f"The load that can be carried if the rotational speed is constant is {w2:.2f} N")

#     • The speed in case the load is constant

# The same relation can be used to determine the speed in case the load is constant
# mu1/mu2 = (u1/u2)**0.5 / (w1/w2)**0.5 => mu1/mu2 = (u1/u2)**0.5 => u2 = u1 * (mu1/mu2)**(-2)

u1 = rpm

u2 = u1 * (mu1/mu2)**(-2)
# converting to rps
u2 = u2 / 60

print(f"The speed in case the load is constant is {u2:.2f} rps")


# 5 - stribeck curve plotting



print("Stop")