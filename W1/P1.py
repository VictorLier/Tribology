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
lam, h, rq = sp.symbols('lam h rq')
eq1 = -lam + h / ((rq**2 + rq**2)**0.5)
eq2 = lam - 4
sol = sp.solve([eq1, eq2], (h, rq))
print(sol)


print("Stop")