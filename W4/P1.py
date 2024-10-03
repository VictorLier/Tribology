
import sympy as sp
import math
import numpy as np
import matplotlib.pyplot as plt

# Define the symbols
l1 = 100e-3 # length of shoe in x direction
l2 = 400e-3 # Width of slider perpindicular to x direction
h0 = 40e-6 # minimum film thickness
sh = 55e-6 # shoulder height
ub = 5 # velocity of the slider
eta = 0.05 # viscosity of the lubricant at operation temperature

# slider is consideres infinitely wide

# 1 - load carrying capacity
H0 = h0/sh # dimensionsless term
print("H0 = ", H0)

Wz = 6*math.log((H0+1)/H0) - 12/(1+2*H0)
wz_mark = Wz*eta*ub*l1**2/sh**2
wz = wz_mark*l1
print(f" Wz = {Wz}")
print(f" wz,mark = {wz_mark}")
print(f"Load carrying capacity: {wz} [N]")

# 2 - maximum pressure and sketch pressure distribution

x = np.linspace(0, l1, 100)

X = x/l1
P = 6*X*(1-X)/((H0+1-X)**2*(1+2*H0))

p = P*eta*ub*l1/sh**2

plt.plot(x, p)
plt.xlabel("x [m]")
plt.ylabel("P [Pa]")
plt.title("Pressure distribution")
plt.show()

# maximum pressure

Pmax = 3*eta*ub*l1*sh/(2*h0*(sh+h0)*(sh+2*h0))

print(f"Maximum pressure: {Pmax} [Pa]")

# 3- Friction factor in the bearing
mu = (2*sh*sp.log(H0/(1+H0)) + 3*sh/(1+2*H0)) / (3*l1*sp.log(H0/(1+H0)) + 6*l1/(1+2*H0))
print(f"Friction factor: {mu}")

# 4 - length l_Fx of the bearing
# pivot point is located such that net torque is zero
# Xcp = 1/Wz*int(P*XX dX)

X = sp.symbols('X')
P = 6*X*(1-X)/((H0+1-X)**2*(1+2*H0))
Xcp = 1/Wz*sp.integrate(sp.simplify(P)*X, (X, 0, 1))
xcp = Xcp*l1

print(f"location of the resulting force from pressure: {xcp} [m]")

# !!!!!!!!!!!!!!