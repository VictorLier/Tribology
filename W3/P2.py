
# A water-lubricated journal bearing in a boiler feed pump has a shaft of radius 0.1m which
# rotates at 10 rev/s. The kinematic viscosity in the full fluid film region may be taken
# as directly proportional to the film thickness and has a value of 4 *10^-7 m^2 s^-1 at a film
# thickness equal to the the radial clearance of 0.1 mm. Determine if the bearing is operating
# under laminar or turbulent conditions. If laminar flow is predicted, what change in these
# operating conditions would produce the onset of vortex flow?

import sympy as sp
import numpy as np

# Given data
r = 0.1 # Radius of shaft in m
n = 10 # Revolutions per second
nu = 4e-7 # Kinematic viscosity in m^2 s^-1
h = 0.1e-3 # Radial clearance in m (and film thickness)

# Angular velocity
omega = 2 * np.pi * n   # Angular velocity in rad/s
v = r * omega # Velocity of the shaft in m/s
l0 = 2 * np.pi * r # Circumference of the shaft

# Reynolds number
Re_x = (v * h**2) / (nu * l0)

print(f"The Reynolds number is {Re_x}")

if Re_x < 1000:
    print(f"The bearing is operating under laminar conditions.")
    # print(f"Change in operating conditions to produce the onset of vortex flow is to increase the speed of the shaft.")
else:
    print(f"The bearing is operating under turbulent conditions")
    # print(f"Change in operating conditions to produce the onset of vortex flow is to increase the viscosity of the fluid, or lower velocity of film thickness")

print("---")

# Taylor number:
Ta = v**2 * h**3 / (nu**2 * r)

print(f"The Taylor number is {Ta}")

if Ta > 1700:
    print("Vorticies are expected to form in the fluid film.")

nu_new = sp.sqrt(v**2 * h**3/(1700*r))

print(f"The new kinematic viscosity required to produce the onset of vortex flow is {nu_new}")