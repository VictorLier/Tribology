
# An oil pump is designed like a journal bearing with a constant oil 
# film thickness according to the sketch. The wrap angle is 320°. 
# The film thickness is 0.2 mm and the oil viscosityis 0.022 N s m−2.  
# The shaft diameter is 60 mm and the length of the pump in the axialdirection is 50 mm.  
# Find the volume pumped per unit time as a function of the resisting pressure at the rotational speed Na= 1500 rpm.

import sympy as sp
import numpy as np

# Given data
alpha = 320 # Wrap angle in degrees
h = 0.2e-3 # Film thickness in m
eta = 0.022 # Viscosity in N s m^-2
d = 60e-3 # Shaft diameter in m
l = 50e-3 # Length of pump in axial direction in m
n = 1500 # Revolutions per minute [rpm]
r = d/2 # Radius of shaft in m
omega = 2 * np.pi * n / 60 # Angular velocity in rad/s

p = sp.symbols('p')
u = omega * r # Velocity of the shaft in m/s

q_mark_x = l*(-h**3/(12 * eta)*p/(alpha/360*sp.pi*d) + u*h/2)

print(f"The volume pumped per unit time as a function of the resisting pressure is {q_mark_x.evalf(3)}")
# print(f"The volume pumped per unit time as a function of the resisting pressure is {q_mark_x.evalf(3)}")

# plotting the function q_mark_x

import matplotlib.pyplot as plt


# plotting pressure from 0 to 26 bar
p_values = np.linspace(0, 26e+5, 100)
q_values = [q_mark_x.evalf(subs={p: p_val}) for p_val in p_values]

plt.plot(p_values, q_values)
plt.xlabel('Resisting pressure (p) [Pa]')
plt.ylabel('Volume pumped per unit time (q) [m^3/s]')
plt.title('Volume pumped per unit time as a function of the resisting pressure')
plt.grid()
plt.show()