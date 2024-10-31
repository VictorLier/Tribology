import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Someya import Table1

S = np.flip(Table1[:,0])
E = np.flip(Table1[:,1])
Phi = np.flip(Table1[:,2])

# variables
D = 100e-3  # [m]
b = 50e-3   # [m]
N = np.linspace(10, 30, 100)    # [Hz]
omega = 2 * np.pi * N   # [rad/s]
c = 100e-6  # [m]
w = 5e3    # [N]
eta =  0.08   # [N*s/m^2]
rho = 860  # [kg/m^3]
nu = eta / rho  # [m^2/s]
r_b = D/2
psi = c/r_b


epsi_short = np.zeros(len(N))
epsi_long = np.zeros(len(N))

def short(epsi,omega):
    return eta*omega*r_b*b**3/(4*c**2) * epsi/((1-epsi**2)**2) * (16*epsi**2 + np.pi**2*(1-epsi**2))**(1/2) - w

def long(epsi,omega):
    return eta*omega*r_b*(r_b/c)**2 * 6*epsi * (np.pi**2 - epsi**2 * (np.pi**2 - 4))**(1/2)/((2+epsi**2)*(1-epsi**2)) - w/b

for i in range(len(N)):
    epsi_short[i] = fsolve(short, 0.9, args = (omega[i]))[0]
    epsi_long[i] = fsolve(long, 0.9, args = (omega[i]))[0]

Phi_short = np.arctan(np.pi*(1-epsi_short**2)**(1/2)/(4*epsi_short))
Phi_long = np.arctan(np.pi/(2*epsi_long)*(1-epsi_long**2)**(1/2))
e_short = epsi_short * c
e_long = epsi_long * c

Sommerfeld = eta*N*D*b/(psi**2*w)
epsi_num = np.interp(Sommerfeld, S, E)
Phi_num = np.interp(Sommerfeld, S, Phi)*np.pi/180

# making a polar plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(Phi_short, e_short)
ax.plot(Phi_long, e_long)
ax.plot(Phi_num, epsi_num*c)
ax.set_theta_zero_location('S')
ax.grid(True)
plt.show()