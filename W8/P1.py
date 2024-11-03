import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Someya import Table1

Table1 = np.flip(Table1, axis=0)
S = Table1[:,0]     # Sommerfeld number
E = Table1[:,1]     # eccentricity ratio
Phi = Table1[:,2]   # attitude angle (deg)
Q = Table1[:,3]     # flow rate ratio
P = Table1[:,4]     
T = Table1[:,5]

# variables
D = 100e-3  # [m]
b = 50e-3   # [m]
N = np.linspace(10, 30, 100)    # [Hz]
omega = 2 * np.pi * N   # [rad/s]
c = 100e-6  # [m]
w = 5e3    # [N]
eta =  0.08   # [N*s/m^2]
rho = 876  # [kg/m^3]
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

# speed 30 Hz is now only considered
# take heat balance into account

iterations = 10

w_air = 1 # [m/s] air velocity of surroundings
# Vogelpohl equation
alpha = 9.807 * (0.7 + 1.2 * w_air**(1/2))
# Reitemeyer equation
A = 9 * D * D**(1/2)
rho = 876 # [kg/m^3]
Cp = 1964 # [J/kg*K]
lamb = 1/3

# log(log(nu+0.8)) = -m*log(T) + k
nu_40 = eta/rho*1e6
t_40 = 40
nu_100 = 0.007/rho*1e6
t_100 = 100
m_lub = (np.log(np.log(nu_100 + 0.8)) - np.log(np.log(nu_40 + 0.8))) / (np.log((t_40+273.15)/(t_100+273.15)))
k_lub = np.log(np.log(nu_40 + 0.8)) + m_lub*np.log(t_40+273.15)
t_1 = 30
t_0 = 20
p_f = 2e5 # [Pa]
L_mark = b/2
damping = 0.5

def eta_i(temp):
    return rho*1e-6*(np.exp(np.exp(-m_lub*np.log(temp+273.15) + k_lub)) - 0.8)

t_mean = np.zeros(iterations)
t = np.zeros(iterations)
t[0] = 40
N = 30
omega = 2*np.pi*N

for i in range(iterations-1):
    S_current = eta*N*b*D/w*(r_b/c)**2
    Q_current = np.interp(S_current, S, Q)
    T_current = np.interp(S_current, S, T)
    epsi_current = np.interp(S_current, S, E)

    eta = eta_i(t[i])
    q_f = np.pi*c**3/(3*eta*L_mark/D) * (1+3/2*epsi_current**2)*p_f
    q = r_b*omega*c*b*Q_current + q_f   # assuming chi = 1
    f_J = psi * T_current

    t_new = ((1-lamb)*(f_J*r_b*w*omega + alpha*A*t_0) + Cp*rho*q*t_1) / (Cp*rho*q + alpha*A*(1-lamb))
    t[i+1] = t[i] + damping * (t_new - t[i])

iterations = np.arange(iterations)
plt.plot(iterations, t)
# plt.plot(iterations, t)
plt.show()