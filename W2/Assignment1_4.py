
# Consider two parallel discs that rotate relative to each other (one steady and one
# rotating) using the aforementioned (see Lecture 1 problems) diameter of 250 mm and
# speed of 1000 rpm. Here some basic investigations regarding friction torque and power
# are conducted.
# We may assume that the discs have a hole i.e. inner radius being variable ri ∈ [0; ro] while
# outer radius ro = 250 mm/2 = 125 mm

# * Derive an expression by the integration of the shear stress τ on the bearing surface
# that relates the torque to the key parameters such as viscosity η, speed ω, thickness
# h and radii. (HINT: You may consider an infinitesimal element that is integrated
# over the circular surface between the discs).

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Define symbols
r, theta, ri, ro, h, eta, omega = sp.symbols('r theta ri ro h eta omega')
tau = eta*omega*r/h

# dtau_glob = tau(r)*dr*r*dtheta
tau_glob = sp.integrate(sp.integrate(tau*r, (r, ri, ro)), (theta, 0, 2*sp.pi))
print(sp.simplify(tau_glob))

# • Extend the derivation to an expression for the frictional power.

# dT = tau_glob*ro*dr => dT = tau(r)*r*dr*dtheta
T = sp.integrate(sp.integrate(tau*r*r, (r, ri, ro)), (theta, 0, 2*sp.pi))
print(sp.simplify(T))

# • For a given film height h = 25 μm and inner radius ri = 95mm plot a power vs.
#   time curve and discuss it. (You should use your oil data from question 3. You may
#   also need heat capacity and density - look in the oil data sheet.)

# Define variables:
h_val = 25e-6 # [m] - film height
ri_val = 95e-3 # [m] - inner radius
ro_val = 125e-3 # [m] - outer radius
omega_val = 1000*2*np.pi/60 # [rad/s] - angular velocity
rho = 850 # [kg/m^3] - density of lubricant
c_p = 2000 # [J/kg/K] - specific heat capacity of lubricant
m_lubricant = np.pi * (ro_val**2 - ri_val**2) * h_val * rho  # [kg] - mass of the lubricant
Temp_data = np.array([27, 37.6, 49.0, 60.7, 76.4]) + 273.15  # Temperature in Kelvin
nu_data = np.array([53.7, 31.4, 21.4, 13.9, 9.05])        # Kinematic viscosity in cSt

# 2 - fit the data
modify_nu = np.log(np.log(nu_data+0.8))
modify_T = np.log(Temp_data)
a=np.polyfit(modify_T, modify_nu, 1)
m_coeff = a[0]
C_coeff = a[1]

step = 100

# Define the time
time = np.linspace(0, 1, step)
Power = np.zeros(step)
Temperature = np.zeros(step)
Temperature[0] = 273.15 + 27 # Initial temperature
nu = np.zeros(step)

# Convert the symbolic expression T to a numerical function
T_func = sp.lambdify((eta, omega, h, ri, ro), sp.simplify(T), 'numpy')

for i in range(len(time)):
    nu[i] = np.exp(np.exp(m_coeff*np.log(Temperature[i]) + C_coeff)) - 0.8
    eta_val = nu[i] * rho * 10e-6 # [N.s/m^2] - dynamic viscosity
    Power[i] = T_func(eta_val, omega_val, h_val, ri_val, ro_val) * omega_val
    if i < len(time)-1:
        Temperature[i+1] = time[1] * Power[i] / (c_p * m_lubricant) + Temperature[i]

# plot the power vs time
plt.figure(1)   
plt.clf()   
plt.plot(time[1:], Power[1:], '-', label="Power vs Time")  
plt.grid(True)
plt.xlabel('Time [s]')
plt.ylabel('Power [W]')
plt.title('Power vs Time')
plt.show()

# plotting temperature vs time
plt.figure(2)
plt.clf()
plt.plot(time[1:], Temperature[1:], '-', label="Temperature vs Time")
# plt.plot(time, Power, '-', label="Power vs Time")
plt.legend()
plt.grid(True)
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.title('Temperature vs Time')
plt.show()

plt.figure(3)
plt.clf()
plt.plot(Temperature, nu, '-', label="Kinematic Viscosity vs Time")
plt.legend()
plt.grid(True)
plt.xlabel('Time [s]')
plt.ylabel('Kinematic Viscosity [m^2/s]')
plt.title('Kinematic Viscosity vs Time')
plt.show()