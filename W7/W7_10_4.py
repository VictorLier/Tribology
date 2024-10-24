
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

r_b = 25e-3   # [m]
b = 100e-3      # [m]
omega = 1000 * 2 * np.pi / 60 # [rad/s]
w_r = 5e3       # [N]
eta_0 = 0.0033 # [Pa.s]
c = 2*r_b/1000

lamb = 2*r_b/b
print(f"Since lambda = {lamb} < 1, the bearing is a short bearing")

# full Sommerfeld Solution - infinitely wide

W_r = w_r/b/(r_b*omega*eta_0)*(c/r_b)**2

eps = sp.symbols('eps', positive = True)

def equation(epsi):
    return 12*np.pi*epsi/((2+epsi**2)*(1-epsi**2)**(1/2)) - W_r

solution = fsolve(equation, 0.8)

print(solution)

def equation2(epsi):
    return 6*epsi*(np.pi - epsi**2 * (np.pi - 4))**(1/2)/((2+epsi**2)*(1-epsi**2)) - W_r

# (b/r_b)**2 * epsi / (4*(1-epsi**2)**2) * (16*epsi**2 + np.pi**2*(1-epsi**2))**(1/2)

solution2 = fsolve(equation2, 0.9)

print(solution2)





