
# Problem 2

import math
from math import exp, log

p_Smax = 1e9 # [Pa] - maximum pressure in the EHL contact
p_Hmin = 1e6 # [Pa] - minimum pressure in the EHL contact
eta_38 = 0.2 # [Pa.s] - viscosity of lubricant at 38 degree celcius
eta_99 = 0.0125 # [Pa.s] - viscosity of lubricant at 99 degree celcius
xi_38 = 1.77e-8 # [m^2/N] - pressure-viscosity coefficient at 38 degree celcius
xi_99 = 1.51e-8 # [m^2/N] - pressure-viscosity coefficient at 99 degree celcius
eta_inf = 6.31e-5 # [Pa.s] - viscosity at elevated pressures and temperatures
c_p = 1.96e8 # [Pa] - pressure-viscosity coefficient

# 1 - viscosity at elevated pressures and temperatures

# Barus equation:
# eta = eta_0 * exp(xi*p)

# at 38 degree celcius
eta_38S = eta_38 * exp(xi_38*p_Smax)
eta_38H = eta_38 * exp(xi_38*p_Hmin)

# at 99 degree celcius
eta_99S = eta_99 * exp(xi_99*p_Smax)
eta_99H = eta_99 * exp(xi_99*p_Hmin)

print(f"The viscosity of the lubricant at 38 degree celcius and maximum pressure is: {eta_38S:.3f} Pa.s")   
print(f"The viscosity of the lubricant at 38 degree celcius and minimum pressure is: {eta_38H:.3f} Pa.s")   
print(f"The viscosity of the lubricant at 99 degree celcius and maximum pressure is: {eta_99S:.3f} Pa.s")   
print(f"The viscosity of the lubricant at 99 degree celcius and minimum pressure is: {eta_99H:.3f} Pa.s")   
print("-----------------------------------")
# Roelands equation:
# eta = eta_0 * (eta_inf/eta_0)^(1 - (1 + p/c_p)^Z_1)
# Where Z_1 = xi/(5.1e-9*(ln(eta_0) + 9.67))

# at 38 degree celcius
Z_1_38 = xi_38/(5.1e-9*(log(eta_38) + 9.67))
eta_38S = eta_38 * (eta_inf/eta_38)**(1 - (1 + p_Smax/c_p)**Z_1_38)
eta_38H = eta_38 * (eta_inf/eta_38)**(1 - (1 + p_Hmin/c_p)**Z_1_38)

# at 99 degree celcius
Z_1_99 = xi_99/(5.1e-9*(log(eta_99) + 9.67))
eta_99S = eta_99 * (eta_inf/eta_99)**(1 - (1 + p_Smax/c_p)**Z_1_99)
eta_99H = eta_99 * (eta_inf/eta_99)**(1 - (1 + p_Hmin/c_p)**Z_1_99)

print(f"The viscosity of the lubricant at 38 degree celcius and maximum pressure is: {eta_38S:.3f} Pa.s")
print(f"The viscosity of the lubricant at 38 degree celcius and minimum pressure is: {eta_38H:.3f} Pa.s")
print(f"The viscosity of the lubricant at 99 degree celcius and maximum pressure is: {eta_99S:.3f} Pa.s")
print(f"The viscosity of the lubricant at 99 degree celcius and minimum pressure is: {eta_99H:.3f} Pa.s")

#  When does it matter to introduce these corrections to the viscosity

# When the pressure is high