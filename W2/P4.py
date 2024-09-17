
# Oil properties: ISO VG 160 oil at 40°C
nu = 160 # [cSt] - kinematic viscosity
rho = 870 # [kg/m^3] - density

# Constants
q_mark = 11 # [1/(min*m)] - inflow per unit width
w_z_mark = 10000 # [N/m] - load per unit width
L1 = 60e-3 # [m] - length of the bearing
L2 = 40e-3 # [m] - length of the cavity

# 2 - compute maximum pressure to sustain the load

# NII - vertical direction:
# L2*p_max + (L1-L2)*p_max/2 - w_z_mark = 0
# => p_max = w_z_mark/(L2 + (L1-L2)/2)
p_max = w_z_mark/(L2 + (L1-L2)/2)

print(f"Maximum pressure to sustain the load is {p_max:.2f} [Pa]")
print(f"Maximum pressure to sustain the load [bar] is {p_max/1e5:.2f} [bar]")

# 3 - Determine fluid film heigh in the two small passages necessary to maintain the load capacity

# The flow is a Couette flow, so the flow rate is given by:
# q_mark/2 = -h**3/(12*eta_0)(dp/dx)
# with dp/dx = p_max/((L1-L2)/2)
# and viscocity
eta_0 = nu*rho*1e-6 # [Pa.s] - dynamic viscocity
Gradiant = -2*p_max/(L1-L2)
q_mark = q_mark/60 # converting so SI units
# For the fluid flow out of ine side of the bearing, the flow rate is:
h = (-q_mark/(2*Gradiant)*12*eta_0)**(1/3)

print(f"Fluid film height in the two small passages is {h:.2e} [m]")

A = (3*q_mark*eta_0*(L1-L2)/p_max)**(1/3)
print(A)

# WHY IS THE RESULT NOT RIGHT!!!!!!!!!!