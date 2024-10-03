
import numpy as np
import matplotlib.pyplot as plt

# Oil properties

# # ISO VG 100 - Harmony® AW Hydraulic Oil
# nu_100_1 = 11.8 # [cSt] - kinematic viscosity at 100 degrees Celsius
# nu_40_1 = 97.6 # [cSt] - kinematic viscosity at 40 degrees Celsius
# nu_index_1 = 110 # viscosity index
# # ISO VG 100 - Harmony® R&O Oil 
# nu_100_2 = 11.4 # [cSt] - kinematic viscosity at 100 degrees Celsius
# nu_40_2 = 97.0 # [cSt] - kinematic viscosity at 40 degrees Celsius
# nu_index_2 = 100 # viscosity index

# ISO VG 32 - Capella® WF 32
nu_40_1 = 29.5
nu_100_1 = 4.37
VI_1 = 7

# ISO VG 32 - PureBlu Hydraulic Oil
nu_40_2 = 32
nu_100_2 = 5.4
VI_2 = 100

# define the temperature range and kinematic viscosity for both oils
T_1 = np.array([40, 100]) + 273.15  # Temperature in Kelvin
v_1 = np.array([nu_40_1, nu_100_1])        # Kinematic viscosity in cSt
T_2 = np.array([40, 100]) + 273.15  # Temperature in Kelvin
v_2 = np.array([nu_40_2, nu_100_2])        # Kinematic viscosity in cSt

# Plot using log(T) and log(log(v))
plt.figure(1)
plt.clf()
# plt.plot(np.log(T_1), np.log(np.log(v_1)), '-o', label="Harmony® AW Hydraulic Oil")
# plt.plot(np.log(T_2), np.log(np.log(v_2)), '-o', label="Harmony® R&O Oil")
plt.grid(True)
labels = np.array([5, 10, 32, 100, 500, 1200])
Ytick = np.log(np.log(labels))
plt.gca().set_yticks(Ytick)
plt.gca().set_yticklabels(labels)
xlabels = np.arange(-20, 101, 10) + 273  # Temperature labels from -20 to 100 with a step of 10
Xtick = np.log(xlabels)
plt.gca().set_xticks(Xtick)
plt.gca().set_xticklabels(xlabels)
plt.axis([np.log(273.15-10), np.log(273.15+100), np.log(np.log(5)), np.log(np.log(1200))])
# plt.axis([np.log(273.15+40), np.log(273.15+100), np.log(np.log(15)), np.log(np.log(2500))])
plt.xlabel('Temperature [K]')
plt.ylabel('Kinematic viscosity [cSt = mm^2/s]')
plt.title('Viscosity index')

# 2 - fit the data
modify_nu_1 = np.log(np.log(v_1+0.8))
modify_T_1 = np.log(T_1)
a_1 = np.polyfit(modify_T_1, modify_nu_1, 1)
m_1 = a_1[0]
C_1 = a_1[1]
print(f"Following the ASTM walther relation, the constants for Capella® WF 32 are: m = {m_1:.3f} and C = {C_1:.3f}")
m_simple_1 = (np.log(np.log(v_1[1]+0.8)) - np.log(np.log(v_1[0]+0.8)))/(np.log(T_1[1])-np.log(T_1[0]))
C_simple_1 = np.log(np.log(v_1[0]+0.8)) - m_simple_1*np.log(T_1[0])
print(f"Following the simple relation, the constants for Capella® WF 32 are: m = {m_simple_1:.3f} and C = {C_simple_1:.3f}")
print("-----")
modify_nu_2 = np.log(np.log(v_2+0.8))
modify_T_2 = np.log(T_2)
a_2 = np.polyfit(modify_T_2, modify_nu_2, 1)
m_2 = a_2[0]
C_2 = a_2[1]
print(f"Following the ASTM walther relation, the constants for PureBlu Hydraulic Oil are: m = {m_2:.3f} and C = {C_2:.3f}")
m_simple_2 = (np.log(np.log(v_2[1]+0.8)) - np.log(np.log(v_2[0]+0.8)))/(np.log(T_2[1])-np.log(T_2[0]))
C_simple_2 = np.log(np.log(v_2[0]+0.8)) - m_simple_2*np.log(T_2[0])
print(f"Following the simple relation, the constants for PureBlu Hydraulic Oil are: m = {m_simple_2:.3f} and C = {C_simple_2:.3f}")

T_line = np.linspace(273.15-10, 273.15+100, 100)
plt.plot(np.log(T_line), m_simple_1*np.log(T_line) + C_simple_1, '--', label="Capella® WF 32 - ASTM Walther - VI = 7")
plt.plot(np.log(T_line), m_simple_2*np.log(T_line) + C_simple_2, '--', label="PureBlu Hydraulic Oil - ASTM Walther - VI = 100")
plt.legend()
plt.show()

# find the viscocity ranges for the temperature range [10;60] degrees Celsius
T = np.array([10, 60]) + 273.15
v_1 = np.exp(np.exp(m_simple_1*np.log(T) + C_simple_1))
v_2 = np.exp(np.exp(m_simple_2*np.log(T) + C_simple_2))
print(f"Capella® WF 32 - viscosity range: {v_1}")
print(f"PureBlu Hydraulic Oil - viscosity range: {v_2}")


