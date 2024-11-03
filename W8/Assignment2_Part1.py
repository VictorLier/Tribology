import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from Someya import Table1, Table2, Table3, Table4, Table5, Table6
from scipy.linalg import eig
# from Someya import Tabletest
Tables = [Table1, Table2, Table3, Table4, Table5, Table6]

# Constants
mass = 460          # Mass           [kg]
R_b = 50.902e-3     # Bearing radius [m]
D = 2*50.8e-3       # Shaft diameter [m]
Cp = 102e-6         # Clearance      [m]
L_norm = 50.8e-3    # Bearing length [m]
W = mass*9.82/2     # Load           [N]
psi = Cp/R_b        # Eccentricity ratio
pf = 2e5            # Feed pressure  [N]

m_p = np.array([0, 0, 0.5, 2/3, 0, 0])  # pre load factor [-]
LD_fac = np.array([0.5, 1.0, 0.5, 1.0, 0.5, 0.5]) # L/D factor [-]
L = D * LD_fac      # Bearing length [m]

# fine milled surfaces
Ra_a = 0.4e-6   
Ra_b = 0.8e-6   
Rq_a = Ra_a * 1.25
Rq_b = Ra_b * 1.25

N = np.linspace(0.1, 10, 50)
Lambda_hydro = np.zeros((len(N), 6))

rho = 860 # [kg/m^3]
cp = 2000 # [J/kg*K]

# ISO VG 32
eta = 0.08
nu_40 = 0.08/rho*1e6
t_40 = 40
nu_100 = 0.007/rho*1e6
t_100 = 100
m_lub = (np.log(np.log(nu_100 + 0.8)) - np.log(np.log(nu_40 + 0.8))) / (np.log((t_40+273.15)/(t_100+273.15)))
k_lub = np.log(np.log(nu_40 + 0.8)) + m_lub*np.log(t_40+273.15)


# 1) find minumum angular velocity - determine the speed at which the clearance minus the eccentricity is equal to the minimum film thickness
for j in range(len(Tables)):
    table = np.flip(Tables[j] , axis=0)
    S = table[:,0]     # Sommerfeld number
    E = table[:,1]     # eccentricity ratio
    Sommerfeld = eta*N*D*L[j]/(psi**2*W)
    for i in range(len(N)):
        epsilon = np.interp(Sommerfeld[i], S, E)
        h_min = Cp * (1 - epsilon)
        Lambda_hydro[i, j] = h_min / (Rq_a**2 + Rq_b**2)**(1/2)

plt.plot(N, Lambda_hydro)
plt.xlabel('Speed [Hz]')
plt.ylabel('Lambda')
plt.legend(['Table1', 'Table2', 'Table3', 'Table4', 'Table5', 'Table6'])
plt.ylim(0, 20)
plt.axhline(y=10, color='r', linestyle='--')
# plt.show()

# 2) Find maximum angular velocity

N = np.linspace(10, 1000, 10000)        # Speed [Hz]
M = np.array([[mass, 0], [0, mass]])    # mass matrix           

for j in range(len(Tables)):
    table = np.flip(Tables[j], axis = 0)
    Sommerfeld = eta*N*D*L[j]/(psi**2*W)
    for i in range(len(N)):
        epsilon = np.interp(Sommerfeld[i], table[:,0], table[:,1])
        Phi     = np.interp(Sommerfeld[i], table[:,0], table[:,2]) * np.pi/180
        kxx     = np.interp(Sommerfeld[i], table[:,0], table[:,6]) * W/Cp
        kxy     = np.interp(Sommerfeld[i], table[:,0], table[:,7]) * W/Cp
        kyx     = np.interp(Sommerfeld[i], table[:,0], table[:,8]) * W/Cp
        kyy     = np.interp(Sommerfeld[i], table[:,0], table[:,9]) * W/Cp
        bxx     = np.interp(Sommerfeld[i], table[:,0], table[:,10]) * W/(N[i]*2*np.pi*Cp)
        bxy     = np.interp(Sommerfeld[i], table[:,0], table[:,11]) * W/(N[i]*2*np.pi*Cp)
        byx     = np.interp(Sommerfeld[i], table[:,0], table[:,12]) * W/(N[i]*2*np.pi*Cp)
        byy     = np.interp(Sommerfeld[i], table[:,0], table[:,13]) * W/(N[i]*2*np.pi*Cp)

        K = np.array([[kxx, kxy], [kyx, kyy]])  # stiffness matrix
        B = np.array([[bxx, bxy], [byx, byy]])  # damping matrix
        A1 = np.block([[M, np.zeros(M.shape)], [np.zeros(M.shape), M]])     
        A2 = np.block([[B, K], [-M, np.zeros(M.shape)]])    

        s, u = eig(-A2, A1)                # eigenvalues s and eigenvectors u
        xi = -np.real(s) / np.abs(s)       # damping factors
        wd_hz = np.imag(s) / (2 * np.pi)   # Damped natural frequencies in Hz

        if np.all(np.real(s) < 0):         # if eigenvalues s have negative real part, the system is stable
            omega_max = np.max(np.abs(np.imag(s))) / (2 * np.pi)

        if np.any(np.real(s) > 0):         # if eigenvalues s have positive real part, the system is unstable
            print(f"bearing {j+1} is unstable at speed {N[i]:.3f} Hz with undamped natural frequency: {omega_max:.3f} in Hz")
            break
    
    if i == len(N)-1:                      # if the loop reaches the end, the system is stable at all speeds
        print(f"bearing {j+1} is stable at all speeds with maximum speed: {N[i]:.3f} Hz")

# 3) maximum lubrications consumptions (flow rate)

N = np.linspace(1, 20, 100)        # Speed [Hz]
for j in range(len(Tables)):
    table = np.flip(Tables[j], axis = 0)
    Sommerfeld = eta*N*D*L[j]/(psi**2*W)

    Qs = np.interp(Sommerfeld, table[:,0], table[:,3]) 
    Qe = np.interp(Sommerfeld, table[:,0], table[:,4])
    Qf = 1 # !!!!!!!!!!!!
    qf = 8*h_min**3/eta * pf * Qf
    chi = 1 # !!!!!!!!!!!!
    q = R_b*N*2*np.pi*Cp*L[j]*(Qs + (1 - chi)*Qe) + qf

    plt.plot(N, q)

plt.xlabel('Speed [Hz]')
plt.ylabel('Flow rate [m^3/s]')
plt.legend(['Table1', 'Table2', 'Table3', 'Table4', 'Table5', 'Table6'])
plt.show()

# 4) maximum friction loss 

N = np.linspace(1, 200, 100) 
for j in range(len(Tables)):
    table = np.flip(Tables[j], axis = 0)
    Sommerfeld = eta*N*D*L[j]/(psi**2*W)

    fJ = np.interp(Sommerfeld, table[:,0], table[:,5]) * psi

    H = fJ * R_b * N * 2 * np.pi * W

    plt.plot(N, H)

plt.xlabel('Speed [Hz]')
plt.ylabel('Friction loss [W]')
plt.legend(['Table1', 'Table2', 'Table3', 'Table4', 'Table5', 'Table6'])
plt.show()





