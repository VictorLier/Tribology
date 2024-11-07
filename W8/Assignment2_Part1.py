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
Cb = (1 - m_p) * Cp # Bearing clearance [m]

# fine milled surfaces
Ra_a = 0.4e-6   
Ra_b = 0.8e-6   
Rq_a = Ra_a * 1.25
Rq_b = Ra_b * 1.25

# ISO VG 32
iterations = 20
w_air = 1 # [m/s] air velocity of surroundings
alpha = 9.807 * (0.7 + 1.2 * w_air**(1/2))
A = 9 * D * D**(1/2)
rho = 876 # [kg/m^3]
cp = 1964 # [J/kg*K]
lamb = 1/3

nu_40 = 0.08/rho*1e6
t_40 = 40
nu_100 = 0.007/rho*1e6
t_100 = 100
m_lub = (np.log(np.log(nu_100 + 0.8)) - np.log(np.log(nu_40 + 0.8))) / (np.log((t_40+273.15)/(t_100+273.15)))
k_lub = np.log(np.log(nu_40 + 0.8)) + m_lub*np.log(t_40+273.15)
t_1 = 30
t_0 = 20
p_f = 2e5 # [Pa]
L_mark = L_norm/2
damping = 0.5

def eta_i(temp):
    return rho*1e-6*(np.exp(np.exp(-m_lub*np.log(temp+273.15) + k_lub)) - 0.8)

t_mean = np.zeros(iterations)
t = np.zeros(iterations)
t[0] = 40
N = 30
omega = 2*np.pi*N
eta = eta_i(t[0])
lub_temp = np.zeros(len(Tables))

# 1. 10*np.pi/180 * R_b = 0.00888: for groove with length L, b/a = 0.18, a/L = 3.5: Qf = 2 * 0.15  
# 2. 10*np.pi/180 * R_b = 0.00888: for groove with length L, b/a = 0.08, a/L = 1.7: Qf = 2 * 0.15
Qf = np.array([2 * 0.15, 2*0.15, 0, 0, 0, 0])

for j in range(len(Tables)):
    for i in range(iterations-1):
        table = np.flip(Tables[j] , axis=0)
        eta = eta_i(t[i])
        S_current = eta*N*L[j]*D/W*(R_b/Cp)**2
        Qs = np.interp(S_current, table[:,0], table[:,3]) 
        Qe = np.interp(S_current, table[:,0], table[:,4])
        f_J = psi * np.interp(S_current, table[:,0], table[:,5])
        epsi_current = np.interp(S_current, table[:,0], table[:,1])

        qf = 8*Cp**3/eta * pf * Qf[j]
        chi = 1 # !!!!!!!!!!!!
        q = R_b*N*2*np.pi*Cp*L[j]*(Qs + (1 - chi)*Qe) + qf

        t_new = ((1-lamb)*(f_J*R_b*W*omega + alpha*A*t_0) + Cp*rho*q*t_1) / (Cp*rho*q + alpha*A*(1-lamb))
        t[i+1] = t[i] + damping * (t_new - t[i])
    
    lub_temp[j] = t[-1]

    plt.plot(np.arange(iterations), t, label=f'{j+1}')

plt.xlabel('Iterations')
plt.ylabel('Temperature [C]')
plt.legend()
plt.show()

print(f"lubricant temperature for bearing 1: {lub_temp[0]:.3f} C and dynamic viscosity: {eta_i(lub_temp[0]):.3f} m^2/s")
print(f"lubricant temperature for bearing 2: {lub_temp[1]:.3f} C and dynamic viscosity: {eta_i(lub_temp[1]):.3f} m^2/s")
print(f"lubricant temperature for bearing 3: {lub_temp[2]:.3f} C and dynamic viscosity: {eta_i(lub_temp[2]):.3f} m^2/s")
print(f"lubricant temperature for bearing 4: {lub_temp[3]:.3f} C and dynamic viscosity: {eta_i(lub_temp[3]):.3f} m^2/s")
print(f"lubricant temperature for bearing 5: {lub_temp[4]:.3f} C and dynamic viscosity: {eta_i(lub_temp[4]):.3f} m^2/s")
print(f"lubricant temperature for bearing 6: {lub_temp[5]:.3f} C and dynamic viscosity: {eta_i(lub_temp[5]):.3f} m^2/s")
print("##############################################################")

eta = np.array([eta_i(lub_temp[i]) for i in range(len(Tables))])

# laminar flow condition
print(f"laminar flow condition, Raynolds number is less than 2400")
print(f" the maximum velocity is thereby: {2400 * eta[0] / (Cp * rho ) / (2 * np.pi * R_b):.3f} Hz for bearing 1")
print(f" the maximum velocity is thereby: {2400 * eta[1] / (Cp * rho ) / (2 * np.pi * R_b):.3f} Hz for bearing 2")
print(f" the maximum velocity is thereby: {2400 * eta[2] / (Cp * rho ) / (2 * np.pi * R_b):.3f} Hz for bearing 3")
print(f" the maximum velocity is thereby: {2400 * eta[3] / (Cp * rho ) / (2 * np.pi * R_b):.3f} Hz for bearing 4")
print(f" the maximum velocity is thereby: {2400 * eta[4] / (Cp * rho ) / (2 * np.pi * R_b):.3f} Hz for bearing 5")
print(f" the maximum velocity is thereby: {2400 * eta[5] / (Cp * rho ) / (2 * np.pi * R_b):.3f} Hz for bearing 6")
print(f" calculated using the formula: Re = 2400 = rho * v * Cp / eta")
print("##############################################################")

# 1) find minumum angular velocity - determine the speed at which the clearance minus the eccentricity is equal to the minimum film thickness
N = np.linspace(0.1, 25, 100)
Lambda_hydro = np.zeros((len(N), 6))
for j in range(len(Tables)):
    table = np.flip(Tables[j] , axis=0)
    Sommerfeld = eta[j]*N*D*L[j]/(psi**2*W)

    epsilon = np.interp(Sommerfeld, table[:,0], table[:,1])
    # h_min = - epsilon * Cb[j] + Cp
    h_min = Cp*(1 - epsilon)
    Lambda_hydro[:, j] = h_min / (Rq_a**2 + Rq_b**2)**(1/2)

    h_min_10 = 10 * (Rq_a**2 + Rq_b**2)**(1/2)
    # epsilon_10 = (Cp - h_min_10) / Cb[j]
    epsilon_10 = 1 - h_min_10 / Cp
    table = Tables[j]
    Sommerfeld_10 = np.interp(epsilon_10, table[:,1], table[:,0])
    N_10 = Sommerfeld_10 * psi**2 * W / (eta[j] * D * L[j])
    print(f"bearing {j+1} has minimum speed: {N_10:.3f} Hz")

plt.plot(N, Lambda_hydro)
plt.xlabel('Speed [Hz]')
plt.ylabel('Lambda')
plt.legend(['1', '2', '3', '4', '5', '6'])
plt.ylim(0, 50)
plt.axhline(y=10, color='r', linestyle='--')
plt.show()

print("##############################################################")

# 2) Find maximum angular velocity

N = np.linspace(10, 1000, 100)        # Speed [Hz]
M = np.array([[mass, 0], [0, mass]])    # mass matrix
eigen_RE = np.zeros((len(N), 6))           

for j in range(len(Tables)):
    table = np.flip(Tables[j], axis = 0)
    Sommerfeld = eta[j]*N*D*L[j]/(psi**2*W)
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
                    
                    
        eigen_RE[i, j] = np.real(s[np.abs(np.imag(s)) > 0][0])


        if np.any(np.real(s) > 0):         # if eigenvalues s have positive real part, the system is unstable
            print(f"bearing {j+1} is unstable at speed {N[i]:.3f} Hz with undamped natural frequency: {omega_max:.3f} in Hz")
            break
    
    if i == len(N)-1:                      # if the loop reaches the end, the system is stable at all speeds
        print(f"bearing {j+1} is stable at all speeds with maximum speed: {N[i]:.3f} Hz")

plt.figure()
plt.plot(N, eigen_RE)
plt.xlabel('Speed [Hz]')
plt.ylabel('Real part of eigenvalues')
plt.legend(['1', '2', '3', '4', '5', '6'])
plt.show()

# 3) maximum lubrications consumptions (flow rate)

N = np.linspace(1, 250, 100)        # Speed [Hz]
for j in range(len(Tables)):
    table = np.flip(Tables[j], axis = 0)
    Sommerfeld = eta[j]*N*D*L[j]/(psi**2*W)

    Qs = np.interp(Sommerfeld, table[:,0], table[:,3]) 
    Qe = np.interp(Sommerfeld, table[:,0], table[:,4])
    qf = 8*h_min**3/eta[j] * pf * Qf[j]
    chi = 1 # !!!!!!!!!!!!
    q = R_b*N*2*np.pi*Cp*L[j]*(Qs + (1 - chi)*Qe) + qf

    plt.plot(N, q)

plt.xlabel('Speed [Hz]')
plt.ylabel('Flow rate [m^3/s]')
plt.legend(['1', '2', '3', '4', '5', '6'])
plt.show()

# 4) maximum friction loss 

N = np.linspace(1, 200, 100) 
for j in range(len(Tables)):
    table = np.flip(Tables[j], axis = 0)
    Sommerfeld = eta[j]*N*D*L[j]/(psi**2*W)

    fJ = np.interp(Sommerfeld, table[:,0], table[:,5]) * psi
    H = fJ * R_b * N * 2 * np.pi * W
    plt.plot(N, H)

plt.xlabel('Speed [Hz]')
plt.ylabel('Friction loss [W]')
plt.legend(['1', '2', '3', '4', '5', '6'])
plt.show()