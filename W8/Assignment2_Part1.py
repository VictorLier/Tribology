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
pf = 1e5            # Feed pressure  [N]

m_p = np.array([0, 0, 0.5, 2/3, 0, 0])  # pre load factor [-]
LD_fac = np.array([0.5,  1.0, 0.5, 1.0, 0.5, 0.5]) # L/D factor [-]
L = D * LD_fac      # Bearing length [m]
Cb = (1 - m_p) * Cp # Bearing clearance [m]
Qf = np.array([2 * 0.15, 2*0.15, 2 * 0.15, 2*0.15, 4*0.15, 4*0.15])

# fine milled surfaces
Rq_a = 0.4e-6 * 1.25
Rq_b = 0.8e-6 *1.25

# ISO VG 32
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
damping = 0.4
chi = 1 # !

# Calculate the lubricant viscosity
def eta_i(temp):
    return rho*1e-6*(np.exp(np.exp(-m_lub*np.log(temp+273.15) + k_lub)) - 0.8)

# Calculate the lubricant temperature
def calculate_lub_temp(j, N, plotswitch):
    t = [50]  # Initial temperature
    omega = 2 * np.pi * N
    table = np.flip(Tables[j], axis=0)
    lub_temp = np.zeros(len(Tables))
    
    i = 0
    while True:
        eta = eta_i(t[i])
        S_current = eta * N * L[j] * D / W * (R_b / Cp) ** 2
        Qs = np.interp(S_current, table[:, 0], table[:, 3])
        # Qe = np.interp(S_current, table[:, 0], table[:, 4])
        f_J = psi * np.interp(S_current, table[:, 0], table[:, 5])
        epsi_current = np.interp(S_current, table[:, 0], table[:, 1])
        h1 = Cp - epsi_current * Cb[j]
        qf = 8 * h1 ** 3 / eta * pf * Qf[j]
        q = R_b * omega * Cp * L[j] * Qs + qf
        t_new = ((1 - lamb) * (f_J * R_b * W * omega + alpha * A * t_0) + cp * rho * q * t_1) / (cp * rho * q + alpha * A * (1 - lamb))
        t.append(t[i] + damping * (t_new - t[i]))

        if i > 50:
            t[-1] = np.mean(t[-2:])
            break
        if np.abs((t[i + 1] - t[i])/t[i]) < 1e-3:
            break
        i += 1
    lub_temp = t[-1]
    if plotswitch == True:
        plt.plot(np.arange(len(t)), t, label=f'{j + 1}')
        np.savetxt('ASSIGN2_AUST/temp' + str(j) + '.txt', np.array([np.arange(len(t)),t]).T)
    return lub_temp

# Example usage: save txt file
plt.figure()
lub_temp = np.zeros(len(Tables))
for j in range(len(Tables)):
    N = 30
    lub_temp[j] = calculate_lub_temp(j, N, True)
plt.xlabel('Iterations')
plt.ylabel('Temperature [C]')
plt.legend()
plt.show()

print(f"lubricant temperature for bearing 1: {lub_temp[0]:.3f} C and dynamic viscosity: {eta_i(lub_temp[0]):.3f} m^2/s for N = 30 Hz")
print(f"lubricant temperature for bearing 2: {lub_temp[1]:.3f} C and dynamic viscosity: {eta_i(lub_temp[1]):.3f} m^2/s for N = 30 Hz")
print(f"lubricant temperature for bearing 3: {lub_temp[2]:.3f} C and dynamic viscosity: {eta_i(lub_temp[2]):.3f} m^2/s for N = 30 Hz")
print(f"lubricant temperature for bearing 4: {lub_temp[3]:.3f} C and dynamic viscosity: {eta_i(lub_temp[3]):.3f} m^2/s for N = 30 Hz")
print(f"lubricant temperature for bearing 5: {lub_temp[4]:.3f} C and dynamic viscosity: {eta_i(lub_temp[4]):.3f} m^2/s for N = 30 Hz")
print(f"lubricant temperature for bearing 6: {lub_temp[5]:.3f} C and dynamic viscosity: {eta_i(lub_temp[5]):.3f} m^2/s for N = 30 Hz")
print("##############################################################")

# converged speeds for Re = 1000
def find_converged_speed(j, target_Re=1000, tol=1e-3, max_iter=100):
    N_low, N_high = 0.1, 1000  # Initial bounds for speed [Hz]
    for _ in range(max_iter):
        N_mid = (N_low + N_high) / 2
        eta_mid = eta_i(calculate_lub_temp(j, N_mid, False))
        Re_mid = rho * N_mid * 2 * np.pi * R_b * Cp / eta_mid

        if abs(Re_mid - target_Re) < tol:
            return N_mid

        if Re_mid < target_Re:
            N_low = N_mid
        else:
            N_high = N_mid

    return N_mid  # Return the best estimate if convergence is not reached

# Find converged speeds
converged_speeds = np.array([find_converged_speed(j) for j in range(len(Tables))])

# laminar flow condition
print(f"laminar flow condition, Raynolds number is less than 1000")
for j in range(len(Tables)):
    print(f" the maximum velocity is thereby: {converged_speeds[j]:.3f} Hz for bearing {j + 1}")
print(f" calculated using the formula: Re = 1000 = rho * v * Cp / eta")
print("##############################################################")

# 1) find minumum angular velocity - determine the speed at which the clearance minus the eccentricity is equal to the minimum film thickness
N = np.linspace(0, 5, 100)
Lambda_hydro = np.zeros((len(N), 6))
for j in range(len(Tables)):
    table = np.flip(Tables[j] , axis=0)
    eta = np.zeros(len(N))
    for i in range(len(N)):
        eta[i] = eta_i(calculate_lub_temp(j, N[i], False))
    Sommerfeld = eta*N*D*L[j]/(psi**2*W)

    # for the plot
    epsilon = np.interp(Sommerfeld, table[:,0], table[:,1])
    # h_min = - epsilon * Cb[j] + Cp
    h_min = Cp*(1 - epsilon)
    Lambda_hydro[:, j] = h_min / (Rq_a**2 + Rq_b**2)**(1/2)

    # finding the exact speed where lambda = 10 using interpolation
    N_10 = np.interp(10, Lambda_hydro[:, j], N)
    # h_min_10 = 10 * (Rq_a**2 + Rq_b**2)**(1/2)
    # # epsilon_10 = (Cp - h_min_10) / Cb[j]
    # epsilon_10 = 1 - h_min_10 / Cp
    # table = Tables[j]
    # Sommerfeld_10 = np.interp(epsilon_10, table[:,1], table[:,0])
    # N_10 = Sommerfeld_10 * psi**2 * W / (eta * D * L[j])
    print(f"bearing {j+1} has minimum speed: {N_10:.3f} Hz")

plt.figure()
plt.plot(N, Lambda_hydro)
plt.xlabel('Speed [Hz]')
plt.ylabel('Lambda')
plt.legend(['1', '2', '3', '4', '5', '6'])
plt.ylim(0, 50)
plt.axhline(y=10, color='r', linestyle='--')
plt.show()

# save the data
np.savetxt('ASSIGN2_AUST/Nmin1.txt', np.array([N, Lambda_hydro[:,0]]).T)
np.savetxt('ASSIGN2_AUST/Nmin2.txt', np.array([N, Lambda_hydro[:,1]]).T)
np.savetxt('ASSIGN2_AUST/Nmin3.txt', np.array([N, Lambda_hydro[:,2]]).T)
np.savetxt('ASSIGN2_AUST/Nmin4.txt', np.array([N, Lambda_hydro[:,3]]).T)
np.savetxt('ASSIGN2_AUST/Nmin5.txt', np.array([N, Lambda_hydro[:,4]]).T)
np.savetxt('ASSIGN2_AUST/Nmin6.txt', np.array([N, Lambda_hydro[:,5]]).T)

print("##############################################################")

# 2) Find maximum angular velocity - determine the speed at which the system becomes unstable
N = np.linspace(10, 2000, 1000)        # Speed [Hz]
M = np.array([[mass/2, 0], [0, mass/2]])    # mass matrix
eigen_RE1 = np.zeros((len(N), 6))
eigen_RE2 = np.zeros((len(N), 6))
eigen_RE3 = np.zeros((len(N), 6))
eigen_RE4 = np.zeros((len(N), 6))
eigen = np.zeros((len(N), 6, 4))
N_stability = np.zeros(6)
for j in range(len(Tables)):
    eigenswitch = 0
    table = np.flip(Tables[j], axis = 0)
    eta = np.zeros(len(N))
    for i in range(len(N)):
        eta = eta_i(calculate_lub_temp(j, N[i], False))
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
            
        # eigen_RE[i, j] = np.real(s[np.abs(np.imag(s)) > 0][0])

        # make a vector eigen_Re with the real part of the first eigenvalue
        s = np.sort(np.real(s))
        eigen_RE1[i, j] = np.real(s[0])
        eigen_RE2[i, j] = np.real(s[1])
        eigen_RE3[i, j] = np.real(s[2])
        eigen_RE4[i, j] = np.real(s[3])
        eigen[i, j, :] = np.real(s)

        if np.any(np.real(s) > 0) and eigenswitch == 0:         # if eigenvalues s have positive real part, the system is unstable
            eigenswitch = 1
            print(f"bearing {j+1} is unstable at speed {N[i]:.3f} Hz with undamped natural frequency: {omega_max:.3f} in Hz")
            N_stability[j] = N[i]
    
    if i == len(N)-1 and eigenswitch == 0:                      # if the loop reaches the end, the system is stable at all speeds
        print(f"bearing {j+1} is stable at all speeds with maximum speed: {N[i]:.3f} Hz")
        N_stability[j] = N[i]

# make a figure with 4 subplots, each subplot contains the real part of the first 4 eigenvalues
plt.figure()
plt.subplot(2, 2, 1)
plt.plot(N, eigen_RE1)
plt.hlines(0, 0, 2000, colors='r', linestyles='--')
plt.xlabel('Speed [Hz]')
plt.ylabel('Real part of eigenvalues')
plt.ylim(-100, 100)
plt.legend(['1', '2', '3', '4', '5', '6'])

plt.subplot(2, 2, 2)
plt.plot(N, eigen_RE2)
plt.hlines(0, 0, 2000, colors='r', linestyles='--')
plt.xlabel('Speed [Hz]')
plt.ylabel('Real part of eigenvalues')
plt.ylim(-100, 100)
plt.legend(['1', '2', '3', '4', '5', '6'])

plt.subplot(2, 2, 3)
plt.plot(N, eigen_RE3)
plt.hlines(0, 0, 2000, colors='r', linestyles='--')
plt.xlabel('Speed [Hz]')
plt.ylabel('Real part of eigenvalues')
plt.ylim(-100, 100)
plt.legend(['1', '2', '3', '4', '5', '6'])

plt.subplot(2, 2, 4)
plt.plot(N, eigen_RE4)
plt.hlines(0, 0, 2000, colors='r', linestyles='--')
plt.xlabel('Speed [Hz]')
plt.ylabel('Real part of eigenvalues')
plt.ylim(-100, 100)
plt.legend(['1', '2', '3', '4', '5', '6'])
plt.show()

plt.figure()
plt.plot(N, eigen[:,0,:], color = 'b', label='1')
plt.plot(N, eigen[:,1,:], color = 'orange', label='2')
plt.plot(N, eigen[:,2,:], color = 'g', label='3')
plt.plot(N, eigen[:,3,:], color = 'r', label='4')
plt.plot(N, eigen[:,4,:], color = 'purple', label='5')
plt.plot(N, eigen[:,5,:], color = 'brown', label='6')
plt.legend()
plt.xlabel('Speed [Hz]')
plt.ylabel('Real part of eigenvalues')
plt.hlines(0, 0, 2000, colors='r', linestyles='--')
plt.ylim(-100, 100)
plt.show()

# save the data
for i in range(6):
    for k in range(4):
        # remove datapoints from eigen that are above 100 and below -100, and make a N that matches the new eigen
        eigenp = eigen[(eigen[:,i,k] < 50) & (eigen[:,i,k] > -100), i, k]
        Np = N[(eigen[:,i,k] < 50) & (eigen[:,i,k] > -100)]
        np.savetxt(f'ASSIGN2_AUST/eigen{i+1}{k+1}.txt', np.array([Np[::2], eigenp[::2]]).T)

# 3) maximum lubrications consumptions (flow rate)
for j in range(len(Tables)):
    N_max = np.min([converged_speeds[j], N_stability[j]])
    if N_max > 500:
        N_max = 500
    N = np.linspace(0, N_max, 100)       # Speed [Hz]
    table = np.flip(Tables[j], axis = 0)
    eta = np.zeros(len(N))
    for i in range(len(N)):
        eta = eta_i(calculate_lub_temp(j, N[i], False))
    Sommerfeld = eta*N*D*L[j]/(psi**2*W)

    Qs = np.interp(Sommerfeld, table[:,0], table[:,3]) 
    Qe = np.interp(Sommerfeld, table[:,0], table[:,4])
    epsilon = np.interp(Sommerfeld, table[:,0], table[:,1])
    h1 = Cp - epsilon * Cb[j]
    qf = 8*h1**3/eta * pf * Qf[j]
    chi = 1 # !!!!!!!!!!!!
    q = R_b*N*2*np.pi*Cp*L[j]*(Qs + (1 - chi)*Qe) + qf

    plt.plot(N, q*60*60*1000)
    np.savetxt('ASSIGN2_AUST/flowrate' + str(j+1) + '.txt', np.array([N, q*60*60*1000]).T)

plt.xlabel('Speed [Hz]')
plt.ylabel('Flow rate [L/h]')
plt.legend(['1', '2', '3', '4', '5', '6'])
plt.show()

# 4) maximum friction loss 
for j in range(len(Tables)):
    N_max = np.min([converged_speeds[j], N_stability[j]])
    if N_max > 500:
        N_max = 500
    N = np.linspace(0, N_max, 100)       # Speed [Hz]
    table = np.flip(Tables[j], axis = 0)
    eta = np.zeros(len(N))
    for i in range(len(N)):
        eta = eta_i(calculate_lub_temp(j, N[i], False))
    Sommerfeld = eta*N*D*L[j]/(psi**2*W)

    fJ = np.interp(Sommerfeld, table[:,0], table[:,5]) * psi
    H = fJ * R_b * N * 2 * np.pi * W
    np.savetxt('ASSIGN2_AUST/loss' + str(j+1) + '.txt', np.array([N, H*1e-3]).T)
    plt.plot(N, H*1e-3)

plt.xlabel('Speed [Hz]')
plt.ylabel('Friction loss [kW]')
plt.legend(['1', '2', '3', '4', '5', '6'])
plt.show()