import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
import scipy.sparse as sps
import copy
from scipy.interpolate import interp1d

class Bearing:
    def __init__(self, bearingtype:int = 0, PadLoes:int=0, ShaftDiameter:float=101.6e-3, RadialClearance:float=102e-6, preLoadFactor:float=0, LengthDiameterFactor:float=0.5, Load:float=460/2, surfaceRoughness:float=1.2e-6, oil:int=0, Revolutions:float=70, feedPressure:float = 2e5) -> None:
        '''
        Creates a bearing object with the given parameters

        Args:
            bearingtype (int): 0 for cylindrical, 1 for two lobes, 2 for pad
            PadLoes (int): 0 for LBP, 1 for LOP
            shaftDiameter (float): The diameter of the shaft in meters. Variable name: D
            radialClearance (float): The radial clearance in meters. Variable name: C_p
            preLoadFactor (float): The preload factor. Either 0, 0.5, or 2/3. Variable name: mp
            lengthDiameterFactor (float): The length diameter factor. Either 0.5 or 1. Variable name: L
            load (float): The load on one bearing in Kg
            surfaceRoughness (float): The surface roughness in meters. Variable name: Ra - See table 3.2
            oil (int): The oil type. 0 for ISO VG 32, 1 for ISO VG 46
            Revolutions (float): Rotational speed in Hz
            feedPressure (float): The feed pressure in Pa

        Raises:
            ValueError: If the preLoadFactor is not 0, 0.5, or 2/3
            ValueError: If the lengthDiameterFactor is not 0.5 or 1
            ValueError: If the oil type is not 0 or 1
        '''
        self.bearingtype = bearingtype
        self.PadLoes = PadLoes
        self.D = ShaftDiameter
        self.r = ShaftDiameter/2
        self.C_p = RadialClearance
        self.C_b = (1 - preLoadFactor) * RadialClearance

        if preLoadFactor not in [0, 0.5, 2/3]:
            raise ValueError("Invalid preload factor")
        else:
            self.mp = preLoadFactor
        
        if LengthDiameterFactor not in [0.5, 1]:
            raise ValueError("Invalid length diameter factor")
        else:
            self.L = LengthDiameterFactor * self.D
            self.LD = LengthDiameterFactor

        self.mass = Load
        self.W = self.mass * 9.82 # Convert to N
        self.R_a = surfaceRoughness
        
        if oil == 0: # ISO VG 32
            self.nu_40 = 32e-6 # [m^2/s] - Kinematic viscosity - ISO 3448
            self.nu_100 = 5.4e-6 # [m^2/s] - Kinematic viscosity - https://www.engineeringtoolbox.com/iso-grade-oil-d_1207.html
            self.rho = 857 # [kg/m^3] - Density -  https://www.engineeringtoolbox.com/iso-grade-oil-d_1207.html
            self.eta_40 = self.nu_40 * self.rho  # [Pa*s] - Absolute viscosity
            self.CP = 1964 # [J/kg*K] - Specific heat capacity
        
        elif oil == 1: # ISO VG 46
            self.nu_40 = 46e-6 # [m^2/s] - Kinematic viscosity - ISO 3448
            self.nu_100 = 6.8e-6 # [m^2/s] - Kinematic viscosity - https://www.engineeringtoolbox.com/iso-grade-oil-d_1207.html
            self.rho = 861 # [kg/m^3] - Density -  https://www.engineeringtoolbox.com/iso-grade-oil-d_1207.html
            self.eta_40 = self.nu_40 * self.rho  # [Pa*s] - Absolute viscosity
            self.CP = 1964 # [J/kg*K] - Specific heat capacity

        else:
            raise ValueError("Invalid oil type")
        
        self.N = Revolutions
        self.omega = 2*np.pi*Revolutions

        self.nu = self.nu_40
        self.T = 40 # starting temperature [c]
        self.pf = feedPressure
        self.psi = self.C_p / (self.D/2)

        self.U = self.N * self.D / 2 * np.pi # [m/s] - Surface speed


    def get_someya(self) -> None:
        '''
        Finds the someya parameters table for the bearing

        Attributes:
            table (np.ndarray): The someya parameters table

        Raises:
            ValueError: If the bearing type is invalid
            Valuerror: If the length diameter factor is invalid
        '''
        if self.bearingtype == 0:
            self.Qf = 2 * 0.15  # Approximation
            if self.LD == 0.5:
                from Someya import Table1 as table
                self.table = table
            elif self.LD == 1:
                from Someya import Table2 as table
                self.table = table
            else:
                raise ValueError("Invalid length diameter factor")
            
        elif self.bearingtype == 1:
            self.Qf = 2 * 0.15  # Approximation
            if self.LD == 0.5 and self.mp == 0.5:
                from Someya import Table3 as table
                self.table = table
            elif self.LD == 1 and self.mp == 2/3:
                from Someya import Table4 as table
                self.table = table
            else:
                raise ValueError("Invalid length diameter factor or preload factor")
            
        elif self.bearingtype == 2:
            self.Qf = 4 * 0.15  # Approximation
            if self.LD == 0.5 and self.mp == 0 and self.PadLoes == 0:
                from Someya import Table5 as table
                self.table = table
            elif self.LD == 0.5 and self.mp == 0 and self.PadLoes == 1:
                from Someya import Table6 as table
                self.table = table
            else:
                raise ValueError("Invalid length diameter factor or preload factor or pad load")
            
        else:
            raise ValueError("Invalid bearing type")


    def someya(self) -> None:
        '''
        Finds the someya parameters for the bearing

        Attributes:
            S (np.ndarray): Sommerfeld number
            E (np.ndarray): eccentricity ratio
            Phi (np.ndarray): attitude angle (deg)
            Q (np.ndarray): flow rate ratio
            P (np.ndarray): 
            T (np.ndarray):
        '''
        # self.table = np.flip(self.table, axis=0)

        self.S_table = self.table[:,0]     # Sommerfeld number
        self.E_table = self.table[:,1]     # eccentricity ratio
        self.Phi_table = self.table[:,2]   # attitude angle (deg)
        self.Q_table = self.table[:,3]     # flow rate ratio
        self.P_table = self.table[:,4]     
        self.T_table = self.table[:,5]
        self.kx_table = self.table[:,6]
        self.kxy_table = self.table[:,7]
        self.kyx_table = self.table[:,8]
        self.kyy_table = self.table[:,9]
        self.Bxx_table = self.table[:,10]
        self.Bxy_table = self.table[:,11]
        self.Byx_table = self.table[:,12]
        self.Byy_table = self.table[:,13]


    def temp_constants(self, printbool: bool = False) -> None:
        '''
        Finds the constants for the temperature calculation

        Args:
            printbool (bool): If True, the constants are printed
        
        Attributes:
            m (float): The constant m
            k (float): The constant k
        '''
        self.m = (np.log(np.log(self.nu_100*1e6+0.8)) - np.log(np.log(self.nu_40*1e6+0.8))) / (np.log(40+273.15) - np.log(100+273.15))
        self.k = np.log(np.log(self.nu_40*1e6+0.8)) + self.m * np.log(40+273.15)

        if printbool:
            print(f"The constants m is {self.m:.3g} and k is {self.k:.3g}")


    def temp_relation(self, printbool: bool = False) -> float:
        '''
        Calculates the absolute viscosity of the oil at the current temperature

        Args:
            T (float): The temperature in Celsius
            printbool (bool): If True, the constants are printed
        
        Returns:
            float: The absolute viscosity of the oil at the given temperature [m^2/s]
        '''
        self.eta = self.rho*1e-6 * (np.exp(np.exp(-self.m * np.log(self.T+273.15) + self.k)) - 0.8) # [Pa*s] - Someya (35)

        if printbool:
            print(f"The absolute viscosity of the oil at {self.T} K is {self.eta:.3g} m^2/s")


    def sommerfeld(self, printbool:bool = False) -> None:
        '''
        Calculates the sommerfeld number for the bearing

        Args:
            printbool (bool): If True, the sommerfeld number is printed

        Attributes:
            S (float): Sommerfeld number
        '''
        self.S = self.eta * self.N * self.L * self.D / self.W * ((self.r)/self.C_p)**2 # Someya (6) - (XIII)

        if printbool:
            print(f"The sommerfeld number is {self.S:.3g}")


    def find_temp_visc(self, plot:bool = False, printbool:bool = False, convergenceTol:float = 1e-6, maxConvergence:int=1000) -> None:
        '''
        Finds the temperature and viscusity of the oil in the bearing by iterating through the temperature and viscosity

        Args:
            plot (bool): If True, the temperature and viscosity are plotted
            printbool (bool): If True, the resulting temperature and viscosity are printed
            convergenceTol (float): The convergence tolerance
            maxConvergence (int): The maximum number of iterations

        Attributes:
            T (float): The temperature of the oil in the bearing
            nu (float): The kinematic viscosity of the oil in the bearing
        
        Raises:
            ValueError: If the temperature does not converge within the maximum number of iterations
        '''
        self.temp_constants() # Find constants
        A = 9 * self.D * np.sqrt(self.D) # [m^2] - Area of the bearing - Someya (26)
        w = 1 # [m/s] sourrounding air - assumed 1
        alpha = 9.807 * (0.7 + 1.2 * np.sqrt(w)) # [W/m^2*K] - Heat transfer coefficient - Someya (25)

        lamb = 1/3 # S8-27 (31)

        t_0 = 20#+273.15 # [K] - sourrounding air temp - assumed 20
        t_1 = 30#+273.15 # [K] - oil inlet temp - assumed 30


        iteration = np.arange(maxConvergence+1)
        T = np.zeros(maxConvergence+1)

        T[0] = self.T # Initial temperature

        for i in range(maxConvergence):
            self.T = T[i]
            self.temp_relation()
            self.sommerfeld()
            Qs = interp1d(self.S_table, self.Q_table, kind='linear', fill_value=(self.Q_table[-1], self.Q_table[0]), bounds_error=False)(self.S) # From table

            f_J = self.psi * interp1d(self.S_table, self.T_table, kind='linear', fill_value=(self.T_table[-1], self.T_table[0]), bounds_error=False)(self.S) # From table
            epsi = interp1d(self.S_table, self.E_table, kind='linear', fill_value=(self.E_table[-1], self.E_table[0]), bounds_error=False)(self.S) # From table

            h1 = self.C_p - epsi * self.C_b
            qf = 8 * h1**3 / self.eta * self.pf * self.Qf
            q = self.r * self.omega * self.C_p * self.L * Qs + qf

            t_new = ((1-lamb) * (f_J * self.r * self.W * self.omega + alpha * A * t_0) + self.CP * self.rho * q * t_1) / (self.CP * self.rho * q + alpha * A * (1-lamb))
            T[i+1] = T[i] + 0.1 * (t_new - T[i])
            self.T = T[i+1]

            if abs(t_new - T[i+1]) < convergenceTol:
                self.temp_relation()
                self.nu = self.eta / self.rho

                if printbool:
                    print(f"The sommerfeld number is {self.S:.3g}")
                    print(f"The temperature of the oil is {self.T:.5g} K and the kinematic viscosity is {self.nu:.5g} m^2/s - Found in {i} iterations")

                break

        else:
            raise ValueError("Temperature did not converge within the maximum number of iterations")

        if plot:
            T = T[:i+1]
            iteration = iteration[:i+1] 
            plt.figure()
            plt.plot(iteration, T)
            plt.xlabel("Iteration")
            plt.ylabel("Temperature [K]")
            plt.title("Temperature convergence")
            plt.show()


    def reynolds_number(self, printbool:bool = False, laminarRE:float = 1500) -> None:
        '''
        Calculates the reynolds number for the bearing and checks if the flow is laminar

        Args:
            printbool (bool): If True, the reynolds number is printed
            laminarRE (float): The maximum reynolds number for the flow to be considered laminar
        
        Attributes:
            Re (float): Reynolds number
            laminar (bool): If True, the flow is laminar
        '''
        # self.Re = self.C_p * self.U / self.nu # Someya (XIII)

        self.Re = self.rho * self.omega * self.C_p * self.r / self.eta


        if self.Re < laminarRE:
            self.laminar = True
        else:
            self.laminar = False
        
        if printbool:
            print(f"The reynolds number is {self.Re:.3g} - The flow is {'laminar' if self.laminar else 'turbulent'}")


    def minimum_film_thickness(self, printbool:bool = False) -> None:
        '''
        Finds the minimum film thickness for the bearing

        Args:
            printbool (bool): If True, the minimum film thickness is printed

        Attributes:
            h_min (float): Minimum film thickness
        '''
        epsilon = interp1d(self.S_table, self.E_table, kind='linear', fill_value='extrapolate')(self.S)

        self.h_min = self.C_p * (1 - epsilon)

        if printbool:
            print(f"The minimum film thickness is {self.h_min:.3g} m")


    def film_parameter(self, printbool:bool = False, Hydro_limit:float = 10) -> None:
        '''
        Finds the non-dimensional film parameter for the bearing. If the film parameter is less than Hydro_limit, the flow is no longer hydrodynamic.

        Args:
            printbool (bool): If True, the film parameter is printed
            Hydro_limit (float): The minimum film parameter for the flow to be considered hydrodynamic
        
        Attributes:
            Lambda (float): Film parameter
            Hydrodynamic (bool): If True, the flow is hydrodynamic
        '''
        self.Lambda = self.h_min / (self.R_a**2 + self.R_a**2)**(1/2) # Bogen (3.22) - Assuming same roughness on both sides

        if self.Lambda > Hydro_limit:
            self.Hydrodynamic = True
        else:
            self.Hydrodynamic = False

        if printbool:
            print(f"The film parameter is {self.Lambda:.3g} - The flow is {'hydrodynamic' if self.Hydrodynamic else 'not hydrodynamic'}")


    def stability(self, printbool:bool = False) -> None:
        '''
        Checks the stability of the bearing

        Args:
            printbool (bool): If True, the stability of the bearing is printed
            
        Attributes:
            stable (bool): If True, the bearing is stable
        '''
        kxx = interp1d(self.S_table, self.kx_table, kind='linear', fill_value=(self.kx_table[-1], self.kx_table[0]), bounds_error=False)(self.S) * self.W / self.C_p
        kxy = interp1d(self.S_table, self.kxy_table, kind='linear', fill_value=(self.kxy_table[-1], self.kxy_table[0]), bounds_error=False)(self.S) * self.W / self.C_p
        kyx = interp1d(self.S_table, self.kyx_table, kind='linear', fill_value=(self.kyx_table[-1], self.kyx_table[0]), bounds_error=False)(self.S) * self.W / self.C_p
        kyy = interp1d(self.S_table, self.kyy_table, kind='linear', fill_value=(self.kyy_table[-1], self.kyy_table[0]), bounds_error=False)(self.S) * self.W / self.C_p
        Bxx = interp1d(self.S_table, self.Bxx_table, kind='linear', fill_value=(self.Bxx_table[-1], self.Bxx_table[0]), bounds_error=False)(self.S) * self.W / (self.omega * self.C_p)
        Bxy = interp1d(self.S_table, self.Bxy_table, kind='linear', fill_value=(self.Bxy_table[-1], self.Bxy_table[0]), bounds_error=False)(self.S) * self.W / (self.omega * self.C_p)
        Byx = interp1d(self.S_table, self.Byx_table, kind='linear', fill_value=(self.Byx_table[-1], self.Byx_table[0]), bounds_error=False)(self.S) * self.W / (self.omega * self.C_p)
        Byy = interp1d(self.S_table, self.Byy_table, kind='linear', fill_value=(self.Byy_table[-1], self.Byy_table[0]), bounds_error=False)(self.S) * self.W / (self.omega * self.C_p)


        M = np.array([[self.mass, 0], [0, self.mass]])
        K = np.array([[kxx, kxy], [kyx, kyy]])
        B = np.array([[Bxx, Bxy], [Byx, Byy]])

        A1 = np.block([[M, np.zeros((2,2))], [np.zeros((2,2)), M]])
        A2 = np.block([[K, B], [-M, np.zeros(M.shape)]])

        s, u = eig(-A2, A1)

        xi = -np.real(s) / np.abs(s)
        wd_hz = np.imag(s) / (2*np.pi)

        if np.any(np.real(s) > 0):
            self.stable = False
        else:
            self.stable = True
        
        if printbool:
            print(f"The bearing is {'stable' if self.stable else 'unstable'} as the real part of the eigenvalues are {np.real(s)}")


    def lubrication_consumption(self, printbool:bool = False) -> None:
        '''
        Finds the lubrication consumption of the bearing

        Args:
            printbool (bool): If True, the lubrication consumption is printed
        
        Attributes:
            q (float): Lubrication consumption [m^3/s]
        '''
        # Qs = np.interp(self.S, self.S_table, self.Q_table)
        Qs = interp1d(self.S_table, self.Q_table, kind='linear', fill_value='extrapolate')(self.S)
        # Qe = np.interp(self.S, self.S_table, self.P_table)
        Qe = interp1d(self.S_table, self.P_table, kind='linear', fill_value='extrapolate')(self.S)

        # epsilon = np.interp(self.S, self.S_table, self.E_table)
        epsilon = interp1d(self.S_table, self.E_table, kind='linear', fill_value='extrapolate')(self.S)

        h1 = self.C_p - epsilon * self.C_b
        qf = 8 * h1**3 / self.eta * self.pf * self.Qf

        chi = 1 # Antagelse !!!!!!!!!!!!!!!!
        self.q = self.r * self.omega * self.C_p * self.L * (Qs + (1 - chi) * Qe) + qf

        if printbool:
            print(f"The lubrication consumption is {self.q:.3g} m^3/s")


    def friction_loss(self, printbool:bool = False) -> None:
        '''
        Finds the friction loss of the bearing

        Args:
            printbool (bool): If True, the friction loss is printed
        
        Attributes:
            H (float): Friction loss [W]
        '''
        # fj = np.interp(self.S, self.S_table, self.T_table) * self.psi
        fj = interp1d(self.S_table, self.T_table, kind='linear', fill_value='extrapolate')(self.S) * self.psi
        self.H = fj * self.r * self.W * self.omega

        if printbool:
            print(f"The friction loss is {self.H:.3g} W")


    def run_all(self) -> None:
        '''
        Runs all the functions for the bearing
        '''
        self.get_someya()
        self.someya()
        self.find_temp_visc()
        self.sommerfeld()
        self.reynolds_number()
        self.minimum_film_thickness()
        self.film_parameter()
        self.stability()
        self.lubrication_consumption()
        self.friction_loss()


    def Analytical(self, n_:int = 10, epsilon:float = 0.2) -> list:
        '''
        Calculates the analytical solution pressure distribution of the bearing

        Args:
            n (int): Resolution of the solution
            epsilon (float): The eccentricity ratio
        
        Returns:
            P (list): Pressure distribution
            y (list): y coordinates
            phi (list): phi coordinates
        ''' 
        y = np.linspace(-self.L/2, self.L/2, n_)
        phi = np.linspace(0, np.pi, n_)
        phi, y = np.meshgrid(phi, y)
        P = 3 * self.eta * self.omega * epsilon / self.C_p**2 * (self.L**2 / 4 - y**2) * np.sin(phi)/(1 + epsilon * np.cos(phi))**3
        return P, y, phi


    def Numerical(self, n_:int = 10, epsilon:float = 0.2) -> list:
        '''
        Calculates the pressure distribution with finite differnce method of the bearing

        Args:
            n (int): Resolution of the solution
            epsilon (float): The eccentricity ratio
        
        Returns:
            P (list): Pressure distribution
            y (list): y coordinates
            phi (list): phi coordinates
        '''
        y = np.linspace(-self.L/2, self.L/2, n_)
        phi = np.linspace(0, np.pi, n_)
        d_phi = abs(phi[1] - phi[0])
        d_y = abs(y[1] - y[0])
        
        M = sps.eye(n_**2, format='csr')
        rhs = np.zeros(n_**2)

        for i in range(1, n_-1):
            for j in range(1, n_-1):
                c = j + i * n_
                n = c + 1   # North
                s = c - 1   # South
                e = j + n_*(i+1)   # East
                w = j + n_*(i-1)   # West

                M[c, c] = 2 / d_phi**2 + 2 / d_y**2
                M[c, n] = -1 / d_y**2
                M[c, s] = -1 / d_y**2
                M[c, e] = -1 / d_phi**2 + 3 * epsilon * np.sin(phi[i]) / (1 + epsilon * np.cos(phi[i])) / (2 * d_phi)
                M[c, w] = -1 / d_phi**2 - 3 * epsilon * np.sin(phi[i]) / (1 + epsilon * np.cos(phi[i])) / (2 * d_phi)

                rhs[c] = 6 * self.omega * self.eta * epsilon * np.sin(phi[i]) / (self.C_p**2 * (1 + epsilon * np.cos(phi[i]))**3)
        
        p = sps.linalg.spsolve(M, rhs)
        P = p.reshape(n_, n_, order='F')

        phi, y = np.meshgrid(phi, y)
        return P, y, phi


    def Compare_sommerfield(self, epsilon_:float = 0.2) -> float:
        '''
        Calculates the sommerfield number with the analytical and numerical solution 

        Args:
            epsilon (float): The eccentricity ratio
        
        Returns:
            S_a (float): Sommerfield number from analytical solution
            S_n (float): Sommerfield number from numerical solution
            phi_analytical (np.ndarray): phi angle from analytical solution
            phi_numerical (np.ndarray): phi angle from numerical solution
        '''
        # Analytical
        w_x_a = self.eta * self.omega * self.r * self.L**3 / (4 * self.C_p**2) * np.pi * epsilon_ / ((1 - epsilon_**2)**(3/2))
        w_y_a = self.eta * self.omega * self.r * self.L**3 / (self.C_p**2) * epsilon_**2 / ((1 - epsilon_**2)**(2))
        self.W = np.sqrt(w_x_a**2 + w_y_a**2)
        self.sommerfeld()
        S_a = copy.copy(self.S)
        phi_analytical = np.arctan(w_y_a / w_x_a)*180/np.pi

        # Numerical
        P_n, y, phi = self.Numerical(epsilon=epsilon_)
        
        d_y = np.abs(y[0,0] - y[1,0])
        d_phi = np.abs(phi[0,0] - phi[0,1])

        w_x_n = 0
        w_y_n = 0

        for i in range(len(phi)):
            for j in range(len(y)):
                w_x_n = w_x_n + P_n[i,j] * self.r *  np.sin(phi[i,j]) * d_phi * d_y
                w_y_n = w_y_n + P_n[i,j] * self.r *  np.cos(phi[i,j]) * d_phi * d_y
        
        self.W = np.sqrt(w_x_n**2 + w_y_n**2)
        self.sommerfeld()
        S_n = copy.copy(self.S)
        phi_numerical = np.abs(np.arctan(w_y_n / w_x_n)*180/np.pi)

        return S_a, S_n, phi_analytical, phi_numerical



if __name__ == "__main__":
    if True: # Test
        print('test')
        cyl = Bearing(surfaceRoughness=1.2e-6, Revolutions=0.5)
        cyl.get_someya()
        cyl.someya()
        cyl.find_temp_visc(printbool=True, plot=True)
        cyl.reynolds_number(printbool=True)
        cyl.minimum_film_thickness(printbool=True)
        cyl.film_parameter(printbool=True)
        cyl.stability(printbool=True)
        # cyl.lubrication_consumption(printbool=True)
        # cyl.friction_loss(printbool=True)


    if False: # Part 1
        print("Friction Loss")

        N = np.linspace(0.01, 500, 100)
        bearings = [[], [], [], [], [], []]
        lillen = [[], [], [], [], [], []]
        friction = [[], [], [], [], [], []]
        consumption = [[], [], [], [], [], []]

        for n in N:
            bearings[0].append(Bearing(bearingtype=0, LengthDiameterFactor=0.5, preLoadFactor=0, Revolutions=n))
            bearings[1].append(Bearing(bearingtype=0, LengthDiameterFactor=1, preLoadFactor=0, Revolutions=n))
            bearings[2].append(Bearing(bearingtype=1, LengthDiameterFactor=0.5, preLoadFactor=0.5, Revolutions=n))
            bearings[3].append(Bearing(bearingtype=1, LengthDiameterFactor=1, preLoadFactor=2/3, Revolutions=n))
            bearings[4].append(Bearing(bearingtype=2, LengthDiameterFactor=0.5, preLoadFactor=0, PadLoes=0, Revolutions=n))
            bearings[5].append(Bearing(bearingtype=2, LengthDiameterFactor=0.5, preLoadFactor=0, PadLoes=1, Revolutions=n))

        for i in range(6):
            for j in range(len(N)):
                bearings[i][j].run_all()

        

        for i in range(6):
            hydroflag = True
            stableflag = True
            laminarflag = True
            for j in range(len(N)):
                if bearings[i][j].Hydrodynamic and bearings[i][j].stable and bearings[i][j].laminar:
                    friction[i].append(bearings[i][j].H)
                    consumption[i].append(bearings[i][j].q)
                    lillen[i].append(bearings[i][j].N)
                elif not bearings[i][j].stable:
                    if stableflag:
                        print(f"Bearing type {i+1} is unstable from {N[j]:.3g} Hz")
                        stableflag = False
                elif not bearings[i][j].laminar:
                    if laminarflag:
                        print(f"Bearing type {i+1} is turbulent from {N[j]:.3g} Hz")
                        laminarflag = False
                elif bearings[i][j].Hydrodynamic:
                    if hydroflag:
                        print(f"Bearing type {i+1} is hydrodynamic from {N[j]:.3g} Hz")
                        hydroflag = False
                    
        plt.figure()
        for i in range(6):
            plt.plot(lillen[i], friction[i], label=f"Bearing type {i+1}")
        plt.xlabel("Revolutions [Hz]")
        plt.ylabel("Friction loss [W]")
        plt.legend()

        plt.figure()
        for i in range(6):
            plt.plot(lillen[i], consumption[i], label=f"Bearing type {i+1}")
        plt.xlabel("Revolutions [Hz]")
        plt.ylabel("Lubrication consumption [m^3/s]")
        plt.legend()
        plt.show()

        # Save data
        for i in range(6):
            np.savetxt(f"Assignment2/1-data/friction_{i+1}.txt", np.array([lillen[i], friction[i]]).T)
            np.savetxt(f"Assignment2/1-data/consumption_{i+1}.txt", np.array([lillen[i], consumption[i]]).T)


    if False: # Part 2-1
        Part2 = Bearing()
        Part2.temp_profile()
        Part2.get_someya()
        Part2.someya()
        Part2.find_temp_visc()

        P_analytical, y, phi = Part2.Analytical()

        plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(y, phi, P_analytical)

        P_nummerical, y, phi = Part2.Numerical()

        plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(y, phi, P_nummerical)

        # Save data
        np.savetxt("Assignment2/2-data/pressure_analytical.txt", np.array([y.flatten('F'), phi.flatten('F'), P_analytical.flatten('F')]).T)
        np.savetxt("Assignment2/2-data/pressure_numerical.txt", np.array([y.flatten('F'), phi.flatten('F'), P_nummerical.flatten('F')]).T)

        # Compare the data
        mse = np.mean((P_analytical - P_nummerical)**2)
        rmse = np.sqrt(mse)
        print(f"The RMSE between the analytical and numerical solution is {rmse:.3g}")

        # Correlation
        corr = np.corrcoef(P_analytical.flatten(), P_nummerical.flatten())[0,1]
        print(f"The correlation between the analytical and numerical solution is {corr:.3g}")

        # procent Max difference
        max_diff = np.max(np.abs(np.max(P_analytical) - np.max(P_nummerical)))
        max_val = np.max(np.abs(P_analytical))
        procent_diff = max_diff / max_val * 100
        print(f"The maximum difference between the analytical and numerical solution is {procent_diff:.3g}%")

        plt.show()


    if False: # Part 2-2
        Part22 = Bearing()
        Part22.temp_profile()
        Part22.get_someya()
        Part22.someya()
        Part22.find_temp_visc()

        Phi = np.flip(Part22.Phi_table)
        Epsilon = Part22.E_table
        Sommerfeld = Part22.S_table
        S_analytical = []
        S_numerical = []
        Phi_analytical = []
        Phi_numerical = []

        for i in range(len(Epsilon)):
            S_a, S_n, phi_a, phi_n = Part22.Compare_sommerfield(epsilon_=Epsilon[i])
            S_analytical.append(S_a)
            S_numerical.append(S_n)
            Phi_analytical.append(phi_a)
            Phi_numerical.append(phi_n)
        
        plt.figure()
        plt.plot(Epsilon, Phi_analytical, label="Analytical")
        plt.plot(Epsilon, Phi_numerical, label="Numerical")
        plt.plot(Epsilon, Phi, label="Someya")
        plt.xlabel("Epsilon")
        plt.ylabel("Phi [deg]")
        plt.legend()

        # Save data
        np.savetxt("Assignment2/2-data/Phi_analytical.txt", np.array([Epsilon, Phi_analytical]).T)
        np.savetxt("Assignment2/2-data/Phi_numerical.txt", np.array([Epsilon, Phi_numerical]).T)
        np.savetxt("Assignment2/2-data/Phi.txt", np.array([Epsilon, Phi]).T)

        plt.figure()
        plt.plot(Epsilon, S_analytical, label="Analytical")
        plt.plot(Epsilon, S_numerical, label="Numerical")
        plt.plot(Epsilon, Sommerfeld, label="Someya")
        plt.xlabel("Epsilon")
        plt.ylabel("Sommerfeld number")
        plt.legend()

        # Save data
        np.savetxt("Assignment2/2-data/S_analytical.txt", np.array([Epsilon, S_analytical]).T)
        np.savetxt("Assignment2/2-data/S_numerical.txt", np.array([Epsilon, S_numerical]).T)
        np.savetxt("Assignment2/2-data/S.txt", np.array([Epsilon, Sommerfeld]).T)

        plt.show()