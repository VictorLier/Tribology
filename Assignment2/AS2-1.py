import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig

class Bearing:
    def __init__(self, bearingtype:int = 0, PadLoes:int=0, ShaftDiameter:float=101.6e-3, RadialClearance:float=102e-6, preLoadFactor:float=0, LengthDiameterFactor:float=0.5, Load:float=460/2, surfaceRoughness:float=0.5e-6, oil:int=0, Revolutions:float=20, StartTemp:float = 40+273.15, feedPressure:float = 2e5) -> None:
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
            StartTemp (float): The starting temperature of the oil in Kelvin
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
        self.T = StartTemp
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
            if self.LD == 0.5:
                from Someya import Table1 as table
                self.table = table
            elif self.LD == 1:
                from Someya import Table2 as table
                self.table = table
            else:
                raise ValueError("Invalid length diameter factor")
        elif self.bearingtype == 1:
            if self.LD == 0.5 and self.mp == 0.5:
                from Someya import Table3 as table
                self.table = table
            elif self.LD == 1 and self.mp == 2/3:
                from Someya import Table4 as table
                self.table = table
            else:
                raise ValueError("Invalid length diameter factor or preload factor")
        elif self.bearingtype == 2:
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
        self.table = np.flip(self.table, axis=0)

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


    def temp_profile(self, printbool: bool = False) -> None:
        '''
        Finds the constant for the 'linear' fit of the temperature-viscosity relationship
        
        Args:
            printbool (bool): If True, the constants are printed

        Attributes:
            A1 (float): Constant for the 'linear' fit
            A2 (float): Constant for the 'linear' fit
        '''
        temp = np.array([40, 100]) + 273.15  # Temperature in Kelvin
        visc = np.array([self.nu_40, self.nu_100])
        visc = np.log(np.log(visc*1e6+0.8))
        temp = np.log(temp)

        a_1 = np.polyfit(temp, visc, 1)
        self.A1 = a_1[0]
        self.A2 = a_1[1]

        if printbool:
            print(f"Following the ASTM walther relation, the constants are: m = {self.A1:.3g} and C = {self.A2:.3g}")


    def temp_relation(self, T:float, printbool: bool = False) -> float:
        '''
        Calculates the kinematic viscosity of the oil at a given temperature

        Args:
            T (float): The temperature in Kelvin
            printbool (bool): If True, the constants are printed
        
        Returns:
            float: The kinematic viscosity of the oil at the given temperature [m^2/s]
        '''
        nu = np.exp(np.exp(self.A1*np.log(T) + self.A2)) - 0.8 #[cSt]
        nu = nu * 1e-6 # Convert to m^2/s

        if printbool:
            print(f"The kinematic viscosity of the oil at {T} K is {nu:.3g} m^2/s")
        
        return nu


    def sommerfeld(self, printbool:bool = False) -> None:
        '''
        Calculates the sommerfeld number for the bearing

        Args:
            printbool (bool): If True, the sommerfeld number is printed

        Attributes:
            S (float): Sommerfeld number
        '''
        self.S = (self.nu*self.rho) * self.N * self.L * self.D / self.W * ((self.D/2)/self.C_p)**2

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
        A = 9 * self.D * self.D**(1/2)
        w = 1 # [m/s] sourrounding air - assumed 1
        alpha = 9.807 * (0.7 + 1.2 * np.sqrt(w)) # [W/m^2*K] - Heat transfer coefficient

        lamb = 1/3 # S8-27 (31)

        Qf = 2 * 0.15 # Aner ikke hvor det kommer fra......

        t_0 = 20#+273.15 # [C] - sourrounding air temp
        t_1 = 30#+273.15 # [C] - oil inlet temp


        iteration = np.arange(maxConvergence+1)
        T = np.zeros(maxConvergence+1)

        T[0] = self.T # Initial temperature

        for i in range(maxConvergence):
            self.nu = self.temp_relation(T[i]) 
            eta = self.nu * self.rho
            self.sommerfeld()
            Qs = np.interp(self.S, self.S_table, self.Q_table)


            f_J = self.psi * np.interp(self.S, self.S_table, self.T_table)
            epsi = np.interp(self.S, self.S_table, self.E_table)

            h1 = self.C_p - epsi * self.C_b
            qf = 8 * h1**3 / eta * self.pf * Qf
            q = self.r * self.omega * self.C_p * self.L * Qs + qf

            t_new = ((1-lamb) * f_J * self.r * self.W * self.omega + alpha * A * t_0) + self.CP * self.rho * q * t_1 / (self.CP * self.rho * q + alpha * A * (1-lamb))
            T[i+1] = T[i] + 0.05 * (t_new - T[i])
            self.T = T[i+1]

            if abs(t_new - T[i+1]) < convergenceTol:
                self.nu = self.temp_relation(self.T)
                self.eta = self.nu * self.rho

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
        epsilon = np.interp(self.S, self.S_table, self.E_table)

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
        self.Lambda = self.h_min / (self.R_a**2 + self.R_a**2)**0.5 # Bogen (3.22) - Assuming same roughness on both sides

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
        epsilon = np.interp(self.S, self.S_table, self.E_table)
        phi = np.interp(self.S, self.S_table, self.Phi_table) * np.pi / 180
        kxx = np.interp(self.S, self.S_table, self.kx_table) * self.W / self.C_p
        kxy = np.interp(self.S, self.S_table, self.kxy_table) * self.W / self.C_p
        kyx = np.interp(self.S, self.S_table, self.kyx_table) * self.W / self.C_p
        kyy = np.interp(self.S, self.S_table, self.kyy_table) * self.W / self.C_p
        Bxx = np.interp(self.S, self.S_table, self.Bxx_table) * self.W / (self.omega * self.C_p)
        Bxy = np.interp(self.S, self.S_table, self.Bxy_table) * self.W / (self.omega * self.C_p)
        Byx = np.interp(self.S, self.S_table, self.Byx_table) * self.W / (self.omega * self.C_p)
        Byy = np.interp(self.S, self.S_table, self.Byy_table) * self.W / (self.omega * self.C_p)

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
        Qf = 2 * 0.15 # Aner ikke hvor det kommer fra......

        Qs = np.interp(self.S, self.S_table, self.Q_table)
        Qe = np.interp(self.S, self.S_table, self.P_table)

        epsilon = np.interp(self.S, self.S_table, self.E_table)

        h1 = self.C_p - epsilon * self.C_b
        qf = 8 * h1**3 / self.eta * self.pf * Qf

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
        fj = np.interp(self.S, self.S_table, self.T_table) * self.psi
        self.H = fj * self.r * self.W * self.omega

        if printbool:
            print(f"The friction loss is {self.H:.3g} W")


    def run_all(self) -> None:
        '''
        Runs all the functions for the bearing
        '''
        self.get_someya()
        self.someya()
        self.temp_profile()
        self.find_temp_visc()
        self.sommerfeld()
        self.reynolds_number()
        self.minimum_film_thickness()
        self.film_parameter()
        self.stability()
        self.lubrication_consumption()
        self.friction_loss()




if __name__ == "__main__":
    if False: # Test
        print('test')
        cyl = Bearing(surfaceRoughness=5e-6, Revolutions=100)
        cyl.temp_profile()
        cyl.get_someya()
        cyl.someya()
        cyl.find_temp_visc(printbool=True, plot=True)
        cyl.reynolds_number(printbool=True)
        cyl.minimum_film_thickness(printbool=True)
        cyl.film_parameter(printbool=True)
        cyl.stability(printbool=True)
        cyl.lubrication_consumption(printbool=True)
        cyl.friction_loss(printbool=True)


    if True:
        print("Friction Loss")

        N = np.linspace(0.1, 120, 500)

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
                    if hydroflag:
                        print(f"Bearing type {i+1} is hydrodynamic from {N[j]:.3g} Hz")
                        hydroflag = False
                    friction[i].append(bearings[i][j].H)
                    consumption[i].append(bearings[i][j].q)
                    lillen[i].append(N[j])
                elif not bearings[i][j].stable:
                    if stableflag:
                        print(f"Bearing type {i+1} is unstable from {N[j]:.3g} Hz")
                        stableflag = False
                elif not bearings[i][j].laminar:
                    if laminarflag:
                        print(f"Bearing type {i+1} is turbulent from {N[j]:.3g} Hz")
                        laminarflag = False
                    
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