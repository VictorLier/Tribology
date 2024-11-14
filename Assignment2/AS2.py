import numpy as np
import matplotlib.pyplot as plt

class Bearing:
    def __init__(self, ShaftDiameter:float=101.6e-3, RadialClearance:float=102e-6, preLoadFactor:float=0, LengthDiameterFactor:float=0.5, Load:float=460/2, surfaceRoughness:float=0.5e-6, oil:int=0, Revolutions:float=20, StartTemp:float = 40+273.15, feedPressure:float = 2e5) -> None:
        '''
        Creates a bearing object with the given parameters

        Args:
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

        self.W = Load * 9.82 # Convert to N
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


    def find_temp_visc(self, plot:bool = False, printbool:bool = False, convergenceTol:float = 1e-6, maxConvergence:int=100) -> None:
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
            T[i+1] = T[i] + 0.1 * (t_new - T[i])
            self.T = T[i+1]

            if abs(t_new - T[i+1]) < convergenceTol:
                self.nu = self.temp_relation(self.T)

                if printbool:
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


      








class Cylindrical(Bearing):
    def __init__(self, ShaftDiameter: float = 0.1016, RadialClearance: float = 0.000102, preLoadFactor: float = 0, LengthDiameterFactor: float = 0.5, Load: float = 460 / 2, surfaceRoughness: float = 5e-7, oil: int = 0) -> None:
        '''
        Creates a Cylindrical bearing object with the given parameters

        Args:
            shaftDiameter (float): The diameter of the shaft in meters. Variable name: D
            radialClearance (float): The radial clearance in meters. Variable name: C_p
            preLoadFactor (float): The preload factor. Either 0, 0.5, or 2/3. Variable name: mp
            lengthDiameterFactor (float): The length diameter factor. Either 0.5 or 1. Variable name: L
            load (float): The load on one bearing in Kg
            surfaceRoughness (float): The surface roughness in meters. Variable name: Ra - See table 3.2
            oil (int): The oil type. 0 for ISO VG 32, 1 for ISO VG 46

            
        Raises:
            ValueError: If the preLoadFactor is not 0, 0.5, or 2/3
            ValueError: If the lengthDiameterFactor is not 0.5 or 1
            ValueError: If the oil type is not 0 or 1
        '''
        super().__init__(ShaftDiameter, RadialClearance, preLoadFactor, LengthDiameterFactor, Load, surfaceRoughness, oil)


    def get_someya(self) -> None:
        '''
        Get the someya parameters for the bearing from tables 
        '''
        if self.LD == 0.5:
            from Someya import Table1 as table
        elif self.LD == 1:
            from Someya import Table2 as table
        else:
            raise ValueError("No someya table found - Invalid length diameter factor")

        self.table = table



class TwoLobes(Bearing):
    def __init__(self, ShaftDiameter: float = 0.1016, RadialClearance: float = 0.000102, preLoadFactor: float = 0, LengthDiameterFactor: float = 0.5, Load: float = 460 / 2, surfaceRoughness: float = 5e-7, oil: int = 0) -> None:
        '''
        Creates a Two-Lobes bearing object with the given parameters

        Args:
            shaftDiameter (float): The diameter of the shaft in meters. Variable name: D
            radialClearance (float): The radial clearance in meters. Variable name: C_p
            preLoadFactor (float): The preload factor. Either 0, 0.5, or 2/3. Variable name: mp
            lengthDiameterFactor (float): The length diameter factor. Either 0.5 or 1. Variable name: L
            load (float): The load on one bearing in Kg
            surfaceRoughness (float): The surface roughness in meters. Variable name: Ra - See table 3.2
            oil (int): The oil type. 0 for ISO VG 32, 1 for ISO VG 46
            
        Raises:
            ValueError: If the preLoadFactor is not 0, 0.5, or 2/3
            ValueError: If the lengthDiameterFactor is not 0.5 or 1
            ValueError: If the oil type is not 0 or 1
        '''
        super().__init__(ShaftDiameter, RadialClearance, preLoadFactor, LengthDiameterFactor, Load, surfaceRoughness, oil)


    def get_someya(self) -> None:
        '''
        Get the someya parameters for the bearing from tables 
        '''
        if self.LD == 0.5 and self.mp == 0.5:
            from Someya import Table3 as table
        elif self.LD == 1 and self.mp == 2/3:
            from Someya import Table4 as table
        else:
            raise ValueError("No someya table found - Invalid length diameter factor or preload factor")

        self.table = table



class Pad(Bearing):
    def __init__(self, ShaftDiameter: float = 0.1016, RadialClearance: float = 0.000102, preLoadFactor: float = 0, LengthDiameterFactor: float = 0.5, Load: float = 460 / 2, surfaceRoughness: float = 5e-7, oil: int = 0, PadLoes: int = 0) -> None:
        '''
        Creates a pad bearing object with the given parameters

        Args:
            shaftDiameter (float): The diameter of the shaft in meters. Variable name: D
            radialClearance (float): The radial clearance in meters. Variable name: C_p
            preLoadFactor (float): The preload factor. Either 0, 0.5, or 2/3. Variable name: mp
            lengthDiameterFactor (float): The length diameter factor. Either 0.5 or 1. Variable name: L
            load (float): The load on one bearing in Kg
            surfaceRoughness (float): The surface roughness in meters. Variable name: Ra - See table 3.2
            oil (int): The oil type. 0 for ISO VG 32, 1 for ISO VG 46
            PadLoes (int): 0 for no pad loss, 1 for pad loss
            
        Raises:
            ValueError: If the preLoadFactor is not 0, 0.5, or 2/3
            ValueError: If the lengthDiameterFactor is not 0.5 or 1
            ValueError: If the oil type is not 0 or 1
            ValueError: If the pad load is not 0 or 1
        '''
        super().__init__(ShaftDiameter, RadialClearance, preLoadFactor, LengthDiameterFactor, Load, surfaceRoughness, oil)
    
        if PadLoes not in [0, 1]:
            raise ValueError("Invalid pad load")
        else:
            self.PadLoes = PadLoes
    

    def get_someya(self) -> None:
        '''
        Get the someya parameters for the bearing from tables 
        '''
        if self.LD == 0.5 and self.mp == 0 and self.PadLoes == 0:
            from Someya import Table5 as table
        elif self.LD == 0.5 and self.mp == 0 and self.PadLoes == 1:
            from Someya import Table6 as table
        else:
            raise ValueError("No someya table found - Invalid length diameter factor or preload factor or pad load")
        
        self.table = table




if __name__ == "__main__":
    if True: # Test
        print('test')
        cyl = Cylindrical()
        cyl.temp_profile()
        cyl.get_someya()
        cyl.someya()
        cyl.find_temp_visc(printbool=True, plot=True)