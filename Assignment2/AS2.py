import numpy as np

class Bearing:
    def __init__(self, ShaftDiameter:float=101.6e-3, RadialClearance:float=102e-6, preLoadFactor:float=0, LengthDiameterFactor:float=0.5, Load:float=460/2, surfaceRoughness:float=0.5e-6, oil:int=0, Revolutions:float=20) -> None:
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

        Raises:
            ValueError: If the preLoadFactor is not 0, 0.5, or 2/3
            ValueError: If the lengthDiameterFactor is not 0.5 or 1
            ValueError: If the oil type is not 0 or 1
        '''
        self.D = ShaftDiameter
        self.C_p = RadialClearance

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
            self.eta_40 = self.nu_40 * self.rho  # [Pa*s] - Dynamic/kinematic viscosity
            self.CP = 1964 # [J/kg*K] - Specific heat capacity
        
        elif oil == 1: # ISO VG 46
            self.nu_40 = 46e-6 # [m^2/s] - Kinematic viscosity - ISO 3448
            self.nu_100 = 6.8e-6 # [m^2/s] - Kinematic viscosity - https://www.engineeringtoolbox.com/iso-grade-oil-d_1207.html
            self.rho = 861 # [kg/m^3] - Density -  https://www.engineeringtoolbox.com/iso-grade-oil-d_1207.html
            self.eta_40 = self.nu_40 * self.rho  # [Pa*s] - Dynamic/kineamtic viscosity
            self.CP = 1964 # [J/kg*K] - Specific heat capacity

        else:
            raise ValueError("Invalid oil type")
        
        self.N = Revolutions
        self.omega = 2*np.pi*Revolutions

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

        self.S = self.table[:,0]     # Sommerfeld number
        self.E = self.table[:,1]     # eccentricity ratio
        self.Phi = self.table[:,2]   # attitude angle (deg)
        self.Q = self.table[:,3]     # flow rate ratio
        self.P = self.table[:,4]     
        self.T = self.table[:,5]
    

    def find_temp_visc(self) -> None:
        '''
        Finds the temperature and viscusity of the oil in the bearing 
        '''
        nu_40 = self.nu_40*1e6
        nu_100 = self.nu_100*1e6
        m_lub = (np.log(np.log(nu_100 + 0.8)) - np.log(np.log(nu_40 + 0.8))) / (np.log((40+273.15)/(100+273.15)))
        k_lub = np.log(np.log(nu_40 + 0.8)) + m_lub*np.log(40+273.15)
 
 


        def eta_i(temp):
            return self.rho*1e-6*(np.exp(np.exp(-m_lub*np.log(temp+273.15) + k_lub)) - 0.8)
        
        t_new = 30
        T_current = 40
        eta = self.eta_40
        N = 30
        L_mark = self.L/2
        p_f = 2e5 # [Pa]
        r_b = self.D/2
        psi = self.C_p/r_b
        lamb = 1/3

        w_air = 1 # [m/s] air velocity of surroundings
        # Vogelpohl equation
        alpha = 9.807 * (0.7 + 1.2 * w_air**(1/2))
        # Reitemeyer equation
        A = 9 * self.D * self.D**(1/2)

        t_1 = 30
        t_0 = 20

        damping = 0.5

        while abs(T_current - t_new) > 0.01:
            S_current = eta*N*self.L*self.D/self.W*(r_b/self.C_p)**2
            Q_current = np.interp(S_current, self.S, self.Q)
            T_current = np.interp(S_current, self.S, self.T)
            epsi_current = np.interp(S_current, self.S, self.E)

            eta = eta_i(T_current)
            T_current = t_new
            q_f = np.pi*self.C_p**3/(3*eta*L_mark/self.D) * ( 1 + 3/2 * epsi_current**2 )*p_f
            q = r_b*self.omega*self.C_p*self.L*Q_current + q_f   # assuming chi = 1
            f_J = psi * T_current

            t_new = ((1-lamb)*(f_J*r_b*self.W*self.omega + alpha*A*t_0) + self.CP*self.rho*q*t_1) / (self.CP*self.rho*q + alpha*A*(1-lamb))
            t_new = T_current + damping * (t_new - T_current)
            print(f"Current: {T_current} - New: {t_new}")








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
        cyl = Cylindrical()
        cyl.get_someya()
        cyl.someya()
        cyl.find_temp_visc()