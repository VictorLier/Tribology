import numpy as np
import matplotlib.pyplot as plt

class Bearing:
    '''
    Cylindrical bearing class for Assignment 3
    '''
    def __init__(self, InnerRaceDia:float=69.5e-3, OuterRaceDia:float=85.5e-3, RollerDia:float=8e-3, RollerLength:float=8e-3, NoOfRollers:int=20, MaxLoad:float=8000, InnerRaceSpeed:float=837, OuterRaceSpeed:float=0, EModulus:float=210e9, PoissonRatio:float=0.3, AbsoluteViscosity:float=0.01, PressViscCoef:float=2e-8)->None:
        '''
        Creates a new instance of the Bearing class with the given parameters.

        Args:
            InnerRaceDia (float) [m] - Inner race diameter
            OuterRaceDia (float) [m] - Outer race diameter
            RollerDia (float) [m] - Diameter of the rollers
            RollerLength (float) [m] - Length of the rollers
            NoOfRollers (int) - Number of rollers
            MaxLoad (float) [N] - Maximum radial load
            InnerRaceSpeed (float) [rad/s] - Angular velocity of inner race
            OuterRaceSpeed (float) [rad/s] - Angular velocity of outer race
            EModulus (float) [Pa] - Young's Modulus
            PoissonRatio (float) - Poisson's ratio
            AbsoluteViscosity (float) [Pa s] - Base absolute viscosity of the lubricant
            PressViscCoef (float) [1 / (Pa s)] - Pressure-viscosity coefficient
        
        Attributes:
            d_i (float) [m] - Inner race diameter
            d_o (float) [m] - Outer race diameter
            d (float) [m] - Diameter of the rollers
            l (float) [m] - Length of the rollers
            i (int) - Number of rollers
            w_z (float) [N] - Maximum radial load
            omega_i (float) [rad/s] - Angular velocity of inner race
            omega_o (float) [rad/s] - Angular velocity of outer race
            E (float) [Pa] - Young's Modulus
            nu (float) - Poisson's ratio
            eta_0 (float) [Pa s] - Base absolute viscosity of the lubricant
            xi (float) [1 / (Pa s)] - Pressure-viscosity coefficient
        '''
        self.d_i = InnerRaceDia
        self.d_o = OuterRaceDia
        self.d = RollerDia
        self.l = RollerLength
        self.n = NoOfRollers
        self.w_z = MaxLoad
        self.omega_i = InnerRaceSpeed
        self.omega_o = OuterRaceSpeed
        self.E = EModulus
        self.nu = PoissonRatio
        self.eta_0 = AbsoluteViscosity
        self.xi = PressViscCoef

        self.r_ax = self.d / 2 # (fig. 17.2)
        self.r_ay = np.inf # (fig. 17.2)
        self.r_bx = self.d_i / 2 # (fig. 17.2)
        self.r_ay = 1e15 # np.inf # (fig. 17.2)
        self.r_by = 1e15 # np.inf # (fig. 17.2)


    def max_load(self, printbool:bool=False)->None:
        '''
        Calculates the load of the roller with the highest load
        
        Args:
            printbool (bool) - Prints the maximum load if True

        Attributes:
            w_zm (float) [N] - Maximum radial load
        '''
        self.w_zm = self.w_z*4 / self.n   # (21.47)

        if printbool:
            print(f'Maximum load on a roller: w_zm = {self.w_zm:.4g} N')


    def min_film_thickness(self, printbool:bool=False, Lambda:float=3)->None:
        '''
        Calculates the minimum film thickness assuming worst case scenario in the Elasto-Hydrodynamic Lubrication (EHL) regime
        
        Args:
            printbool (bool) - Prints the minimum film thickness if True
            Lambda (float) - Non-dimensional film parameter (default: 3) - Lowest in the EHL regime - p. 57

        Attributes:
            h_min (float) [m] - Minimum film thickness
            Lambda (float) - Non-dimensional parameter
            self.R_qt (float) [m] - RMS roughness of the track
            self.R_qr (float) [m] - RMS roughness of the roller
        '''
        self.Lambda = Lambda
        R_at = 0.3e-6 # Worst case from figure 3.8
        R_ar = 0.12e-6 # Worst case from figure 3.8
        self.R_qt = R_at * 1.11 # (3.5)
        self.R_qr = R_ar * 1.11 # (3.5)

        self.h_min = self.Lambda * (self.R_qt**2 + self.R_qr**2)**0.5 # (3.22)

        if printbool:
            print(f'Minimum film thickness: h_min = {self.h_min:.4g} m')


    def effective_radius(self, printbool:bool=False)->None:
        '''
        Calculates the effective radius in the x and y direction
        
        Args:
            printbool (bool) - Prints the effective radius if True
        
        Attributes:
            R_x (float) [m] - Effective radius in the x direction
            R_y (float) [m] - Effective radius in the y direction
        '''
        self.R_x = (self.r_ax * self.r_bx) / (self.r_ax + self.r_bx) # (17.4)
        self.R_y = (self.r_ay * self.r_by) / (self.r_ay + self.r_by) # (17.5)
        if np.isnan(self.R_y):
            self.R_y = np.inf

        if printbool:
            print(f'Effective radius in x: R_x = {self.R_x:.4g} m')
            print(f'Effective radius in y: R_y = {self.R_y:.4g} m')


    def curvature(self, printbool:bool=False)->None:
        '''
        Calculates the curvature sum, curvature difference and radius ratio
        
        Args:
            printbool (bool) - Prints the curvature difference if True

        Attributes:
            R (float) - Curvature sum
            R_d (float) - Curvature difference
            alpha_r (float) - Radius ratio
        
        Raises:
            ValueError - The radius ratio smaller than one - eq.(17.1) is not satisfied
        '''
        self.R = (self.R_x * self.R_y) / (self.R_x + self.R_y) # (17.2)
        if np.isnan(self.R):
            self.R = np.inf
        
        self.R_d = self.R * (1 / self.R_x - 1 / self.R_y) # (17.3)
        self.alpha_r = self.R_y / self.R_x # (16.57)

        if self.alpha_r < 1:
            raise ValueError('Radius ratio smaller than one - eq.(17.1) is not satisfied')

        if printbool:
            print(f'Curvature sum: R = {self.R:.4g}')
            print(f'Curvature difference: R_d = {self.R_d:.4g}')
            print(f'Radius ratio: alpha_r = {self.alpha_r:.4g}')


    def ellipticity_parameter(self, k_init:float=0.5, convergence_tol:float=1e-7, max_iter:int=10000, phi_no:int=1000, printbool:bool=False, plotbool:bool=False)->None:
        '''
        Calculates the ellipticity parameter of the roller bearing with a iterative method

        Args:
            k_init (float) - Initial guess for the ellipticity parameter (default: 0.5)
            convergence_tol (float) - Convergence tolerance for the iterative method (default: 1e-7 - p.437)
            max_iter (int) - Maximum number of iterations (default: 1000)
            phi_no (int) - Number of points to evaluate the elliptic integrals (default: 1000)
            printbool (bool) - Prints the ellipticity parameter if True (default: False)
            plotbool (bool) - Plots the convergence of the ellipticity parameter if True (default: False)
        
        Attributes:
            k (float) - Ellipticity parameter
            epsilon (float) - complete elliptic integral of second kind
            F (float) - complete elliptic integral of first kind

        Raises:
            ValueError - If the iterative method does not converge
        '''
        k = np.zeros(max_iter+1)
        k[0] = k_init

        phi = np.linspace(0, np.pi/2, phi_no)

        for i, k_current in enumerate(k[:-1]):
            F = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**(-0.5)) # (17.10)
            epsilon = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**0.5, phi) # (17.11)
            k[i+1] = ( (2*F - epsilon*(1+self.R_d)) /(epsilon * (1-self.R_d)) )**(1/2) # (17.10)

            if abs(k[i] - k[i+1]) < convergence_tol:
                break

        np.trim_zeros(k)
        self.k = k[-1]
        self.F = F
        self.epsilon = epsilon
        iter = np.arange(len(k))

        if plotbool:
            plt.plot(iter, k)
            plt.title('Convergence of the ellipticity parameter')
            plt.xlabel('Iteration')
            plt.ylabel('Ellipticity parameter')
            plt.title('Convergence of the ellipticity parameter')
        
        if len(iter) == max_iter:
            raise ValueError('Iterative method did not converge')
        
        if printbool:
            print(f'Ellipticity parameter: k = {self.k:.4g}')
            print(f'Complete elliptic integral of second kind: Epsilon = {self.epsilon:.4g}')
            print(f'Complete elliptic integral of first kind: F = {self.F:.4g}')


    def effective_elastic_modulus(self, printbool:bool=False)->None:
        '''
        Calculates the effective elastic modulus assuming same material for the roller and the track

        Args:
            printbool (bool) - Prints the effective elastic modulus if True
        
        Attributes:
            E_prime (float) [Pa] - Effective elastic modulus
        '''
        self.E_prime = 1 / ( (1-self.nu**2) / self.E ) # (17.17)

        if printbool:
            print(f'Effective elastic modulus: E_prime = {self.E_prime:.4g} Pa')


    def diameter_contact(self, printbool:bool=False)->None:
        '''
        Calculates the diameter of the contact ellipse in both directions

        Args:
            printbool (bool) - Prints the diameter of the contact ellipse if True
        
        Attributes:
            D_x (float) [m] - Diameter of the contact ellipse in the x direction
            D_y (float) [m] - Diameter of the contact ellipse in the y direction
        '''
        self.D_x = 2 * (6 * self.epsilon * self.w_zm * self.R/(np.pi * self.k * self.E_prime))**(1/3) # (17.15)
        self.D_y = 2 * (6 * self.k**2 * self.epsilon * self.w_zm * self.R/(np.pi * self.E_prime))**(1/3) # (17.16)

        if printbool:
            print(f'Diameter of the contact ellipse in x: D_x = {self.D_x:.4g} m')
            print(f'Diameter of the contact ellipse in y: D_y = {self.D_y:.4g} m')


    def maximum_pressure(self, printbool:bool=False)->None:
        '''
        Calculates the maximum contact pressure in the bearing.

        Args:
            printbool (bool) - Prints the maximum contact pressure if True
        
        Attributes:
            p_m (float) [Pa] - Maximum contact pressure
        '''
        self.p_m = 6 * self.w_zm / (np.pi * self.D_x * self.D_y) # (17.7)

        if printbool:
            print(f'Maximum contact pressure: p_m = {self.p_m:.4g} Pa')


    def maximum_deform(self, printbool:bool=False)->None:
        '''
        Calculates the maximum elastic deformation
        
        Args:
            printbool (bool) - Prints the maximum elastic deformation if True
        
        Attributes:
            delta_m (float) [m] - Maximum elastic deformation
        '''
        self.delta_m = self.F * (9/(2 * self.epsilon * self.R) * (self.w_zm / (np.pi * self.k * self.E_prime))**2)**(1/3) # (17.16)

        if printbool:
            print(f'Maximum elastic deformation: delta_m = {self.delta_m:.4g} m')


    def rectangular_dimensionless_load(self, printbool:bool=False)->None:
        '''
        Calculates the dimensionless load and load per unit width of rectangular contact
        
        Args:
            printbool (bool) - Prints the contact semi-width if True
        
        Attributes:
            W_prime (float) - dimensionless load
            w_x_prime (float) [N/m] - load per unit width
        '''
        self.w_x_prime = self.w_zm / self.l # (p. 448)
        self.W_prime = self.w_x_prime / (self.E_prime * self.R_x)   # (17.38)

        if printbool:
            print(f'Dimensionless load: W_prime = {self.W_prime:.4g}')
            print(f'Load per unit width: w_x_prime = {self.w_x_prime:.4g} N/m')
    
    def rectangular_max_deformation(self, printbool:bool=False)->None:
        '''
        Calculates the maximum elastic deformation of the rectangular contact
        
        Args:
            printbool (bool) - Prints the maximum elastic deformation if True
        
        Attributes:
            delta_m (float) [m] - Maximum elastic deformation
        '''
        self.delta_m = 2 * self.W_prime * self.R_x / np.pi * (np.log(2 * np.pi/self.W_prime) - 1) # (17.39)

        if printbool:
            print(f'Maximum elastic deformation: delta_m = {self.delta_m:.4g} m')


    def rectangular_max_pressure(self, printbool:bool=False)->None:
        '''
        Calculates the maximum contact pressure of the rectangular contact
        
        Args:
            printbool (bool) - Prints the maximum contact pressure if True
        
        Attributes:
            p_m (float) [Pa] - Maximum contact pressure
        '''
        self.p_m = self.E_prime * (self.W_prime/(2 * np.pi))**(1/2) # (17.40)

        if printbool:
            print(f'Maximum contact pressure: p_m = {self.p_m:.4g} Pa')


if __name__ == '__main__':
    if False: # Test
        bearing = Bearing()


    if True: # Question 1
        print('Question 1')
        Q1 = Bearing()
        Q1.max_load(printbool=True)
        Q1.min_film_thickness(printbool=True)


    if True: # Question 2
        print('Question 2')
        Q2 = Bearing()
        Q2.max_load()
        Q2.effective_radius(printbool=True)
        Q2.effective_elastic_modulus(printbool=True)
        Q2.rectangular_dimensionless_load(printbool=True)
        Q2.rectangular_max_deformation(printbool=True)
        Q2.rectangular_max_pressure(printbool=True)    





