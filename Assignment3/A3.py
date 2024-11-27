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
            r_ax (float) [m] - Curvature radius of the roller in the x direction
            r_ay (float) [m] - Curvature radius in the roller y direction
            r_bxi (float) [m] - Curvature radius of the inner track in the x direction
            r_bxo (float) [m] - Curvature radius of the outer track in the x direction
            r_byi (float) [m] - Curvature radius of the inner track in the y direction
            r_byo (float) [m] - Curvature radius of the outer track in the y direction
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

        self.r_ax = self.d / 2      # (fig. 17.2)
        self.r_ay = np.inf          # (fig. 17.2)
        self.r_bxi = self.d_i / 2   # (fig. 17.2)
        self.r_bxo = -self.d_o / 2  # (fig. 17.2)
        self.r_byi = -np.inf        # (fig. 17.2)
        self.r_byo = np.inf         # (fig. 17.2)


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
            R_xi (float) [m] - Effective radius in the x direction for the inner track
            R_xo (float) [m] - Effective radius in the x direction for the outer track
            R_yi (float) [m] - Effective radius in the y direction for the inner track
            R_yo (float) [m] - Effective radius in the y direction for the outer track
        '''
        self.R_xi = (self.r_ax * self.r_bxi) / (self.r_ax + self.r_bxi) # (17.4)
        self.R_xo = (self.r_ax * self.r_bxo) / (self.r_ax + self.r_bxo) # (17.4)
        self.R_yi = (self.r_ay * self.r_byi) / (self.r_ay + self.r_byi) # (17.5)
        self.R_yo = (self.r_ay * self.r_byo) / (self.r_ay + self.r_byo) # (17.5)
        if np.isnan(self.R_yi):
            self.R_yi = np.inf
        if np.isnan(self.R_yo):
            self.R_yo = np.inf


        if printbool:
            print(f'Effective radius for inner track in x: R_xi = {self.R_xi:.4g} m')
            print(f'Effective radius for outer track in x: R_xo = {self.R_xo:.4g} m')
            print(f'Effective radius for inner track in y: R_yi = {self.R_yi:.4g} m')
            print(f'Effective radius for outer track in y: R_yo = {self.R_yo:.4g} m')


    def curvature(self, printbool:bool=False)->None:
        '''
        Calculates the curvature sum, curvature difference and radius ratio
        
        Args:
            printbool (bool) - Prints the curvature difference if True

        Attributes:
            R_i (float) - Curvature sum for the inner track
            R_o (float) - Curvature sum for the outer track
            R_di (float) - Curvature difference for the inner track
            R_do (float) - Curvature difference for the outer track
            alpha_ri (float) - Radius ratio for the inner track
            alpha_ro (float) - Radius ratio for the outer track
            
        Raises:
            ValueError - The radius ratio smaller than one - eq.(17.1) is not satisfied
        '''
        self.R_i = (self.R_xi * self.R_yi) / (self.R_xi + self.R_yi) # (17.2)
        self.R_o = (self.R_xo * self.R_yo) / (self.R_xo + self.R_yo)
        if np.isnan(self.R_i):
            self.R_i = np.inf
        if np.isnan(self.R_o):
            self.R_o = np.inf
        
        self.R_di = self.R_i * (1 / self.R_xi - 1 / self.R_yi) # (17.3)
        self.R_do = self.R_o * (1 / self.R_xo - 1 / self.R_yo) # (17.3)

        self.alpha_ri = self.R_yi / self.R_xi # (16.57)
        self.alpha_ro = self.R_yo / self.R_xo

        if self.alpha_ri < 1:
            raise ValueError('Radius ratio for inner track is smaller than one - eq.(17.1) is not satisfied')
        if self.alpha_ro < 1:
            raise ValueError('Radius ratio for outer track is smaller than one - eq.(17.1) is not satisfied')

        if printbool:
            print(f'Curvature sum for the inner track: R_i = {self.R_i:.4g} m and for the outer track: R_o = {self.R_o:.4g} m')
            print(f'Curvature difference for the inner track: R_di = {self.R_di:.4g} m and for the outer track: R_do = {self.R_do:.4g} m')
            print(f'Radius ratio: alpha_ri = {self.alpha_ri:.4g} and alpha_ro = {self.alpha_ro:.4g}')


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
            k_i (float) - Ellipticity parameter for the inner track
            k_o (float) - Ellipticity parameter for the outer track
            epsilon_i (float) - complete elliptic integral of second kind for the inner track
            epsilon_o (float) - complete elliptic integral of second kind for the outer track
            F_i (float) - complete elliptic integral of first kind for the inner track
            F_o (float) - complete elliptic integral of first kind for the outer track
            
        Raises:
            ValueError - If the iterative method does not converge
        '''
        # Inner track
        k_i = np.zeros(max_iter+1)
        k_i[0] = k_init

        phi = np.linspace(0, np.pi/2, phi_no)

        for i, k_current in enumerate(k_i[:-1]):
            F_i = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**(-0.5)) # (17.10)
            epsilon_i = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**0.5, phi) # (17.11)
            k_i[i+1] = ( (2*F_i - epsilon_i*(1+self.R_di)) /(epsilon_i * (1-self.R_di)) )**(1/2) # (17.10)

            if abs(k_i[i] - k_i[i+1]) < convergence_tol:
                break

        np.trim_zeros(k_i)
        self.k_i = k_i[-1]
        self.F_i = F_i
        self.epsilon_i = epsilon_i
        iter_i = np.arange(len(k_i))

        # Outer track
        k_o = np.zeros(max_iter+1)
        k_o[0] = k_init

        for i, k_current in enumerate(k_o[:-1]):
            F_o = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**(-0.5)) # (17.10)
            epsilon_o = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**0.5, phi) # (17.11)
            k_o[i+1] = ( (2*F_o - epsilon_o*(1+self.R_do)) /(epsilon_o * (1-self.R_do)) )**(1/2) # (17.10)

            if abs(k_o[i] - k_o[i+1]) < convergence_tol:
                break
        
        np.trim_zeros(k_o)
        self.k_o = k_o[-1]
        self.F_o = F_o
        self.epsilon_o = epsilon_o
        iter_o = np.arange(len(k_o))

        if plotbool:
            plt.plot(iter_i, k_i, label='Inner track')
            plt.plot(iter_o, k_o, label='Outer track')
            plt.title('Convergence of the ellipticity parameter')
            plt.xlabel('Iteration')
            plt.ylabel('Ellipticity parameter')
            plt.title('Convergence of the ellipticity parameter')
        
        if len(iter_i) == max_iter or len(iter_o) == max_iter:
            raise ValueError('Iterative method did not converge')
        
        if printbool:
            print(f'Ellipticity parameter for inner track: k_i = {self.k_i:.4g}')
            print(f'Ellipticity parameter for outer track: k_o = {self.k_o:.4g}')
            print(f'Complete elliptic integral for inner track: Epsilon_i = {self.epsilon_i:.4g}')
            print(f'Complete elliptic integral for outer track: Epsilon_o = {self.epsilon_o:.4g}')
            print(f'Complete elliptic integral for inner track: F_i = {self.F_i:.4g}')
            print(f'Complete elliptic integral for outer track: F_o = {self.F_o:.4g}')


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
            D_xi (float) [m] - Diameter of the contact ellipse in the x direction for the inner track
            D_xo (float) [m] - Diameter of the contact ellipse in the x direction for the outer track
            D_yi (float) [m] - Diameter of the contact ellipse in the y direction for the inner track
            D_yo (float) [m] - Diameter of the contact ellipse in the y direction for the outer track
        '''
        self.D_xi = 2 * (6 * self.epsilon_i * self.w_zm * self.R_i/(np.pi * self.k_i * self.E_prime))**(1/3) # (17.15)
        self.D_xo = 2 * (6 * self.epsilon_o * self.w_zm * self.R_o/(np.pi * self.k_o * self.E_prime))**(1/3) # (17.15)

        self.D_yi = 2 * (6 * self.k_i**2 * self.epsilon_i * self.w_zm * self.R_i/(np.pi * self.E_prime))**(1/3) # (17.16)
        self.D_yo = 2 * (6 * self.k_o**2 * self.epsilon_o * self.w_zm * self.R_o/(np.pi * self.E_prime))**(1/3) # (17.16)

        if printbool:
            print(f'Diameter of the contact ellipse in x: D_xi = {self.D_xi:.4g} m for inner track and D_xo = {self.D_xo:.4g} m for outer track')
            print(f'Diameter of the contact ellipse in y: D_yi = {self.D_yi:.4g} m for inner track and D_yo = {self.D_yo:.4g} m for outer track')


    def maximum_pressure(self, printbool:bool=False)->None:
        '''
        Calculates the maximum contact pressure in the bearing.

        Args:
            printbool (bool) - Prints the maximum contact pressure if True
        
        Attributes:
            p_mi (float) [Pa] - Maximum contact pressure for the inner track
            p_mo (float) [Pa] - Maximum contact pressure for the outer track
        '''
        self.p_mi = 6 * self.w_zm / (np.pi * self.D_xi * self.D_yi) # (17.7)
        self.p_mo = 6 * self.w_zm / (np.pi * self.D_xo * self.D_yo) # (17.7)

        if printbool:
            print(f'Maximum contact pressure for inner track: p_mi = {self.p_mi:.4g} Pa')
            print(f'Maximum contact pressure for outer track: p_mo = {self.p_mo:.4g} Pa')


    def maximum_deform(self, printbool:bool=False)->None:
        '''
        Calculates the maximum elastic deformation
        
        Args:
            printbool (bool) - Prints the maximum elastic deformation if True
        
        Attributes:
            delta_mi (float) [m] - Maximum elastic deformation for the inner track
            delta_mo (float) [m] - Maximum elastic deformation for the outer track
        '''
        self.delta_mi = self.F_i * (9/(2 * self.epsilon_i * self.R_i) * (self.w_zm / (np.pi * self.k_i * self.E_prime))**2)**(1/3) # (17.16)
        self.delta_mo = self.F_o * (9/(2 * self.epsilon_o * self.R_o) * (self.w_zm / (np.pi * self.k_o * self.E_prime))**2)**(1/3) # (17.16)

        if printbool:
            print(f'Maximum elastic deformation for inner track: delta_mi = {self.delta_mi:.4g} m')
            print(f'Maximum elastic deformation for outer track: delta_mo = {self.delta_mo:.4g} m')


    def rectangular_dimensionless_load(self, printbool:bool=False)->None:
        '''
        Calculates the dimensionless load and load per unit width of rectangular contact
        
        Args:
            printbool (bool) - Prints the contact semi-width if True
        
        Attributes:
            W_primei (float) - dimensionless load for the inner track
            W_primeo (float) - dimensionless load for the outer track
            w_x_prime (array) [N/m] - load per unit width
        '''
        self.w_x_prime = self.w_zm / self.l # (p. 448)
        
        self.W_primei = self.w_x_prime / (self.E_prime * self.R_xi) # (17.38)
        self.W_primeo = self.w_x_prime / (self.E_prime * self.R_xo) # (17.38)

        if printbool:
            print(f'Dimensionless load: for inner track W_prime = {self.W_primei:.4g}, for outer track W_prime = {self.W_primeo:.4g}')
            print(f'Load per unit width: w_x_prime = {self.w_x_prime:.4g} N/m')


    def rectangular_max_deformation(self, printbool:bool=False)->None:
        '''
        Calculates the maximum elastic deformation of the rectangular contact
        
        Args:
            printbool (bool) - Prints the maximum elastic deformation if True
        
        Attributes:
            delta_mi (float) [m] - Maximum elastic deformation for the inner track
            delta_mo (float) [m] - Maximum elastic deformation for the outer track
        '''
        self.delta_mi = 2 * self.W_primei * self.R_xi / np.pi * (np.log(2 * np.pi/self.W_primei) - 1) # (17.39)
        self.delta_mo = 2 * self.W_primeo * self.R_xo / np.pi * (np.log(2 * np.pi/self.W_primeo) - 1) # (17.39)

        if printbool:
            print(f'Maximum elastic deformation: for inner track: delta_mi = {self.delta_mi:.4g} m, for outer track: delta_mo = {self.delta_mo:.4g} m')


    def rectangular_max_pressure(self, printbool:bool=False)->None:
        '''
        Calculates the maximum contact pressure of the rectangular contact
        
        Args:
            printbool (bool) - Prints the maximum contact pressure if True
        
        Attributes:
            p_mi (float) [Pa] - Maximum contact pressure for the inner track
            p_mo (float) [Pa] - Maximum contact pressure for the outer track
        '''
        self.p_mi = self.E_prime * (self.W_primei/(2 * np.pi))**(1/2) # (17.40)
        self.p_mo = self.E_prime * (self.W_primeo/(2 * np.pi))**(1/2) # (17.40)

        if printbool:
            print(f'Maximum contact pressure: for inner track p_mi = {self.p_mi:.4g} Pa, for outer track p_mo = {self.p_mo:.4g} Pa')

if __name__ == '__main__':
    if False: # Test
        bearing = Bearing()
        bearing.max_load(printbool=True)
        bearing.min_film_thickness(printbool=True)
        bearing.effective_radius(printbool=True)

    if False: # Question 1
        print('Question 1')
        Q1 = Bearing()
        Q1.max_load(printbool=True)
        Q1.min_film_thickness(printbool=True)


    if True: # Question 2
        print('Question 2')
        Q2 = Bearing()
        Q2.max_load()
        Q2.effective_radius()
        Q2.effective_elastic_modulus()
        Q2.rectangular_dimensionless_load()
        Q2.rectangular_max_deformation(printbool=True)
        Q2.rectangular_max_pressure(printbool=True)







