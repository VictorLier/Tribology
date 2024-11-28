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
            r_bx (array) [m] - Curvature radius of the inner, outer track in the x direction
            r_by (array) [m] - Curvature radius of the inner, outer track in the y direction
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
        self.r_bx = np.array([self.d_i/2,  -self.d_o/2])    # (fig. 17.2)
        self.r_by = np.array([-np.inf, np.inf])            # (fig. 17.2)


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
            print(f'Maximum load on a roller: w_zm = {self.w_zm:.3g} N')


    def min_film_thickness(self, printbool:bool=False, Lambda:float=3, R_at:float=0.3e-6, R_ar:float=0.12e-6)->None:
        '''
        Calculates the minimum film thickness assuming worst case scenario in the Elasto-Hydrodynamic Lubrication (EHL) regime
        
        Args:
            printbool (bool) - Prints the minimum film thickness if True
            Lambda (float) - Non-dimensional film parameter (default: 3) - Lowest in the EHL regime - p. 57
            R_at (float) [m] - Roughness of the track (default: 0.3e-6) - Worst case from figure 3.8
            R_ar (float) [m] - Roughness of the roller (default: 0.12e-6) - Worst case from figure 3.
        '''
        Lambda
        R_a = np.array([R_at, R_ar])
        R_q = R_a * 1.11 # (3.5)

        h_min = Lambda * (R_q[0]**2 + R_q[1]**2)**0.5 # (3.22)

        if printbool:
            print(f'Minimum film thickness: h_min = {h_min:.3g} m')


    def effective_radius(self, printbool:bool=False)->None:
        '''
        Calculates the effective radius in the x and y direction. Will handle infinite values
        
        Args:
            printbool (bool) - Prints the effective radius if True
        
        Attributes:
            R_x (array) [m] - Effective radius of inner, outer track in the x direction
            R_y (array) [m] - Effective radius of inner, outer track in the y direction
        '''
        self.R_x = (self.r_ax * self.r_bx) / (self.r_ax + self.r_bx) # (17.4)
        self.R_y = (self.r_ay * self.r_by) / (self.r_ay + self.r_by) # (17.5)

        # self.R_x = np.nan_to_num(self.R_x, nan=np.inf)  # Replace nan with inf
        # self.R_y = np.nan_to_num(self.R_y, nan=np.inf)  # Replace nan with inf

        if printbool:
            print(f'Effective radius for inner track in x: R_xi = {self.R_x[0]:.3g} m')
            print(f'Effective radius for outer track in x: R_xo = {self.R_x[1]:.3g} m')
            print(f'Effective radius for inner track in y: R_yi = {self.R_y[0]:.3g} m')
            print(f'Effective radius for outer track in y: R_yo = {self.R_y[1]:.3g} m')


    def curvature(self, printbool:bool=False)->None:
        '''
        Calculates the curvature sum, curvature difference and radius ratio
        
        Args:
            printbool (bool) - Prints the curvature difference if True

        Attributes:
            R (array) [m] - Curvature sum for the inner, outer track
            R_d (array) [m] - Curvature difference for the inner, outer track
            alpha_r (array) - Radius ratio for the inner, outer track
            
        Raises:
            ValueError - The radius ratio smaller than one - eq.(17.1) is not satisfied
        '''
        self.R = (self.R_x * self.R_y) / (self.R_x + self.R_y) # (17.2)
        self.R_d = self.R * (1 / self.R_x - 1 / self.R_y) # (17.3)
        self.alpha_r = self.R_y / self.R_x # (16.57)

        if self.alpha_r[0] < 1:
            raise ValueError('Radius ratio for inner track is smaller than one - eq.(17.1) is not satisfied')
        if self.alpha_r[1] < 1:
            raise ValueError('Radius ratio for outer track is smaller than one - eq.(17.1) is not satisfied')

        if printbool:
            print(f'Curvature sum for the inner track: R_i = {self.R[0]:.3g} m')
            print(f'Curvature sum for the outer track: R_o = {self.R[1]:.3g} m')
            print(f'Curvature difference for the inner track: R_di = {self.R_d[0]:.3g} m')
            print(f'Curvature difference for the outer track: R_do = {self.R_d[1]:.3g} m')
            print(f'Radius ratio for the inner track: alpha_ri = {self.alpha_r[0]:.3g}')
            print(f'Radius ratio for the outer track: alpha_ro = {self.alpha_r[1]:.3g}')


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
            k (array) - Ellipticity parameter for the inner, outer track
            epsilon (array) - complete elliptic integral of second kind for the inner, outer track
            F (array) - complete elliptic integral of first kind for the inner, outer track
            
        Raises:
            ValueError - If the iterative method does not converge
        '''
        self.k = np.zeros(2)
        self.epsilon = np.zeros(2)
        self.F = np.zeros(2)

        # Inner track
        k_i = np.zeros(max_iter+1)
        k_i[0] = k_init

        phi = np.linspace(0, np.pi/2, phi_no)

        for i, k_current in enumerate(k_i[:-1]):
            F_i = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**(-0.5)) # (17.10)
            epsilon_i = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**0.5, phi) # (17.11)
            k_i[i+1] = ( (2*F_i - epsilon_i*(1+self.R_d[0])) /(epsilon_i * (1-self.R_d[0])) )**(1/2) # (17.10)

            if abs(k_i[i] - k_i[i+1]) < convergence_tol:
                break

        np.trim_zeros(k_i)
        self.k[0] = k_i[-1]
        self.F[0] = F_i
        self.epsilon[0] = epsilon_i
        iter_i = np.arange(len(k_i))

        # Outer track
        k_o = np.zeros(max_iter+1)
        k_o[0] = k_init

        for i, k_current in enumerate(k_o[:-1]):
            F_o = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**(-0.5)) # (17.10)
            epsilon_o = np.trapezoid((1 - (1 - 1/k_current**2) * np.sin(phi)**2)**0.5, phi) # (17.11)
            k_o[i+1] = ( (2*F_o - epsilon_o*(1+self.R_d[1])) /(epsilon_o * (1-self.R_d[1])) )**(1/2) # (17.10)

            if abs(k_o[i] - k_o[i+1]) < convergence_tol:
                break
        
        np.trim_zeros(k_o)
        self.k[1] = k_o[-1]
        self.F[1] = F_o
        self.epsilon[1] = epsilon_o
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
            print(f'Ellipticity parameter for inner track: k_i = {self.k_i:.3g}')
            print(f'Ellipticity parameter for outer track: k_o = {self.k_o:.3g}')
            print(f'Complete elliptic integral for inner track: Epsilon_i = {self.epsilon_i:.3g}')
            print(f'Complete elliptic integral for outer track: Epsilon_o = {self.epsilon_o:.3g}')
            print(f'Complete elliptic integral for inner track: F_i = {self.F_i:.3g}')
            print(f'Complete elliptic integral for outer track: F_o = {self.F_o:.3g}')


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
            print(f'Effective elastic modulus: E_prime = {self.E_prime:.3g} Pa')


    def rectangular_dimensionless_load(self, printbool:bool=False)->None:
        '''
        Calculates the dimensionless load and load per unit width of rectangular contact
        
        Args:
            printbool (bool) - Prints the contact semi-width if True
        
        Attributes:
            W_prime (array) - Dimensionless load for the inner, outer track
            w_x_prime (array) [N/m] - load per unit width
        '''
        self.w_x_prime = self.w_zm / self.l # (p. 448)        
        self.W_prime = self.w_x_prime / (self.E_prime * self.R_x) # (17.38)

        if printbool:
            print(f'Dimensionless load for inner track: W_primei = {self.W_prime[0]:.3g}')
            print(f'Dimensionless load for outer track: W_primeo = {self.W_prime[1]:.3g}')
            print(f'Load per unit width: w_x_prime = {self.w_x_prime:.3g} N/m')


    def rectangular_contact_width(self, printbool:bool=False)->None:
        '''
        Calculates the contact width of the rectangular contact

        attributes:
            D_x (array) [m] - Contact width for the inner, outer track
        '''
        self.D_x = 2 * self.R_x * (8 * self.W_prime / np.pi)**0.5 # (17.37)

        if printbool:
            print(f'Contact width for inner track: b_i = {self.D_x[0]:.3g} m')
            print(f'Contact width for outer track: b_o = {self.D_x[1]:.3g} m')


    def rectangular_max_deformation(self, printbool:bool=False)->None:
        '''
        Calculates the maximum elastic deformation of the rectangular contact
        
        Args:
            printbool (bool) - Prints the maximum elastic deformation if True
        
        Attributes:
            delta_m (array) [m] - Maximum elastic deformation for the inner, outer track
        '''
        self.delta_m = 2 * self.W_prime * self.R_x / np.pi * (np.log(2 * np.pi/self.W_prime) - 1) # (17.39)

        if printbool:
            print(f'Maximum elastic deformation: for inner track: delta_mi = {self.delta_m[0]:.3g} m')
            print(f'Maximum elastic deformation: for outer track: delta_mo = {self.delta_m[1]:.3g} m')


    def rectangular_max_pressure(self, printbool:bool=False)->None:
        '''
        Calculates the maximum contact pressure of the rectangular contact
        
        Args:
            printbool (bool) - Prints the maximum contact pressure if True
        
        Attributes:
            p_m (array) [Pa] - Maximum contact pressure for the inner, outer track
        '''
        self.p_m = self.E_prime * (self.W_prime/(2 * np.pi))**(1/2) # (17.40)

        if printbool:
            print(f'Maximum contact pressure: for inner track p_mi = {self.p_m[0]:.3g} Pa')
            print(f'Maximum contact pressure: for outer track p_mo = {self.p_m[1]:.3g} Pa')


    def rectangular_dimensionless_deforamtion(self, n_x:int=100, plotbool:bool=False)->None:
        '''
        Calculates the dimensionless eleastic deformation of the rectangular contact

        Args:
            n_x (int) - Number of points to evaluate the deformation (default: 100)
        
        Attributes:
            X (array) - Dimensionless x coordinate
            delta_bar (array) - Dimensionless elastic deformation
        '''
        self.X = np.linspace(-1,1,n_x)
        P = -1 * self.X**2 + 1 # Dimensionless parabolic pressure distribution
        delta = 0
        self.delta_bar = np.zeros(n_x)
        for i in range(1, n_x-1):
            for j in range(n_x):
                delta = delta + P[j] * np.log(np.abs((self.X[i+1] + self.X[i])/2 - self.X[j]) * np.abs((self.X[i-1] + self.X[i])/2 - self.X[j])) # (18.31)
            self.delta_bar[i] = - (self.X[i+1] - self.X[i]) / (2*np.pi) * delta
            delta = 0
        
        if plotbool:
            plt.figure()
            plt.title('Dimensionless elastic deformation')
            plt.plot(self.X, self.delta_bar)
            plt.xlabel('X')
            plt.ylabel('delta_bar')


    def rectangular_deformation(self, plotbool:bool=False)->None:
        '''
        Calculates the deformation of the rectangular contact from the dimensionless deformation
        
        Args:
            plotbool (bool) - Plots the deformation if True

        Attributes:
            delta (array) [m] - Elastic deformation of the inner, outer rectangular contact
        '''
        self.delta = np.zeros((2, len(self.X)))
        self.delta[0] = self.D_x[0]**2 * self.delta_bar / self.R_x[0] # (18.23)
        self.delta[1] = self.D_x[1]**2 * self.delta_bar / self.R_x[1] # (18.23)

        if plotbool:
            x1 = np.linspace(-self.D_x[0], self.D_x[0], len(self.X))
            x2 = np.linspace(-self.D_x[1], self.D_x[1], len(self.X))
            plt.figure()
            plt.title('Elastic deformation of the rectangular contact')
            plt.plot(x1, self.delta[0], label='Inner track')
            plt.plot(x2, self.delta[1], label='Outer track')
            plt.xlabel('X')
            plt.ylabel('delta [m]')
            plt.legend()


    def rectangular_pressure_distribution(self, normalized:bool=False, num_points:int=100, saveplot:bool=False)->None:
        '''
        Plots the pressure distribution of the rectangular contact

        Args:
            normalized (bool) - Normalizes the pressure distribution if True (default: False)
            num_points (int) - Number of points to evaluate the pressure distribution (default: 100)
            saveplot (bool) - Saves the plot if True (default: False)
        '''
        x = np.linspace(-self.l/2, self.l/2, num_points)
        p_xi = self.p_mi * np.sqrt(1 - (2*x/self.l)**2) # (17.46)
        p_xo = self.p_mo * np.sqrt(1 - (2*x/self.l)**2) # (17.46)

        if normalized:
            x /= self.l
            p_xi /= self.p_mi
            p_xo /= self.p_mo
            plt.xlabel('x/l')
            plt.ylabel('p/p_mi')    
        else:
            plt.xlabel('x [m]')
            plt.ylabel('p [Pa]')
    
        plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
        plt.plot(x, p_xi, label='Inner track')
        plt.plot(x, p_xo, label='Outer track')

        if saveplot:
            np.savetxt('Assignment3/data/inner_track_pressure.txt', np.array([x, p_xi]).T)
            np.savetxt('Assignment3/data/outer_track_pressure.txt', np.array([x, p_xo]).T)


    def rectangular_deformation_distribution(self, normalized:bool=False, num_points:int=100, saveplot:bool=False)->None:
        '''
        Plot the rectangular deformation distribution
        
        Args:
            normalized (bool) - Normalizes the deformation distribution if True (default: False)
            num_points (int) - Number of points to evaluate the deformation distribution (default: 100)
            saveplot (bool) - Saves the plot if True (default: False)
        '''
        X = np.linspace(-0.5, 0.5, num_points)
        Delta = X[1] - X[0]
        delta_ii = np.zeros(num_points)
        delta_io = np.zeros(num_points)
        P = np.sqrt(1 - (2*X)**2) # (17.46)
        




if __name__ == '__main__':
    if True: # Test
        bearing = Bearing()
        bearing.max_load(printbool=True)
        bearing.effective_elastic_modulus(printbool=True)
        bearing.effective_radius(printbool=True)
        bearing.rectangular_dimensionless_load(printbool=True)
        bearing.rectangular_contact_width(printbool=True)
        bearing.rectangular_dimensionless_deforamtion(n_x=5, plotbool=True)
        bearing.rectangular_deformation(plotbool=True)
        
        
        
        plt.show()



    if False: # Question 1
        print('Question 1')
        Q1 = Bearing()
        Q1.max_load(printbool=True)
        Q1.min_film_thickness(printbool=True)


    if False: # Question 2
        print('Question 2')
        Q2 = Bearing()
        Q2.max_load()
        Q2.effective_radius(printbool=True)
        Q2.effective_elastic_modulus(printbool=True)
        Q2.rectangular_dimensionless_load(printbool=True)
        Q2.rectangular_max_deformation(printbool=True)
        Q2.rectangular_max_pressure(printbool=True)
        Q2.rectangular_pressure_distribution(normalized=True)







