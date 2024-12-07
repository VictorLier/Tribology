import numpy as np
import matplotlib.pyplot as plt
import copy
import scipy.sparse as sps

class Bearing:
    '''
    Cylindrical bearing class for Assignment 3
    '''
    def __init__(self, InnerRaceDia:float=69.5e-3, OuterRaceDia:float=85.5e-3, RollerDia:float=8e-3, RollerLength:float=8e-3, NoOfRollers:int=20, MaxLoad:float=8000, InnerRaceSpeed:float=837, OuterRaceSpeed:float=0, EModulus:float=210e9, PoissonRatio:float=0.3, AbsoluteViscosity:float=0.01, PressViscCoef:float=2e-8, N:int=100)->None:
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
            N (int) - Number of points to evaluate the elliptic integrals (default: 100)
        
        Attributes:
            D (array) [m] - Inner and outer race diameter
            d (float) [m] - Diameter of the rollers
            l (float) [m] - Length of the rollers
            i (int) - Number of rollers
            w_z (float) [N] - Maximum radial load
            omega (array) [rad/s] - Angular velocity of inner and outer race
            E (float) [Pa] - Young's Modulus
            nu (float) - Poisson's ratio
            eta_0 (float) [Pa s] - Base absolute viscosity of the lubricant
            xi (float) [1 / (Pa s)] - Pressure-viscosity coefficient
            N (int) - Number of points to evaluate the elliptic integrals
            r_ax (float) [m] - Curvature radius of the roller in the x direction
            r_ay (float) [m] - Curvature radius in the roller y direction
            r_bx (array) [m] - Curvature radius of the inner, outer track in the x direction
            r_by (array) [m] - Curvature radius of the inner, outer track in the y direction
        '''
        self.D = np.array([InnerRaceDia, OuterRaceDia])
        self.d = RollerDia
        self.l = RollerLength
        self.n = NoOfRollers
        self.w_z = MaxLoad
        self.omega = np.array([InnerRaceSpeed, OuterRaceSpeed])
        self.E = EModulus
        self.nu = PoissonRatio
        self.eta_0 = AbsoluteViscosity
        self.xi = PressViscCoef
        self.N = N

        self.r_ax = self.d / 2      # (fig. 17.2)
        self.r_ay = np.inf          # (fig. 17.2)
        self.r_bx = np.array([self.D[0]/2,  -self.D[1]/2])    # (fig. 17.2)
        self.r_by = np.array([-np.inf, np.inf])            # (fig. 17.2)

        self.clearance()


    def clearance(self, printbool:bool=False)->None:
        '''
        Calculates the clearance of the bearing

        Args:
            printbool (bool) - Prints the clearance if True
        
        Attributes:
            c_d (float) [m] - Clearance of the bearing
            clear (bool) - True if the bearing has clearance
        '''
        self.c_d = self.D[1] - self.D[0] - 2 * self.d # (21.2)

        if self.c_d > 0:
            self.clear = True
        elif self.c_d == 0:
            self.clear = False
        else:
            raise ValueError('Negative clearance')

        if printbool:
            print(f'Bearing has clearance: {self.clear}')
            print(f'Clearance of the bearing: c_d = {self.c_d:.4g} m')


    def max_load(self, printbool:bool=False)->None:
        '''
        Calculates the load of the roller with the highest load
        Accounts for clearance
        
        Args:
            printbool (bool) - Prints the maximum load if True

        Attributes:
            w_zm (float) [N] - Maximum radial load
        '''
        if not self.clear:
            self.w_zm = self.w_z*4 / self.n   # (21.47)

        elif self.clear:
            Z_up = np.pi * (1 - self.c_d/2 * self.delta_m)**(3/2)
            Z_down = 2.491 * ((1 + ( (1 - self.c_d/2 * self.delta_m)/1.23 )**2 )**(1/2) - 1)
            Z_w = Z_up / Z_down # (21.46)
            self.w_zm = self.w_z * Z_w / self.n # (21.45)
            self.w_zm = self.w_zm[0] # Both sides should be the same

        if printbool:
            print(f'Load of the roller with the highest load: w_zm = {self.w_zm:.3g} N')
        

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

        self.h_min = Lambda * (R_q[0]**2 + R_q[1]**2)**0.5 # (3.22)

        if printbool:
            print(f'Minimum film thickness: h_min = {self.h_min:.3g} m')


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
        if self.r_by[0] == -np.inf:
            self.R_y = np.inf
        else:
            self.R_y = (self.r_ay * self.r_by) / (self.r_ay + self.r_by) # (17.5)

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
        self.w_zm_prime = self.w_zm / self.l # (p. 448)        
        self.W_prime = self.w_zm_prime / (self.E_prime * self.R_x) # (17.38)

        if printbool:
            print(f'Dimensionless load for inner track: W_primei = {self.W_prime[0]:.3g}')
            print(f'Dimensionless load for outer track: W_primeo = {self.W_prime[1]:.3g}')
            print(f'Load per unit width: w_x_prime = {self.w_zm_prime:.3g} N/m')


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


    def rectangular_pressure_distribution(self, plotbool:bool=False)->None:
        '''
        Finds the pressure distribution of the rectangular contact

        Args:
            plotbool (bool) - Plots the pressure distribution if True
        
        Attributes:
            x (array) [m] - x coordinate for inner, outer track
            p (array) [Pa] - Pressure distribution for inner, outer track
        '''
        self.x = np.linspace(-self.D_x, self.D_x, self.N)
        self.p = self.p_m * np.sqrt(1 - (2*self.x/self.D_x)**2) # (17.6) - dy inf
        self.p = np.nan_to_num(self.p) # Replace nan with zero

        if plotbool:
            plt.figure()
            plt.title('Pressure distribution')
            plt.plot(self.x[:,0], self.p[:,0], label='Inner track')
            plt.plot(self.x[:,1], self.p[:,1], label='Outer track')
            plt.xlabel('x [m]')
            plt.ylabel('p [Pa]')
            plt.grid()
            plt.legend()

            np.savetxt('Assignment3/data/pressure_distribution_inner.txt', np.array([self.x[:,0], self.p[:,0]]).T)
            np.savetxt('Assignment3/data/pressure_distribution_outer.txt', np.array([self.x[:,1], self.p[:,1]]).T) 


    def rectangular_dimensionless_deforamtion(self, plotbool:bool=False)->None:
        '''
        Calculates the dimensionless eleastic deformation of the rectangular contact

        Args:
            plotbool (bool) - Plots the dimensionless deformation and pressure if True
        
        Attributes:
            X (array) - Dimensionless x coordinate
            delta_bar (array) - Dimensionless elastic deformation
        '''
        X = np.linspace(-2,2,self.N)
        extra_points = np.array([X[0] - (X[1] - X[0]), X[-1] + (X[1] - X[0])])
        self.X = np.concatenate(([extra_points[0]], X, [extra_points[1]]))
        P = np.sqrt(1-self.X**2) # Dimensionless parabolic pressure distribution - Hertz distribution - (17.6)
        P = np.nan_to_num(P, 0) # Replace nan with zero
        delta = 0
        self.delta_bar = np.zeros(self.N+2)
        for i in range(1, self.N+1):
            for j in range(self.N+2):
                term1 = np.abs((self.X[i+1] + self.X[i])/2 - self.X[j])
                term2 = np.abs((self.X[i-1] + self.X[i])/2 - self.X[j])
                delta_new = delta + P[j] * np.log(term1 * term2) # (18.31)
                delta = delta_new
            Delta = self.X[i+1] - self.X[i]
            self.delta_bar[i] = - Delta / (2*np.pi) * delta
            delta = 0
        self.X = self.X[1:-1]
        self.delta_bar = self.delta_bar[1:-1]
        P = P[1:-1]

        if plotbool:
            plt.figure()
            plt.title('Dimensionless elastic deformation')
            plt.plot(self.X, self.delta_bar)
            plt.xlabel('X')
            plt.ylabel('delta_bar')
            plt.grid()

            plt.figure()
            plt.title('Dimensionless pressure distribution')
            plt.plot(self.X, P)
            plt.xlabel('X')
            plt.ylabel('P')
            plt.grid()

            np.savetxt('Assignment3/data/dimensionless_deformation.txt', np.array([self.X, self.delta_bar]).T)
            np.savetxt('Assignment3/data/dimensionless_pressure.txt', np.array([self.X, P]).T)


    def rectangular_deformation(self, plotbool:bool=False)->None:
        '''
        Calculates the deformation of the rectangular contact from the dimensionless deformation
        
        Args:
            plotbool (bool) - Plots the deformation if True

        Attributes:
            delta (array) [m] - Elastic deformation of the inner, outer rectangular contact
        '''
        self.delta = np.zeros((self.N,2))
        self.delta[:,0] = (self.D_x[0]/2)**2 * self.delta_bar / self.R_x[0] # (17.23)
        self.delta[:,1] = (self.D_x[1]/2)**2 * self.delta_bar / self.R_x[1] # (17.23)

        if plotbool:
            plt.figure()
            plt.title('Elastic deformation of the rectangular contact')
            plt.plot(self.x[:,0], self.delta[:,0], label='Inner track')
            plt.plot(self.x[:,1], self.delta[:,1], label='Outer track')
            plt.xlabel('x [m]')
            plt.ylabel('delta [m]')
            plt.grid()
            plt.legend()

            np.savetxt('Assignment3/data/deformation_inner.txt', np.array([self.x[:,0], self.delta[:,0]]).T)
            np.savetxt('Assignment3/data/deformation_outer.txt', np.array([self.x[:,1], self.delta[:,1]]).T)
        

    def rectangular_roller_geometry(self, plotbool:bool=False)->None:
        '''
        Plots the geometry of the roller with the rectangular contact

        Args:
            plotbool (bool) - Plots the roller geometry if True
        '''
        S = self.x**2 / (2 * self.R_x)
        h = S + self.delta

        if plotbool:
            plt.figure()
            plt.title('Inner track roller geometry')
            plt.plot(self.x[:,0], S[:,0], label='Undeformed roller', linestyle='--')
            plt.plot(self.x[:,0], h[:,0], label='Deformed roller', linestyle='-')

            plt.xlabel('x [m]')
            plt.ylabel('y [m]')
            plt.legend()
            plt.grid()

            plt.figure()
            plt.title('Outer track roller geometry')
            plt.plot(self.x[:,1], S[:,1], label='Undeformed roller', linestyle='--')
            plt.plot(self.x[:,1], h[:,1], label='Deformed roller', linestyle='-')

            plt.xlabel('x [m]')
            plt.ylabel('y [m]')
            plt.legend()
            plt.grid()

            np.savetxt('Assignment3/data/undeformed_inner.txt', np.array([self.x[:,0], S[:,0]]).T)
            np.savetxt('Assignment3/data/deformed_inner.txt', np.array([self.x[:,0], h[:,0]]).T)

            np.savetxt('Assignment3/data/undeformed_outer.txt', np.array([self.x[:,1], S[:,1]]).T)
            np.savetxt('Assignment3/data/deformed_outer.txt', np.array([self.x[:,1], h[:,1]]).T)


    def velocity(self, printbool:bool=False)->None:
        '''
        Calculates the velocity of the roller and the mean surface velocity

        Args:
            printbool (bool) - Prints the velocity of the roller if True
        
        Attributes:
            u (array) [m/s] - Velocity of the inner and outer track
            u_bar (float) [m/s] - Mean surface velocity
        '''
        self.u = np.array([9,9])
        self.u_bar = np.mean(self.u)

        if printbool:
            print(f'Velocity of the inner track: u_i = {self.u[0]:.3g} m/s')
            print(f'Velocity of the outer track: u_o = {self.u[1]:.3g} m/s')
            print(f'Mean surface velocity: u_bar = {self.u_bar:.3g} m/s')


    def finite_difference(self, Nx:int=10, plotbool:bool=False, printbool:bool=False)->None:
        '''
        Solves the Reynolds equation with a finite difference method
        
        Args:
            Nx (int) - Number of points in each direction (default: 10)
            plotbool (bool) - Plots the solution if True
            printbool (bool) - Prints the solution if True
        '''
        X, Y = np.meshgrid(np.linspace(0, 1, Nx), np.linspace(0, 1, Nx))

        x_ = X.flatten('F')        
        
        dX = X[1] - X[0]
        dY = Y[1] - Y[0]

        h0 = 0.03e-6    # Assumed minimum film thickness [m]

        M = sps.eye(Nx**2, format='csr') # Matrix
        rhs = np.zeros(Nx**2) # Right hand side

        H = 1 + (self.d/2)/(2*h0) * (x_ - 1)**2
        lam = self.l / self.d

        r = self.d/2

        for i in range(1, Nx - 1):
            for j in range(1, Nx - 1):
                c = j + i * Nx
                n = c + 1
                s = c - 1
                e = j + Nx * (i + 1)
                w = j + Nx * (i - 1)

                M[c, c] = (-2 / dX**2 - 2 / (lam**2 * dY**2) - 3 * (((np.sqrt(H[e]) - np.sqrt(H[w])) * (H[e] - H[w])) / 4 + np.sqrt(H[c]) * (H[w] - 2 * H[c] + H[e])) / (2 * H[c]**(3/2) * dX**2))
                M[c, n] = 1 / (lam**2 * dY**2)
                M[c, s] = 1 / (lam**2 * dY**2)
                M[c, e] = 1 / dX**2
                M[c, w] = 1 / dX**2

                rhs[c] = 1 / (H[c]**(3/2)) * (H[e] - H[w]) / (2 * dX)

        g = sps.linalg.spsolve(M, rhs)
        G = np.reshape(g, (Nx, Nx))
        P = G * (H**(-3/2))
        p = (6 * self.eta_0 * (self.u[0] + self.u[1]) * r * P) / (h0**2)
        h = h0 * H

        # Visc
        p_new = -(1/self.xi) * np.log(np.abs(1 - self.xi * p))
        eta_new = self.eta_0 * np.exp(self.xi * p_new)

        if plotbool:
            plt.figure()
            plt.contourf(X, Y, h)
            plt.colorbar()
            plt.title('Film thickness')
            plt.xlabel('X')
            plt.ylabel('Y')

            plt.figure()
            plt.contourf(X, Y, p)
            plt.colorbar()
            plt.title('Pressure')
            plt.xlabel('X')
            plt.ylabel('Y')

            plt.figure()
            plt.contourf(X, Y, eta_new)
            plt.colorbar()
            plt.title('Viscosity')
            plt.xlabel('X')
            plt.ylabel('Y')


    def dimenionless_speed(self, printbool:bool=False)->None:
        '''
        Calculates the dimensionless speed parameter for the outer track

        Args:
            printbool (bool) - Prints the dimensionless speed if True
        
        Attributes:
            U (float) - Dimensionless speed parameter
        '''
        self.U = self.eta_0 * self.u_bar / (self.E_prime * self.R_x[1]**2) # 18.10

        if printbool:
            print(f'Dimensionless speed parameter: U = {self.U:.3g}')


    def dimensionless_load(self, printbool:bool=False)->None:
        '''
        Calculates the dimensionless load parameter for the outer track

        Args:
            printbool (bool) - Prints the dimensionless load if True
        
        Attributes:
            W (float) - Dimensionless load parameter
        '''
        self.W = self.w_zm / (self.E_prime * self.R_x[1])   # 18.11

        if printbool:
            print(f'Dimensionless load parameter: W = {self.W:.3g}')


    def dimensionless_material(self, printbool:bool=False)->None:
        '''
        Calculates the dimensionless material parameter for the outer track

        Args:
            printbool (bool) - Prints the dimensionless material parameter if True
        
        Attributes:
            G (float) - Dimensionless material parameter
        '''
        self.G = self.xi * self.E_prime  # 18.14

        if printbool:
            print(f'Dimensionless material parameter: G = {self.G:.3g}')


    def pressure_spike(self, printbool:bool=False)->None:
        '''
        Calculates the pressure spike amplitude and location from the operating parameters

        Args:
            printbool (bool) - Prints the pressure spike if True
        
        Attributes:
            p_sk (float) [Pa] - Pressure spike amplitude
            x_sk (float) [m] - Pressure spike location
        '''
        self.p_sk = 0.648 * self.W**(0.185) * self.U**(0.275) * self.G**(0.391) * self.E_prime # 18.70
        self.x_sk = 1.111 * self.W**(0.606) * self.U**(-0.021) * self.G**(0.077) * self.R_x[1] # 18.71

        if printbool:
            print(f'Pressure spike amplitude: p_sk = {self.p_sk:.3g} Pa')
            print(f'Pressure spike location: x_sk = {self.x_sk:.3g} m')


    def minimum_film_thickness(self, printbool:bool=False)->None:
        '''
        Calculates the minimum film thickness and location from the operating parameters

        Args:
            printbool (bool) - Prints the minimum film thickness if True
        
        Attributes:
            h_min (float) [m] - Minimum film thickness
            x_min (float) [m] - Minimum film thickness location
        '''
        self.h_min = 1.714 * self.W**(-0.128) * self.U**(0.694) * self.G**(0.568) * self.R_x[1] # 18.72
        self.x_min = 1.439 * self.W**(0.548) * self.U**(-0.011) * self.G**(0.026) * self.R_x[1] # 18.75

        if printbool:
            print(f'Minimum film thickness: h_min = {self.h_min:.3g} m')
            print(f'Minimum film thickness location: x_min = {self.x_min:.3g} m')


    def center_of_pressure(self, printbool:bool=False)->None:
        '''
        Calculates the center of pressure from the operating parameters

        Args:
            printbool (bool) - Prints the center of pressure if True
        
        Attributes:
            x_cp (float) - Center of pressure
        '''
        self.x_cp = -3.595 * self.W**(-1.019) * self.U**(0.638) * self.G**(-0.358) * self.R_x[1] # 18.76

        if printbool:
            print(f'Center of pressure: x_cp = {self.x_cp:.3g} m')


    def center_film_thickness(self, printbool:bool=False)->None:
        '''
        Calculates the center film thickness from the operating parameters

        Args:
            printbool (bool) - Prints the center film thickness if True
        
        Attributes:
            h_c (float) - Center film thickness
        '''
        self.h_c = 2.922 * self.W**(-0.166) * self.U**(0.692) * self.G**(0.470) * self.R_x[1] # 18.74

        if printbool:
            print(f'Center film thickness: h_c = {self.h_c:.3g} m')






if __name__ == '__main__':
    if False: # Test
        print('Test')


    if False: # Question 1
        print('Question 1')
        Q1 = Bearing()
        Q1.max_load(printbool=True)
        Q1.min_film_thickness(printbool=True)


    if False: # Question 2
        print('Question 2')
        Q2 = Bearing(N=100)
        Q2.max_load()
        Q2.min_film_thickness()
        Q2.effective_radius()
        Q2.effective_elastic_modulus()
        Q2.rectangular_dimensionless_load()
        
        print("Part 1")
        Q2.rectangular_max_deformation(printbool=True)
        Q2.rectangular_max_pressure(printbool=True)

        print("Part 2")
        print(f"For the inner track the deformation is {Q2.delta_m[0]:.3g} m and the minimum film thickness is {Q2.h_min:.3g} m")

        print("Part 3")
        Q2.rectangular_dimensionless_deforamtion()
        Q2.rectangular_contact_width()
        Q2.rectangular_pressure_distribution(plotbool=True)
        Q2.rectangular_deformation(plotbool=True)
        Q2.rectangular_roller_geometry(plotbool=True)

        print("Part 5")
        Q2.D[1] = 85.51e-3    # [m]
        Q2.clearance(printbool=True)
        Q2.max_load(printbool=True)
        
        plt.show()


    if False: # Question 3
        print('Question 3')
        Q3 = Bearing()
        Q3.velocity(printbool=True)
        Q3.finite_difference(Nx=10, plotbool=True)


    if False: # Question 4
        print('Question 4')
        Q4 = Bearing()
        Q4.max_load()
        Q4.effective_elastic_modulus()
        Q4.effective_radius()
        Q4.velocity()
        Q4.dimenionless_speed()
        Q4.dimensionless_load()
        Q4.dimensionless_material()

        print("Part 1")
        Q4.pressure_spike(printbool=True)
        Q4.minimum_film_thickness(printbool=True)
        Q4.center_of_pressure(printbool=True)
        Q4.center_film_thickness(printbool=True)

        print("Part 2")
        plt.figure()
        plt.plot()
        plt.scatter(Q4.x_sk, Q4.p_sk, label='Pressure spike')
        plt.scatter(Q4.x_min, Q4.h_min, label='Minimum film thickness')
        plt.scatter(Q4.x_cp, 0, label='Center of pressure')
        plt.scatter(Q4.x_cp, Q4.h_c, label='Center film thickness')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.legend()
        plt.grid()
        plt.show()


