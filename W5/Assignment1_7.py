import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt

class AS7:
    def __init__(self, number_pads: int = 20, type: int = 0, notch_size: float = 0.1, outer_radius: float = 0.250/2, inner_radius: float = 0.095 , max_load: float = 4e5, kinematic_viscosity: float = 32e-6, minimum_height: float = 60e-6, density: float = 860, specific_heat: float = 2000, running_speed: float = 1000, shroud: bool = False) -> None:
        '''
        number_pads (int): number of pads of the bearing

        type (int): type of the bearing, 0 for fixed inclined pads, 1 for parralel-step

        notch_size (float): The proportional size of the notch (0 < notch_size < 1)

        outer_radius (float): The outer radius of the bearing [m]

        inner_radius (float): The inner radius of the bearing [m]

        max_load (float): The maximum load of the bearing [N]

        kinematic_viscosity (float): The kinematic viscosity of the oil at p = 0 [m^2/s]

        minimum_height (float): The minimum height of the bearing [m]

        density (float): The density of the oil [kg/m^3]

        specific_heat (float): The specific heat of the oil [J/kg K]

        running_speed (float): The running speed of the bearing [rpm]

        shroud (bool): If True, the bearing has a shroud
        '''
        if type != 0 and type != 1:
            raise ValueError('Type must be 0 or 1')
        if notch_size < 0 or notch_size > 1:
            raise ValueError('Notch size must be between 0 and 1')
        self.number_pads = number_pads
        self.type = type
        self.notch_size = notch_size
        self.r_o = outer_radius
        self.r_i = inner_radius
        self.r_m = (self.r_o + self.r_i) / 2    # The mean radius of the bearing
        self.b = self.r_o - self.r_i            # The width of the bearing
        self.max_load = max_load
        self.nu = kinematic_viscosity
        self.h_0 = minimum_height
        self.rho_0 = density
        self.eta_0 = self.nu * self.rho_0       # The dynamic viscosity
        self.Cp = specific_heat
        self.rpm = running_speed
        self.shroud = shroud
        self.calc_length()
        self.calc_velocity()
        self.load_per_unit_width()


    def calc_length(self) -> None:
        '''
        Calculates the length of the pad
        '''
        self.l = 2*np.pi * self.r_m/self.number_pads * (1 - self.notch_size) # The length of the pad


    def calc_velocity(self) -> None:
        '''
        Calculates the velocity of the bearing at the mid point of the pad
        '''
        # Er lidt i tvil om ub er korrekt
        self.u_b = 2 * np.pi * self.r_m * self.rpm / 60


    def load_per_unit_width(self) -> None:
        '''
        Calculates the load per unit width
        '''
        self.wm_z = self.max_load / self.number_pads / self.b # (8.1)


    def geometry_parameters(self, plot: bool = False, print_bol: bool = False):
        '''
        Prints and calcualets the expressions the "optimal" geometry parameters of the bearing

        plot (bool): If True, plots the shoulder height
        '''
        if self.type == 0: # Fixed inclined pads
            self.x = np.linspace(0, self.l, 25)
            self.H_0 = np.sqrt(2)/2 # 8.28
            self.W_z = 6*np.log((self.H_0 + 1)/self.H_0) - 12/(1+2*self.H_0) #8.30
            self.s_h = np.sqrt(self.wm_z * self.W_z * self.eta_0 * self.u_b) * self.l / self.wm_z
            self.h_0 = self.s_h * self.H_0
            self.h = self.h_0 + self.s_h * (1 - self.x/self.l)
            if print_bol:
                print(f"Fixed inclined pads: The optimal being: s_h = {self.s_h:.3g} and h_0 = {self.h_0:.3g}")
            if plot:
                plt.figure()
                plt.plot(self.x, self.h)
                plt.xlabel('x [m]')
                plt.ylabel('h [m]')
                plt.ylim(0, np.max(self.h)*1.1)
                plt.title('Shoulder height - Fixed inclined pads')
                plt.show()

                np.savetxt('W4/data/fixedProfile.txt', np.array([self.x, self.h]).T)

        if self.type == 1: # Parallel-step
            self.n_s = 0.7182
            self.H_0 = 1.155
            self.W_z = 3 * self.n_s * (1 - self.n_s) / ( (1 - self.n_s) * (self.H_0 + 1)**3 + self.n_s * self.H_0**3)
            self.s_h = np.sqrt(self.wm_z*self.W_z*self.eta_0*self.u_b)*self.l / self.wm_z
            self.h_0 = self.s_h * self.H_0


            self.x = np.linspace(0, self.l, 25)
            self.h = np.ones(len(self.x)) * self.h_0
            index = int(self.n_s * len(self.x))
            self.h[:index] = self.h_0 + self.s_h
            
            
            # self.x = np.array([0, self.n_s*self.l, self.n_s*self.l, self.l])
            # self.h = np.array([self.h_0+self.s_h, self.h_0+self.s_h, self.h_0, self.h_0])
            if print_bol:
                print(f"Parallel-step: h_0 = {self.h_0:.3g} and s_h = {self.s_h:.3g}")
            if plot:
                plt.figure()
                plt.plot(self.x, self.h)
                plt.xlabel('x [m]')
                plt.ylabel('h [m]')
                plt.ylim(0, np.max(self.h)*1.1)
                plt.title('Shoulder height - Parallel-step')
                plt.show()

                np.savetxt('W4/data/parallelProfile.txt', np.array([self.x, self.h]).T)


    def pressure_distrubution(self, plot = False, printbol = False) -> None:
        '''
        Calculates the pressure distribution with a finite difference method

        plot (bool): If True, plots the pressure distribution

        printbol (bool): If True, prints the load capacity
        '''
        dx = abs(self.x[0] - self.x[1])

        # Create the diagonals
        north = np.zeros(len(self.h))
        for i in range(len(self.h)-2):
            i = i+1
            north[i+1] = 3 * self.h[i]**2 * (self.h[i+1] - self.h[i-1]) / (2 * dx) * 1/(2*dx) + self.h[i]**3 * 1/(dx**2)
        
        middle = np.zeros(len(self.h))
        for i in range(len(self.h)):
            middle[i] = self.h[i]**3 * -2/(dx**2)
        
        south = np.zeros(len(self.h))
        for i in range(len(self.h)-2):
            i = i+1
            south[i-1] = 3 * self.h[i]**2 * (self.h[i+1] - self.h[i-1]) / (2 * dx) * -1/(2*dx) + self.h[i]**3 * 1/(dx**2)
        
        # Bondary conditions
        middle[0] = 1
        middle[-1] = 1

        # Create the matrix
        data = np.array([north, middle, south])
        A = sps.spdiags(data, [1, 0, -1], len(self.h), len(self.h), format = 'csr').T

        # rhs
        rhs = np.zeros(len(self.h))
        for i in range(len(self.h)-2):
            i=i+1
            rhs[i] = 6 * self.eta_0 * self.u_b * (self.h[i+1] - self.h[i-1]) / (2 * dx)

        # Solve the system
        self.p = sps.linalg.spsolve(A, rhs)

        # Calculate the load capacity
        self.F = np.trapezoid(self.p, self.x) 
        self.F = self.F * self.b
        if printbol:
            print(f"The load capacity is: {self.F:.3g} N")
            print(f"The total load capacity is: {self.F*self.number_pads:.3g} N")


        if plot:
            plt.figure()
            plt.plot(self.x, self.p)
            plt.xlabel('x [m]')
            plt.ylabel('p [Pa]')
            plt.title('Pressure distribution')
            plt.show()

            np.savetxt(f'W4/data/pressureDistribution{self.type}.txt', np.array([self.x, self.p]).T)


    def pressure_distrubution_2d(self, nx: int = 25, plot = False, printbol = False) -> None:
        '''
        Calculates the pressure distribution with a finite difference method in 2D
        
        n (int): Number of points in the x and y direction

        plot (bool): If True, plots the pressure distribution

        printbol (bool): If True, prints the load capacity
        '''

        self.X = np.linspace(0, self.l, nx)
        self.Y = np.linspace(0, self.b, nx)
        dx = abs(self.X[0] - self.X[1])
        dy = abs(self.Y[0] - self.Y[1])



        h = np.zeros((nx, nx))

        if self.type == 0: # Linear incline
            for i in range(nx):
                h[i,:] = self.h_0 + self.s_h * (1 - self.X/self.l)
        elif self.type == 1: # Step pads
            for i in range(nx):
                if self.Y[i] < self.l*self.n_s:
                    h[:,i] = self.h_0 + self.s_h
                else:
                    h[:,i] = self.h_0

        
        if self.shroud:
            for i in range(nx):
                for j in range(nx):
                    if self.Y[j] > self.X[i]*3: #or self.Y[j] > self.X[i]*1 - self.b:
                        h[i,j] = self.h_0
                    
                    if self.Y[j] > self.l*3 - self.X[i]*3:
                        h[i,j] = self.h_0


        M = sps.eye(nx**2)
        M = M.tocsr()

        rhs = np.zeros(nx**2)

        for i in range(1, nx-1):
            for j in range(1, nx-1):
                c = j + nx * (i)  # Current index
                n = c + 1  # North neighbor
                s = c - 1  # South neighbor
                e = j + nx * (i+1)  # East neighbor
                w = j + nx * (i - 1)  # West neighbor

                # Filling the M matrix
                M[c, c] = -2 * h[j, i]**3 * (1/(dx**2) + 1/(dy**2))
                M[c, n] = 3 * h.flatten('F')[c]**2 * (h.flatten('F')[n] - h.flatten('F')[s]) / (2 * dy**2) + h.flatten('F')[c]**3 / (dy**2)
                M[c, s] = -3 * h.flatten('F')[c]**2 * (h.flatten('F')[n] - h.flatten('F')[s]) / (2 * dy**2) + h.flatten('F')[c]**3 / (dy**2)
                M[c, e] = 3 * h.flatten('F')[c]**2 * (h.flatten('F')[e] - h.flatten('F')[w]) / (2 * dx**2) + h.flatten('F')[c]**3 / (dx**2)
                M[c, w] = -3 * h.flatten('F')[c]**2 * (h.flatten('F')[e] - h.flatten('F')[w]) / (2 * dx**2) + h.flatten('F')[c]**3 / (dx**2)

                # Filling the rhs vector
                rhs[c] = 6 * self.u_b * self.eta_0 * (h.flatten('F')[e] - h.flatten('F')[w]) / (2 * dx)
        
        p2 = sps.linalg.spsolve(M, rhs)

        self.p2 = p2.reshape(nx, nx, order = 'F')

        if plot:
            X, Y = np.meshgrid(self.X, self.Y)

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(X, Y, self.p2)
            
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111, projection='3d')
            ax1.plot_surface(X, Y, h)
            plt.show()
            plt.show()

            np.savetxt('W5/data/pressureDistributionShroudHeight.txt', np.array([X.flatten('F'), Y.flatten('F'), h.flatten('F')]).T)


        if printbol:
            self.F2 = np.trapezoid(np.trapezoid(self.p2, self.X), self.Y)
            print(f"The load capacity is: {self.F2:.3g} N")


    def film_thick(self, plot: bool = False, printbol: bool = False) -> None:
        '''
        Calculates the film thicknes of a finite width bearing with the intended loading

        plot (bool): If True, plots the film thickness

        printbol (bool): If True, prints the film thickness
        '''

        if self.type == 0: # Fixed inclined pads
            h0 = 5.5e-6
            H0 = h0 / self.s_h
            print("H0 er : ", H0)
            wsome = float(input("Wz: "))
            Wz = 1 / wsome

            wz = Wz * self.eta_0 * self.u_b * self.b * H0**2 * self.l**2 / h0**2

            force_error = wz - self.max_load/self.number_pads

        if self.type == 1: # Parallel-step
            h0 = 5.5e-6
            H0 = h0 / self.s_h
            print("H0 er : ", H0)
            Wz =  float(input("Wz (omkring 0.15): "))
            wz = Wz * self.eta_0 * self.u_b * self.b * self.l**2 / self.s_h**2

            force_error = wz - self.max_load/self.number_pads

        if printbol:
            print(f"The film thickness is: {h0:.3g} m")
            print(f"The force is: {wz:.3g} N")
            print(f"The force error is:", force_error)


    

if __name__ == '__main__':
    if False: # Part a
        print('Part a')
        print('Fixed inclined pads')
        inc = AS7(type = 0)
        inc.geometry_parameters()
        inc.pressure_distrubution(printbol = True)
        inc.pressure_distrubution_2d(printbol=True)

        print('Parallel-step')
        par = AS7(type = 1)
        par.geometry_parameters()
        par.pressure_distrubution(printbol = True)
        par.pressure_distrubution_2d(printbol=True)

    
    if False: # Part b
        inc = AS7(type = 0)
        inc.geometry_parameters(plot = False)
        inc.pressure_distrubution(plot = False)
        inc.pressure_distrubution_2d(plot = False)

        X, Y = np.meshgrid(inc.X, inc.Y)

        np.savetxt('W5/data/pressureDistribution02d.txt', np.array([X.flatten('F'), Y.flatten('F'), inc.p2.flatten('F')]).T)
        Pres_1d0 = np.tile(inc.p, len(inc.Y))
        np.savetxt('W5/data/pressureDistribution01d.txt', np.array([Y.flatten('F'), X.flatten('F'), Pres_1d0]).T)

        par = AS7(type = 1)
        par.geometry_parameters(plot = False)
        par.pressure_distrubution(plot = False)
        par.pressure_distrubution_2d(plot = False)

        X, Y = np.meshgrid(par.X, par.Y)

        np.savetxt('W5/data/pressureDistribution12d.txt', np.array([X.flatten('F'), Y.flatten('F'), par.p2.flatten('F')]).T)
        Pres_1d1 = np.tile(par.p, len(inc.Y))
        np.savetxt('W5/data/pressureDistribution11d.txt', np.array([Y.flatten('F'), X.flatten('F'), Pres_1d1]).T)


    if False: # Part c
        print('Part c')
        inc = AS7(type = 0)
        inc.geometry_parameters(plot = False)
        inc.film_thick(printbol = True)


    if True: # Part e
        print('Part e')
        shroud = AS7(type = 1, shroud = True)
        shroud.geometry_parameters(plot = False)
        shroud.pressure_distrubution_2d(plot = True, printbol = True)

        # Save height and pressure distribution
        X, Y = np.meshgrid(shroud.X, shroud.Y)
        np.savetxt('W5/data/pressureDistributionShroud2d.txt', np.array([X.flatten('F'), Y.flatten('F'), shroud.p2.flatten('F')]).T)
