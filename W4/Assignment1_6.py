import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps


class AS6:
    def __init__(self, number_pads: int = 20, type: int = 0, notch_size: float = 0.1, outer_radius: float = 0.250/2, inner_radius: float = 0.095 , max_load: float = 4e5, kinematic_viscosity: float = 32e-6, minimum_height: float = 60e-6, density: float = 860, specific_heat: float = 2000, running_speed: float = 1000) -> None:
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
            self.x = np.linspace(0, self.l, 100)
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


            self.x = np.linspace(0, self.l, 1000)
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


    def friction(self):
        '''
        Calculates the friction coefficient
        '''
        if self.type == 0:
            up = 2 * self.s_h * np.log(self.H_0/(self.H_0+1)) + 3*self.s_h/(1+2*self.H_0)
            ned = 3 * self.l * np.log(self.H_0/(self.H_0+1)) + 6* self.l/(1+2*self.H_0)
            self.mu = up/ned # (8.34)
        
        if self.type == 1:
            self.P_m = 6 * self.n_s * (1 - self.n_s) / ( (1 - self.n_s) * (self.H_0 + 1)**3 + self.n_s * self.H_0**3) # (8.65
            self.mu = self.s_h / self.l * (1 + (2 * (self.H_0 + 1 - self.n_s)) / (self.P_m * self.H_0*(1+self.H_0))) # (8.80)


    def volume_flow(self):
        '''
        Calculates the dimenslionless volume flow rate per unite length
        '''
        if self.type == 0:
            self.Q = 2 * self.H_0 * (1 + self.H_0) / (1 + 2* self.H_0) # (8.36)
        
        if self.type == 1:
            self.Q = - self.P_m * (self.H_0 + 1)**3 / (6 * self.n_s) + self.H_0 + 1  # (8.81)


    def power_loss(self):
        '''
        Calculates the dimensionless power loss
        '''
        if self.type == 0:
            self.H_p = -4 * np.log(self.H_0 / (self.H_0 + 1)) - 6 / (1 + 2 * self.H_0) # (8.37)
        
        if self.type == 1:
            F_b = - self.P_m/2 - (self.H_0 + 1 - self.n_s) / (self.H_0 * (1 + self.H_0)) # (8.79)
            self.H_p = - F_b


    def temp_rise(self):
        '''
        Calculates the temperature rise
        '''
        self.Dt_m = 2 * self.u_b * self.l * self.eta_0 / (self.rho_0 * self.Cp * self.s_h**2) * self.H_p/self.Q # (8.13) - Adiabatic temperature rise


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
        if printbol:
            print(f"The load capacity is: {self.F:.3g} N")


        if plot:
            plt.figure()
            plt.plot(self.x, self.p)
            plt.xlabel('x [m]')
            plt.ylabel('p [Pa]')
            plt.title('Pressure distribution')
            plt.show()

            np.savetxt(f'W4/data/pressureDistribution{self.type}.txt', np.array([self.x, self.p]).T)



if __name__ == '__main__':
    if False: # Test
        print("Test")
        test = AS6()


    if False: # Part a
        print('Part a')
        inc = AS6(type = 0)
        inc.geometry_parameters(plot = True, print_bol = True)

        par = AS6(type = 1)
        par.geometry_parameters(plot = True, print_bol = True)


    if False: # Part b
        print("Part b:")
        inc = AS6(type = 0)
        inc.geometry_parameters()
        inc.friction()

        par = AS6(type = 1)
        par.geometry_parameters()
        par.friction()

        if inc.mu < par.mu:
            print(f"The friction coefficient for the fixed inclined pads is lower than the parallel-step: {inc.mu:.4g} < {par.mu:.4g}")
        else:
            print(f"The friction coefficient for the parallel-step is lower than the fixed inclined pads: {par.mu:.4g} < {inc.mu:.4g}")


    if False: # Part c
        print("Part c:")
        inc = AS6(type = 0)
        inc.geometry_parameters()
        inc.friction()
        inc.volume_flow()
        inc.power_loss()
        inc.temp_rise()

        par = AS6(type = 1)
        par.geometry_parameters()
        par.friction()
        par.volume_flow()
        par.power_loss()
        par.temp_rise()

        print(f"The temperature rise for the fixed inclined pads is: {inc.Dt_m} K")
        print(f"The temperature rise for the parallel-step is: {par.Dt_m} K")

        # Skal lige sammenlignes med Assignment1_3.py for at fine den viscosity ændring


    if False: # Part d
        print("Part d:")
        inc = AS6(type = 0)
        inc.geometry_parameters(plot=False)
        inc.pressure_distrubution(plot = True, printbol=True)
        print(inc.wm_z)

        par = AS6(type = 1)
        par.geometry_parameters(plot=False)
        par.pressure_distrubution(plot = True, printbol=True)
        print(par.wm_z)

    if True: # Part e
        print("Part e:")
        loads = np.linspace(1e5, 10e5, 100)
        film_thickness_0 = np.zeros(len(loads))
        for i, load in enumerate(loads):
            inc = AS6(max_load = load, type = 0)
            inc.geometry_parameters()
            film_thickness_0[i] = inc.h_0
        
        film_thickness_1 = np.zeros(len(loads))
        for i, load in enumerate(loads):
            par = AS6(max_load = load, type = 1)
            par.geometry_parameters()
            film_thickness_1[i] = par.h_0
        
        plt.figure()
        plt.plot(loads, film_thickness_0, label = 'Fixed inclined pads')
        plt.plot(loads, film_thickness_1, label = 'Parallel-step')
        plt.xlabel('Load [N]')
        plt.ylabel('Film thickness [m]')
        plt.title('Film thickness vs load')
        plt.legend()
        plt.show()

        np.savetxt('W4/data/filmThickness0.txt', np.array([loads, film_thickness_0]).T)
        np.savetxt('W4/data/filmThickness1.txt', np.array([loads, film_thickness_1]).T)

        inc_normal = AS6(max_load = 4e5, type = 0)
        inc_plus = AS6(max_load = 4e5*1.1, type = 0)
        inc_minus = AS6(max_load = 4e5*0.9, type = 0)

        inc_normal.geometry_parameters()
        inc_plus.geometry_parameters()
        inc_minus.geometry_parameters()

        print(f"For the incline pads, the film thickness is: {inc_normal.h_0:.3g} m")
        print(f"For a 10% increase in load, the film thickness changed by: {inc_plus.h_0-inc_normal.h_0:.3g} m")
        print(f"For a 10% decrease in load, the film thickness changed by: {inc_minus.h_0-inc_normal.h_0:.3g} m")

        par_normal = AS6(max_load = 4e5, type = 1)
        par_plus = AS6(max_load = 4e5*1.1, type = 1)
        par_minus = AS6(max_load = 4e5*0.9, type = 1)

        par_normal.geometry_parameters()
        par_plus.geometry_parameters()
        par_minus.geometry_parameters()

        print(f"For the parallel-step, the film thickness is: {par_normal.h_0:.3g} m")
        print(f"For a 10% increase in load, the film thickness changed by: {par_plus.h_0-par_normal.h_0:.3g} m")
        print(f"For a 10% decrease in load, the film thickness changed by: {par_minus.h_0-par_normal.h_0:.3g} m")