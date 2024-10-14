import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps
from P3 import finite


class AS6:
    def __init__(self, number_pads: int = 20, type: int = 0, notch_size: float = 0.1, outer_radius: float = 0.250, inner_radius: float = 0.095 , max_load: float = 10e5, kinematic_viscosity: float = 32e-6, minimum_height: float = 60e-6, density: float = 860, specific_heat: float = 2000, running_speed: float = 1000) -> None:
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
        self.r_m = (self.r_o + self.r_i) / 2
        self.max_load = max_load
        self.nu = kinematic_viscosity
        self.h_0 = minimum_height
        self.rho_0 = density
        self.eta_0 = self.nu * self.rho_0
        self.Cp = specific_heat
        self.rpm = running_speed
        self.calc_length()
        self.calc_velocity()


    def calc_length(self) -> None:
        '''
        Calculates the length of the pad
        '''
        pad_angle = 360 / self.number_pads * (1 - self.notch_size) # The angle of the pad
        self.l = 2 * np.pi * self.r_m * pad_angle / 360 # The length of the pad


    def calc_velocity(self) -> None:
        '''
        Calculates the velocity of the bearing at the mid point of the pad
        '''
        # Er lidt i tvil om ub er korrekt
        self.u_b = 2 * np.pi * self.r_m * self.rpm / 60


    def geometry_parameters(self, plot: bool = False, print_bol: bool = False):
        '''
        Prints and calcualets the expressions the "optimal" geometry parameters of the bearing

        plot (bool): If True, plots the shoulder height
        '''
        if self.type == 0: # Fixed inclined pads
            if print_bol:
                print("Fixed inclined pads: h = h_0 + s_h * (1 - x/l) with the optimal being: s_h = sqrt(2) * h_0") # (8.16) (8.28)
            x = np.linspace(0, self.l, 100)
            self.s_h = np.exp(2) * self.h_0
            self.H_0 = self.h_0 / self.s_h
            self.P_m = 3 / (2* self.H_0 * (1 + self.H_0) * (1 + 2 * self.H_0)) # (8.26)
            h = self.h_0 + self.s_h * (1 - x/self.l)
            if plot:
                plt.figure()
                plt.plot(x, h)
                plt.xlabel('x [m]')
                plt.ylabel('h [m]')
                plt.ylim(0, np.max(h)*1.1)
                plt.title('Shoulder height - Fixed inclined pads')
                plt.show()

        if self.type == 1: # Parallel-step
            if print_bol:
                print("Parallel-step: H_0 = h_0 / s_h = 1.155 and n_s = 0.7182") # (8.68) (8.69)
            self.s_h = self.h_0 / 1.155
            self.H_0 = 1.155
            self.n_s = 0.7182
            self.P_m = 6*self.n_s / ((1-self.n_s) * (self.H_0 + 1)**3 + self.n_s * self.H_0**3) # (8.65)
            x = np.array([0, self.n_s*self.l, self.n_s*self.l, self.l])
            h = np.array([self.h_0+self.s_h, self.h_0+self.s_h, self.h_0, self.h_0])
            if plot:
                plt.figure()
                plt.plot(x, h)
                plt.xlabel('x [m]')
                plt.ylabel('h [m]')
                plt.ylim(0, np.max(h)*1.1)
                plt.title('Shoulder height - Parallel-step')
                plt.show()


    def friction(self):
        '''
        Calculates the friction coefficient
        '''
        if self.type == 0:
            up = 2 * self.s_h * np.log(self.H_0/(self.H_0+1)) + 3*self.s_h/(1+2*self.H_0)
            ned = 3 * self.l * np.log(self.H_0/(self.H_0+1)) + 6* self.l/(1+2*self.H_0)
            self.mu = up/ned # (8.34)
        
        if self.type == 1:
            self.mu = self.s_h / self.l * (1 + (2 * (self.H_0 + 1 - self.n_s)) / (self.P_m * self.H_0*(1+self.H_0))) # (8.80)


    def volume_flow(self):
        '''
        Calculates the dimenslionless volume flow rate per unite length
        '''
        if self.type == 0:
            self.Q = 2 * self.H_0 * (1 + self.H_0) / (1 + 2* self.H_0) # (8.36)
        
        if self.type == 1:
            self.Q = - self.P_m * (self.H_0 + 1)**3 / (6 * self.n_s) + self.H_0 + 1


    def power_loss(self):
        '''
        Calculates the dimensionless power loss
        '''
        if self.type == 0:
            self.H_p = -4 * np.log(self.H_0 / (self.H_0 + 1)) - 6 / (1 + 2 * self.H_0) # (8.37)
        
        if self.type == 1:
            F_b = - self.P_m/2 + (self.H_0 + 1 - self.n_s) / (self.H_0 * (1 + self.H_0)) # (8.79
            self.H_p = - F_b


    def temp_rise(self):
        '''
        Calculates the temperature rise
        '''
        # u_b = 10 # Fatter ikke lige hvorfor den er med?? Se side 189
        # Tror måske det bare er hastigheden
        self.Dt_m = 2 * self.u_b * self.l * self.eta_0 / (self.rho_0 * self.Cp * self.s_h**2) * self.H_p/self.Q # (8.13) - Adiabatic temperature rise


    def pressure_distrubution(self) -> None:
        '''
        Calculates the pressure distribution
        '''
        print("Pressure distribution")
        # if self.type == 0: # Fixed inclined pads




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
            print(f"The friction coefficient for the fixed inclined pads is lower than the parallel-step: {inc.mu:.4f} < {par.mu:.4f}")
        else:
            print(f"The friction coefficient for the parallel-step is lower than the fixed inclined pads: {par.mu:.4f} < {inc.mu:.4f}")


    if True: # Part c
        print("Part c:")
        inc = AS6(type = 0)
        inc.geometry_parameters()
        inc.volume_flow()
        inc.power_loss()
        inc.temp_rise()

        par = AS6(type = 1)
        par.geometry_parameters()
        par.volume_flow()
        par.power_loss()
        par.temp_rise()

        print(f"The temperature rise for the fixed inclined pads is: {inc.Dt_m} K")
        print(f"The temperature rise for the parallel-step is: {par.Dt_m} K")

        # Skal lige sammenlignes med Assignment1_3.py for at fine den viscosity ændring


    if False: # Part d
        print("Part d:")
        inc = AS6(type = 0)





        x_fin, p_fin = finite(100, inc.u_b, inc.h_0, inc.l, inc.eta_0)


        plt.figure()
        plt.plot(x_fin, p_fin, label='Finite difference')
        plt.xlabel('x [m]')
        plt.ylabel('p [Pa]')
        plt.title('Pressure distribution')
        plt.legend()
        plt.show()


