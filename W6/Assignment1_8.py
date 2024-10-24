import numpy as np
import matplotlib.pyplot as plt


class AS8:
    def __init__(self, outer_radius: float = 0.125, film_thickness: float = 10e-6, kinematic_viscosity: float = 32e-6, density: float = 860, running_speed: float = 1000, load: float = 4e5) -> None:
        '''
        outer_radius (float): The outer radius of the bearing [m]
        
        film_thickness (float): The film thickness of the bearing [m]

        kinematic_viscosity (float): The kinematic viscosity of the oil at p = 0 [m^2/s]

        density (float): The density of the oil [kg/m^3]

        running_speed (float): The running speed of the bearing [rpm]

        load (float): The load of the bearing [N]
        '''
        self.r_o = outer_radius
        self.h_0 = film_thickness
        self.nu = kinematic_viscosity
        self.rho_0 = density
        self.eta_0 = self.nu * self.rho_0
        self.rpm = running_speed
        self.omega = self.rpm * 2 * np.pi / 60
        self.w_z = load

    def calc_geometry(self, print_bol: bool = False) -> None:
        '''
        Calculates the 'optimal' geometry of the bearing

        print_bi (bool): If True, the function will print the results
        '''
        # self.r_i = self.r_o * 0.53 # p.327


        self.r_i = self.r_o * 0.53


        self.Delta = 10 * self.h_0 # p.323


        if print_bol:
            print("The inner radius should be: ", self.r_i)
            print("The notch depth should be at least: ", self.Delta)


    def vicous_power_dissipation(self) -> None:
        '''
        Calculates the viscous power dissipation of the bearing (H_v)
        '''
        self.H_v = np.pi * self.eta_0 * self.omega**2 / (2 * self.h_0) * (self.r_o**4 - self.r_i**4) # (13.8)


    def pumping_power(self) -> None:
        '''
        Calculates the pumping power of the bearing (H_p)
        '''
        self.p_r = 2 * self.w_z * np.log(self.r_o/self.r_i) / (np.pi * (self.r_o**2 - self.r_i**2)) # (13.6)

        self.H_p = np.pi * self.h_0**3 * self.p_r**2 /(6 * self.eta_0 * np.log(self.r_o/self.r_i)) # (13.9)


    def total_power(self, print_bool: bool = False) -> None:
        '''
        Calculates the total power of the bearing (H)
        
        print_bool (bool): If True, the function will print the results
        '''

        self.calc_geometry()
        self.vicous_power_dissipation()
        self.pumping_power()
        self.H = self.H_v + self.H_p # (13.10)

        if print_bool:
            print("The viscous power dissipation is: ", self.H_v, "[W]")
            print("The pumping power is: ", self.H_p, "[W]")
            print("The total power is: ", self.H, "[W]")

if __name__ == '__main__':
    if False: # Part a
        print("Part a")
        bearing = AS8()
        bearing.calc_geometry(print_bol = True)
    
    if True: # Part b
        print("Part b")
        bearing = AS8()
        bearing.total_power(print_bool=True)
        # Shit bror. Det er jo sindsyg vicous friction fordi film thickness er så lille

        bearing_meretyk = AS8(film_thickness=50e-6)
        bearing_meretyk.total_power(print_bool=True)
        # Det var jo bedre

        bearing_totalfed = AS8(film_thickness=10e-5)
        bearing_totalfed.total_power(print_bool=True)

        # Vi kan lige lave et plot

        film_tykkelse = np.linspace(1e-6, 150e-6, 100)
        
        bearings = np.zeros(len(film_tykkelse), dtype=object)

        totalPower = np.zeros(len(film_tykkelse))
        viscPower = np.zeros(len(film_tykkelse))
        pumpPower = np.zeros(len(film_tykkelse))

        for i, height in enumerate(film_tykkelse):
            bearing = AS8(film_thickness=height)
            bearing.total_power()
            totalPower[i] = bearing.H
            viscPower[i] = bearing.H_v
            pumpPower[i] = bearing.H_p

        # Convering units
        film_tykkelse = film_tykkelse * 1e6
        totalPower = totalPower / 1e3
        viscPower = viscPower / 1e3
        pumpPower = pumpPower / 1e3

        plt.figure()
        plt.plot(film_tykkelse, totalPower, 'r', label="Total")
        plt.plot(film_tykkelse, viscPower, 'b', label="Viscous")
        plt.plot(film_tykkelse, pumpPower, 'g', label="Pumping")
        plt.axis([0, 150, 0, 20])
        plt.legend()
        plt.xlabel('Film thickness [$\mu $m]')
        plt.ylabel('Power [kW]')
        plt.grid(True)
        plt.show()

        np.savetxt(f'W6/data/totalPower.txt', np.array([film_tykkelse, totalPower]).T)
        np.savetxt(f'W6/data/viscPower.txt', np.array([film_tykkelse, viscPower]).T)
        np.savetxt(f'W6/data/pumpPower.txt', np.array([film_tykkelse, pumpPower]).T)

        print(min(totalPower))
        # print film thickness for minimum power
        print(film_tykkelse[np.argmin(totalPower)])
