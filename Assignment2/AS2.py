class Bearing:
    def __init__(self, ShaftDiameter:float=101.6e-3, RadialClearance:float=102e-6, preLoadFactor:float=0, LengthDiameterFactor:float=0.5, Load:float=460/2, surfaceRoughness:float=0.5e-6) -> None:
        '''
        Creates a bearing object with the given parameters

        Args:
            shaftDiameter (float): The diameter of the shaft in meters. Variable name: D
            radialClearance (float): The radial clearance in meters. Variable name: Cp
            preLoadFactor (float): The preload factor. Variable name: mp
            lengthDiameterFactor (float): The length diameter factor
            load (float): The load in Kg
            surfaceRoughness (float): The surface roughness in meters. Variable name: Ra - See table 3.2
        '''
        self.D = ShaftDiameter
        self.Cp = RadialClearance
        self.mp = preLoadFactor
        self.L = LengthDiameterFactor * self.D
        self.W = Load * 9.82 # Convert to N
        self.Ra = surfaceRoughness


        

class Cylindrical(Bearing):
    pass



class TwoLobes(Bearing):
    pass


class Pad(Bearing):
    pass


if __name__ == "__main__":
    if True: # Test
        Lobes = TwoLobes()
        print(Lobes.W)