import numpy as np


class Flow:
    def __init__(self, V_Cruise, ISA, Wing):
        self.V = V_Cruise
        self.a = ISA.SpeedOfSound()
        self.Mu = ISA.DynamicViscosity()
        self.MAC = Wing.MAC()
        self.Rho = ISA.Density()

    def Mach(self):
        return self.V / self.a

    def Reynolds(self):
        return self.V * self.MAC * self.Rho / self.Mu

    def Beta(self):
        return np.sqrt(1 - self.Mach() ** 2)
