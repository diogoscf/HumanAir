import numpy as np


class ISA:
    def __init__(self, Height, TempOffset, TemperatureGradient):
        self.Height = Height
        self.TempOffset = TempOffset
        self.TemperatureGradient = TemperatureGradient
        self.T_0 = 288.15  # K
        self.Rho_0 = 1.225  #
        self.g = 9.80665
        self.R = 287.05
        self.Mu_0 = 1.716e-5
        self.S_mu = 111  # Sutherland's constant air
        self.gamma = 1.4

    def Temperature(self):
        return self.T_0 + self.TemperatureGradient * self.Height + self.TempOffset

    def Density(self):
        Temp_no_offset = self.T_0 + self.TemperatureGradient * self.Height
        Correction = Temp_no_offset / (Temp_no_offset + self.TempOffset)
        return (
            Correction
            * self.Rho_0
            * (1 + self.TemperatureGradient * self.Height / self.T_0)
            ** (-(self.g / (self.R * self.TemperatureGradient) + 1))
        )

    def DynamicViscosity(self):
        return (
            self.Mu_0
            * (self.Temperature() / (self.T_0 - 15.15)) ** (3 / 2)
            * (self.T_0 - 15.15 + self.S_mu)
            / (self.Temperature() + self.S_mu)
        )

    def SpeedOfSound(self):
        return np.sqrt(self.gamma * self.R * self.Temperature())
