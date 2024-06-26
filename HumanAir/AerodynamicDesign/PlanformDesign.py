import numpy as np
import matplotlib.pyplot as plt


class Planform:
    def __init__(self, AR, Taper, QuarterChordSweep, tc, MTOW=None, WS=None, S=None):
        self.MTOW = MTOW
        self.WS = WS
        self.AR = AR
        self.Taper = Taper
        self.QuarterChordSweep = QuarterChordSweep
        self.g = 9.81
        self.S = S
        self.tc = tc
        if MTOW is None and WS is None and S is None:
            raise ValueError("Please provide a value for either MTOW and WS or S!")

    def WingSurfaceArea(self):
        if self.S is not None:
            return self.S
        else:
            return self.MTOW / self.WS

    def WingSpan(self):
        return np.sqrt(self.WingSurfaceArea() * self.AR)

    def RootChord(self):
        return 2 * self.WingSurfaceArea() / ((1 + self.Taper) * self.WingSpan())

    def TipChord(self):
        return self.Taper * self.RootChord()

    def MAC(self):
        return (2 * self.RootChord() / 3) * (1 + self.Taper + self.Taper**2) / (1 + self.Taper)

    def MAC_y(self):
        return (self.WingSpan() / 6) * (1 + 2 * self.Taper) / (1 + self.Taper)

    def HalfChordSweep(self):
        x_c4_root = self.RootChord() * 3 / 4
        x_c4_tip = x_c4_root - self.WingSpan() / 2 * np.tan(np.deg2rad(self.QuarterChordSweep))

        x_c2_root = self.RootChord() * 1 / 2
        x_c2_tip = x_c4_tip - self.TipChord() * 1 / 4
        print(x_c2_tip)
        print(x_c2_root)
        angle = np.rad2deg(np.arctan((np.abs(x_c2_root - x_c2_tip)) / (self.WingSpan() / 2)))
        return angle

    def t_root_max(self):
        return self.tc * self.RootChord()

    def PlotWingPlanform(self):  # pragma: no cover
        plt.figure()
        plt.plot([0, 0], [self.RootChord(), 0], color="black")  # Plotting Root Chord
        tip_x_qc = (
            self.RootChord() - self.RootChord() / 4 - np.tan(np.deg2rad(self.QuarterChordSweep)) * self.WingSpan() / 2
        )
        plt.plot(
            [self.WingSpan() / 2, self.WingSpan() / 2],
            [tip_x_qc + self.TipChord() / 4, tip_x_qc - 3 * self.TipChord() / 4],
            color="black",
        )  # Plotting Tip Chord
        plt.plot(
            [0, self.WingSpan() / 2], [self.RootChord(), tip_x_qc + self.TipChord() / 4], color="black"
        )  # Plot Leading Edge
        plt.plot([0, self.WingSpan() / 2], [0, tip_x_qc - 3 * self.TipChord() / 4], color="black")  # Plot Trailing Edge
        plt.plot([0, self.WingSpan() / 2], [self.RootChord() - self.RootChord() / 4, tip_x_qc])  # Plot Quarter Chord
        plt.axis("equal")
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    WingConventional = Planform(7.44, 0.4, 0, 0.1367, 25753, 618)
    print("S = ", WingConventional.WingSurfaceArea())
    print("b = ", WingConventional.WingSpan())
    print("Root Chord Length = ", WingConventional.RootChord())
    print("Tip Chord Length = ", WingConventional.TipChord())
    print("Max root thickness = ", WingConventional.t_root_max())

    WingConventional.PlotWingPlanform()
