import numpy as np
import matplotlib.pyplot as plt
import json


class LongitudinalStability:
    def __init__(
        self, CLh, CLah, Xcgh, XLEMAC, CgAft, CgFwd, SM, deda, VhV, FuselageLength, Wing, airfoil_wing, airfoil_h
    ):
        self.CLh = CLh
        self.CLah = CLah
        self.FuselageLength = FuselageLength
        self.lh = Xcgh * FuselageLength - (XLEMAC + 0.4 * Wing.MAC())
        self.VhV = VhV
        self.c = Wing.MAC()
        self.AR = Wing.AR
        self.QuarterChordSweep = Wing.QuarterChordSweep
        self.Xac = 0.25
        self.dxcg = 0.01
        self.Xcg = np.arange(0, 2 + self.dxcg, self.dxcg)
        self.SM = SM
        self.deda = deda
        self.CgFwd = CgFwd
        self.CgAft = CgAft

        with open(airfoil_wing) as WingAirfoil:
            data = json.load(WingAirfoil)
            self.Clalphaah = data["C_L_Alpha"]
            self.Cm_0_Wing = data["Cm_0"]

        with open(airfoil_h) as HSAirfoil:
            data = json.load(HSAirfoil)
            self.Clalphah = data["C_L_Alpha"]
            self.Cm_0_HS = data["Cm_0"]
        # print("Longitudinal Stability Initialized")

    def CMac_Wing(self):
        return self.Cm_0_Wing * (
            self.AR
            * (np.cos(np.deg2rad(self.QuarterChordSweep))) ** 2
            / (self.AR + 2 * np.cos(np.deg2rad(self.QuarterChordSweep)))
        )

    def Stability(self):
        return (self.Xcg - self.Xac + self.SM) / (
            (self.Clalphah / 1.1) / (self.Clalphaah / 1.1) * (1 - self.deda) * self.lh / self.c * self.VhV**2
        )

    def Stability_NoMargin(self):
        return (self.Xcg - self.Xac) / (
            (self.Clalphah / 1.1) / (self.Clalphaah / 1.1) * (1 - self.deda) * self.lh / self.c * self.VhV**2
        )

    def Controllability(self):
        return (self.Xcg + self.CMac_Wing() / self.CLah - self.Xac) / (
            self.CLh / self.CLah * self.lh / self.c * self.VhV**2
        )

    def ShS(self):
        Cgmin = self.CgFwd
        Cgmax = self.CgAft

        ShSStability = np.interp(Cgmax, self.Xcg, self.Stability())
        ShSControllability = np.interp(Cgmin, self.Xcg, self.Controllability())

        return np.max([ShSStability, ShSControllability])

    def Plotting(self):  # pragma: no cover
        stability = self.Stability()
        stability_no_margin = self.Stability_NoMargin()
        controllability = self.Controllability()

        positive_stability_mask = stability > 0
        positive_stability_no_margin_mask = stability_no_margin > 0
        positive_controllability_mask = controllability > 0

        plt.figure()

        plt.plot(
            self.Xcg[positive_stability_mask],
            stability[positive_stability_mask],
            label="Stability",
            color="green",
            linestyle="solid",
        )
        plt.plot(
            self.Xcg[positive_stability_no_margin_mask],
            stability_no_margin[positive_stability_no_margin_mask],
            label="Stability (No Margin)",
            color="green",
            linestyle="dashed",
        )
        plt.plot(
            self.Xcg[positive_controllability_mask],
            controllability[positive_controllability_mask],
            label="Controllability",
            color="blue",
            linestyle="solid",
        )
        plt.plot(
            [self.CgFwd, self.CgAft], [self.ShS(), self.ShS()], label="Cg Excursion", color="red", linestyle="solid"
        )

        plt.legend()
        plt.xlabel("Xcg")
        plt.ylabel("Value")
        plt.title("Longitudinal Stability Analysis")
        plt.grid(True)
        plt.show()
