import numpy as np

# import matplotlib.pyplot as plt
# import json
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data


class WeightEstm:
    def __init__(self, dict):
        self.dict = dict

    def OEW_prime(self):
        return (
            self.dict["Iterations Class I"]["A"] * self.dict["Iterations Class I"]["MTOW_kg"]
            + self.dict["Iterations Class I"]["B"]
        )

    def PowertrainWeight(self, bat):
        return (
            9.81
            * self.dict["Iterations Class I"]["MTOW_kg"]
            / 1000
            / self.dict["Performance"]["W/P_N/W"]
            / self.dict["Power_prop"]["eta_powertrain"]
            * (1 - bat)
            / self.dict["Power_prop"]["P_ptr_kW/kg"]
        )

    def BatteryWeight(self, bat):
        return (
            9.81
            * self.dict["Iterations Class I"]["MTOW_kg"]
            / self.dict["Performance"]["W/P_N/W"]
            * self.dict["Performance"]["endurance"]
            * bat
            / self.dict["Power_prop"]["E_bat_Wh/kg"]
            / self.dict["Power_prop"]["eta_bat"]
            / self.dict["Power_prop"]["DoD_bat"]
            / self.dict["Power_prop"]["eta_electricmotor"]
        )

    def FuelWeight(self, bat):
        return (
            1.15
            * 9.81
            * self.dict["Iterations Class I"]["MTOW_kg"]
            / self.dict["Performance"]["W/P_N/W"]
            * (1 - bat)
            * self.dict["Performance"]["endurance"]
            / self.dict["Power_prop"]["E_fuel_Wh/kg"]
            / self.dict["Power_prop"]["eta_generator"]
        )

    def WingWeight(self):
        return (
            self.dict["Iterations Class I"]["Aw"] * self.dict["Iterations Class I"]["MTOW_kg"]
            + self.dict["Iterations Class I"]["Bw"]
        )

    def Iterations(self, bat):
        MTOW_new = 0
        MTOW_old = self.dict["Iterations Class I"]["MTOW_kg"]
        ok = False

        while (
            np.abs((MTOW_new - self.dict["Iterations Class I"]["MTOW_kg"]) / self.dict["Iterations Class I"]["MTOW_kg"])
            > 0.02
        ):
            if ok:
                self.dict["Iterations Class I"]["MTOW_kg"] = MTOW_new

            OEW_prime = self.OEW_prime()
            PowertrainWeight = self.PowertrainWeight(bat)
            BatteryWeight = self.BatteryWeight(bat)
            FuelWeight = self.FuelWeight(bat)
            WingWeight = self.WingWeight()

            MTOW_new = (
                OEW_prime
                + PowertrainWeight
                + BatteryWeight
                + FuelWeight
                + WingWeight
                + self.dict["Iterations Class I"]["Wpl_des_kg"]
            )

            if MTOW_new > 8000:
                break

            ok = True

        if MTOW_new < 4000:
            self.dict["Iterations Class I"]["MTOW_kg"] = MTOW_old
            return (
                MTOW_new,
                self.dict["Contingency"] * MTOW_new,
                self.dict["Contingency"] * OEW_prime,
                self.dict["Contingency"] * PowertrainWeight,
                self.dict["Contingency"] * BatteryWeight,
                self.dict["Contingency"] * FuelWeight,
                self.dict["Contingency"] * WingWeight,
                self.dict["Contingency"] * self.dict["Iterations Class I"]["Wpl_des_kg"],
            )
        else:
            self.dict["Iterations Class I"]["MTOW_kg"] = MTOW_old
            return (
                0,
                self.dict["Contingency"] * MTOW_new,
                self.dict["Contingency"] * OEW_prime,
                self.dict["Contingency"] * PowertrainWeight,
                self.dict["Contingency"] * BatteryWeight,
                self.dict["Contingency"] * FuelWeight,
                self.dict["Contingency"] * WingWeight,
                self.dict["Contingency"] * self.dict["Iterations Class I"]["Wpl_des_kg"],
            )

    def PolynomialRegression(self, bat):
        lst_P = []
        lst_bat = []

        for pbat in bat:
            row = self.Iterations(pbat)
            if row[0] != 0:
                lst_P.append(9.81 * row[1] / self.dict["Performance"]["W/P_N/W"])
                lst_bat.append(pbat)

        lst_bat = np.array(lst_bat)
        lst_P = np.array(lst_P)

        # Filter out non-positive values in lst_P
        valid_indices = lst_P > 0
        lst_bat = lst_bat[valid_indices]
        lst_P = lst_P[valid_indices]

        if len(lst_P) == 0:
            return np.array([20, 20]), np.array([20, 20])

        coeff_exp = np.polyfit(lst_bat, np.log(lst_P), 1)
        coeff_pol = np.polyfit(lst_bat, lst_P, 2)

        # y_pol = coeff_pol[0] * lst_bat**2 + coeff_pol[1] * lst_bat + coeff_pol[2]
        # y_exp = np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * lst_bat)

        return coeff_exp, coeff_pol


if __name__ == "__main__":  # pragma: no cover
    data = WeightEstm(aircraft_data)

    bat = 0.11
    row = data.Iterations(bat)

    bat_lst = np.arange(0, 0.15, 0.001)
    coeff_exp, coeff_pol = data.PolynomialRegression(bat_lst)
    print(coeff_exp)
