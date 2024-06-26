import numpy as np
import matplotlib.pyplot as plt


import sys
import os
import math

# Integration in progress
# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
sys.path.append(project_root)

from HumanAir import aircraft_data

from HumanAir.LoadingDiagram.Parameters import Parameters_CurrentDesign as p
import HumanAir.LoadingDiagram.Equations as eq
import HumanAir.LoadingDiagram.Plotting as plot


class WP_WS:
    # ========== 0: Initialise WS ==========
    def __init__(self):
        self.step = 0.1
        self.final_value_WS = 2000
        self.WS = np.arange(0.0001, self.final_value_WS, self.step)
        self.ylst = np.linspace(0, 1, int(self.final_value_WS / self.step))

        self.ReferenceWS = [
            449.0307667,
            602.1288757,
            643.1831172,
            725.2916003,
            725.2916003,
            578.0625962,
            789.2055899,
            780.3381114,
            930.8329016,
            925.2755192,
            1274.01538,
            1286.972761,
            1236.320571,
            987.5590745,
            1321.946577,
            1662.527136,
            1437.958113,
            1871.550807,
            1837.537433,
        ]

        self.ReferenceWP = [
            0.089508018,
            0.082049016,
            0.080130987,
            0.068752535,
            0.068752535,
            0.058831561,
            0.066557244,
            0.069219534,
            0.088135561,
            0.078032631,
            0.06465638,
            0.064313168,
            0.062743365,
            0.066302235,
            0.055426119,
            0.066832653,
            0.065302355,
            0.093120238,
            0.042435741,
        ]

        #       ========== 1: Calculate Values ==========
        self.Stallspeed_x = eq.Stallspeedx(
            0, 0, p.V_stall, p.Clmax_Land
        )  # Choose if we want to use Cl with clean config or landing config
        self.Takeoff_y = eq.Takeoff(p.TOP, self.WS, p.h_TO, p.temp_offset, p.Clmax_TO)
        self.Landing_x = eq.Landingx(p.Clmax_Land, p.h_Land, p.temp_offset, p.s_land, p.f)
        self.Cruise_y = eq.Cruise(
            p.eta_p, p.h_Cruise, p.temp_offset, p.Cdo, p.V_cruise, self.WS, p.A, p.e, p.CruisePower, p.CruiseWeight
        )
        self.Climbrate_y = eq.Climbrate(p.eta_p, p.A, p.e, 0, 0, p.Cdo, p.climbrate, self.WS)
        self.Climbgradient_y = eq.Climbgradient(
            p.eta_p, self.WS, p.climbgradient, p.V_climb, p.A, p.e, p.Clmax_TO, p.Cl_SafetyFactor, p.Cdo, 0, 0
        )

        # Convert the x-values to y-values
        self.Landing_y = self.Landing_x * self.ylst
        self.Stallspeed_y = self.Stallspeed_x * self.ylst

    def calculate_optimal_point(self):
        # Create a matrix where each column is one of the y-values
        y_matrix = np.vstack(
            (self.Stallspeed_y, self.Landing_y, self.Takeoff_y, self.Cruise_y, self.Climbrate_y, self.Climbgradient_y)
        )

        # Find the minimum values across all curves for each WS
        min_envelope = np.min(y_matrix, axis=0)

        # Find the optimal W/P and W/S
        optimal_WP = np.max(min_envelope)
        if optimal_WP is None:
            optimal_WP = 0.0001
        optimal_WS = self.WS[0]

        if all(math.isnan(value) for value in min_envelope):
            idx = 1
        else:
            idx = self.WS[np.where(min_envelope == optimal_WP)[0][0]]

        # change this value to tune it but just god knows which is the optimal one
        WP_tolerance = 1e-4  # Define the tolerance for W/P change
        WS_tolerance = 0.05 * idx  # Define the tolerance for W/S change

        for y in range(len(min_envelope) - 1):
            if (
                np.abs(min_envelope[y] - optimal_WP) < WP_tolerance
            ):  # check if the numbers are increasing and update the optimal values
                if np.abs(self.WS[y] - optimal_WS) > WS_tolerance:
                    optimal_WP = min_envelope[y]
                    optimal_WS = self.WS[y]

            elif (
                min_envelope[y] - optimal_WP > WP_tolerance
            ):  # check if the values for W/P are decreasing and update the optimal values
                optimal_WP = min_envelope[y]
                optimal_WS = self.WS[y]

        return optimal_WP, optimal_WS

    #  ========== 2: Plot Lines =========="""
    def plot(self, saving=None):
        # not sure if altitude correction should be applied to power,
        # assumed here and in FlightPerformance.takeoff_landing.py that it is not used
        # power_alt_correction = 1.132 * (eq.Density(p.h_TO, p.temp_offset)/1.225) - 0.132
        # Gagg and Ferrar model, Gudmundsen eq 7-16
        acf = aircraft_data.aircraft_data
        W = acf["CL2Weight"]["MTOW_N"]
        S_clean = acf["Aero"]["S_Wing"]
        S_to = acf["Aero"]["S_Wing"] * acf["Flaps"]["Sprime_S_takeoff"]
        S_ld = acf["Aero"]["S_Wing"] * acf["Flaps"]["Sprime_S_landing"]
        P_to = (
            acf["Power_prop"]["Fuel_engine_P_TO_W"] * acf["Power_prop"]["eta_powertrain"]
        )  # * power_alt_correction # prop eff should not be included
        P_cont = (
            acf["Power_prop"]["Fuel_engine_P_max_cont_W"] * acf["Power_prop"]["eta_powertrain"]
        )  # prop eff taken into account in formulas

        our_design_to = [W / S_to, W / P_to]
        our_design_ld = [W / S_ld, W / P_cont]
        our_design_clean = [W / S_clean, W / P_cont]

        optimal_point = self.calculate_optimal_point()

        plt.figure(figsize=(10, 7))

        plot.Plotx(self.Stallspeed_x, self.ylst, "Stall Speed")
        plot.Plotx(self.Landing_x, self.ylst, "Landing Distance")
        plot.Ploty(self.WS, self.Takeoff_y, "Take-off Distance")
        plot.Ploty(self.WS, self.Cruise_y, "Cruise")
        plot.Ploty(self.WS, self.Climbrate_y, "Climb Rate")
        plot.Ploty(self.WS, self.Climbgradient_y, "Climb Gradient")

        plt.fill_between(self.WS, self.Takeoff_y, 1, color="red", alpha=0.1)
        plt.fill_between(self.WS, self.Cruise_y, 1, color="red", alpha=0.1)
        plt.fill_between(self.WS, self.Climbrate_y, 1, color="red", alpha=0.1)
        plt.fill_between(self.WS, self.Climbgradient_y, 1, color="red", alpha=0.1)
        plt.fill_betweenx(self.ylst, self.Stallspeed_x, 2000, color="red", alpha=0.1)
        plt.fill_betweenx(self.ylst, self.Landing_x, 2000, color="red", alpha=0.1)
        plt.scatter(self.ReferenceWS, self.ReferenceWP, color="orange", alpha=0.5, label="Reference Aircraft")
        # plt.scatter(WS[index],Cruise_y[index], label="Chosen Design Point", s=100, color='red')
        plt.scatter(1040.95, 0.076651773, label="Cessna 206", s=100, color="blue")
        plt.scatter(optimal_point[1], optimal_point[0], label="Design Point", s=100, color="green")
        plt.scatter(our_design_to[0], our_design_to[1], label="Current Design, Takeoff", s=100, color="orange")
        plt.scatter(our_design_ld[0], our_design_ld[1], label="Current Design, Landing", s=100, color="black")
        plt.scatter(our_design_clean[0], our_design_clean[1], label="Current Design, clean", s=100, color="brown")
        plt.xlabel(r"W/S (N/$m^2$)")
        plt.ylabel("W/P (N/W)")
        plt.ylim(0, 0.25)
        plt.xlim(0, self.final_value_WS)
        plt.subplots_adjust(right=0.75)
        plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
        plt.title(p.name)

        if saving:
            plt.savefig("Performance_FlyingWing.svg")

        plt.show()


if __name__ == "__main__":  # pragma: no cover
    wp = WP_WS()
    print(wp.calculate_optimal_point())
    wp.plot()
