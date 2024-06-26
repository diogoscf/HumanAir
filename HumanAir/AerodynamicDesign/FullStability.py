# import pandas as pd
import numpy as np
import sys
from math import tan, sqrt, pi, atan

# import time
import os

# import json
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data
from HumanAir.Weights_and_CG.weight_fractions_C2 import (
    calculate_lh,
    optimised_xlemac_landing_gears,
    # cg_excursion,
    # component_mass,
)
from HumanAir.isa import isa


def Geometry(acd=aircraft_data):
    taper = acd["Aero"]["Taper_Wing"]
    # Root chord
    Cr = 3 / 2 * acd["Aero"]["MAC_wing"] * ((1 + taper) / (1 + taper + taper**2))

    # Net wing area, outside of fuselage
    Cfus = Cr * (1 - (1 - taper) * acd["Geometry"]["fus_width_m"] / acd["Aero"]["b_Wing"])
    Snet = acd["Aero"]["S_Wing"] - (Cr + Cfus) * acd["Geometry"]["fus_width_m"] / 2

    # Length to net root chord, outside of fuselage
    Ymgc = acd["Aero"]["b_Wing"] / (2 * (1 - taper)) * (1 - 2 / 3 * ((1 + taper + taper**2) / (1 + taper)))
    tanLE = -4 / acd["Aero"]["AR"] * (-0.25 * (1 - taper) / (1 + taper))
    l_fn = acd["Geometry"]["XLEMAC_m"] - tanLE * (Ymgc - acd["Geometry"]["fus_width_m"] / 2)

    # Stabiliser affected by propwash
    taper_h = acd["Aero"]["Taper_HS"]
    Cr_h = 3 / 2 * acd["Aero"]["MAC_HS"] * ((1 + taper_h) / (1 + taper_h * taper_h**2))
    Cpw = Cr_h * (1 - (1 - taper_h) * acd["Power_prop"]["Dp_m"] / acd["Aero"]["b_h"])
    Shslip = (Cr_h + Cpw) * acd["Power_prop"]["Dp_m"] / 2

    # Distance between quarter chord MAC of wing and stabiliser

    return Cr, Snet, l_fn, Shslip


def TailAero(l_H, acd=aircraft_data):
    # Get and convert values to imp*rial
    Shslip = Geometry(acd)[3] * 10.7639104
    Sh = acd["Aero"]["S_h"] * 10.7639104
    U1 = acd["Performance"]["Vc_m/s"] * 3.2808399
    Pav = 0.98 * 400 * acd["Power_prop"]["eta_p"]
    Dp = acd["Power_prop"]["Dp_m"] / 0.3048
    qbar = (
        0.5
        * isa(acd["Performance"]["Altitude_Cruise_m"], acd["Performance"]["Temp_offset_TO_Land_cruise"])[2]
        * acd["Performance"]["Vc_m/s"] ** 2
    ) * 0.020885

    # Magic Roskam equation, only god knows how this works
    eta_H = 1 + (Shslip / Sh) * (2200 * Pav / (qbar * U1 * pi * Dp**2))
    VhVcorr = acd["Aero"]["VhV"] * eta_H

    # Downwash Gradient supposedly from Slingerland whoever that may be
    r = float(l_H * 2 / acd["Aero"]["b_Wing"])
    mtv = float(0 * 2 / acd["Aero"]["b_Wing"])
    deda = (
        (
            (r / (r**2 + mtv**2)) * 0.4876 / np.sqrt(r**2 + 0.6319 + mtv**2)
            + (1 + (r**2 / (r**2 + 0.7915 + 5.0734 * mtv**2)) ** 0.3113) * (1 - np.sqrt(mtv**2 / (1 + mtv**2)))
        )
        * acd["Aero"]["CLalpha"]
        / (pi * acd["Aero"]["AR"])
    )

    return VhVcorr, deda


def Liftrate(l_H, acd=aircraft_data):
    # Tail lift rate
    SweepHS_05 = np.tan(np.radians(acd["Aero"]["QuarterChordSweep_HS_deg"])) - 4 / acd["Aero"]["AR_HS"] * (
        0.25 * (1 - acd["Aero"]["Taper_HS"]) / (1 + acd["Aero"]["Taper_HS"])
    )
    ClaH = (
        2 * pi * acd["Aero"]["AR_HS"] / (2 + sqrt(4 + (acd["Aero"]["AR_HS"] / 0.95) ** 2 * (1 + tan(SweepHS_05) ** 2)))
    )

    # Aircraft less tail lift rate
    Snet = Geometry(acd)[1]
    CLaAH = (
        acd["Aero"]["CLalpha"]
        * (1 + 2.15 * acd["Geometry"]["fus_width_m"] / acd["Aero"]["b_Wing"])
        * Snet
        / acd["Aero"]["S_Wing"]
        + pi / 2 * acd["Geometry"]["fus_width_m"] ** 2 / acd["Aero"]["S_Wing"]
    )

    # Controllability tail lift
    CLH = -0.35 * acd["Aero"]["AR_HS"] ** (1 / 3)

    # Controllability aircraft-less-tail lift
    CLAH = acd["Aero"]["CLmax_Land"] - CLH * TailAero(l_H, acd)[0] ** 2 * acd["Aero"]["S_h"] / acd["Aero"]["S_Wing"]

    return ClaH, CLaAH, CLH, CLAH


def Xacplane(l_H, acd=aircraft_data):
    # From ADSEE data for no wing sweep
    Xac_wing = 0.205 + 0.005 * acd["Aero"]["AR"]

    # Fuselage contribution
    CLaAH = Liftrate(l_H, acd)[1]
    l_fn = Geometry(acd)[2]
    Xac_fus = (
        -(1.8 / CLaAH)
        * (acd["Geometry"]["fus_width_m"] * acd["Geometry"]["fus_height_m"] * l_fn)
        / (acd["Aero"]["S_Wing"] * acd["Aero"]["MAC_wing"])
    )

    # Total Xac
    Xac = Xac_fus + Xac_wing

    return Xac


def MomentCoefficient(l_H, acd=aircraft_data):
    # Wing contribution

    Cm_0_Wing = acd["Aero"]["Cm_0_wing"]

    Cmacw = Cm_0_Wing * (
        acd["Aero"]["AR"]
        * np.cos(acd["Aero"]["QuarterChordSweep_Wing_deg"]) ** 2
        / (acd["Aero"]["AR"] + 2 * np.cos(acd["Aero"]["QuarterChordSweep_Wing_deg"]))
    )

    # Flap Contribution
    ClaH, CLaAH, CLH, CLAH = Liftrate(l_H, acd)
    Xac = Xacplane(l_H, acd)
    Cmacflap = acd["Flaps"]["mu2_land"] * (
        -acd["Flaps"]["mu1_land"] * acd["Flaps"]["deltaCLmax_land"] * acd["Flaps"]["cprime_c_landing"]
        - (CLAH + acd["Flaps"]["deltaCLmax_land"] * (1 - acd["Flaps"]["Swf"] / acd["Aero"]["S_Wing"]))
        * 1
        / 8
        * acd["Flaps"]["cprime_c_landing"]
        * (acd["Flaps"]["cprime_c_landing"] - 1)
    ) - acd["Aero"]["CLmax_Land"] * (0.25 - Xac / acd["Aero"]["MAC_wing"])

    # Fuselage contribution
    Cmacfus = (
        -1.8
        * (1 - 2.5 * acd["Geometry"]["fus_width_m"] / acd["Geometry"]["fus_length_m"])
        * (pi * acd["Geometry"]["fus_width_m"] * acd["Geometry"]["fus_height_m"] * acd["Geometry"]["fus_length_m"])
        / (4 * acd["Aero"]["S_Wing"] * acd["Aero"]["MAC_wing"])
        * acd["Flaps"]["CL0"]
        / CLaAH
    )

    return Cmacw + Cmacflap + Cmacfus


def StabControl(acd=aircraft_data):
    # Initiate xcg/MAC
    Xcg = np.arange(-0.2, 1.20, 0.001)

    # l_H iteration/calculation, XLEMAC placement
    # TODO: add this
    # TODO: done
    lh_sv = calculate_lh(ac_data=aircraft_data)
    l_H = acd["Stability"]["QCW_to_QCh"]

    # Get necessary values from functions
    Xac = Xacplane(l_H, acd)
    ClaH, CLaAH, CLH, CLAH = Liftrate(l_H, acd)
    VhVcorr, deda = TailAero(l_H, acd)
    Cmac = MomentCoefficient(l_H, acd)

    # Stability Sh/S
    StabSM = (Xcg - Xac + 0.05) / ((ClaH / CLaAH) * (1 - deda) * l_H / acd["Aero"]["MAC_wing"] * VhVcorr**2)
    StabNeutral = (Xcg - Xac) / ((ClaH / CLaAH) * (1 - deda) * l_H / acd["Aero"]["MAC_wing"] * VhVcorr**2)

    # Controllability Sh/S
    Control = (Xcg - Xac + Cmac / CLAH) / (CLH / CLAH * l_H / acd["Aero"]["MAC_wing"] * VhVcorr**2)

    return StabSM, StabNeutral, Control, Xcg, lh_sv


# negative values does not work because the nose gear will colapse, the begin value has to be larger than 0.05
def TailIteration(ac_datafile=aircraft_data, begin_value=0.2, end_value=0.6, step=10):
    # define the minimum value for Sh/S
    Sh_S_min = 10000

    for x_percentage in range(int(begin_value * 100), int(end_value * 100), step):
        # iterate such that the batteries and landing gear fits
        # print(x_percentage/100)

        # get the xlemac such that the batteries and landing gear fits
        optimised_xlemac_landing_gears(ac_data=ac_datafile, percentage=x_percentage / 100, bat_xcg_init=0.2)

        lh_converged = False
        aft_cg = round(ac_datafile["Stability"]["Cg_Aft"], 2)
        front_cg = round(ac_datafile["Stability"]["Cg_Front"], 2)

        Xcg_excursion = np.arange(front_cg, aft_cg, 0.001)

        # initialise the surfaces areas
        S_h_old = 10000
        S = ac_datafile["Aero"]["S_Wing"]

        # loop to iterate for lh and find the optimal value with Sh
        while not lh_converged:
            # get the stability and control lines with the updated values
            StabSM, _, Control, Xcg, lh_sv = StabControl(acd=ac_datafile)

            # round the values to 4 decimals for the Xcg range
            Xcg = np.array([round(i, 4) for i in Xcg])

            # get the index for the control and stability lines such that the xcg fits
            Sh_S_control_idx = np.max(np.where(np.isclose(front_cg, Xcg, atol=0.0001)))
            Sh_S_stab_idx = np.min(np.where(np.isclose(aft_cg, Xcg, atol=0.0001)))

            # update the horizontal tail surface area and Sh_S depending on the limiting case
            if StabSM[Sh_S_stab_idx] <= Control[Sh_S_control_idx]:
                Sh_S = Control[Sh_S_control_idx]
                S_h = Sh_S * S

                # update the dictionary such that the lh gets updated
                ac_datafile["Aero"]["S_h"] = S_h
                ac_datafile["Aero"]["b_h"] = sqrt(ac_datafile["Aero"]["AR_HS"] * ac_datafile["Aero"]["S_h"])
                ac_datafile["Aero"]["c_root_HS"] = (
                    2
                    * ac_datafile["Aero"]["S_h"]
                    / (ac_datafile["Aero"]["b_h"] * (1 + ac_datafile["Aero"]["Taper_HS"]))
                )
                ac_datafile["Aero"]["c_tip_HS"] = ac_datafile["Aero"]["Taper_HS"] * ac_datafile["Aero"]["c_root_HS"]
                ac_datafile["Aero"]["MAC_HS"] = (
                    (2 / 3 * ac_datafile["Aero"]["c_root_HS"])
                    * (1 + ac_datafile["Aero"]["Taper_HS"] + ac_datafile["Aero"]["Taper_HS"] ** 2)
                    / (1 + ac_datafile["Aero"]["Taper_HS"])
                )

                # save the best configuration
                if Sh_S < Sh_S_min:
                    save_value_xcg_wing = x_percentage
                    Sh_S_min = Sh_S
                    Xcg_excursion_min = Xcg_excursion

                    dummy_dict = {}
                    dummy_dict["Aero"] = {}
                    dummy_dict["Geometry"] = {}
                    dummy_dict["Geometry"]["XLEMAC_m"] = ac_datafile["Geometry"]["XLEMAC_m"]
                    dummy_dict["Stability"] = {}
                    dummy_dict["Stability"]["QCW_to_QCh"] = float(lh_sv)
                    dummy_dict["Aero"]["S_h"] = S_h
                    dummy_dict["Aero"]["b_h"] = sqrt(ac_datafile["Aero"]["AR_HS"] * ac_datafile["Aero"]["S_h"])
                    dummy_dict["Aero"]["c_root_HS"] = (
                        2
                        * ac_datafile["Aero"]["S_h"]
                        / (ac_datafile["Aero"]["b_h"] * (1 + ac_datafile["Aero"]["Taper_HS"]))
                    )
                    dummy_dict["Aero"]["c_tip_HS"] = ac_datafile["Aero"]["Taper_HS"] * ac_datafile["Aero"]["c_root_HS"]
                    dummy_dict["Aero"]["MAC_HS"] = (
                        (2 / 3 * ac_datafile["Aero"]["c_root_HS"])
                        * (1 + ac_datafile["Aero"]["Taper_HS"] + ac_datafile["Aero"]["Taper_HS"] ** 2)
                        / (1 + ac_datafile["Aero"]["Taper_HS"])
                    )

            else:
                Sh_S = StabSM[Sh_S_stab_idx]
                S_h = Sh_S * S

                # update the dictionary such that the lh gets updated
                ac_datafile["Aero"]["S_h"] = S_h
                ac_datafile["Aero"]["b_h"] = sqrt(ac_datafile["Aero"]["AR_HS"] * ac_datafile["Aero"]["S_h"])
                ac_datafile["Aero"]["c_root_HS"] = (
                    2
                    * ac_datafile["Aero"]["S_h"]
                    / (ac_datafile["Aero"]["b_h"] * (1 + ac_datafile["Aero"]["Taper_HS"]))
                )
                ac_datafile["Aero"]["c_tip_HS"] = ac_datafile["Aero"]["Taper_HS"] * ac_datafile["Aero"]["c_root_HS"]
                ac_datafile["Aero"]["MAC_HS"] = (
                    (2 / 3 * ac_datafile["Aero"]["c_root_HS"])
                    * (1 + ac_datafile["Aero"]["Taper_HS"] + ac_datafile["Aero"]["Taper_HS"] ** 2)
                    / (1 + ac_datafile["Aero"]["Taper_HS"])
                )

                # save the best configuration
                if Sh_S < Sh_S_min:
                    save_value_xcg_wing = x_percentage
                    Sh_S_min = Sh_S
                    Xcg_excursion_min = Xcg_excursion

                    dummy_dict = {}
                    dummy_dict["Aero"] = {}
                    dummy_dict["Stability"] = {}
                    dummy_dict["Geometry"] = {}
                    dummy_dict["Geometry"]["XLEMAC_m"] = ac_datafile["Geometry"]["XLEMAC_m"]
                    dummy_dict["Stability"]["QCW_to_QCh"] = float(lh_sv)
                    dummy_dict["Aero"]["S_h"] = S_h
                    dummy_dict["Aero"]["b_h"] = sqrt(ac_datafile["Aero"]["AR_HS"] * ac_datafile["Aero"]["S_h"])
                    dummy_dict["Aero"]["c_root_HS"] = (
                        2
                        * ac_datafile["Aero"]["S_h"]
                        / (ac_datafile["Aero"]["b_h"] * (1 + ac_datafile["Aero"]["Taper_HS"]))
                    )
                    dummy_dict["Aero"]["c_tip_HS"] = ac_datafile["Aero"]["Taper_HS"] * ac_datafile["Aero"]["c_root_HS"]
                    dummy_dict["Aero"]["MAC_HS"] = (
                        (2 / 3 * ac_datafile["Aero"]["c_root_HS"])
                        * (1 + ac_datafile["Aero"]["Taper_HS"] + ac_datafile["Aero"]["Taper_HS"] ** 2)
                        / (1 + ac_datafile["Aero"]["Taper_HS"])
                    )

            # check for convergence
            if np.abs(S_h_old - S_h) / S_h_old < 0.0001:
                lh_converged = True

            # save the old value of Sh for the lh iteration
            S_h_old = S_h

    # save the best configuration
    ac_datafile["Geometry"]["XLEMAC_m"] = dummy_dict["Geometry"]["XLEMAC_m"]
    ac_datafile["Aero"]["S_h"] = dummy_dict["Aero"]["S_h"]
    ac_datafile["Aero"]["b_h"] = dummy_dict["Aero"]["b_h"]
    ac_datafile["Aero"]["c_root_HS"] = dummy_dict["Aero"]["c_root_HS"]
    ac_datafile["Aero"]["c_tip_HS"] = dummy_dict["Aero"]["c_tip_HS"]
    ac_datafile["Aero"]["MAC_HS"] = dummy_dict["Aero"]["MAC_HS"]
    ac_datafile["Stability"]["QCW_to_QCh"] = float(dummy_dict["Stability"]["QCW_to_QCh"])
    ac_datafile["Stability"]["Xcg_oew_wing_mac"] = save_value_xcg_wing

    tan_LE_sweep = tan(0) - 4 / ac_datafile["Aero"]["AR_HS"] * (
        (-3 / 4) * (1 - ac_datafile["Aero"]["Taper_HS"]) / (1 + ac_datafile["Aero"]["Taper_HS"])
    )
    quarter_chord_sweep = atan(
        tan_LE_sweep
        + ac_datafile["Aero"]["c_root_HS"] / (2 * ac_datafile["Aero"]["b_h"]) * (ac_datafile["Aero"]["Taper_HS"] - 1)
    )
    ac_datafile["Aero"]["QuarterChordSweep_HS_deg"] = quarter_chord_sweep * 180 / np.pi

    # plot the optimised xcg for the landing gear
    Plotting(acd=ac_datafile, show=False)

    Sh_S_list = np.ones(np.shape(Xcg_excursion_min)[0]) * Sh_S_min
    plt.plot(Xcg_excursion_min, Sh_S_list, label="Optimised Xcg for Landing Gear")
    plt.savefig("ScissorPlot.svg")
    plt.show()

    print(f"The aircraft has a horizontal tail with a surface area of {round(ac_datafile['Aero']['S_h'], 2)} [m^2]")


def Plotting(acd=aircraft_data, show=True):  # pragma: no cover
    # Get data to plot from previous functions
    StabSM, StabNeutral, Control, Xcg, lh_sv = StabControl(acd)
    StabSM = np.array(StabSM, dtype=np.float64)
    StabNeutral = np.array(StabNeutral, dtype=np.float64)
    Control = np.array(Control, dtype=np.float64)
    Xcg = np.array(Xcg, dtype=np.float64)
    # TailSizing(ac_datafile=aircraft_data)
    # Actual plotting
    plt.plot(Xcg, StabSM, label="Stability with Safety Margin", color="limegreen", linewidth=2.2)
    plt.plot(Xcg, StabNeutral, label="Neutral Stability", linestyle="--", color="red")
    plt.plot(Xcg, Control, label="Controllability", color="dodgerblue", linewidth=2.2)
    plt.fill_between(Xcg, 0, StabSM, color="crimson", alpha=0.2)
    plt.fill_between(Xcg, 0, Control, color="crimson", alpha=0.2)
    plt.xlim(-0.2, 1.2)
    plt.ylim(0, 0.6)
    plt.legend()

    if show:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    # Plotting()
    TailIteration(ac_datafile=aircraft_data)
