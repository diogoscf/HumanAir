# import pandas as pd
import numpy as np
import sys

# from math import tan, sqrt, pi
# import time
import os

# import json
# import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data

from HumanAir.isa import isa


def chord(data=aircraft_data, location=1.9):
    result = (
        2
        * data["Aero"]["S_Wing"]
        / (1 + data["Aero"]["Taper_Wing"])
        / data["Aero"]["b_Wing"]
        * (1 - (1 - data["Aero"]["Taper_Wing"]) / data["Aero"]["b_Wing"] * 2 * location)
    )
    return result


def cy_y(acd=aircraft_data, y=1):
    return (
        2
        * acd["Aero"]["S_Wing"]
        / (1 + acd["Aero"]["Taper_Wing"])
        / acd["Aero"]["b_Wing"]
        * (y**2 / 2 - (1 - acd["Aero"]["Taper_Wing"]) / acd["Aero"]["b_Wing"] * 2 * y**3 / 3)
    )


def cy_y2(acd=aircraft_data, y=1):
    return (
        2
        * acd["Aero"]["S_Wing"]
        / (1 + acd["Aero"]["Taper_Wing"])
        / acd["Aero"]["b_Wing"]
        * (y**3 / 3 - (1 - acd["Aero"]["Taper_Wing"]) / acd["Aero"]["b_Wing"] * 2 * y**4 / 4)
    )


def AileronDerivatives(acd=aircraft_data):  # sizing method starting from p 466 Roskam Part VI
    Ch_0 = 0.035

    # Step 1: find the control derivative w.r.t AoA

    # set parameters based on Roskam VI for 2D AoA derivative
    print("Check the parameter set in AileronStickForce if the geometry of the aileron changes")
    # RN = 7 * 10**6  # use the 10^7 line from Roskam
    # Y_90 = 1.272 / 100 * 100  # from CATIA
    # Y_99 = 0.199 / 100 * 100  # from CATIA
    # tan_half_TE = (Y_90 / 2 - Y_99 / 2) / 9

    PHI_doubleprime_te = 5.237 / 180.0 * np.pi  # taken from CATIA and defined in Roskam VI page 467 in radians

    cl_alpha_cl_alpha_theory = 0.9  # taken from graph based on the value of tan_half_TE Roskam VI page 469

    cl_alpha_theory = 7.0  # [/rad] taken from graph Roskam VI page 469

    c_h_alpha_theory = -0.5  # taken from graph Roskam VI page 468

    c_prime_h_alpha_c_h_alpha_theory = 0.75  # taken from graph Roskam VI page 468

    c_prime_h_alpha = c_prime_h_alpha_c_h_alpha_theory * c_h_alpha_theory

    c_doubleprime_h_alpha = c_prime_h_alpha + 2 * cl_alpha_theory * (1 - cl_alpha_cl_alpha_theory) * (
        np.tan(PHI_doubleprime_te / 2) - 0.137
    )  # NOT SURE IF IT IS MAX or LOCAL t/c

    # balance_ratio = np.sqrt(
    #     (0.05 / (acd["Aileron"]["hinge_position"] - 0.05)) ** 2
    #     - (6.31 / 2 / (acd["Aileron"]["hinge_position"] - 0.05) / 100) ** 2
    # )

    c_h_alpha_bal_c_doubleprime_h_alpha = 0.96  # select from Roskam VI page 471 using CATIA and assuming ELLIPTIC nose

    c_h_alpha_bal = c_h_alpha_bal_c_doubleprime_h_alpha * c_doubleprime_h_alpha

    c_h_alpha_M = c_h_alpha_bal / np.sqrt((1 - 0.182**2))

    eta_i = acd["Aileron"]["start"]
    eta_o = acd["Aileron"]["end"]
    K_alpha_eta_i = 2.625  # taken from Roskam VI page 483
    K_alpha_eta_o = 4.0  # taken from Roskam VI page 483
    K_alpha = (K_alpha_eta_i * (1 - eta_i) - K_alpha_eta_o * (1 - eta_o)) / (eta_o - eta_i)

    B2 = 1.0  # taken from Roskam VI page 483 based on cf/c and cb/cf

    deltaC_h_alpha = (
        0.0025 * cl_alpha_cl_alpha_theory * cl_alpha_theory * B2 * K_alpha * 1
    )  # 1 stand for the cosine of the quarter chord sweep of the wing

    C_h_alpha = (acd["Aero"]["AR"] * 1) / (acd["Aero"]["AR"] + 2) * c_h_alpha_M + deltaC_h_alpha

    # Step 2: find the control derivative w.r.t Aileron deflection

    # RN = 7 * 10**6  # use the 10^7 line from Roskam
    # Y_90 = 1.272 / 100 * 100  # from CATIA
    # Y_99 = 0.199 / 100 * 100  # from CATIA
    # tan_half_TE = (Y_90 / 2 - Y_99 / 2) / 9

    PHI_doubleprime_te = 5.237 / 180.0 * np.pi  # taken from CATIA and defined in Roskam VI page 467 in radians

    c_prime_h_delta_c_h_delta_theory = 0.77  # taken from graph based on the value of tan_half_TE Roskam VI page 475
    c_h_delta_theory = -0.85

    cl_delta_theory = 4.56  # taken from Roskam VI page 228
    cl_delta_cl_delta_theory = 0.83  # taken from Roskam VI page 230

    c_prime_h_delta = c_prime_h_delta_c_h_delta_theory * c_h_delta_theory

    c_doubleprime_h_delta = c_prime_h_delta + 2 * cl_delta_theory * (1 - cl_delta_cl_delta_theory) * (
        np.tan(PHI_doubleprime_te / 2) - 0.137
    )  # NOT SURE IF IT IS MAX or LOCAL t/c

    c_h_delta_bal_c_doubleprime_h_delta = 0.89  # select from Roskam VI page 471 using CATIA and assuming ELLIPTIC nose

    c_h_delta_bal = c_h_delta_bal_c_doubleprime_h_delta * c_doubleprime_h_delta

    c_h_delta_M = c_h_delta_bal / np.sqrt((1 - 0.182**2))

    # eta_i = acd["Aileron"]["start"]
    # eta_o = acd["Aileron"]["end"]
    # K_delta_eta_i = 2.43  # taken from Roskam VI page 485
    # K_delta_eta_o = 4.0  # taken from Roskam VI page 485
    # K_delta = (K_delta_eta_i * (1 - eta_i) - K_delta_eta_o * (1 - eta_o)) / (eta_o - eta_i)

    B2 = 1.0  # taken from Roskam VI page 483 based on cf/c and cb/cf
    alpha_sigma = 0.37  # taken from Roskam VI page 230 based on the aileron deflection angle

    deltaC_h_delta = (
        0.0025 * cl_delta_cl_delta_theory * cl_delta_theory * B2 * K_alpha * 1
    )  # 1 stand for the cosine of the quarter chord sweep of the wing

    C_h_delta = (c_h_delta_M + alpha_sigma * c_h_alpha_M * 2 / (acd["Aero"]["AR"] + 2)) + deltaC_h_delta

    # Step 3: update dictionary
    acd["Aileron"]["Ch_0"] = Ch_0
    acd["Aileron"]["C_h_alpha"] = C_h_alpha
    acd["Aileron"]["C_h_delta"] = C_h_delta


def StickArm(acd=aircraft_data, alpha=0.0, delta=0.0, h=3000, V=60.0):
    alpha_rad = alpha * np.pi / 180.0
    delta_rad = delta * np.pi / 180.0

    rho = isa(h, acd["Performance"]["Temp_offset_TO_Land_cruise"])[2]

    Ch = acd["Aileron"]["Ch_0"] + acd["Aileron"]["C_h_alpha"] * alpha_rad + acd["Aileron"]["C_h_delta"] * delta_rad

    chord_a = chord(
        data=acd, location=(acd["Aileron"]["start"] + acd["Aileron"]["end"]) / 2 * (acd["Aero"]["b_Wing"] / 2)
    )

    surface_a = (
        0.5
        * acd["Aileron"]["hinge_position"]
        * (
            chord(data=acd, location=acd["Aileron"]["start"] * (acd["Aero"]["b_Wing"] / 2))
            + chord(data=acd, location=acd["Aileron"]["end"] * (acd["Aero"]["b_Wing"] / 2))
        )
        * (acd["Aileron"]["end"] - acd["Aileron"]["start"])
        * (acd["Aero"]["b_Wing"] / 2)
    )

    HM = 0.5 * rho * V**2 * acd["Aileron"]["hinge_position"] * chord_a * surface_a * Ch

    d_aileron = 0.95 * 0.0063 / 0.1 * chord_a

    ratio = abs(HM) / d_aileron / 133.0

    stick_arm = d_aileron * ratio

    acd["Aileron"]["chord_a"] = chord_a
    acd["Aileron"]["surface_a"] = surface_a
    acd["Aileron"]["d_aileron"] = d_aileron
    acd["Aileron"]["stick_arm"] = stick_arm


def AileronSizing(acd=aircraft_data):
    # define slope of aileron effectiveness graph taken from ADSEE slides
    start_x = 0.2
    start_y = 0.4
    end_x = 0.4
    end_y = 0.6
    slope = (end_y - start_y) / (end_x - start_x)

    # positon of aileron hinge can be varied but always needs to be placed behind
    pos_lst = np.arange(0.25, 0.76, 0.05)

    design_found = False

    for pos in pos_lst:
        tau = slope * (pos - start_x) + start_y

        CL_delta_a = (
            2
            * acd["Aero"]["cl_alpha_airfoil_deg"]
            / (180 / 3.14)
            * tau
            / acd["Aero"]["S_Wing"]
            / acd["Aero"]["b_Wing"]
            * (
                cy_y(acd, 0.95 * acd["Aero"]["b_Wing"] / 2)
                - cy_y(acd, (acd["Flaps"]["flap_end"] + 0.05) * acd["Aero"]["b_Wing"] / 2)
            )
        )

        CL_P = (
            -4
            * (acd["Aero"]["cl_alpha_airfoil_deg"] / (180 / 3.14) + acd["Aero"]["cd_0_airfoil"])
            / acd["Aero"]["S_Wing"]
            / acd["Aero"]["b_Wing"]
            * (cy_y2(acd, acd["Aero"]["b_Wing"] / 2) - cy_y2(acd, 0))
        )

        P_rad = -CL_delta_a / CL_P * 14 * (2 * acd["Performance"]["Vc_m/s"]) / acd["Aero"]["b_Wing"]

        P_deg = P_rad * (180 / 3.14)

        turn_time = round(60 / P_deg, 2)

        if turn_time < (acd["CL2Weight"]["MTOW_N"] / 9.81 + 200) / 590:
            design_found = True
            acd["Aileron"]["start"] = acd["Flaps"]["flap_end"] + 0.05
            acd["Aileron"]["end"] = 0.95
            acd["Aileron"]["roll_rate_deg"] = P_deg
            acd["Aileron"]["roll_rate_rad"] = P_rad
            acd["Aileron"]["CL_delta_a"] = CL_delta_a
            acd["Aileron"]["CL_P"] = CL_P
            acd["Aileron"]["turn_time"] = turn_time
            acd["Aileron"]["hinge_position"] = pos
            # print(acd["Aileron"])
            break

    if design_found:
        return acd
    else:
        raise Exception("No suitable design found.")


if __name__ == "__main__":  # pragma: no cover
    # AileronDerivatives()
    # StickArm(acd=aircraft_data, alpha=0.0, delta=14.0, h=3000.0, V=60.0)
    AileronSizing()
