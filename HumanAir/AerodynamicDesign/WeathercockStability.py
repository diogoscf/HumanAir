# import pandas as pd
import numpy as np
import sys
from math import tan, sqrt  # , pi

# import time
import os

# import json
# import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data

# from HumanAir.isa import isa


def VerticalTailSizing(acd=aircraft_data):
    # change this using CATIA
    print("DONT FORGET TO REVISE PARAMETER FROM CATIA")
    S_fus = 19
    hf1 = 2.22
    hf2 = 1.8
    bf1 = 1.82
    bf2 = 1.82

    k_beta = (
        0.3
        * (
            acd["Stability"]["Cg_Aft"] * acd["Aero"]["MAC_wing"]
            + acd["Stability"]["Cg_Front"] * acd["Aero"]["MAC_wing"]
            + 2 * acd["Geometry"]["XLEMAC_m"]
        )
        / 2
        / acd["Geometry"]["fus_length_m"]
        + 0.75 * acd["Geometry"]["fus_height_m"] / acd["Geometry"]["fus_length_m"]
        - 0.105
    )

    Cn_beta_f = (
        -1
        * k_beta
        * S_fus
        * acd["Geometry"]["fus_length_m"]
        / acd["Aero"]["S_Wing"]
        / acd["Aero"]["b_Wing"]
        * sqrt(hf1 / hf2)
        * (bf2 / bf1) ** (1 / 3)
    )

    Cn_beta_p = (
        -0.053
        * 4
        * (
            (
                acd["Stability"]["Cg_Aft"] * acd["Aero"]["MAC_wing"]
                + acd["Stability"]["Cg_Front"] * acd["Aero"]["MAC_wing"]
                + 2 * acd["Geometry"]["XLEMAC_m"]
            )
            / 2
            + acd["Stability"]["Xcg_prop_m"]
        )
        * acd["Power_prop"]["Dp_m"]
        / acd["Aero"]["S_Wing"]
        / acd["Aero"]["b_Wing"]
    )

    Cn_beta_i = -0.017

    Cn_beta_AH = Cn_beta_f + Cn_beta_p + Cn_beta_i

    end_point_x = -0.12
    end_point_y = 0.09
    start_point_x = -0.02
    start_point_y = 0.04
    slope = (end_point_y - start_point_y) / (end_point_x - start_point_x)

    value = (Cn_beta_AH - start_point_x) * slope + start_point_y

    AR_v_lst = np.arange(1.5, 5, 0.1)

    found = False

    for AR_v in AR_v_lst:
        acd["Aero"]["AR_v"] = AR_v
        acd["Aero"]["S_v"] = (
            value
            * acd["Aero"]["S_Wing"]
            * acd["Aero"]["b_Wing"]
            / (
                (
                    0.9 * acd["Geometry"]["fus_length_m"]
                    - (
                        acd["Stability"]["Cg_Aft"] * acd["Aero"]["MAC_wing"]
                        + acd["Stability"]["Cg_Front"] * acd["Aero"]["MAC_wing"]
                        + 2 * acd["Geometry"]["XLEMAC_m"]
                    )
                    / 2
                )
            )
        )
        acd["Aero"]["b_v"] = sqrt(acd["Aero"]["AR_v"] * acd["Aero"]["S_v"])
        acd["Aero"]["c_root_v"] = 2 * acd["Aero"]["S_v"] / ((1 + acd["Aero"]["Taper_v"]) * acd["Aero"]["b_v"])
        acd["Aero"]["c_tip_v"] = acd["Aero"]["Taper_v"] * acd["Aero"]["c_root_v"]
        acd["Aero"]["MAC_v"] = (
            2
            * acd["Aero"]["c_root_v"]
            / 3
            * (1 + acd["Aero"]["Taper_v"] + acd["Aero"]["Taper_v"] ** 2)
            / (1 + acd["Aero"]["Taper_v"])
        )
        acd["Aero"]["MAC_y"] = acd["Aero"]["b_v"] / 6 * (1 + 2 * acd["Aero"]["Taper_v"]) / (1 + acd["Aero"]["Taper_v"])

        if acd["Aero"]["c_root_v"] * 0.75 < 0.1 * acd["Geometry"]["fus_length_m"] + acd["Aero"]["MAC_y"] / tan(0.628):
            found = True
            break

    if not found:
        raise Exception("Not valid configuration found. Revise/update design parameters!")
    # else:
    #     print(
    #         acd["Aero"]["c_root_v"], acd["Aero"]["c_tip_v"],
    # acd["Aero"]["S_v"], acd["Aero"]["b_v"], acd["Aero"]["AR_v"]
    #     )
    else:
        return acd


if __name__ == "__main__":  # pragma: no cover
    VerticalTailSizing()
