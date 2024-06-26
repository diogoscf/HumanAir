# import pandas as pd
# import numpy as np
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data

# from isa import isa


def BrakeCheck(acd=aircraft_data, P_max=0.0, A_calliper=0.0, r_in=0.0, r_out=0.0):
    KE_req = 0.5 * acd["CL2Weight"]["MTOW_N"] / 9.81 * 29.2**2 / 2

    d_ground = 508  # without taking into account reverse thrust

    F_retarded = KE_req / d_ground

    F_calliper = F_retarded * (acd["Landing_gear"]["Dwm_m"] / 2) / ((r_in + r_out) / 2)

    mu = 0.15

    # print(2 * P_max * 10**6 * A_calliper * mu, F_calliper)

    return 2 * P_max * 10**6 * A_calliper * mu > F_calliper


if __name__ == "__main__":  # pragma: no cover
    # AileronDerivatives()
    print(BrakeCheck(acd=aircraft_data, P_max=1.2, A_calliper=0.01, r_in=0.1, r_out=0.25))
    # AileronSizing()
