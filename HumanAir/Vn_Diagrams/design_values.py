import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.Vn_Diagrams.loading_diagram import calculate_manoeuvre_velocities, calc_nmax_nmin_manoeuvre
from HumanAir.Vn_Diagrams.gust_diagram import calculate_gust_diagram_loads
from HumanAir.aircraft_data import aircraft_data
from HumanAir.isa import isa


def calculate_load_design_values(aircraft_data):
    h_cruise = aircraft_data["Performance"]["Altitude_Cruise_m"]
    h_land = aircraft_data["Performance"]["Altitude_Land_m"]
    temp_offset = aircraft_data["Performance"]["Temp_offset_TO_Land_cruise"]

    Vc_ms, Vd_ms, _, V_S1, _, _ = calculate_manoeuvre_velocities(aircraft_data)
    n_max_manoeuvre, _ = calc_nmax_nmin_manoeuvre(aircraft_data["Weights"]["MTOW_N"])
    n_max_gust_cruise, *_ = calculate_gust_diagram_loads(
        aircraft_data, Vc_ms, Vd_ms, V_S1, h=h_cruise, temp_offset=temp_offset
    )

    n_max_gust_land, *_ = calculate_gust_diagram_loads(
        aircraft_data, Vc_ms, Vd_ms, V_S1, h=h_land, temp_offset=temp_offset
    )

    n_max_cruise = max(n_max_manoeuvre, n_max_gust_cruise)
    n_ult_cruise = 1.5 * n_max_cruise

    n_max_land = max(n_max_manoeuvre, n_max_gust_land)
    n_ult_land = 1.5 * n_max_land

    V_H = Vc_ms / 0.9

    gamma, R = 1.4, 287
    T = isa(h_cruise)[0]
    M_D = Vd_ms / np.sqrt(gamma * R * T)

    return M_D, V_H, n_ult_cruise, n_ult_land


def save_to_acdata_dict(aircraft_data, M_D, V_H, n_ult_c, n_ult_l):  # pragma: no cover
    aircraft_data["Performance"]["n_ult"] = n_ult_c
    aircraft_data["Performance"]["n_ult_l"] = n_ult_l
    aircraft_data["Performance"]["M_D"] = M_D
    aircraft_data["Performance"]["Vh_m/s"] = V_H


if __name__ == "__main__":  # pragma: no cover
    M_D, V_H, n_ult_c, n_ult_l = calculate_load_design_values(aircraft_data)
    print(f"M_D = {M_D}")
    print(f"V_H = {V_H}")
    print(f"n_ult = {n_ult_c}")
    print(f"n_ult.l = {n_ult_l}")
