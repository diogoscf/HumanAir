import os
import sys
from math import isclose


sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.Vn_Diagrams.design_values import calculate_load_design_values

# Mock aircraft data for testing
mock_aircraft_data = {
    "Performance": {
        "Altitude_Cruise_m": 10000.0,
        "Altitude_Land_m": 0.0,
        "Temp_offset_TO_Land_cruise": 15.0,
        "Vc_m/s": 100.0,
        "W/S_N/m2": 40,
    },
    "Weights": {
        "MTOW_N": 45350.92,
    },
    "Aero": {"CLmax_clean": 1.8, "CLmax_Land": 2.0, "MAC_wing": 2, "CLalpha": 1.2},
}


def test_calculate_load_design_values():
    M_D, V_H, n_ult_cruise, n_ult_land = calculate_load_design_values(mock_aircraft_data)

    assert isclose(M_D, 0.465, rel_tol=0.01)
    assert isclose(V_H, 111.11, rel_tol=0.01)
    assert isclose(n_ult_cruise, 24.55, rel_tol=0.01)
    assert isclose(n_ult_land, 14.62, rel_tol=0.01)
