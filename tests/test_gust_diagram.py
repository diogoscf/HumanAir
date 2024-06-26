import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.Vn_Diagrams.gust_diagram import calculate_gust_diagram_loads


def test_calculate_gust_diagram_loads():
    # Mock aircraft data
    aircraft_data = {
        "Performance": {
            "W/S_N/m2": 300,
        },
        "Aero": {
            "MAC_wing": 5.0,
            "CLalpha": 0.1,
            "CLmax_clean": 1.8,
        },
    }

    # Test case 1: Normal scenario
    Vc_ms = 80.0
    Vd_ms = 100.0
    V_S1 = 60.0
    (
        n_max,
        n_min,
        V_B,
        n_cruise_pve,
        n_cruise_nve,
        n_dive_pve,
        n_dive_nve,
        n_B_pve,
        n_B_nve,
    ) = calculate_gust_diagram_loads(aircraft_data, Vc_ms, Vd_ms, V_S1)

    assert n_max > n_min  # Basic validation
    assert V_B > 0.0  # V_B should be positive
    assert n_cruise_pve > n_cruise_nve  # Load factors should be consistent

    # Test case 2: Edge case with commuter aircraft
    Vc_ms = 100.0
    Vd_ms = 120.0
    V_S1 = 25.0
    commuter_ac = True
    (
        n_max,
        n_min,
        V_B,
        n_cruise_pve,
        n_cruise_nve,
        n_dive_pve,
        n_dive_nve,
        n_B_pve,
        n_B_nve,
    ) = calculate_gust_diagram_loads(aircraft_data, Vc_ms, Vd_ms, V_S1, commuter_ac=commuter_ac)

    assert n_max > n_min  # Basic validation
    assert V_B > 0.0  # V_B should be positive
    assert n_cruise_pve > n_cruise_nve  # Load factors should be consistent
