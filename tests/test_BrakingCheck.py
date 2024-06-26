import os
import sys

# from math import isclose

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from HumanAir.Other.BrakingCheck import BrakeCheck

# Mock aircraft data for testing
mock_aircraft_data = {
    "CL2Weight": {
        "MTOW_N": 98000.0,  # Maximum takeoff weight in Newtons
    },
    "Landing_gear": {
        "Dwm_m": 1.0,  # Distance from wheel center to caliper mounting point in meters
    },
}


def test_BrakeCheck():
    # Define test parameters
    P_max = 5.0  # Maximum pressure in MPa
    A_calliper = 0.05  # Calliper area in square meters
    r_in = 0.2  # Inner radius in meters
    r_out = 0.4  # Outer radius in meters

    # Expected result based on calculations
    expected_result = True

    result = BrakeCheck(mock_aircraft_data, P_max, A_calliper, r_in, r_out)

    assert result == expected_result
