from math import isclose
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from HumanAir.AerodynamicDesign.FlapsDesign import chord, flaps_design


# Test data
test_data = {
    "Aero": {"S_Wing": 40.0, "Taper_Wing": 0.5, "b_Wing": 10.0, "CLmax_clean": 1.4, "CLmax_Land": 2, "CLmax_TO": 1.8},
    "Geometry": {"fus_width_m": 3.0},
    "Flaps": {"deflection_landing": 15, "deflection_takeoff": 15},
}


def test_chord():
    location = 1.9
    expected_result = 2 * 40 / (1 + 0.5) / 10 * (1 - (1 - 0.5) / 10 * 2 * location)
    result = chord(data=test_data, location=1.9)
    assert isclose(result, expected_result, rel_tol=1e-9)


def test_flaps_design():
    ac_data = test_data.copy()
    flaps_design(ac_data)

    assert "flap_start" in ac_data["Flaps"]
    assert "flap_end" in ac_data["Flaps"]
    assert "Swf" in ac_data["Flaps"]
    assert "cf_c" in ac_data["Flaps"]
    assert "deltaCLmax_land" in ac_data["Flaps"]
    assert "deltaCLmax_takeoff" in ac_data["Flaps"]
    assert "cprime_c_landing" in ac_data["Flaps"]
    assert "cprime_c_takeoff" in ac_data["Flaps"]
    assert "AoA_landing" in ac_data["Flaps"]
    assert "AoA_takeoff" in ac_data["Flaps"]
    assert "CL_AoA0_landing" in ac_data["Flaps"]
    assert "CL_AoA0_takeoff" in ac_data["Flaps"]
    assert "Sprime_S_landing" in ac_data["Flaps"]
    assert "Sprime_S_takeoff" in ac_data["Flaps"]

    # Additional checks can be performed for specific expected values if known
