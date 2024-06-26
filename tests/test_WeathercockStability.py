import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from HumanAir.AerodynamicDesign.WeathercockStability import VerticalTailSizing

# Test data
test_data = {
    "Stability": {"Cg_Aft": 0.4, "Cg_Front": 0.3, "Xcg_prop_m": 2.0},
    "Aero": {
        "MAC_wing": 1.5,
        "S_Wing": 30.0,
        "b_Wing": 10.0,
        "Taper_v": 0.5,
    },
    "Geometry": {"XLEMAC_m": 2.0, "fus_length_m": 30.0, "fus_height_m": 5.0},
    "Power_prop": {"Dp_m": 2.0},
}


def test_vertical_tail_sizing():
    ac_data = test_data.copy()
    result = VerticalTailSizing(ac_data)

    assert "AR_v" in result["Aero"]
    assert "S_v" in result["Aero"]
    assert "b_v" in result["Aero"]
    assert "c_root_v" in result["Aero"]
    assert "c_tip_v" in result["Aero"]
    assert "MAC_v" in result["Aero"]
    assert "MAC_y" in result["Aero"]
