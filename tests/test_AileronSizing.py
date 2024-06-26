import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from HumanAir.AerodynamicDesign.AileronSizing import chord, cy_y, cy_y2, AileronDerivatives, StickArm, AileronSizing


test_data = {
    "Aero": {
        "S_Wing": 30.0,
        "Taper_Wing": 0.5,
        "b_Wing": 10.0,
        "cl_alpha_airfoil_deg": 5.73,
        "cd_0_airfoil": 0.015,
        "AR": 7.5,
    },
    "Aileron": {
        "start": 0.3,
        "end": 0.7,
        "hinge_position": 0.5,
    },
    "Flaps": {
        "flap_end": 0.3,
    },
    "Performance": {
        "Vc_m/s": 70,
        "Temp_offset_TO_Land_cruise": 0.0,
    },
    "CL2Weight": {
        "MTOW_N": 200000,
    },
}


def test_chord():
    result = chord(test_data, location=1.9)
    assert isinstance(result, float)


def test_cy_y():
    result = cy_y(test_data, y=1)
    assert isinstance(result, float)


def test_cy_y2():
    result = cy_y2(test_data, y=1)
    assert isinstance(result, float)


def test_AileronDerivatives():
    ac_data = test_data.copy()
    AileronDerivatives(ac_data)
    assert "Ch_0" in ac_data["Aileron"]
    assert "C_h_alpha" in ac_data["Aileron"]
    assert "C_h_delta" in ac_data["Aileron"]


def test_StickArm():
    ac_data = test_data.copy()
    StickArm(ac_data, alpha=0.0, delta=14.0, h=3000.0, V=60.0)
    assert "chord_a" in ac_data["Aileron"]
    assert "surface_a" in ac_data["Aileron"]
    assert "d_aileron" in ac_data["Aileron"]
    assert "stick_arm" in ac_data["Aileron"]


def test_AileronSizing():
    ac_data = test_data.copy()
    result = AileronSizing(ac_data)
    assert "start" in result["Aileron"]
    assert "end" in result["Aileron"]
    assert "roll_rate_deg" in result["Aileron"]
    assert "roll_rate_rad" in result["Aileron"]
    assert "CL_delta_a" in result["Aileron"]
    assert "CL_P" in result["Aileron"]
    assert "turn_time" in result["Aileron"]
    assert "hinge_position" in result["Aileron"]
