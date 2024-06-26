import sys
import os
from math import isclose
import tempfile
import json

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.AerodynamicDesign.LongitudinalStability import LongitudinalStability
from HumanAir.AerodynamicDesign.PlanformDesign import Planform

AR = 8.0
Taper = 0.3
QuarterChordSweep = 25.0
tc = 0.1
MTOW = 5000.0
WS = 100.0
S = 1

wing = Planform(AR, Taper, QuarterChordSweep, tc, MTOW=MTOW, WS=WS, S=S)


def test_CMac_Wing():
    # Define sample inputs for the LongitudinalStability class
    CLh = 0.5
    CLah = 5.7
    Xcgh = 0.5
    XLEMAC = 2.0
    CgAft = 0.6
    CgFwd = 0.4
    SM = 0.1
    deda = 0.1
    VhV = 1.0
    FuselageLength = 10.0

    airfoil_wing_data = {"C_L_Alpha": 5.7, "Cm_0": -0.1}
    airfoil_h_data = {"C_L_Alpha": 4.5, "Cm_0": -0.05}

    with (
        tempfile.NamedTemporaryFile(delete=False, mode="w") as wing_file,
        tempfile.NamedTemporaryFile(delete=False, mode="w") as h_file,
    ):
        json.dump(airfoil_wing_data, wing_file)
        json.dump(airfoil_h_data, h_file)
        wing_file_name = wing_file.name
        h_file_name = h_file.name

    stability = LongitudinalStability(
        CLh, CLah, Xcgh, XLEMAC, CgAft, CgFwd, SM, deda, VhV, FuselageLength, wing, wing_file_name, h_file_name
    )
    expected_CMac_Wing = -0.0669
    assert isclose(stability.CMac_Wing(), expected_CMac_Wing, rel_tol=1e-3)


def test_Stability():
    # Define sample inputs for the LongitudinalStability class
    CLh = 0.5
    CLah = 5.7
    Xcgh = 0.5
    XLEMAC = 2.0
    CgAft = 0.6
    CgFwd = 0.4
    SM = 0.1
    deda = 0.1
    VhV = 1.0
    FuselageLength = 10.0

    airfoil_wing_data = {"C_L_Alpha": 5.7, "Cm_0": -0.1}
    airfoil_h_data = {"C_L_Alpha": 4.5, "Cm_0": -0.05}

    with (
        tempfile.NamedTemporaryFile(delete=False, mode="w") as wing_file,
        tempfile.NamedTemporaryFile(delete=False, mode="w") as h_file,
    ):
        json.dump(airfoil_wing_data, wing_file)
        json.dump(airfoil_h_data, h_file)
        wing_file_name = wing_file.name
        h_file_name = h_file.name

    stability = LongitudinalStability(
        CLh, CLah, Xcgh, XLEMAC, CgAft, CgFwd, SM, deda, VhV, FuselageLength, wing, wing_file_name, h_file_name
    )
    expected_stability = 0.067
    assert isclose(stability.Stability()[50], expected_stability, rel_tol=1e-2)  # Checking a specific index


def test_Controllability():
    # Define sample inputs for the LongitudinalStability class
    CLh = 0.5
    CLah = 5.7
    Xcgh = 0.5
    XLEMAC = 2.0
    CgAft = 0.6
    CgFwd = 0.4
    SM = 0.1
    deda = 0.1
    VhV = 1.0
    FuselageLength = 10.0

    airfoil_wing_data = {"C_L_Alpha": 5.7, "Cm_0": -0.1}
    airfoil_h_data = {"C_L_Alpha": 4.5, "Cm_0": -0.05}

    with (
        tempfile.NamedTemporaryFile(delete=False, mode="w") as wing_file,
        tempfile.NamedTemporaryFile(delete=False, mode="w") as h_file,
    ):
        json.dump(airfoil_wing_data, wing_file)
        json.dump(airfoil_h_data, h_file)
        wing_file_name = wing_file.name
        h_file_name = h_file.name

    stability = LongitudinalStability(
        CLh, CLah, Xcgh, XLEMAC, CgAft, CgFwd, SM, deda, VhV, FuselageLength, wing, wing_file_name, h_file_name
    )
    expected_controllability = 0.37
    assert isclose(stability.Controllability()[50], expected_controllability, rel_tol=1e-2)  # Checking a specific index


def test_ShS():
    # Define sample inputs for the LongitudinalStability class
    CLh = 0.5
    CLah = 5.7
    Xcgh = 0.5
    XLEMAC = 2.0
    CgAft = 0.6
    CgFwd = 0.4
    SM = 0.1
    deda = 0.1
    VhV = 1.0
    FuselageLength = 10.0

    airfoil_wing_data = {"C_L_Alpha": 5.7, "Cm_0": -0.1}
    airfoil_h_data = {"C_L_Alpha": 4.5, "Cm_0": -0.05}

    with (
        tempfile.NamedTemporaryFile(delete=False, mode="w") as wing_file,
        tempfile.NamedTemporaryFile(delete=False, mode="w") as h_file,
    ):
        json.dump(airfoil_wing_data, wing_file)
        json.dump(airfoil_h_data, h_file)
        wing_file_name = wing_file.name
        h_file_name = h_file.name

    stability = LongitudinalStability(
        CLh, CLah, Xcgh, XLEMAC, CgAft, CgFwd, SM, deda, VhV, FuselageLength, wing, wing_file_name, h_file_name
    )
    expected_ShS = 0.214
    assert isclose(stability.ShS(), expected_ShS, rel_tol=1e-2)
