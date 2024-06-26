import os
import sys
from math import isclose, tan, sqrt

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from HumanAir.Weights_and_CG.weight_fractions import find_lg, component_mass, iterate_cg_lg


aircraft_data = {
    "Weights": {
        "MTOW_N": 5000,  # Maximum Takeoff Weight [N]
        "OEW_N": 300,  # Operating Empty Weight [N]
        "Wbat_N": 200,  # Weight of Batteries [N]
        "Ww_N": 100,  # Weight of Wings [N]
        "Wfuel_N": 500,  # Weight of Fuel [N]
        "Wpl_des_kg": 150,  # Design Payload Weight [kg]
    },
    "Geometry": {
        "fus_height_m": 2.5,  # Fuselage height [m]
        "fus_length_m": 10.0,  # Fuselage length [m]
        "tail_length_m": 2.0,  # Length of tail [m]
        "XLEMAC_m": 2.0,  # Distance from nose to Leading Edge of MAC [m]
    },
    "Aero": {
        "MAC_wing": 1.5,  # Mean Aerodynamic Chord of the wing [m]
    },
    "Landing_gear": {
        # Placeholder for calculated values
        "lm_m": None,  # Main gear length [m]
        "ln_m": None,  # Nose gear length [m]
        "ymin_m": None,  # Minimum ground clearance [m]
        "Hs_m": None,  # Strike height [m]
        "Dwm_m": None,  # Main wheel diameter [m]
        "Dwn_m": None,  # Nose wheel diameter [m]
        "Wtm_m": None,  # Main wheel load capacity [N]
        "Wtn_m": None,  # Nose wheel load capacity [N]
        "PMW_N": None,  # Main wheel load at MTOW [N]
        "PNW_N": None,  # Nose wheel load at MTOW [N]
        "Xmw_m": None,  # Main gear CG location [m]
        "Xnw_m": None,  # Nose gear CG location [m]
    },
}


def test_find_lg():
    # initialise components
    nose_loading = 0.08
    aftcg = 5
    l_m, l_n, Pmg, Pnw, H_s, Pmw, ymin, Dw_m, Wt_m, Dw_n, Wt_n = find_lg(nose_loading, aftcg, ac_datafile=aircraft_data)

    # expected values
    expected_l_m = tan(0.279253) * (1.25 + 0.771 + 0.5 * 0.329565)
    expected_l_n = expected_l_m * 468.9092 / 40.774
    expected_ymin = (expected_l_m + expected_l_n) / (
        sqrt(expected_l_n**2 * tan(0.959931) ** 2 / (1.25 + 0.771 + 0.5 * 0.329565) ** 2 - 1)
    )

    # Add actual expected values for your test cases
    assert isclose(l_m, expected_l_m, rel_tol=1e-2)
    assert isclose(l_n, expected_l_n, rel_tol=1e-2)
    assert isclose(Pmg, 468.9092, rel_tol=1e-2)
    assert isclose(Pnw, 40.774, rel_tol=1e-2)
    assert isclose(H_s, 0.77, rel_tol=1e-2)
    assert isclose(Pmw, 234.4546, rel_tol=1e-1)
    assert isclose(ymin, expected_ymin, rel_tol=1e-2)
    assert isclose(Dw_m, 0.329565, rel_tol=1e-2)
    assert isclose(Wt_m, 0.127, rel_tol=1e-2)
    assert isclose(Dw_n, 0.329565, rel_tol=1e-2)
    assert isclose(Wt_n, 0.127, rel_tol=1e-2)


def test_component_mass():
    wcg = component_mass(aircraft_data)
    # expected values
    OEW_frac = 300 / 5000
    Wbat_frac = 200 / 5000
    Ww_frac = 100 / 5000

    # test the weight fractions
    assert isclose(wcg[0, 0], Ww_frac, rel_tol=1e-2)
    assert isclose(wcg[0, -2], Wbat_frac, rel_tol=1e-2)
    assert isclose(wcg[0, -1], OEW_frac, rel_tol=1e-2)

    # test the mass values
    assert isclose(wcg[1, 0], (Ww_frac * 1000) / 1.9627, rel_tol=1e-2)
    assert isclose(wcg[1, -2], (Wbat_frac * 2000) / 3.924, rel_tol=1e-2)
    assert isclose(wcg[1, -1], (OEW_frac * 300) / 0.5886, rel_tol=1e-2)

    # test the cg values
    assert isclose(wcg[2, 0], 0, rel_tol=1e-2)
    assert isclose(wcg[2, -2], 0, rel_tol=1e-2)
    assert isclose(wcg[2, -1], 0, rel_tol=1e-2)


def test_iterate_cg_lg():
    wcg, CGlist, xlemac = iterate_cg_lg(aircraft_data, PERCENTAGE=0.5)

    Ww_frac = 100 / 5000

    powertrain_location = 0.05 * 10

    # test the xcg locations
    assert isclose(wcg[0, 0], Ww_frac, rel_tol=1e-2)
    assert isclose(wcg[2, 2], powertrain_location, rel_tol=1e-2)  # powertrain location

    # test the first position of the cg loadout
    assert isclose(CGlist[0], 3.925, rel_tol=1e-2)

    # test the xlemac
    assert isclose(xlemac, 3.2, rel_tol=1e-2)
