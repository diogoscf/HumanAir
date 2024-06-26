import os
import sys
from math import isclose
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.Vn_Diagrams.loading_diagram import (
    nmax_manoeuvre,
    nmin_manoeuvre,
    calc_nmax_nmin_manoeuvre,
    calculate_manoeuvre_velocities,
)


def test_nmax_manoeuvre():
    MTOW_lbs = 10000

    # test both cases
    assert isclose(nmax_manoeuvre(MTOW_lbs), 3.3, rel_tol=0.01)
    assert isclose(nmax_manoeuvre(320), 3.8, rel_tol=0.01)


def test_nmin_manoeuvre():
    MTOW_lbs = 10000

    # test both cases
    assert isclose(nmin_manoeuvre(MTOW_lbs), -1.32, rel_tol=0.01)
    assert isclose(nmin_manoeuvre(320), -1.52, rel_tol=0.01)


def test_calc_nmax_nmin_manoeuvre():
    MTOW_N = 45350.92
    nmax, nmin = calc_nmax_nmin_manoeuvre(MTOW_N)

    assert isclose(nmax, 3.3, rel_tol=0.01)
    assert isclose(nmin, -1.32, rel_tol=0.01)


mock_aircraft_data = {
    "Weights": {"MTOW_N": 45350.92},
    "Performance": {
        "W/S_N/m2": 500,
        "Vc_m/s": 100,
    },
    "Aero": {
        "CLmax_clean": 1.5,
        "CLmax_Land": 2.0,
    },
}


def test_calculate_manoeuvre_velocities():
    Vc_ms, Vd_ms, V_A, V_S1, V_HH, V_S0 = calculate_manoeuvre_velocities(mock_aircraft_data)

    V_A_expected = np.sqrt(2 * 3.3 * 500 / (1.5 * 1.225))
    V_S1_expected = np.sqrt(2 * 500 / (1.5 * 1.225))
    V_HH_expected = np.sqrt(2 * 1.32 * 500 / (1.5 * 1.225))
    V_SO_expected = np.sqrt(2 * 500 / (2 * 1.225))

    assert isclose(Vc_ms, 100, rel_tol=0.01)
    assert isclose(Vd_ms, 1.4 * 100, rel_tol=0.01)
    assert isclose(V_A, V_A_expected, rel_tol=0.01)
    assert isclose(V_S1, V_S1_expected, rel_tol=0.01)
    assert isclose(V_HH, V_HH_expected, rel_tol=0.01)
    assert isclose(V_S0, V_SO_expected, rel_tol=0.01)
