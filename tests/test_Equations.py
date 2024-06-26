import os
import sys
from math import isclose
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.LoadingDiagram.Equations import (
    Density,
    Cd,
    Stallspeedx,
    Takeoff,
    Landingx,
    Cruise,
    Climbrate,
    Climbgradient,
)


def test_Density():
    assert isclose(Density(0, 0), 1.225, rel_tol=1e-4)


def test_Cd():
    # expected value
    expected_cd = 0.02 + 0.5**2 / (np.pi * 8 * 0.8)

    assert isclose(Cd(0.02, 8, 0.8, 0.5), expected_cd, rel_tol=1e-4)


def test_Stallspeedx():
    # expected value
    expected_stallspeed = 0.5 * 1.225 * 30**2 * 1.2
    assert isclose(Stallspeedx(0, 0, 30, 1.2), expected_stallspeed, rel_tol=1e-4)


def test_Takeoff():
    # expected value
    expected_takeoff = 100 * 1.5 / 500  # / 1.21

    assert isclose(Takeoff(100, 500, 0, 0, 1.5), expected_takeoff, rel_tol=1e-4)


def test_Landingx():
    # expected value
    expected_landing = (2 * 1.225 * 500 / 0.305199384478051392) / (2 * 2)

    assert isclose(Landingx(2.0, 0, 0, 500, 2), expected_landing, rel_tol=1e-4)


def test_Cruise():
    # expected value
    expected_cruise = (
        (0.9 / 0.95)
        * 0.85
        * (1.225 / 1.225) ** 0.75
        * (0.02 * 0.5 * 1.225 * 100**3 / (0.9 * 500) + 500 * 0.95 / (np.pi * 8 * 0.8 * 0.5 * 1.225 * 100)) ** (-1)
    )

    assert isclose(Cruise(0.85, 0, 0, 0.02, 100, 500, 8, 0.8, 0.9, 0.95), expected_cruise, rel_tol=1e-1)


def test_Climbrate():
    # expected value
    expected_climbrate = 0.85 / (5 + (np.sqrt(500) * np.sqrt(2 / 1.225)) / (1.345 * (8 * 0.8) ** 0.75 / (0.02**0.25)))

    assert isclose(Climbrate(0.85, 8, 0.8, 0, 0, 0.02, 5, 500), expected_climbrate, rel_tol=1e-4)


def test_Climbgradient():
    # expected value
    expected_climbgradient = 2.41 * 0.85 / (np.sqrt(500) * (0.02 + 0.02 / 0.5) / (0.5) * np.sqrt(2 / (1.225 * 0.5)))

    assert isclose(
        Climbgradient(0.85, 500, 0.02, 100, 8, 0.8, 1.2, 1.1, 0, 0, temp_offset=0), expected_climbgradient, rel_tol=1e-1
    )
