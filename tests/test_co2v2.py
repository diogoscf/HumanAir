import os
import sys
import numpy as np
from math import isclose

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.CO2_Calculator.co2v2 import (
    calculate_co2_reduction_flightdist,
    calculate_new_co2_flightdist,
    calculate_standard_co2_flightdist,
)


def test_calculate_standard_co2_flightdist():
    expected_co2 = 0.12
    actual_co2 = calculate_standard_co2_flightdist(flightdist=np.array([1, 1, 1]))

    assert isclose(actual_co2[1], expected_co2, rel_tol=1e-1)


def test_calculate_new_co2_flightdist():
    expected_co2 = 0.135
    actual_co2 = calculate_new_co2_flightdist(
        flightdist=np.array([1, 1, 1]), mission_freqs=np.array([[1, 1, 1], [1, 1, 1]])
    )

    assert isclose(actual_co2[1], expected_co2, rel_tol=1e-1)


def test_calculate_co2_reduction_flightdist():
    expected_co2 = 0.87
    actual_co2 = calculate_co2_reduction_flightdist(flightdist=np.array([1, 1, 1]))

    assert isclose(actual_co2, expected_co2, rel_tol=1e-1)
