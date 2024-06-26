import sys
import os
import numpy as np
from math import isclose

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from HumanAir.CO2_Calculator.conceptual_co2 import (
    _calculate_standard_co2_flight_lengths,
    _calculate_new_co2_flight_lengths,
    _calculate_co2_reduction_flight_lengths,
    calculate_standard_co2_average_flight,
    calculate_new_co2_average_flight,
    calculate_co2_reduction_average_flight,
)


def test__calculate_standard_co2_flight_lengths():
    maintanance_co2 = 0.48
    co2 = _calculate_standard_co2_flight_lengths(mission_freqs=np.array([[1, 1, 1], [1, 1, 1]]))

    assert isclose(maintanance_co2, co2[1], rel_tol=0.1)


def test__calculate_new_co2_flight_lengths():
    maintanance_co2 = 511
    co2 = _calculate_new_co2_flight_lengths(mission_freqs=np.array([[1, 1, 1], [1, 1, 1]]))

    assert isclose(maintanance_co2, co2[1], rel_tol=0.01)


def test__calculate_co2_reduction_flight_lengths():
    expected_co2 = 25 / 100
    co2 = _calculate_co2_reduction_flight_lengths()

    assert isclose(expected_co2, co2, rel_tol=0.1)


def test_calculate_standard_co2_average_flight():
    expected_co2 = 12 / 100
    co2 = calculate_standard_co2_average_flight(avg_prefuel_flight_length_nm=1)

    assert isclose(expected_co2, co2[1], rel_tol=0.01)


def test_calculate_new_co2_average_flight():
    expected_co2 = 13.5 / 100
    co2 = calculate_new_co2_average_flight(avg_prefuel_flight_length_nm=1, average_legs=1)

    assert isclose(expected_co2, co2, rel_tol=0.01)


def test_calculate_co2_reduction_average_flight():
    expected_co2 = 48 / 100
    co2 = calculate_co2_reduction_average_flight()

    assert isclose(expected_co2, co2, rel_tol=0.01)
