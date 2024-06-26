from math import isclose
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from HumanAir.CO2_Calculator.get_flight_distribution import fit_log_normal


def test_fit_log_normal():
    # Test case 1: median = 1, mean > median
    mean = 3
    median = 1
    result = fit_log_normal(mean, median)
    assert isclose(result.mean(), mean, rel_tol=1e-9), f"Test case 1 failed: expected mean {mean}, got {result.mean()}"

    # Test case 2: median = 1, mean = median
    mean = 2
    median = 1
    result = fit_log_normal(mean, median)
    assert isclose(result.mean(), mean, rel_tol=1e-9), f"Test case 2 failed: expected mean {mean}, got {result.mean()}"
