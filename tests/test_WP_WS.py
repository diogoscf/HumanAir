import numpy as np
from math import isclose
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.LoadingDiagram.Main import WP_WS


def test_initialise_WP_WS():
    wp_ws = WP_WS()

    # Verify WS is correctly initialized
    np.testing.assert_array_almost_equal(wp_ws.WS[:5], np.array([0.0001, 0.1001, 0.2001, 0.3001, 0.4001]), decimal=4)
    assert len(wp_ws.WS) == int(2000 / 0.1)

    # Verify ReferenceWS
    expected_ReferenceWS = [
        449.0307667,
        602.1288757,
        643.1831172,
        725.2916003,
        725.2916003,
        578.0625962,
        789.2055899,
        780.3381114,
        930.8329016,
        925.2755192,
        1274.01538,
        1286.972761,
        1236.320571,
        987.5590745,
        1321.946577,
        1662.527136,
        1437.958113,
        1871.550807,
        1837.537433,
    ]
    np.testing.assert_array_almost_equal(wp_ws.ReferenceWS, expected_ReferenceWS, decimal=4)

    # Verify ReferenceWP
    expected_ReferenceWP = [
        0.089508018,
        0.082049016,
        0.080130987,
        0.068752535,
        0.068752535,
        0.058831561,
        0.066557244,
        0.069219534,
        0.088135561,
        0.078032631,
        0.06465638,
        0.064313168,
        0.062743365,
        0.066302235,
        0.055426119,
        0.066832653,
        0.065302355,
        0.093120238,
        0.042435741,
    ]
    np.testing.assert_array_almost_equal(wp_ws.ReferenceWP, expected_ReferenceWP, decimal=4)


def test_calculate_optimal_point():
    wp_ws = WP_WS()
    optimal_WP, optimal_WS = wp_ws.calculate_optimal_point()

    # Add the expected values for optimal_WP and optimal_WS
    expected_optimal_WP = 0.1195
    expected_optimal_WS = 796.4

    assert isclose(optimal_WP, expected_optimal_WP, rel_tol=1e-2)
    assert isclose(optimal_WS, expected_optimal_WS, rel_tol=1e-2)
