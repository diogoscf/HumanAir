from math import isclose
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.FuselageSizing.FuselageSizing import FuselageSizing

mock_aircraft_data = {
    "General": {"N_pax": 150},
    "Geometry": {
        "w_engine": 2.5,
        "h_engine": 1.2,
        "l_engine": 3.0,
        "s_engine": 0.2,
        "h_battery": 0.3,
        "volume_battery": 0.5,
        "tail_length_m": 2.5,
        "h_tail": 1,
        "l_empty": 0.5,
        "w_aisle": 1,
        "h_aisle": 1.5,
        "w_motor": 0.8,
        "l_motor": 2.0,
    },
    "Landing_gear": {
        "Dwn_m": 0.8,
        "Wtn_m": 0.3,
        "Dwm_m": 1.0,
        "lm_m": 5.0,
        "ln_m": 3.0,
        "tail_length_m": 2.5,
        "h_tail": 1.0,
        "ymin_m": 2.0,
        "Xnw_m": 4.0,
        "Xmw_m": 6.0,
        "Wtm_m": 0.4,
        "Dtm_m": 0.5,
        "l_s_m": 0.5,
        "l_s_n": 0.2,
    },
}

# Creating an instance of FuselageSizing
fuselage = FuselageSizing(mock_aircraft_data, 0.35)


def test_fuselage_sizing_init():
    # Verify instance variables are initialized correctly
    assert fuselage.bat_density == 1000
    assert fuselage.n_seat == 151
    assert fuselage.w_engine == 2.5
    assert fuselage.h_engine == 1.2
    assert fuselage.l_engine == 3.0
    assert fuselage.s_engine == 0.2
    assert fuselage.D_nose == 0.8
    assert fuselage.h_nose == 0.3
    assert fuselage.D_main == 1.0
    assert fuselage.lm == 5.0
    assert fuselage.ln == 3.0
    assert fuselage.h_battery == 0.3
    assert fuselage.l_tailcone == 2.5
    assert fuselage.h_tail == 1.0
    assert fuselage.V_battery == 0.5
    assert fuselage.l_empty == 0.5
    assert fuselage.w_aisle == 1.0
    assert fuselage.h_aisle == 1.5
    assert fuselage.w_motor == 0.8
    assert fuselage.l_motor == 2.0
    assert fuselage.l_empty == 0.5


def test_n_row():
    assert fuselage.n_row() == 76  # Expected number of rows


def test_l_nosecone():
    assert fuselage.l_nosecone() == 5.15  # Expected length of nose cone


def test_l_cabin():
    assert fuselage.l_cabin() == 82.84  # Expected length of cabin


def test_length_main_strut():
    assert isclose(fuselage.length_main_strut(0.1), 2.15, rel_tol=1e-2)  # Expected length of main strut


def test_l_battery():
    assert isclose(fuselage.l_battery(0.1), 16.67, rel_tol=1e-2)  # Expected length of battery


def test_w_battery():
    assert isclose(fuselage.w_battery(0.1), 16.67, rel_tol=1e-2)  # Expected width of battery


def test_l_end_nose_land():
    assert isclose(fuselage.l_end_nose_land(), 5.75, rel_tol=1e-2)  # Expected length of nose landing gear


def test_l_end_main_land():
    assert isclose(fuselage.l_end_main_land(0.1), 6 - 2.15 - 0.5, rel_tol=1e-2)  # Expected length of main landing gear


def test_check_back():
    assert fuselage.check_back(0.1)


def test_battery_dim():
    l_battery, w_battery, s_gear = fuselage.battery_dim(0.1)

    # get the expected values for the dimensions

    expected_s_gear = 0.1 + 0.01 * 13
    expected_w_battery = expected_s_gear * 2 + 0.5

    assert isclose(l_battery, fuselage.l_battery(expected_w_battery), rel_tol=1e-2)  # Expected length of battery
    assert isclose(w_battery, expected_w_battery, rel_tol=1e-2)  # Expected width of battery
    assert s_gear == expected_s_gear  # Expected gear separation


def test_top_width():
    assert isclose(fuselage.top_width(), (0.56 + 0.03 + 0.04) * 2 + 1, rel_tol=1e-2)  # Expected top width


def test_bigger_mag():
    assert fuselage.bigger_mag(0.1, 0.2) == 0.2  # Expected bigger magnitude


def test_bottom_width():
    assert isclose(fuselage.bottom_width(0.1), 2.56, rel_tol=1e-2)  # Expected bottom width


def test_height():
    assert isclose(fuselage.height(), 1.4 + 2 * 0.04 + 0.1 + 0.05 + 0.03 + 0.4, rel_tol=1e-2)


def test_y_floorheight():
    assert isclose(fuselage.y_floorheight(), 0.05 + 0.4 + 0.04, rel_tol=1e-2)


def test_length_fus():
    assert isclose(fuselage.length_fus(), 5.15 + 82.84 + 2.5 + 1.42, rel_tol=1e-2)


def test_fuselage_wetted():
    assert isclose(fuselage.fuselage_wetted(0.1), 832, rel_tol=1e-2)


def test_maximum_perimeter():
    assert isclose(fuselage.maximum_perimeter(0.1), 2.56 + (0.56 + 0.03 + 0.04) * 2 + 1 + 4.13, rel_tol=1e-2)


def test_above_position():
    expected_dict = {
        "frontwall": (0, 0.04),
        "empty_space": (0.04, 0.54),
        "motor": (0.54, 2.54),
        "engine": (2.54, 5.54),
        "firewall": (5.54, 5.69),
        "cockpit": (5.69, 7.11),
        "row1": (7.11, 8.2),
        "row2": (8.2, 9.29),
        "row3": (9.2, 10.38),
        "tailcone": (10.38, 12.88),
        "backwall": (12.88, 12.92),
    }

    actual_dict = fuselage.above_position()

    assert actual_dict["frontwall"][0] == expected_dict["frontwall"][0]
    assert actual_dict["backwall"][1] == expected_dict["backwall"][1]


def test_below_position():
    expected_dict = {
        "frontwall": (0, 0.04),
        "nose landing gear": (4.0, 5.65),
        "main landing gear": (3.34, 6.0),
        "battery": (31.3, 33.04),
    }

    actual_dict = fuselage.below_position(0.1)

    assert actual_dict["frontwall"][0] == expected_dict["frontwall"][0]
    assert actual_dict["nose landing gear"][0] == expected_dict["nose landing gear"][0]


def test_calculate_battery_split():
    xcg1, l1, xcg2, l2 = fuselage.calculate_battery_split(0.6)

    assert isclose(xcg1, 4.5, rel_tol=1e-1)
    assert isclose(l1, -1.39, rel_tol=1e-1)
    assert isclose(xcg2, 14.10, rel_tol=1e-1)
    assert isclose(l2, 2.15, rel_tol=1e-1)
