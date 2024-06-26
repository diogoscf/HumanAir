import sys
import os
import numpy as np
import pandas as pd
from io import StringIO

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.StructuralAnalysis.WingStructure import WingStructure, find_nearest, interp

mock_data = {
    "Aero": {"S_Wing": 34.56, "Taper_Wing": 0.4, "b_Wing": 19.93, "c_root_wing": 2.5},
    "Geometry": {
        "t_spar_tip": 0.010,
        "t_spar_root": 0.025,
        "t_skin_wing": 0.007,
        "spar_pos": [0.15, 0.5],
        "wing_stringer_area_m2": 0.02,
        "wing_stringer_number": [5],
        "wing_stringer_sections": [1],
        "wingbox_material": "Aluminium",
        "fus_width_m": 5.5,
    },
    "Materials": {"Aluminium": {"rho": 2700, "E": 70000, "G": 27000}},
    "Structures": {"structural_wing_weight": 1500},
}

simplified_airfoil_data = """\
1.000000  0.000000
0.750000  0.073047
0.500000  0.115624
0.250000  0.114919
0.000000  0.002411
0.250000 -0.020431
0.500000  0.003409
0.750000  0.026291
1.000000  0.000000
"""


# Function to read airfoil data from the multi-line string
def read_airfoil_data():
    return pd.read_csv(StringIO(simplified_airfoil_data), sep="\\s+", header=None, names=["x", "y"])


# Use the read_airfoil_data function to create airfoil_data_df
airfoil_data_df = read_airfoil_data()


def test_calc_spar_dist():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    spar_dist = ws.calc_spar_dist()

    assert len(spar_dist) == ws.nodes
    assert np.all(spar_dist >= 0)


def test_calc_hmax_dist():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    hmax_dist = ws.calc_hmax_dist()
    assert len(hmax_dist) == ws.nodes
    assert np.all(hmax_dist >= 0)


def test_calc_stringer_dist():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    stringer_dist = ws.calc_stringer_dist()
    assert len(stringer_dist) == ws.nodes
    assert np.all(stringer_dist >= 0)


def test_calc_airfoil_division():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    df_1, df_2 = ws.calc_airfoil_division()
    assert len(df_1) + len(df_2) == len(airfoil_data_df["x"])


def test_calculate_chord_distribution():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    chord_dist = ws.calculate_chord_distribution()
    assert len(chord_dist) == ws.nodes
    assert np.all(chord_dist >= 0)


def test_y_spar():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    airfoil_division = ws.calc_airfoil_division()
    df = pd.DataFrame(airfoil_division[0], columns=["x", "y"])
    y_spar = ws.y_spar(df.values)
    assert len(y_spar) > 0
    assert np.all(np.array(y_spar) >= 0)


def test_y_updown():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    y_updown = ws.y_updown()
    assert len(y_updown) == 2
    assert np.any(np.array(y_updown) >= 0)


def test_diff():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    diff = ws.diff([[1, 1], [2, 2], [3, 3]])
    assert len(diff) == 3


def test_y():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    y = ws.y([[1, 1], [2, 2], [3, 3]])
    assert len(y) == 3


def test_d_s1s2():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    h_s1s2 = ws.h_s1s2()
    assert len(h_s1s2) == 2
    assert np.all(np.array(h_s1s2) >= 0)


def test_h_s1s2():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    h_s1s2 = ws.h_s1s2()
    assert len(h_s1s2) == 2
    assert np.all(np.array(h_s1s2) >= 0)


def test_centroid():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    centroid = ws.centroid()
    assert len(centroid) == ws.nodes
    assert np.all(np.array(centroid) >= 0)


def test_ixx():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    ixx = ws.Ixx()
    assert ixx.all() >= 0


def test_shear_centre():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    shear_centre = ws.shear_centre()
    assert len(shear_centre) == ws.nodes
    assert np.all(np.array(shear_centre) >= 0)


def test_torsional_constant():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    torsional_constant = ws.torsional_constant()
    assert torsional_constant.all() >= 0


def test_weight_dist():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    weight_dist = ws.weight_dist()
    assert len(weight_dist) == ws.nodes
    assert np.all(weight_dist >= 0)


def test_calc_enclosed_area_dist():
    ws = WingStructure(mock_data, airfoil_data_df, nodes=501)
    enclosed_area_dist = ws.calc_enclosed_area_dist()
    assert len(enclosed_area_dist) == ws.nodes
    assert np.all(enclosed_area_dist >= 0)


def test_find_nearest():
    mock_array = [1, 2, 3]
    value = 1.5

    assert find_nearest(mock_array, value) == 0


def test_interp():
    x = 1.5
    x_given = [1, 2, 3]
    y_given = [10, 20, 30]
    assert interp(x, x_given, y_given) == 10
