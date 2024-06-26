import os
import sys
import numpy as np
import pandas as pd
from io import StringIO
from math import isclose

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from HumanAir.StructuralAnalysis.LoadDistributions import (
    interpolate_Cl_Cd_Cm,
    get_deflection,
    get_twist,
    force_distribution,
    moment_distribution,
    weight_distribution,
    axial_distribution_ground,
    axial_distribution_flight,
    read_points_from_load_dist,
    strut_error_calculation,
    strut_force,
    InternalLoads,
)
from HumanAir.isa import isa
from HumanAir.StructuralAnalysis.WingStructure import WingStructure


def test_get_deflection():
    MOI = np.array([1, 1, 1, 1, 1])
    y = np.array([0, 1, 2, 3, 4])
    M = np.array([0, 0, 0, 0, 0])
    E = 210e9
    w_fuselage = 1
    deflection = get_deflection(MOI, y, M, E, w_fuselage)
    assert np.allclose(deflection, np.zeros_like(y))


def test_get_twist():
    J = np.array([1, 1, 1, 1, 1])
    y = np.array([0, 1, 2, 3, 4])
    T = np.array([0, 0, 0, 0, 0])
    G = 80e9
    twist = get_twist(J, y, T, G)
    assert np.allclose(twist, np.zeros_like(y))


def test_force_distribution():
    AoA = 0
    altitude = 10000
    V = 100
    chord_dist = np.array([1, 1, 1, 1, 1])
    Cl_DATA = {0: {"y_span": [0, 1, 2, 3, 4], "coefficient": [1, 1, 1, 1, 1]}}
    Cdi_DATA = {0: {"y_span": [0, 1, 2, 3, 4], "coefficient": [0.01, 0.01, 0.01, 0.01, 0.01]}}
    L, D = force_distribution(AoA, altitude, V, chord_dist, Cl_DATA, Cdi_DATA)
    assert np.allclose(L, np.ones_like(chord_dist) * 0.5 * isa(altitude)[2] * V**2)
    assert np.allclose(D, np.ones_like(chord_dist) * 0.005 * isa(altitude)[2] * V**2)


def test_moment_distribution():
    AoA = 0
    altitude = 10000
    V = 100
    chord_dist = np.array([1, 1, 1, 1, 1])
    Cm_DATA = {0: {"y_span": [0, 1, 2, 3, 4], "coefficient": [0.1, 0.1, 0.1, 0.1, 0.1]}}
    ac_data = {"Aero": {"MAC_wing": 2}}
    M = moment_distribution(AoA, altitude, V, chord_dist, Cm_DATA, ac_data)
    assert np.allclose(M, np.ones_like(chord_dist) * 0.1 * 0.5 * isa(altitude)[2] * V**2 * 2)


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


def test_weight_distribution():
    chord_dist = np.array([1, 1, 1, 1, 1])
    wing_structure = WingStructure(mock_data, airfoil_data_df, nodes=5)
    ac_data = {"CL2Weight": {"Wfuel_N": 10}, "Geometry": {"strut_loc_b/2": 0.5}, "Aero": {"b_Wing": 10}}

    W, W_fuel, idxs, c_between_struts = weight_distribution(chord_dist, wing_structure, ac_data)

    assert idxs == (1, 4)
    assert np.allclose(c_between_struts, np.ones(3))


def test_axial_distribution_ground():
    nodes = 5
    Vy_strut = 100
    ac_data = {"Geometry": {"strut_loc_b/2": 0.5}}
    A_ground = axial_distribution_ground(nodes, Vy_strut, ac_data)
    assert np.allclose(A_ground, np.array([100, 100, 0]))


def test_axial_distribution_flight():
    nodes = 5
    Vy_strut = 100
    ac_data = {"Geometry": {"strut_loc_b/2": 0.5}}
    A_cruise = axial_distribution_flight(nodes, Vy_strut, ac_data)
    assert np.allclose(A_cruise, np.array([-100, -100, 0]))


def test_read_points_from_load_dist():
    L_cruise = np.array([1, 2, 3, 4, 5])
    W_cruise = np.array([5, 4, 3, 2, 1])
    W_fuel = np.array([1, 1, 1, 1, 1])
    idxs = (1, 3)
    a, b, c, d, e, f = read_points_from_load_dist(L_cruise, W_cruise, W_fuel, idxs)
    assert a == 5
    assert b == 1
    assert c == 1
    assert d == 2
    assert e == 3
    assert f == 5


def test_strut_error_calculation():
    P = 1
    Vz = np.array([1, 2, 3, 4, 5])
    Vz = Vz.astype(np.float64)
    n_before_strut = 2
    theta_strut = 0.5
    y_points_halfspan = np.array([0, 1, 2, 3, 4])
    MOI = np.array([1, 1, 1, 1, 1])
    E = 1000
    l_strut = 1
    A_strut = 1
    w_fuselage = 1000
    error = strut_error_calculation(
        P, Vz, n_before_strut, theta_strut, y_points_halfspan, MOI, E, l_strut, A_strut, w_fuselage
    )
    assert isclose(error, -0.1, rel_tol=1)


def test_strut_force():
    MOI = np.array([1, 1, 1, 1, 1])
    y_points_halfspan = np.array([0, 1, 2, 3, 4])
    Vz = np.array([1, 2, 3, 4, 5])
    Vz = Vz.astype(np.float64)
    w_fuselage = 1
    ac_data = {
        "Geometry": {
            "strut_loc_b/2": 0.5,
            "fus_height_m": 1,
            "strut_section_area_m2": 1,
            "wingbox_material": "aluminum",
        },
        "Aero": {"b_Wing": 10},
        "Materials": {"aluminum": {"E": 10}},
    }

    Vz_strut, Vy_strut, V_strut = strut_force(
        MOI, y_points_halfspan, Vz, w_fuselage, max_iter=0, tol=1e-6, ac_data=ac_data
    )

    assert np.isclose(Vz_strut, 0)
    assert np.isclose(Vy_strut, 0)
    assert np.isclose(V_strut, 0)


def test_internal_loads():
    L = np.array([1, 1, 1, 1, 1])
    D = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
    M = np.array([0.01, 0.01, 0.01, 0.01, 0.01])
    wing_structure = WingStructure(mock_data, airfoil_data_df, nodes=5)

    ac_data = {
        "Aero": {"b_Wing": 10, "QuarterChordSweep_Wing_deg": 0},
        "CL2Weight": {"Wfuel_N": 100},
        "Geometry": {
            "strut_loc_b/2": 0.5,
            "fus_height_m": 1,
            "strut_section_area_m2": 1,
            "wingbox_material": "aluminum",
        },
        "Materials": {"aluminum": {"E": 210e9}},
    }

    Vx, _, _, _, _, _ = InternalLoads(L, D, M, wing_structure, ac_data)

    assert np.allclose(Vx, np.array([-1, -1 / 2, 0]), rtol=1e-1)


def test_interpolate_cl_cd_cm():
    Cl_data = {0: {"y_span": [0.1, 0.2, 0.3, 0.4, 0.5], "coefficient": [0.11, 0.21, 0.31, 0.41, 0.51]}}

    Cdi_data = {0: {"y_span": [0.01, 0.02, 0.03, 0.04, 0.05], "coefficient": [0.011, 0.021, 0.031, 0.041, 0.051]}}

    Cm_data = {
        0: {"y_span": [0.001, 0.002, 0.003, 0.004, 0.005], "coefficient": [0.0011, 0.0021, 0.0031, 0.0041, 0.0051]}
    }

    interpolated_values = interpolate_Cl_Cd_Cm(Cl_data, Cdi_data, Cm_data, 1)
    assert interpolated_values == (
        {0: {"coefficient": 0.51, "y_span": 1}},
        {0: {"coefficient": 0.051, "y_span": 1}},
        {0: {"coefficient": 0.0051, "y_span": 1}},
    )
