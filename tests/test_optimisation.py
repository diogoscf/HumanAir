import os
import sys
import numpy as np
import pandas as pd
from io import StringIO
from math import isclose

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from HumanAir.StructuralAnalysis.optimisation import (
    get_stringers_at_nodes,
    full_skin_weight,
    get_weight,
    get_axial_forces,
    get_shear_stress,
    get_max_axial_stress,
    get_I,
    get_A,
    get_Q,
    compare_approximate,
    get_force_distributions,
    # stiffened_skin_buckling_constraint,
    get_critical_skin_buckling_stress,
    # shear_buckling_constraint,
    # tensile_failure_constraint,
    # stringer_buckling_constraint,
    # deflection_constraint,
    # constraint_data_from_variables,
    # get_loads_from_acdata,
    unstack_variables,
    stack_variables,
)


def test_get_stringers_at_nodes():
    stringer_sections = [1, 2, 3]
    no_stringers = [1, 1, 1]
    nodes_halfspan = 4

    stringers_at_nodes = get_stringers_at_nodes(stringer_sections, no_stringers, nodes_halfspan)

    assert np.allclose(stringers_at_nodes, np.ones(nodes_halfspan))


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


def test_full_skin_weight():
    # Define test parameters
    t_skin = 0.1
    rho_material = 2700.0
    chord_dist = np.array([1, 1, 1, 1, 1])
    rho = 1000
    non_fus_idx = 4

    # Calculate full skin weight
    weight = full_skin_weight(t_skin, chord_dist, rho, non_fus_idx, airfoil_shape=airfoil_data_df)

    # Define expected weight
    expected_weight = t_skin * rho_material * 0.76

    assert np.isclose(weight, expected_weight, rtol=1)


def test_get_weight():
    # Define test parameters
    variables = np.array([1.0, 2.0, 3.0, 4.0])
    htot = 1.0
    rho_material = 2700.0
    len_nodes = 3
    A_stringer = 1
    stringer_sections_half_span = np.ones(1)
    n_half_span = 1
    chord_dist = np.array([1, 1, 1])
    non_fus_dist = 2

    # Calculate weight using get_weight function
    weight = get_weight(
        variables,
        htot,
        rho_material,
        len_nodes,
        A_stringer,
        stringer_sections_half_span,
        n_half_span,
        chord_dist,
        non_fus_dist,
    )

    # Define expected weight based on the formula used in get_weight
    expected_weight = variables[0] * variables[1] * variables[2] * variables[3] * rho_material * 2

    assert np.isclose(weight, expected_weight, rtol=1)


def test_get_axial_forces():
    Mx = np.array([1, 1, 1])
    Vy = np.array([0, 0, 0])
    MOI = 10
    h_max = 1
    stringers_at_nodes = np.array([1, 1, 1])
    A_stringer = 1

    P = get_axial_forces(Mx, Vy, MOI, h_max, stringers_at_nodes, A_stringer)

    assert np.allclose(P, np.array([1, 1, 1]) / 10)


def test_get_shear_stress():
    Vz = np.array([1, 1, 1])
    My = np.array([1, 1, 1])
    Q = np.array([1, 1, 1])
    MOI = 10
    enclosed_area = 1
    t_spar = 1

    shear_stress = get_shear_stress(Vz, My, Q, MOI, enclosed_area, t_spar)

    assert np.allclose(shear_stress, np.ones(len(My)) * (0.1 + 0.5), rtol=1e-2)


def test_get_max_axial_stress():
    Mx = np.array([1, 1, 1])
    Vy = np.array([1, 1, 1])
    MOI = 10
    hmax = 1
    area = 10

    axial_stress = get_max_axial_stress(Mx, Vy, MOI, hmax, area)

    assert np.allclose(axial_stress, np.ones(len(Mx)) * (0.1 + 0.1), rtol=1e-2)


def test_get_I():
    t_spar = 1
    t_skin = 1
    no_stringers = 1
    A_stringers = 10
    w_top = 1
    w_bottom = 1
    h_avemax = 1
    h_frontspar = 10
    h_rearspar = 1

    MOI = get_I(t_spar, t_skin, no_stringers, A_stringers, w_top, w_bottom, h_avemax, h_frontspar, h_rearspar)

    assert isclose(MOI, 2 + 83.417 + 10, rel_tol=1e-2)


def test_get_A():
    t_skin = 1
    t_spar = 1
    no_stringers = 1
    A_stringer = 10
    spar_pos = [1, 1]
    h_frontspar = 10
    h_rearspar = 1
    chord_dist = 1

    A = get_A(t_skin, t_spar, no_stringers, A_stringer, spar_pos, h_frontspar, h_rearspar, chord_dist)

    assert isclose(A, 10 + 10 + 1 + 0, rel_tol=1e-2)


def test_get_Q():
    t_spar = 10
    h_frontspar = 10
    h_rearspar = 1
    Q = get_Q(t_spar, h_frontspar, h_rearspar)

    assert isclose(Q, 100 / 8 * 10 + 1 / 8 * 10, rel_tol=1e-2)


def test_compare_approximate():
    first = {"a": np.array([1.0, 2.0, 3.0]), "b": np.array([4.0, 5.0, 6.0])}
    second = {"a": np.array([1.0, 2.0, 3.0]), "b": np.array([4.0, 5.0, 6.0])}

    assert compare_approximate(first, second)


def test_get_force_distributions():
    AoA = 0.0
    altitude = 10000.0
    Vc = 250.0

    Cl_data = {0: {"y_span": [0.1, 0.1, 0.1], "coefficient": [0.1, 0.1, 0.1]}}
    Cdi_data = {0: {"y_span": [0.1, 0.1, 0.1], "coefficient": [0.1, 0.1, 0.1]}}
    Cm_data = {0: {"y_span": [0.1, 0.1, 0.1], "coefficient": [0.1, 0.1, 0.1]}}

    chord_dist = np.array([1.0, 1.0, 1.0])

    L_cruise, _, _ = get_force_distributions(AoA, altitude, Vc, Cl_data, Cdi_data, Cm_data, chord_dist)

    assert np.allclose(L_cruise, np.ones(3) * 1290, rtol=1)


def test_get_critical_skin_buckling_stress():
    # Example input data for testing
    spar_pos = np.array([0, 1])
    stringers_at_nodes = 5
    A_stringer = 1
    E = 10
    nu = 0.5
    sigma_yield = np.array([100])
    t_skin = 1
    t_stiffener = 1
    b_stiffener = 1

    # Call the function
    critical_stress = get_critical_skin_buckling_stress(
        spar_pos, stringers_at_nodes, A_stringer, E, nu, sigma_yield, t_skin, t_stiffener, b_stiffener
    )

    # Assertions
    assert isinstance(critical_stress, np.ndarray)  # Ensure the output is a NumPy array
    assert critical_stress.shape == (len(spar_pos) - 1,)  # Ensure the shape matches the number of spar positions
    assert np.all(critical_stress >= 0)  # Ensure all critical stresses are non-negative


# def test_stiffened_skin_buckling_constraint():
#     variables = np.array([1.0, 1.0, 1.0, 1.0])
#     spar_pos = [1.0, 1.0]
#     A_stringer = 10.0
#     hmax = 1.0
#     t_stiffener = 1.0
#     b_stiffener = 1.0
#     E = 1.0
#     nu = 0.3
#     sigma_yield = 1.0
#     MOI_args = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#
#     constraint_value = stiffened_skin_buckling_constraint(
#         variables, spar_pos, A_stringer, hmax, t_stiffener, b_stiffener, E, nu, sigma_yield, MOI_args
#     )
#
#     assert constraint_value >= 0


# def test_shear_buckling_constraint():
#     assert False
#
#
# def test_tensile_failure_constraint():
#     assert False
#
#
# def test_stringer_buckling_constraint():
#     assert False
#
#
# def test_deflection_constraint():
#     assert False
#
#


def test_unstack_variables():
    # Define test variables
    variables = np.array([1.0, 2.0, 3.0, 4.0])
    n_skin, n_spar, n_rib = 2, 3, 4

    # Call unstack_variables function
    t_tip, t_root, t_skin, no_stringers = unstack_variables(variables)

    expected_rib = variables[n_skin + n_spar : n_skin + n_spar + n_rib]

    # Assert equality
    assert np.allclose(t_skin, variables[2] / 1000)
    assert np.allclose(t_tip, variables[0] / 1000)
    assert np.allclose(t_root, expected_rib)


def test_stack_variables():
    # Define test variables
    t_tip = 1
    t_root = 1
    t_skin = 1
    no_stringers = np.array([1, 1, 1, 1])

    # Call stack_variables function
    variables = stack_variables(t_tip, t_root, t_skin, no_stringers)

    # Assert equality
    assert np.allclose(len(variables), 7)
