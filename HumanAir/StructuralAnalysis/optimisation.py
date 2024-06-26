import numpy as np
from scipy.optimize import minimize  # type: ignore[import-untyped]

# from scipy.integrate import cumulative_trapezoid  # type: ignore[import-untyped]
import sys
import os
import matplotlib.pyplot as plt  # noqa: F401
import time
import copy

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.StructuralAnalysis.WingStructure import WingStructure
from HumanAir.StructuralAnalysis.LoadDistributions import (
    InternalLoads,
    force_distribution,
    moment_distribution,
    # weight_distribution,
    interpolate_Cl_Cd_Cm,
    get_deflection,
)
from HumanAir.aircraft_data import (
    aircraft_data,
    airfoil_shape,
    Cl_data_wing,
    Cdi_data_wing,
    Cm_data_wing,
    save_ac_data_to_json,
)

# Caching for SPEEEEEEEEEEEEEEEEEEEED
# force_distributions = None
# force_dist_params = None # Check whether parameters to get the distribution have changed


def get_stringers_at_nodes(stringer_sections, no_stringers, nodes_halfspan):
    if len(stringer_sections) != len(no_stringers):
        raise ValueError("The number of stringer sections should be equal to the number of stringers per section")

    stringers = np.ones(nodes_halfspan)
    stringer_section_ends = np.rint(np.cumsum(stringer_sections) * len(stringers)).astype(int)
    stringer_section_starts = np.insert(stringer_section_ends[:-1], 0, 0)
    for i in range(len(stringer_sections)):
        start = stringer_section_starts[i]
        end = stringer_section_ends[i]
        # print(start, end)
        stringers[start:end] = no_stringers[i]

    return stringers


def full_skin_weight(t_skin, chord_dist, rho, non_fus_idx, airfoil_shape=airfoil_shape):
    """
    Calculate the weight of the skin of the wing

    Parameters
    ----------
    t_skin : float
        Thickness of the skin [m]
    chord_dist : np.array
        Distribution of the chord length along the wing [m]
    rho : float
        Density of the material [kg/m^3]
    non_fus_idx : int
        Index of the first non-fuselage node
    airfoil_shape : pd.DataFrame
        DataFrame containing the coordinates of the airfoil shape

    Returns
    -------
    W_skin : np.array
        Weight of the skin at each node [N]
    """

    x_vals = airfoil_shape["x"].values
    y_vals = airfoil_shape["y"].values
    l_skin = np.sum(np.sqrt(np.diff(x_vals) ** 2 + np.diff(y_vals) ** 2))
    W_skin = rho * l_skin * t_skin * chord_dist[non_fus_idx:]
    return W_skin


def get_weight(
    variables,
    htot,
    rho,
    len_nodes,
    A_stringer,
    stringer_sections_halfspan,
    n_halfspan,
    chord_dist,
    non_fus_idx,
    print_skin=False,
):
    t_spar_tip, t_spar_root, t_skin, no_stringers = unstack_variables(variables)
    stringers_at_nodes = get_stringers_at_nodes(stringer_sections_halfspan, no_stringers, n_halfspan)

    t_spar = np.linspace(t_spar_root, t_spar_tip, n_halfspan)

    W_spar = np.sum(rho * t_spar * htot * len_nodes)
    # W_skin = rho * t_skin * Sw # Two sides, but half a wing
    W_skin = np.sum(full_skin_weight(t_skin, chord_dist, rho, non_fus_idx, airfoil_shape=airfoil_shape) * len_nodes)
    W_stringers = np.sum(rho * A_stringer * stringers_at_nodes * len_nodes)
    weight = W_skin + W_spar + W_stringers

    if print_skin:
        print(W_skin)
    return weight


n_iter = 0
time_at_last_iter = None


def objective(
    variables,
    htot,
    rho,
    len_nodes,
    A_stringer,
    stringer_sections_halfspan,
    n_halfspan,
    chord_dist,
    non_fus_idx,
    verbose=False,
):  # pragma: no cover
    weight = get_weight(
        variables, htot, rho, len_nodes, A_stringer, stringer_sections_halfspan, n_halfspan, chord_dist, non_fus_idx
    )
    obj_val = (weight - 40) / 100  # Lowest possible is around 60 kg
    # print(obj_val, flush=True)

    if verbose:
        global n_iter, time_at_last_iter
        n_iter += 1
        if time_at_last_iter is None:
            time_at_last_iter = time.time()

        if n_iter % 10 == 0:
            now = time.time()
            if time_at_last_iter is not None:
                print(n_iter, f"{weight:.2f}", f"{(now - time_at_last_iter):.2f}s (10 iter)", flush=True)
                time_at_last_iter = now

    return obj_val


def get_axial_forces(Mx, Vy, MOI, h_max, stringers_at_nodes, A_stringer):
    sigma_bending = Mx * h_max / MOI
    sigma_axial = Vy / (A_stringer * stringers_at_nodes)
    sigma = sigma_bending + sigma_axial

    P = sigma * A_stringer
    return P


def get_shear_stress(Vz, My, Q, MOI, enclosed_area, t_spar):
    shear_from_transverse = Vz * Q / (MOI * t_spar)
    shear_from_torsion = My / (2 * enclosed_area * t_spar)
    return shear_from_transverse + shear_from_torsion


def get_max_axial_stress(Mx, Vy, MOI, hmax, area):
    stress_bending = Mx * hmax / MOI
    stress_axial = Vy / area
    stress = np.abs(stress_axial) + np.abs(stress_bending)
    return stress


# def get_axial_stresses(Mx, MOI, h_max):
#     sigma = Mx * h_max / MOI
#     return sigma


def get_I(t_spar, t_skin, no_stringers, A_stringer, w_top, w_bottom, h_avemax, h_frontspar, h_rearspar):
    I_skin = t_skin * (w_top + w_bottom) * h_avemax**2
    I_spar = 1 / 12 * t_spar * (h_frontspar**3 + h_rearspar**3)
    I_stringers = A_stringer * no_stringers * h_avemax**2
    MOI = I_skin + I_spar + I_stringers
    return MOI


def get_A(t_skin, t_spar, no_stringers, A_stringer, spar_pos, h_frontspar, h_rearspar, chord_dist):
    A = (
        h_frontspar * t_spar
        + h_rearspar * t_spar
        + no_stringers * A_stringer
        + 2 * t_skin * (spar_pos[1] - spar_pos[0]) * chord_dist
    )
    return A


def get_Q(t_spar, h_frontspar, h_rearspar):
    Q = h_frontspar**2 / 8 * t_spar + h_rearspar**2 / 8 * t_spar
    return Q


lengths_vars = None
CONVERSION_FACTOR_THICKNESSES = 1000  # Convert to mm and back, for numerical reasons
CONVERSION_FACTOR_STRINGERS = 5  # Convert to something around order 1 for numerical reasons


def stack_variables(t_spar_tip, t_spar_root, t_skin, no_stringers):
    return np.hstack(
        (
            [
                t_spar_tip * CONVERSION_FACTOR_THICKNESSES,
                t_spar_root * CONVERSION_FACTOR_THICKNESSES,
                t_skin * CONVERSION_FACTOR_THICKNESSES,
            ],
            [x / CONVERSION_FACTOR_STRINGERS for x in no_stringers],
        )
    ).flatten()


def unstack_variables(variables):
    t_spar_tip, t_spar_root, t_skin = variables[:3].flatten()
    no_stringers = variables[3:].flatten()
    return (
        t_spar_tip / CONVERSION_FACTOR_THICKNESSES,
        t_spar_root / CONVERSION_FACTOR_THICKNESSES,
        t_skin / CONVERSION_FACTOR_THICKNESSES,
        no_stringers * CONVERSION_FACTOR_STRINGERS,
    )


def compare_approximate(first, second):
    """Return whether two dicts of arrays are roughly equal"""
    if not isinstance(first, type(second)):
        # print("different types", type(first), type(second), flush=True)
        return False

    if not isinstance(first, dict):
        if isinstance(first, np.ndarray):
            return np.allclose(first, second)

        return first == second

    if first.keys() != second.keys():
        # print("different keys", first.keys(), second.keys(), flush=True)
        return False

    for key in first:
        if not isinstance(first[key], type(second[key])):
            # print(key, first[key], second[key], flush=True)
            return False

        if isinstance(first[key], dict):
            if not compare_approximate(first[key], second[key]):
                return False
            continue

        if isinstance(first[key], (np.ndarray, float)):
            if not np.allclose(first[key], second[key]):
                # print(key, first[key], second[key], flush=True)
                return False
            continue

        # print(type(first[key]), type(second[key]), isinstance(first[key], np.ndarray), flush=True)
        if first[key] != second[key]:
            # print(key, first[key], second[key], flush=True)
            return False

    return True


# Caching for SPEEEEEEEEEEEEEEEEEEEED
CURRENT_AC_DATA = None
CURRENT_LOADS = None


def get_loads_from_acdata(ac_data, nlim, nodes, L_cruise, D_cruise, M_cruise):
    global CURRENT_AC_DATA, CURRENT_LOADS
    if compare_approximate(ac_data, CURRENT_AC_DATA):
        return CURRENT_LOADS

    wing_structure = WingStructure(ac_data, airfoil_shape, nodes)

    Vx, Vy, Vz, Mx, My, Mz = InternalLoads(
        L_cruise, D_cruise, M_cruise, wing_structure, ac_data=ac_data, load_factor=nlim
    )

    # print(nlim, np.max(np.abs(L_cruise)), np.max(np.abs(Vz)), np.max(np.abs(Mx)), flush=True)

    CURRENT_LOADS = Vx, Vy, Vz, Mx, My, Mz
    CURRENT_AC_DATA = copy.deepcopy(ac_data)

    return Vx, Vy, Vz, Mx, My, Mz


# def MOI_from_variables(*args):
#     return constraint_data_from_variables(*args)[0]


def constraint_data_from_variables(
    variables,
    A_stringer,
    w_top,
    w_bottom,
    h_avemax,
    h_frontspar,
    h_rearspar,
    n_halfspan,
    n_fullspan,
    stringer_sections_halfspan,
    spar_pos,
    chord_dist,
    ac_data,
    nlim,
    L_cruise,
    D_cruise,
    M_cruise,
):
    t_spar_tip, t_spar_root, t_skin, no_stringers = unstack_variables(variables)

    stringers_at_nodes = get_stringers_at_nodes(stringer_sections_halfspan, no_stringers, n_halfspan)

    t_spar = np.linspace(t_spar_root, t_spar_tip, n_halfspan)

    MOI = get_I(t_spar, t_skin, stringers_at_nodes, A_stringer, w_top, w_bottom, h_avemax, h_frontspar, h_rearspar)

    Q = get_Q(t_spar, h_frontspar, h_rearspar)
    area = get_A(t_skin, t_spar, stringers_at_nodes, A_stringer, spar_pos, h_frontspar, h_rearspar, chord_dist)

    enclosed_area = (h_frontspar + h_rearspar) / 2 * (spar_pos[1] - spar_pos[0]) * chord_dist

    ac_data["Geometry"]["t_spar_root"] = t_spar_root
    ac_data["Geometry"]["t_spar_tip"] = t_spar_tip
    ac_data["Geometry"]["t_skin_wing"] = t_skin
    ac_data["Geometry"]["wing_stringer_number"] = no_stringers

    loads = get_loads_from_acdata(ac_data, nlim, n_fullspan, L_cruise, D_cruise, M_cruise)

    return MOI, stringers_at_nodes, t_spar, t_skin, Q, area, enclosed_area, loads


def deflection_constraint(variables, y, max_deflect, E, w_fuselage, MOI_args):
    MOI, *_, loads = constraint_data_from_variables(variables, *MOI_args)
    Mx = loads[3]
    deflection = get_deflection(MOI, y, Mx, E, w_fuselage)
    # print(np.max(deflection), flush=True)
    return (max_deflect - np.max(deflection)) / max_deflect


def stringer_buckling_constraint(variables, A_stringer, hmax, I_stringer, L_eff, E, MOI_args):
    P_cr = (np.pi) ** 2 * E * I_stringer / L_eff**2  # critical force at which stringer will buckle
    MOI, stringers, *_, loads = constraint_data_from_variables(variables, *MOI_args)[:2]
    Vy, Mx = loads[1], loads[3]

    # compressive stress
    P = get_axial_forces(Mx, Vy, MOI, hmax, stringers, A_stringer)
    # print(P_cr - np.max(P), flush=True)
    return (P_cr - np.max(P)) / P_cr


def tensile_failure_constraint(variables, hmax, sigma_yield, MOI_args):
    MOI, *_, area, _, loads = constraint_data_from_variables(variables, *MOI_args)
    Vy, Mx = loads[1], loads[3]
    sigma = get_max_axial_stress(Mx, Vy, MOI, hmax, area)
    # print((sigma_yield - np.max(sigma)) / sigma_yield, variables, flush=True)
    return (sigma_yield - np.max(sigma)) / sigma_yield


def shear_buckling_constraint(variables, E, nu, h_rearspar, MOI_args, ks=9.5):
    # inertia
    MOI, _, t_spar, _, Q, _, enclosed_area, loads = constraint_data_from_variables(variables, *MOI_args)
    Vz, My = loads[2], loads[4]

    # critical shear stress
    # critical height between the two spars is shorter spar
    shear_critical = (np.pi) ** 2 * ks * E / (12 * (1 - nu) ** 2) * (t_spar / h_rearspar) ** 2

    # shear
    shear = np.abs(get_shear_stress(Vz, My, Q, MOI, enclosed_area, t_spar))

    idx = np.argmin((shear_critical - shear))
    # print(np.pi, ks, E, nu, t_spar[idx], h_rearspar[idx], shear_critical[idx])

    # print(shear_critical[idx]/1e6, shear[idx]/1e6, ((shear_critical - shear) / shear_critical)[idx], flush=True)
    return ((shear_critical - shear) / shear_critical)[idx]


def get_critical_skin_buckling_stress(
    spar_pos,
    stringers_at_nodes,
    A_stringer,
    E,
    nu,
    sigma_yield,
    t_skin,
    t_stiffener,
    b_stiffener,
    C_skin=4,
    C_stiffener=0.425,
    C_we=4,
):
    b = 2 * (spar_pos[1] - spar_pos[0]) / stringers_at_nodes

    # crippling stress of the stiffener
    alpha = 0.8  # empirical values
    n = 0.6  # empirical values
    stress_crip_ratio = alpha * (
        C_stiffener / sigma_yield * np.pi**2 * E / (12 * (1 - nu**2)) * (t_stiffener / b_stiffener) ** 2
    ) ** (1 - n)

    if stress_crip_ratio < 1:  # TODO: stress_crippling underfined if stress_crip_ratio >= 1
        stress_crippling = sigma_yield * stress_crip_ratio  # check it cripples before it yields
    else:
        stress_crippling = sigma_yield

    # calculate effective width of the stiffener (accounting for attachement between the stringer and the skin)
    we = t_skin / 2 * np.sqrt(C_we * np.pi**2 / (12 * (1 - nu**2))) * np.sqrt(E / stress_crippling)

    # Calculate the initial skin buckling taking into account the 2we
    stress_critical_skin = C_skin * np.pi**2 * E / (12 * (1 - nu**2)) * (t_skin / (b - 2 * we)) ** 2

    # final critical stress of the entire panel
    stress_critical = (
        stress_critical_skin * t_skin * (b - 2 * we) + stress_crippling * (A_stringer + 2 * we * t_skin)
    ) / (t_skin * (b - 2 * we) + (A_stringer + 2 * we * t_skin))

    stress_critical[(stress_critical_skin >= stress_crippling) | ((b - 2 * we) <= 0)] = stress_crippling
    # print(np.min(stress_critical)/1e6)
    # if np.min(stress_critical / 1e6) < 60:
    #     idx = np.argmin(stress_critical)
    #     # print(
    #     #     "LOW VALUE:",
    #     #     np.min(stress_critical / 1e6),
    #     #     stress_critical_skin[idx] / 1e6,
    #     #     stress_crippling / 1e6,
    #     #     we,
    #     #     b[idx],
    #     #     (b - 2 * we)[idx],
    #     #     (A_stringer + 2 * we * t_skin),
    #     #     flush=True,
    #     # )

    return stress_critical


def stiffened_skin_buckling_constraint(
    variables,
    spar_pos,
    A_stringer,
    hmax,
    t_stiffener,
    b_stiffener,
    E,
    nu,
    sigma_yield,
    MOI_args,
    C_skin=4,
    C_stiffener=0.425,
    C_we=4,
):
    MOI, stringers_at_nodes, _, t_skin, _, area, _, loads = constraint_data_from_variables(variables, *MOI_args)
    stress_critical = get_critical_skin_buckling_stress(
        spar_pos,
        stringers_at_nodes,
        A_stringer,
        E,
        nu,
        sigma_yield,
        t_skin,
        t_stiffener,
        b_stiffener,
        C_skin=C_skin,
        C_stiffener=C_stiffener,
        C_we=C_we,
    )

    Vy, Mx = loads[1], loads[3]

    stress_applied = get_max_axial_stress(Mx, Vy, MOI, hmax, area)

    # print(np.max(stress_applied), stress_critical[(np.argmax(np.abs(stress_critical)))], flush=True)

    idx = np.argmin(stress_critical - stress_applied)

    return ((stress_critical - stress_applied) / stress_critical)[idx]


def get_force_distributions(AoA, altitude, Vc, Cl_data, Cdi_data, Cm_data, chord_dist, ac_data=aircraft_data):
    L_cruise, D_cruise = force_distribution(AoA, altitude, Vc, chord_dist, Cl_DATA=Cl_data, Cdi_DATA=Cdi_data)

    M_cruise = moment_distribution(AoA, altitude, Vc, chord_dist, Cm_data, ac_data=ac_data)

    # W, W_fuel, idxs, _ = weight_distribution(chord_dist, wing_structure, ac_data=ac_data)

    return L_cruise, D_cruise, M_cruise  # , W, W_fuel, idxs


last_time = time.time()
start = time.time()
times_called = 0


def print_time():  # pragma: no cover
    global last_time, times_called
    now = time.time()
    print(f"{times_called} - Time taken:", now - last_time, flush=True)
    last_time = now
    times_called += 1


# TODO: Make the AoA taken from aircraft_data, need to implement it in the json
def run_optimiser(
    ac_data=aircraft_data,
    AoA=0.0,
    altitude=None,
    nlim=None,
    max_deflect=1.7,
    nodes=401,
    maxiter=None,
    full_return=False,
    verbose=False,
):  # pragma: no cover
    """
    Function to run the optimization of the wing structure

    Parameters
    ----------
    ac_data : dict
        Dictionary containing the aircraft data
    AoA : float, default -5.0
        Angle of attack of the aircraft [deg]
    altitude : float or None, default None
        Altitude of the aircraft [m].
        If None, the altitude is taken from the aircraft data (cruise altitude).
    nlim : float or None, default None
        Load factor of the aircraft
        If None, the load factor is taken from the aircraft data (nult/1.5).
    max_deflect : float, default 1.7
        Maximum deflection of the wing [m]
    nodes : int, default 401
        Number of nodes to discretize the wing structure
    maxiter : int or None, default None
        Maximum number of iterations for the optimization.
        If None, the default value of the optimizer is used.
    full_return : bool, default False
        If True, the function returns additional information about the optimization process

    Returns
    -------
    optimized_t_spar_tip, optimized_t_spar_root, optimized_t_skin : np.array
        Optimized thickness of the spar (tip and root) and skin
    optimized_no_stringers : np.array
        Optimized number of stringers
    weight : float
        Minimum weight of the wing structure (i.e. for the optimised structure)
    """

    altitude = ac_data["Performance"]["Altitude_Cruise_m"] if altitude is None else altitude
    nlim = ac_data["Performance"]["n_max"] if nlim is None else nlim

    extra = 1 if nodes % 2 == 1 else 0
    n_halfspan = (nodes // 2) + extra

    # Initialize Wing Structure Class
    wing_structure = WingStructure(ac_data, airfoil_shape, nodes)
    h_mid, h_s1s2 = wing_structure.h_s1s2()
    l_box_up, l_box_down = wing_structure.d_s1s2()

    h_mid, h_s1s2 = h_mid[(h_mid.shape[0] // 2) :], h_s1s2[(h_s1s2.shape[0] // 2) :]
    l_box_up, l_box_down = l_box_up[(l_box_up.shape[0] // 2) :], l_box_down[(l_box_down.shape[0] // 2) :]
    chord_dist = wing_structure.chord_distribution
    y_points = wing_structure.ypts
    chord_dist_halfspan, y_points_halfspan = (
        chord_dist[(chord_dist.shape[0] // 2) :],
        y_points[(y_points.shape[0] // 2) :],
    )

    hmax = wing_structure.hmax_dist.flatten()  # maximum height of the wingbox, as a function of y
    hmax = hmax[(hmax.shape[0] // 2) :]

    # Parameters
    # rho = 2710  # density of aluminium (kg/m^3)

    h_frontspar = h_s1s2[:, 0]  # height of the spar at 15% of the chord, as a function of y
    h_rearspar = h_s1s2[:, 1]  # height of the spar at 50% of the chord, as a function of y
    htot = (h_frontspar + h_rearspar).flatten()  # total height, just for calculation ease, as a function of y
    h_avemax = htot / 4  # averaged "max" height from the central line, as a function of y
    w_top = l_box_up.flatten()  # width of top "straight" skin, as a function of y
    w_bottom = l_box_down.flatten()  # width of bottom "straight" skin, as a function of y
    # wtot = w_bottom + w_top  # total width, just for calculation ease

    w_fuselage = ac_data["Geometry"]["fus_width_m"] / 2  # width of the fuselage (we are interested in half)
    non_fus_idx = np.argmin(np.abs(y_points_halfspan - w_fuselage))  # index of the first non-fuselage node

    # E = 68e9  # Young's Modulus for aluminium (Pa)
    # sigma_yield = 40e6  # Yield strength for aluminium (Pa)

    # stringer_sections_halfspan = [0.4, 0.3, 0.3]  # Should add up to 1
    # no_string = [50, 30, 20]

    # Initial guess for thickness of spar
    # t_spar = wing_structure.t_spar_dist.flatten()
    t_spar_tip0, t_spar_root0 = wing_structure.t1_spar, wing_structure.t2_spar
    t_skin0 = wing_structure.t_skin
    stringer_sections_halfspan = wing_structure.stringer_sections
    no_string = wing_structure.stringer_number
    # t_spar0 = t_spar0[t_spar0.shape[0] // 2 :]
    # Constant skin thickness [m]
    # t_skin0 = wing_structure.t_skin * np.ones(t_spar0.shape)
    # no_stringers0 = get_stringers_at_nodes(
    #     stringer_sections_halfspan, no_string, n_halfspan
    # )
    # Initial guess
    # Spar thickness uniformly decreases from root to tip, skin thickness is constant
    t0 = stack_variables(t_spar_tip0, t_spar_root0, t_skin0, no_string)

    # Define the bounds for each variable
    # extra = 1 if nodes % 2 == 1 else 0
    bounds_t_spar = [(0.0005 * CONVERSION_FACTOR_THICKNESSES, 0.01 * CONVERSION_FACTOR_THICKNESSES)] * 2
    bounds_t_skin = [
        (0.0005 * CONVERSION_FACTOR_THICKNESSES, 0.002 * CONVERSION_FACTOR_THICKNESSES)
    ]  # [(0.0005, 0.01)]
    bounds_no_stringers = [(1 / CONVERSION_FACTOR_STRINGERS, 60 / CONVERSION_FACTOR_STRINGERS)] * (len(no_string))
    bounds = np.array(bounds_t_spar + bounds_t_skin + bounds_no_stringers)

    Cl_DATA, Cdi_DATA, Cm_DATA = interpolate_Cl_Cd_Cm(Cl_data_wing, Cdi_data_wing, Cm_data_wing, y_points)

    material = ac_data["Geometry"]["wingbox_material"]

    L_cruise, D_cruise, M_cruise = get_force_distributions(
        AoA, altitude, ac_data["Performance"]["Vc_m/s"], Cl_DATA, Cdi_DATA, Cm_DATA, chord_dist, ac_data=ac_data
    )

    # Vx, Vy, Vz, Mx, My, Mz = InternalLoads(
    #     L_cruise, D_cruise, M_cruise, wing_structure, ac_data=ac_data, load_factor=nlim
    # )

    rib_dists = np.diff(ac_data["Geometry"]["wing_rib_pos"], prepend=0.0, append=1.0)
    max_rib_dist_m = np.max(rib_dists) * wing_structure.b / 2

    # Constraints dictionary
    MOI_args = [
        wing_structure.stringer_area,
        w_top,
        w_bottom,
        h_avemax,
        h_frontspar,
        h_rearspar,
        n_halfspan,
        nodes,
        stringer_sections_halfspan,
        wing_structure.spar_pos,
        chord_dist_halfspan,
        ac_data,
        nlim,
        L_cruise,
        D_cruise,
        M_cruise,
    ]

    constraints = [
        {
            "type": "ineq",
            "fun": deflection_constraint,
            "args": (y_points_halfspan, max_deflect, ac_data["Materials"][material]["E"], w_fuselage, MOI_args),
        },
        # Rene says no stringer (column) buckling
        # {
        #     "type": "ineq",
        #     "fun": stringer_buckling_constraint,
        #     "args": (
        #         Mx,
        #         Vy,
        #         wing_structure.stringer_area,
        #         hmax,
        #         ac_data["Geometry"]["wing_stringer_MOI_m4"],
        #         max_rib_dist_m,
        #         ac_data["Materials"][material]["E"],
        #         MOI_args,
        #     ),
        # },
        {
            "type": "ineq",
            "fun": tensile_failure_constraint,
            "args": (hmax, ac_data["Materials"][material]["sigma_y"], MOI_args),
        },
        {
            "type": "ineq",
            "fun": stiffened_skin_buckling_constraint,
            "args": (
                wing_structure.spar_pos,
                wing_structure.stringer_area,
                hmax,
                ac_data["Geometry"]["wing_stringer_thickness_m"],
                ac_data["Geometry"]["wing_stringer_small_length_m"],
                ac_data["Materials"][material]["E"],
                ac_data["Materials"][material]["poisson_ratio"],
                ac_data["Materials"][material]["sigma_y"],
                MOI_args,
            ),
        },
        {
            "type": "ineq",
            "fun": shear_buckling_constraint,
            "args": (
                ac_data["Materials"][material]["E"],
                ac_data["Materials"][material]["poisson_ratio"],
                h_rearspar,
                MOI_args,
            ),
        },
    ]

    # Perform the optimization
    len_node_segment = wing_structure.b / (2 * n_halfspan)

    solver_options = {"maxiter": maxiter, "disp": True} if maxiter is not None else {"disp": True}

    objective_args = (
        htot,
        ac_data["Materials"][material]["rho"],
        len_node_segment,
        wing_structure.stringer_area,
        stringer_sections_halfspan,
        n_halfspan,
        chord_dist_halfspan,
        non_fus_idx,
        verbose,
    )

    solution = minimize(
        objective,
        t0,
        args=objective_args,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options=solver_options,
    )

    # Optimized thickness values
    optimized_t_spar_tip, optimized_t_spar_root, optimized_t_skin, optimized_no_stringers = unstack_variables(
        solution.x
    )
    weight = get_weight(solution.x, *objective_args[:-1], print_skin=True)
    # print(solution.message)
    # I = get_I(optimized_t_spar, optimized_t_skin, optimized_no_stringers, wing_structure.stringer_area)
    # deflection = get_deflection(I, y_points_halfspan, Mx, aircraft_data["Materials"][material]["E"])

    if full_return:
        return (
            optimized_t_spar_tip,
            optimized_t_spar_root,
            optimized_t_skin,
            optimized_no_stringers,
            weight,
            y_points_halfspan,
            MOI_args,
            wing_structure.stringer_area,
            wing_structure.spar_pos,
            h_frontspar,
            h_rearspar,
            max_deflect,
            max_rib_dist_m,
            hmax,
            w_fuselage,
        )

    return optimized_t_spar_tip, optimized_t_spar_root, optimized_t_skin, optimized_no_stringers, weight


if __name__ == "__main__":  # pragma: no cover
    (
        optimized_t_spar_tip,
        optimized_t_spar_root,
        optimized_t_skin,
        optimized_no_stringers,
        weight,
        y_points_halfspan,
        MOI_args,
        A_stringer,
        spar_pos,
        h_frontspar,
        h_rearspar,
        max_deflect,
        max_rib_dist,
        hmax,
        w_fuselage,
    ) = run_optimiser(aircraft_data, nodes=501, maxiter=1000, max_deflect=1.7, full_return=True, verbose=True)
    # print(max_rib_dist)

    # MOI, y, M, E
    material = aircraft_data["Geometry"]["wingbox_material"]

    MOI, stringers_at_nodes, t_spar, t_skin, Q, area, enclosed_area, loads = constraint_data_from_variables(
        stack_variables(optimized_t_spar_tip, optimized_t_spar_root, optimized_t_skin, optimized_no_stringers),
        *MOI_args,
    )

    Vx, Vy, Vz, Mx, My, Mz = loads

    deflection = get_deflection(MOI, y_points_halfspan, Mx, aircraft_data["Materials"][material]["E"], w_fuselage)
    axial_forces = get_axial_forces(Mx, Vy, MOI, h_frontspar / 2, stringers_at_nodes, A_stringer)
    axial_stresses = get_max_axial_stress(Mx, Vy, MOI, hmax, area)
    shear = np.abs(get_shear_stress(Vz, My, Q, MOI, enclosed_area, t_spar))

    max_force_allowed = (
        (np.pi) ** 2
        * aircraft_data["Materials"][material]["E"]
        * aircraft_data["Geometry"]["wing_stringer_MOI_m4"]
        / max_rib_dist**2
    )
    sigma_yield = aircraft_data["Materials"][material]["sigma_y"]
    sigma_critical = get_critical_skin_buckling_stress(
        spar_pos,
        stringers_at_nodes,
        A_stringer,
        aircraft_data["Materials"][material]["E"],
        aircraft_data["Materials"][material]["poisson_ratio"],
        sigma_yield,
        t_skin,
        aircraft_data["Geometry"]["wing_stringer_thickness_m"],
        aircraft_data["Geometry"]["wing_stringer_small_length_m"],
    )

    ks = 9.5
    nu = aircraft_data["Materials"][material]["poisson_ratio"]
    E = aircraft_data["Materials"][material]["E"]
    shear_critical = (np.pi) ** 2 * ks * E / (12 * (1 - nu) ** 2) * (t_spar / h_rearspar) ** 2

    # idx = np.argmin((shear_critical - shear))
    # print("------ FINAL ---------")
    # print(np.pi, ks, E, nu, t_spar[idx], h_rearspar[idx], shear_critical[idx])
    # print(shear_critical[idx] / 1e6, shear[idx] / 1e6, ((shear_critical - shear) / shear_critical)[idx], flush=True)

    print("Total time taken:", time.time() - start, "s")

    plt.figure()
    plt.plot(y_points_halfspan, MOI_args[-3][-len(y_points_halfspan) :], "r-", label="Lift")
    plt.plot(y_points_halfspan, Vz, "b-", label="Vz")
    plt.grid()
    plt.legend()

    plt.figure()
    plt.plot(y_points_halfspan, deflection, "g-", label="Deflection")
    plt.plot(y_points_halfspan[[0, -1]], (max_deflect, max_deflect), "k--", label="$v_{max}$")
    plt.grid()

    # plt.figure()
    # plt.plot(y_points_halfspan, axial_forces, "r-", label="Axial Forces")
    # # plt.plot(y_points_halfspan, Mx, "y-", label="Bending Moment")
    # plt.plot(y_points_halfspan[[0, -1]], (max_force_allowed, max_force_allowed), "k--", label="$P_{cr}$")
    # plt.legend()

    plt.figure()
    plt.plot(y_points_halfspan, axial_stresses / 1e6, "b-", label="Axial Stresses")
    plt.plot(y_points_halfspan[[0, -1]], (sigma_yield / 1e6, sigma_yield / 1e6), "k--", label="$\\sigma_{y}$")
    plt.plot(y_points_halfspan, sigma_critical / 1e6, "r--", label="$\\sigma_{cr}$")
    plt.legend()
    plt.grid()

    plt.figure()
    plt.plot(y_points_halfspan, shear / 1e6, "b-", label="Shear Stresses")
    plt.plot(y_points_halfspan, shear_critical / 1e6, "k--", label="Critical")
    plt.legend()
    plt.grid()

    plt.show()

    # Print the results
    print("Optimized thickness for spar:", optimized_t_spar_root, " to ", optimized_t_spar_tip)
    print("Optimized thickness for skin:", optimized_t_skin)
    print("Optimized number of stringers:", optimized_no_stringers)
    print("Minimum weight:", weight)

    if input("Save optimised data to design.json? [y/N]").lower() == "y":
        aircraft_data["Geometry"]["t_spar_root"] = optimized_t_spar_root
        aircraft_data["Geometry"]["t_spar_tip"] = optimized_t_spar_tip
        aircraft_data["Geometry"]["t_skin_wing"] = optimized_t_skin
        aircraft_data["Geometry"]["wing_stringer_number"] = np.ceil(optimized_no_stringers).astype(int).tolist()

        aircraft_data["Structures"]["structural_wing_weight"] = weight * 2  # Multiply by 2 for both wings

        save_ac_data_to_json(aircraft_data)
        print("Data saved to design.json")
