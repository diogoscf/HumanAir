import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate  # type: ignore[import-untyped]
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data, airfoil_shape, Cl_data_wing, Cm_data_wing, Cdi_data_wing
from HumanAir.StructuralAnalysis.WingStructure import WingStructure

# from HumanAir.StructuralAnalysis.optimisation import get_deflection
from HumanAir.isa import isa

# from HumanAir.unit_conversions import G


# Define the forces along half span
# def chord(Sw, taper_ratio, Cl_DATA, AoA, nodes):
#     b = Cl_DATA[AoA]["y_span"][-1] * 2
#     # Generate spanwise coordinate points
#     y = np.linspace(Cl_DATA[AoA]["y_span"][0], Cl_DATA[AoA]["y_span"][-1], nodes)  # n is the number of nodes
#     # Calculate the chord distribution
#     chord_length = 2 * Sw / (1 + taper_ratio) / b * (1 - (1 - taper_ratio) * np.abs(2 * y / b))
#     return chord_length, y


# NOTE: w_fuselage is the full width, it will be divided by 2
def get_deflection(MOI, y, M, E, w_fuselage):
    non_fus_idx = np.argmin(np.abs(y - (w_fuselage / 2)))
    integrand = M[non_fus_idx:] / MOI[non_fus_idx:]
    # print(np.sum(I), np.sum(M), np.sum(integrand))
    dvdy = -1 / E * integrate.cumulative_trapezoid(integrand, y[non_fus_idx:], initial=0)
    v = -integrate.cumulative_trapezoid(dvdy, y[non_fus_idx:], initial=0)
    v = np.concatenate((np.zeros(non_fus_idx), v))
    return v


def get_twist(J, y, T, G):
    integrand = T / J
    # print(np.sum(I), np.sum(M), np.sum(integrand))
    theta = 1 / G * integrate.cumulative_trapezoid(integrand, y, initial=0)
    return theta


def force_distribution(AoA, altitude, V, chord_dist, Cl_DATA, Cdi_DATA):
    rho = isa(altitude)[2]
    Cl = np.array(Cl_DATA[AoA]["coefficient"])
    Cdi = np.array(Cdi_DATA[AoA]["coefficient"])
    L = Cl * 0.5 * rho * V**2 * chord_dist  # [N/m]
    D = Cdi * 0.5 * rho * V**2 * chord_dist  # [N/m]
    # plt.plot(Cl_DATA[AoA]["y_span"], L, label="Lift")
    # plt.ylim(0, np.max(L) * 1.1)
    # plt.show()
    return L, D


def moment_distribution(AoA, altitude, V, chord_dist, Cm_DATA, ac_data=aircraft_data):
    MAC = ac_data["Aero"]["MAC_wing"]
    M = np.array(Cm_DATA[AoA]["coefficient"]) * 0.5 * isa(altitude)[2] * V**2 * chord_dist * MAC  # [Nm/m]
    return M


def weight_distribution(c, wing_structure, ac_data=aircraft_data):
    fuelweight = ac_data["CL2Weight"]["Wfuel_N"]

    # Fuel is located between the struts
    strut_loc = ac_data["Geometry"]["strut_loc_b/2"]
    span = ac_data["Aero"]["b_Wing"]
    nodes_half_wing = c.shape[0] // 2
    n_before_strut = np.rint(nodes_half_wing * strut_loc).astype(int)
    extra = 0 if c.shape[0] % 2 == 0 else 1
    idx_between_struts = (nodes_half_wing - n_before_strut, nodes_half_wing + n_before_strut + extra)
    c_between_struts = c[idx_between_struts[0] : idx_between_struts[1]]
    enclosed_area = wing_structure.enclosed_area_dist[idx_between_struts[0] : idx_between_struts[1]]

    W_structure = wing_structure.weight_dist()

    # Fuel weight distributed with wingbox enclosed area
    W_avg_fuel = fuelweight / (span * strut_loc)  # [N/m]
    W_fuel = W_avg_fuel * enclosed_area / (np.mean(enclosed_area))
    diff = c.shape[0] - c_between_struts.shape[0]
    W_fuel = np.pad(W_fuel, pad_width=diff // 2, mode="constant")

    # [N/m]
    return (
        W_structure + W_fuel,
        W_fuel,
        (nodes_half_wing - n_before_strut, nodes_half_wing + n_before_strut + extra),
        c_between_struts,
    )


# axial forces distribution - during ground - so only the ends (0.6L)
# Vy_Strut = 6670 / np.tan(0.546)  # deg in radian


def axial_distribution_ground(nodes, Vy_strut, ac_data=aircraft_data):
    strut_loc = ac_data["Geometry"]["strut_loc_b/2"]
    nodes_half_wing = np.ceil(nodes / 2).astype(int)
    n_before_strut = np.rint(nodes_half_wing * strut_loc).astype(int)
    A_ground = np.zeros(nodes_half_wing)
    A_ground[:n_before_strut] = Vy_strut
    return A_ground


# axial forces distribution - during flight - so only on the inner (inside 0.4L))
def axial_distribution_flight(nodes, Vy_strut, ac_data=aircraft_data):
    strut_loc = ac_data["Geometry"]["strut_loc_b/2"]
    nodes_half_wing = np.ceil(nodes / 2).astype(int)
    n_before_strut = np.rint(nodes_half_wing * strut_loc).astype(int)
    A_cruise = np.zeros(nodes_half_wing)
    A_cruise[:n_before_strut] = -Vy_strut
    return A_cruise


def read_points_from_load_dist(L_cruise, W_cruise, W_fuel, idxs):
    a = np.max(L_cruise)
    b = np.min(L_cruise)
    c = np.min(W_cruise)
    d = W_cruise[idxs[1]]
    e = d + W_fuel[idxs[0]]
    f = np.max(W_cruise)

    return a, b, c, d, e, f


def strut_error_calculation(
    P, Vz, n_before_strut, theta_strut, y_points_halfspan, MOI, E, l_strut, A_strut, w_fuselage
):
    Vz_strut = P * np.sin(theta_strut)
    Vz = Vz.copy()
    Vz[:n_before_strut] += Vz_strut  # [N]
    Mx = -integrate.cumulative_trapezoid(np.flip(Vz * y_points_halfspan), y_points_halfspan, initial=0)
    Mx = Mx[::-1]
    v_wing = get_deflection(MOI, y_points_halfspan, Mx, E, w_fuselage)[n_before_strut]
    v_strut = Vz_strut * l_strut / (A_strut * E)

    return v_wing - v_strut


def strut_force(MOI, y_points_halfspan, Vz_orig, w_fuselage, max_iter=1000, tol=1e-6, ac_data=aircraft_data):
    strut_loc = ac_data["Geometry"]["strut_loc_b/2"]
    nodes_half_wing = len(y_points_halfspan)
    n_before_strut = np.rint(nodes_half_wing * strut_loc).astype(int)

    h_fus = ac_data["Geometry"]["fus_height_m"]
    halfspan = ac_data["Aero"]["b_Wing"] / 2
    l_strut = np.sqrt((halfspan * strut_loc) ** 2 + h_fus**2)
    theta_strut = np.arctan(h_fus / (halfspan * strut_loc))

    E = ac_data["Materials"][ac_data["Geometry"]["wingbox_material"]]["E"]
    A_strut = ac_data["Geometry"]["strut_section_area_m2"]
    # P = -Vz_orig[n_before_strut]*3.5  # [N], initial guess

    # Mx = -integrate.cumulative_trapezoid(np.flip(Vz_orig * y_points_halfspan), y_points_halfspan, initial=0)[::-1]
    # v_40_orig = get_deflection(MOI, y_points_halfspan, Mx, E)[n_before_strut]
    # P = strut_force_zero_deflection(MOI, v_40_orig, y_points_halfspan[:n_before_strut], E)

    def strut_error_calc_w_args(p):
        return strut_error_calculation(
            p, Vz_orig, n_before_strut, theta_strut, y_points_halfspan, MOI, E, l_strut, A_strut, w_fuselage
        )

    P_curr = 0
    err_curr = strut_error_calc_w_args(P_curr)

    step = 10000

    for _ in range(max_iter):
        # deflection at 40% of halfspan
        P_next = P_curr + step
        err_next = strut_error_calc_w_args(P_next)

        if np.abs(err_next) < tol:
            P_curr = P_next
            err_curr = err_next
            break

        if err_curr * err_next < 0:  # found the bracket
            step /= 2
            if abs(err_curr) < abs(err_next):
                P_next = P_curr
                err_next = err_curr
            continue

        if np.abs(err_curr) < np.abs(err_next):
            step = -step

        P_curr = P_next
        err_curr = err_next

        # V_strut = P
        # Vz_strut = P * np.cos(theta_strut)
        # Vy_strut = P * np.sin(theta_strut)

        # Vz = Vz_orig.copy()
        # Vz[:n_before_strut] += Vz_strut  # [N]

        # Mx = -integrate.cumulative_trapezoid(np.flip(Vz * y_points_halfspan), y_points_halfspan, initial=0)
        # Mx = Mx[::-1]

        # v = get_deflection(MOI, y_points_halfspan, Mx, E)
        # v_40 = v[n_before_strut]
        # if v_40 > 0.5:
        #     v_40 = 0.5
        # if v_40 < 0:
        #     v_40 = 0

        # deflection at an angle
        # delta = v_40 / np.sin(theta_strut)

        # # strut force
        # P_new = delta * A_strut * E / l_strut

        # # print(P, v_40, delta, P_new, flush=True)

        # # plt.plot(y_points_halfspan, Vz, "r-")
        # # plt.plot(y_points_halfspan, Mx, "g-")
        # # plt.ylim(-0.2e5, 1.2e5)
        # # plt.twinx()
        # # plt.plot(y_points_halfspan, v, "b-")
        # # plt.ylim(-0.2, 1.2)
        # # plt.grid()
        # # plt.show()

        # if abs(P_new - P) < tol:
        #     break

        # # update
        # P = P_new

    V_strut = P_curr
    Vz_strut = P_curr * np.sin(theta_strut)
    Vy_strut = P_curr * np.cos(theta_strut)

    # print(V_strut)
    # print(l_strut)

    return Vz_strut, Vy_strut, V_strut


def InternalLoads(L, D, M, wing_structure, ac_data=aircraft_data, load_factor=1):
    y_points = wing_structure.ypts
    nodes = wing_structure.nodes
    chord_dist = wing_structure.chord_distribution

    W, *_ = weight_distribution(chord_dist, wing_structure, ac_data=ac_data)
    nodes_half_wing = nodes // 2
    # b_half = ac_data["Aero"]["b_Wing"] / 2
    sweep = ac_data["Aero"]["QuarterChordSweep_Wing_deg"]

    # X_forces = -(D * (load_factor**2)) * b_half / (nodes_half_wing)  # drag and thrust act on the x axis (thrust = 0)
    X_force_distribution = -D * (load_factor**2)
    Vx = integrate.cumulative_trapezoid(
        np.flip(X_force_distribution[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0
    )
    Vx = Vx[::-1]
    # Z_force_distribution = W - (L * (b_half / (nodes_half_wing))) * load_factor
    Z_force_distribution = W - L * load_factor
    Vz = integrate.cumulative_trapezoid(
        np.flip(Z_force_distribution[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0
    )
    Vz = Vz[::-1]

    strut_Vz, strut_Vy, _ = strut_force(
        wing_structure.Ixx()[nodes_half_wing:],
        y_points[nodes_half_wing:],
        Vz,
        wing_structure.w_fuselage,
        max_iter=100,
        tol=1e-6,
        ac_data=ac_data,
    )

    # Vy - axial force diagram because of the strut
    Ax_total_flight = axial_distribution_flight(
        nodes, strut_Vy, ac_data
    )  # axial force in the y span, during flight -- this is what we want now cause of cruise
    # Ax_total_ground = axial_distribution_ground(nodes, strut_Vy, ac_data)  # axial force in the y span, on the ground
    # Vy = np.append(Ax_total_flight, [0])
    Vy = Ax_total_flight

    # Vz - shear force because of the strut
    strut_loc = ac_data["Geometry"]["strut_loc_b/2"]
    n_before_strut = np.rint(nodes_half_wing * strut_loc).astype(int)
    Vz[:n_before_strut] += strut_Vz  # [N]

    # add the moment about x
    # print(b, Vz.shape, n, y_points.shape, M.shape)
    Mx = -integrate.cumulative_trapezoid(
        np.flip(Vz * y_points[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0
    )
    Mx = Mx[::-1]
    Mz = -integrate.cumulative_trapezoid(
        np.flip(Vx * y_points[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0
    )
    Mz = Mz[::-1]

    # plt.plot(y_points[nodes_half_wing:], Z_forces[nodes_half_wing:], label="W-L")
    # plt.plot(y_points[nodes_half_wing:], Vz, label="Vz")
    # plt.plot(y_points[nodes_half_wing:], Mx, label="Mx")
    # plt.legend()
    # plt.show()

    # due to lift and weight moment arm (assume both are at c/4)
    My_lw = integrate.cumulative_trapezoid(
        np.flip(Z_force_distribution[nodes_half_wing:] * np.tan(sweep)), y_points[nodes_half_wing:], initial=0
    )
    My_m = integrate.cumulative_trapezoid(
        np.flip((M * load_factor)[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0
    )

    My = (My_lw + My_m)[::-1]

    return Vx, Vy, Vz, Mx, My, Mz


def interpolate_Cl_Cd_Cm(Cl_data, Cdi_data, Cm_data, y_points):
    nodes_data = len(Cl_data[0]["coefficient"])
    ypts_orig = np.linspace(Cl_data[0]["y_span"][0], Cl_data[0]["y_span"][-1], nodes_data)

    # plt.plot(ypts_orig, Cl_data[0.0]["coefficient"])
    # plt.show()

    for angle in Cl_data.keys():
        Cl_data[angle]["coefficient"] = np.interp(y_points, ypts_orig, Cl_data[angle]["coefficient"])
        Cm_data[angle]["coefficient"] = np.interp(y_points, ypts_orig, Cm_data[angle]["coefficient"])
        Cdi_data[angle]["coefficient"] = np.interp(y_points, ypts_orig, Cdi_data[angle]["coefficient"])

        Cl_data[angle]["y_span"] = y_points
        Cm_data[angle]["y_span"] = y_points
        Cdi_data[angle]["y_span"] = y_points

    return Cl_data, Cdi_data, Cm_data


if __name__ == "__main__":  # pragma: no cover
    # Data
    AoA = 0
    # Sw = 39  # [m2]
    # taper_ratio = 0.4
    # Vcruise = 60  # [m/s]
    # rho = 0.9  # [kg/m3]
    # structuralmass = 5250 / 9.81
    # batterymass_w = 0
    nl = 3.8  # Load Factor
    nl2 = -1.52
    # sweep = 0.157
    altitude = 3000  # [m]
    nodes = 401

    wing_structure_data = WingStructure(aircraft_data, airfoil_shape, nodes)
    chord_dist = wing_structure_data.chord_distribution
    y_points = wing_structure_data.ypts
    # print(Cl_DATA[-10.0].keys())

    Cl_DATA, Cdi_DATA, Cm_DATA = interpolate_Cl_Cd_Cm(Cl_data_wing, Cdi_data_wing, Cm_data_wing, y_points)

    # Cm_DATA = np.interp(y_points, ypts_orig, Cm_DATA[AoA]["coefficient"])
    # Cdi_DATA = np.interp(y_points, ypts_orig, Cdi_DATA[AoA]["coefficient"])

    L_cruise, D_cruise = force_distribution(
        AoA, altitude, aircraft_data["Performance"]["Vc_m/s"], chord_dist, Cl_DATA=Cl_DATA, Cdi_DATA=Cdi_DATA
    )

    M_cruise = moment_distribution(
        AoA, altitude, aircraft_data["Performance"]["Vc_m/s"], chord_dist, Cm_DATA, ac_data=aircraft_data
    )

    # print(np.sum(L_cruise))

    # nl = origin
    Vx1, Vy1, Vz1, Mx1, My1, Mz1 = InternalLoads(
        L_cruise, D_cruise, M_cruise, wing_structure_data, ac_data=aircraft_data, load_factor=1
    )

    # nl = 3.8
    Vx, Vy, Vz, Mx, My, Mz = InternalLoads(
        L_cruise, D_cruise, M_cruise, wing_structure_data, ac_data=aircraft_data, load_factor=nl
    )

    # nl = -1
    Vx2, Vy2, Vz2, Mx2, My2, Mz2 = InternalLoads(
        L_cruise, D_cruise, M_cruise, wing_structure_data, ac_data=aircraft_data, load_factor=nl2
    )

    # extra = 1 if nodes % 2 == 0 else 0
    y_points_plot = y_points[(nodes // 2) :]
    plt.rcParams.update({"font.size": 30})
    plt.axhline(linewidth=2, color="k")
    plt.axvline(linewidth=2, color="k")

    # plt.subplot(2, 2, 1)
    # plt.plot(y_points, nl*L_cruise, label='Lift')
    # plt.plot(y_points, -W_cruise, label = 'Weight')
    # plt.title('Lift and weight distribution along Span')
    # plt.xlabel('Spanwise Position [m]')
    # plt.ylabel('Force [N]')
    # plt.legend()
    # plt.xlim(left=0)
    # plt.tight_layout()
    # plt.grid()

    # plt.subplot(2, 2, 1)
    # plt.plot(y_points_plot, Vx, label="n=3.8")
    # plt.plot(y_points_plot, Vx1, label="cruise")
    # plt.plot(y_points_plot, Vx2, label="n=-1")
    # plt.title("Shear Force Vx along Span")
    # plt.xlabel("Spanwise Position [m]")
    # plt.ylabel("Vx [N]")
    # # plt.xlim(left=0)
    # plt.grid()
    # plt.legend()

    # plt.subplot(2, 2, 1)
    # plt.plot(y_points_plot, Vy, "r-", label="n = 3.51", lw=3)
    # plt.plot(y_points_plot, Vy1, "g-", label="Cruise", lw=3)
    # plt.plot(y_points_plot, Vy2, "b-", label="n = -1.4", lw=3)
    # plt.title("Axial")
    # plt.xlabel("Spanwise Position [m]")
    # plt.ylabel("Vy [N]")
    # plt.xlim(left=0)
    # plt.grid()
    # plt.legend()

    # plt.subplot(2, 2, 2)
    # plt.plot(y_points_plot, Vz, "r-", label="n = 3.51", lw=3)
    # plt.plot(y_points_plot, Vz1, "g-", label="Cruise", lw=3)
    # plt.plot(y_points_plot, Vz2, "b-", label="n = -1.4", lw=3)
    # plt.title("Shear")
    # plt.xlabel("Spanwise Position [m]")
    # plt.ylabel("Vz [N]")
    # plt.xlim(left=0)
    # plt.legend()
    # plt.grid()

    # plt.subplot(2, 2, 3)
    plt.plot(y_points_plot, Mx, "r-", label="n = 3.51", lw=3)
    plt.plot(y_points_plot, Mx1, "g-", label="Cruise", lw=3)
    plt.plot(y_points_plot, Mx2, "b-", label="n = -1.4", lw=3)
    plt.title("Bending Moment")
    plt.xlabel("Spanwise Position [m]")
    plt.ylabel("Mx [Nm]")
    # plt.xlim(left=0)
    # plt.legend()
    # plt.grid()

    # plt.subplot(2, 2, 4)
    # plt.plot(y_points_plot, Mz, label="n=3.8")
    # plt.plot(y_points_plot, Mz1, label="cruise")
    # plt.plot(y_points_plot, Mz2, label="n=-1")
    # plt.title("Torque Mz along Span")
    # plt.xlabel("Spanwise Position [m]")
    # plt.ylabel("Mz [Nm]")
    # # plt.xlim(left=0)
    # plt.legend()
    # plt.grid()

    # plt.subplot(2, 2, 4)
    # plt.plot(y_points_plot, My, "r-", label="n = 3.51")
    # plt.plot(y_points_plot, My1, "g-", label="Cruise")
    # plt.plot(y_points_plot, My2, "b-", label="n = -1.4")
    # plt.title("Torque My along Span")
    # plt.xlabel("Spanwise Position [m]")
    # plt.ylabel("My [Nm]")
    # plt.legend()
    # plt.grid()

    plt.legend(loc="upper center", ncol=3, fancybox=True, shadow=True)
    plt.grid()
    plt.gcf().set_size_inches(18, 16)
    plt.tight_layout()
    plt.savefig(
        os.path.join(os.path.dirname(__file__), "..", "..", "Figures", "shear-force.png"), bbox_inches="tight", dpi=200
    )
    plt.show()

    # twist = (
    #     get_twist(
    #         wing_structure_data.torsional_constant()[nodes // 2 :],
    #         wing_structure_data.ypts[nodes // 2 :],
    #         My,
    #         aircraft_data["Materials"][aircraft_data["Geometry"]["wingbox_material"]]["G"],
    #     )
    #     * 180
    #     / np.pi
    # )

    # deflection = get_deflection(
    #     wing_structure_data.Ixx()[nodes // 2 :],
    #     y_points_plot,
    #     Mx,
    #     aircraft_data["Materials"][aircraft_data["Geometry"]["wingbox_material"]]["E"],
    # )
    # plt.plot(y_points_plot, deflection)
    # plt.show()
