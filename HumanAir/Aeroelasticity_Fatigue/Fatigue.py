import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from HumanAir.aircraft_data import aircraft_data, airfoil_shape, Cl_data_wing, Cdi_data_wing, Cm_data_wing
from HumanAir.StructuralAnalysis.WingStructure import WingStructure
from HumanAir.StructuralAnalysis.optimisation import get_max_axial_stress, get_A, get_force_distributions
from HumanAir.StructuralAnalysis.LoadDistributions import InternalLoads, interpolate_Cl_Cd_Cm


def fatigue_life(Sult=None, alpha=None, Smax=None, K_t=None, verification=False, Experimental_SN=False):
    """
    Generate high-level estimate for the fatigue life of a material using the S-N curve.

    Parameters
    ----------
    Sult : float
        The ultimate tensile strength of the material in [MPa].
    alpha : float
        Correction factor on the ultimate tensile strength, necessary for the final asymptote of the S-N curve.
    Smax : float
        The maximum stress the aircraft should carry without plastic deformation, according to CS23 in [MPa].
    K_t : float
        The fatigue stress concentration factor.
    verification: bool
        If True, the verification on the function will be ran only
    experimental_SN: bool
        If True, the experimental data will be used to extract the S-N curve
    Returns
    -------
    Nf : float
        The fatigue life of the material.

    """
    if not verification:
        if not Experimental_SN:
            # General parameters
            Sult = Sult * 1e6
            if Smax is not None:
                Smax = Smax * 1e6
            N_min = 100
            N_max = 1e6

            if Smax is not None:
                if Smax >= Sult or Smax <= alpha * Sult:
                    raise ValueError("The maximum stress should be between 0 and the ultimate tensile strength.")

            if K_t is None:
                # From N = 0 to N = N_min
                x1 = np.linspace(0, np.log(100), 2)
                S1 = np.ones(len(x1)) * np.log(Sult)

                # From N = N_min to N = N_max
                x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
                S2 = np.log(Sult) + (np.log(alpha * Sult) - np.log(Sult)) / (np.log(N_max) - np.log(N_min)) * (
                    x2 - np.log(N_min)
                )

                # From N = N_max to N = infinity
                x3 = np.linspace(np.log(N_max), 20, 2)
                S3 = np.ones(len(x2)) * np.log(alpha * Sult)

                # Merge all coordinates
                x = np.concatenate((x1, x2, x3))
                S = np.concatenate((S1, S2, S3))

                # Compute intersection point
                if Smax is not None:
                    i = (np.log(Smax) - np.log(Sult)) / (
                        (np.log(alpha * Sult) - np.log(Sult)) / (np.log(N_max) - np.log(N_min))
                    ) + (np.log(N_min))

                    # Lifetime in number of cycles
                    Nf = int(np.exp(i) / 4)  # divide by 4 to account for the coarseness of the approach

                # Plot the S-N curve
                plt.figure()
                plt.plot(np.exp(x), np.exp(S), "b")
                if Smax is not None:
                    # Plot intersection point
                    plt.plot(np.exp(i), Smax, "ro")
                    plt.plot([np.exp(i), np.exp(i)], [0, Smax], "r--")
                    plt.plot([0, np.exp(i)], [Smax, Smax], "r--")
                    # Plot the coodinate of the intersection point next to the point,
                    # with a slight offset using numerical values
                    plt.text(
                        np.exp(i),
                        Smax,
                        "({:.2e}, {:.2e})".format(np.exp(i), Smax),
                        fontsize=12,
                        verticalalignment="bottom",
                    )
                    # put the lifetime in a textbox on the top right corner
                    plt.text(
                        0.95,
                        0.95,
                        f"lifetime: {Nf} flights",
                        transform=plt.gca().transAxes,
                        fontsize=14,
                        verticalalignment="top",
                        horizontalalignment="right",
                        bbox=dict(facecolor="white", alpha=0.5),
                    )
                plt.xscale("log")
                plt.yscale("log")
                plt.xlabel("Number of cycles")
                plt.ylabel("Stress")
                plt.title("S-N curve")
                plt.show()

            if K_t is not None:
                # From N = 0 to N = N_min
                x1 = np.linspace(0, np.log(100), 2)
                S1 = np.ones(len(x1)) * np.log(Sult)

                # From N = N_min to N = N_max
                x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
                S2 = np.log(Sult) + (np.log(alpha * Sult / K_t) - np.log(Sult)) / (np.log(N_max) - np.log(N_min)) * (
                    x2 - np.log(N_min)
                )

                # From N = N_max to N = infinity
                x3 = np.linspace(np.log(N_max), 20, 2)
                S3 = np.ones(len(x2)) * np.log(alpha * Sult / K_t)

                # Merge all coordinates
                x = np.concatenate((x1, x2, x3))
                S = np.concatenate((S1, S2, S3))

                # Compute intersection point
                if Smax is not None:
                    i = (np.log(Smax) - np.log(Sult)) / (
                        (np.log(alpha * Sult) - np.log(Sult)) / (np.log(N_max) - np.log(N_min))
                    ) + (np.log(N_min))

                    # Lifetime in number of cycles
                    Nf = int(np.exp(i) / 4)  # divide by 4 to account for the coarseness of the approach

                # Plot the S-N curve
                plt.figure()
                plt.plot(np.exp(x), np.exp(S), "b")
                if Smax is not None:
                    # Plot intersection point
                    plt.plot(np.exp(i), Smax, "ro")
                    plt.plot([np.exp(i), np.exp(i)], [0, Smax], "r--")
                    plt.plot([0, np.exp(i)], [Smax, Smax], "r--")
                    # Plot the coodinate of the intersection point next to the point,
                    # with a slight offset using numerical values
                    plt.text(
                        np.exp(i),
                        Smax,
                        "({:.2e}, {:.2e})".format(np.exp(i), Smax),
                        fontsize=12,
                        verticalalignment="bottom",
                    )
                    # put the lifetime in a textbox on the top right corner
                    plt.text(
                        0.95,
                        0.95,
                        f"lifetime: {Nf} flights",
                        transform=plt.gca().transAxes,
                        fontsize=14,
                        verticalalignment="top",
                        horizontalalignment="right",
                        bbox=dict(facecolor="white", alpha=0.5),
                    )
                plt.xscale("log")
                plt.yscale("log")
                plt.xlabel("Number of cycles")
                plt.ylabel("Stress")
                plt.title("S-N curve")
                plt.show()

        if Experimental_SN:
            Smax = Smax * 1e6
            data = np.array(
                [
                    (5.00e4, 303.3693),
                    (6.70e4, 296.4746),
                    (2.40e4, 282.685),
                    (8.70e4, 282.685),
                    (6.70e4, 275.7903),
                    (1.01e5, 262.0008),
                    (3.60e4, 255.106),
                    (1.90e5, 241.3165),
                    (4.00e5, 241.3165),
                    (4.00e5, 206.8427),
                    (1.20e6, 206.8427),
                    (2.00e6, 199.948),
                    (4.30e6, 193.0532),
                    (1.10e7, 193.0532),
                    (2.10e6, 186.1584),
                    (2.70e7, 172.3689),
                    (2.10e7, 165.4742),
                    (4.00e7, 165.4742),
                    (6.10e7, 158.5794),
                ]
            )
            data[:, 1] = data[:, 1] * 1e6
            if Smax is not None:
                if Smax >= data[0, 1]:
                    print("Instant failure: Max stress exceeds material ultimate stress, will fail instantly")
                    return None
                if Smax <= data[-1, 1]:
                    print("No fatigue: Max stress below material endurance limit, will not fail due to fatigue")
                    return None
            # Linearly interpolate between curve points to get the Nf for the given Smax
            if Smax is not None:
                for i in range(len(data)):
                    if Smax <= data[i, 1] and Smax >= data[i + 1, 1]:
                        Nf = int(
                            np.exp(
                                np.log(data[i, 0])
                                + (np.log(data[i + 1, 0]) - np.log(data[i, 0]))
                                / (np.log(data[i + 1, 1]) - np.log(data[i, 1]))
                                * np.log((Smax / data[i, 1]))
                            )
                        )
                        break
            Nf_actual = Nf / 4

            # plot the experimental data curve
            plt.figure()
            plt.plot(data[:, 0], data[:, 1])
            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("Number of cycles")
            plt.ylabel("Stress")
            plt.title("S-N curve")
            if Smax is not None:
                plt.plot(Nf, Smax, "ro")
                plt.plot([Nf, Nf], [0, Smax], "r--")
                plt.plot([0, Nf], [Smax, Smax], "r--")
                # Plot the coodinate of the intersection point next to the point,
                # with a slight offset using numerical values
                plt.text(Nf, Smax, "({:.2e}, {:.2e})".format(Nf, Smax), fontsize=12, verticalalignment="bottom")
                # put the lifetime in a textbox on the top right corner
                plt.text(
                    0.95,
                    0.95,
                    f"lifetime: {Nf_actual} flights",
                    transform=plt.gca().transAxes,
                    fontsize=14,
                    verticalalignment="top",
                    horizontalalignment="right",
                    bbox=dict(facecolor="white", alpha=0.5),
                )
            plt.show()

    # test = False
    # if test:
    #     # From N = 0 to N = N_min
    #     x1 = np.linspace(0, np.log(100), 2)
    #     S1 = np.ones(len(x1))*np.log(Sult)

    #     # From N = N_min to N = N_max
    #     x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
    #     S2 = np.log(Sult) +
    #       (np.log(alpha * Sult) - np.log(Sult))/(np.log(N_max) - np.log(N_min)) * (x2- np.log(N_min))

    #     # From N = N_max to N = infinity
    #     x3 = np.linspace(np.log(N_max), 20, 2)
    #     S3 = np.ones(len(x2))* np.log(alpha * Sult)

    #     # Compute intersection point
    #     x = (np.log(Smax) - np.log(Sult)) /
    #           ((np.log(alpha * Sult) - np.log(Sult))/(np.log(N_max) - np.log(N_min))) + (np.log(N_min))

    #     # Plot the Basquin curve
    #     plt.figure()
    #     plt.plot(np.exp(x1), np.exp(S1), 'b')
    #     plt.plot(np.exp(x2), np.exp(S2), 'b')
    #     plt.plot(np.exp(x3), np.exp(S3), 'b')
    #     # Draw horizontal and vertical lines to intersection point, and note the values on the axes
    #     # Plot intersection point
    #     plt.plot(np.exp(x), Smax, 'ro')
    #     plt.plot([np.exp(x), np.exp(x)], [0, Smax], 'r--')
    #     plt.plot([0, np.exp(x)], [Smax, Smax], 'r--')
    #     # Plot the coodinate of the intersection point next to the point, with a slight offset using numerical values
    #     plt.text(np.exp(x), Smax, '({:.2e}, {:.2e})'.format(np.exp(x), Smax), fontsize=12, verticalalignment='bottom')

    #     plt.xscale('log')
    #     plt.yscale('log')
    #     plt.xlabel('Number of cycles')
    #     plt.ylabel('Stress')
    #     plt.title('Basquin curve')

    #     # From N = 0 to N = N_min
    #     x1 = np.linspace(0, np.log(100), 2)
    #     S1 = np.ones(len(x1))*np.log(Sult)

    #     # From N = N_min to N = N_max
    #     x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
    #     S2 = np.log(Sult) +
    #       (np.log(alpha * Sult / K_t) - np.log(Sult))/(np.log(N_max) - np.log(N_min)) * (x2- np.log(N_min))

    #     # From N = N_max to N = infinity
    #     x3 = np.linspace(np.log(N_max), 20, 2)
    #     S3 = np.ones(len(x2))* np.log(alpha * Sult / K_t)

    #     # Compute intersection point
    #     x = (np.log(Smax) - np.log(Sult)) /
    #       ((np.log(alpha * Sult / K_t) - np.log(Sult))/(np.log(N_max) - np.log(N_min))) + (np.log(N_min))

    #     # Plot the Basquin curve

    #     plt.plot(np.exp(x1), np.exp(S1), 'b')
    #     plt.plot(np.exp(x2), np.exp(S2), 'b')
    #     plt.plot(np.exp(x3), np.exp(S3), 'b')
    #     # Draw horizontal and vertical lines to intersection point, and note the values on the axes
    #     # Plot intersection point
    #     plt.plot(np.exp(x), Smax, 'ro')
    #     plt.plot([np.exp(x), np.exp(x)], [0, Smax], 'r--')
    #     plt.plot([0, np.exp(x)], [Smax, Smax], 'r--')
    #     # Plot the coodinate of the intersection point next to the point, with a slight offset using numerical values
    #     plt.text(np.exp(x), Smax, '({:.2e}, {:.2e})'.format(np.exp(x), Smax), fontsize=12, verticalalignment='bottom')

    #     plt.xscale('log')
    #     plt.yscale('log')
    #     plt.xlabel('Number of cycles')
    #     plt.ylabel('Stress')
    #     plt.title('Basquin curve')
    #     plt.show()

    if verification:
        data = np.array(
            [
                (5.00e4, 303.3693),
                (6.70e4, 296.4746),
                (2.40e4, 282.685),
                (8.70e4, 282.685),
                (6.70e4, 275.7903),
                (1.01e5, 262.0008),
                (3.60e4, 255.106),
                (1.90e5, 241.3165),
                (4.00e5, 241.3165),
                (4.00e5, 206.8427),
                (1.20e6, 206.8427),
                (2.00e6, 199.948),
                (4.30e6, 193.0532),
                (1.10e7, 193.0532),
                (2.10e6, 186.1584),
                (2.70e7, 172.3689),
                (2.10e7, 165.4742),
                (4.00e7, 165.4742),
                (6.10e7, 158.5794),
            ]
        )
        data[:, 1] = data[:, 1] * 1e6
        Sult = Sult * 1e6
        if Smax is not None:
            Smax = Smax * 1e6
        alpha = alpha
        N_min = 100
        N_max = 1e6
        # From N = 0 to N = N_min
        x1 = np.linspace(0, np.log(100), 2)
        S1 = np.ones(len(x1)) * np.log(Sult)

        # From N = N_min to N = N_max
        x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
        S2 = np.log(Sult) + (np.log(alpha * Sult) - np.log(Sult)) / (np.log(N_max) - np.log(N_min)) * (
            x2 - np.log(N_min)
        )

        # From N = N_max to N = infinity
        x3 = np.linspace(np.log(N_max), 20, 2)
        S3 = np.ones(len(x2)) * np.log(alpha * Sult)

        # Merge all coordinates
        x = np.concatenate((x1, x2, x3))
        S = np.concatenate((S1, S2, S3))

        # Compute intersection point
        if Smax is not None:
            i = (np.log(Smax) - np.log(Sult)) / (
                (np.log(alpha * Sult) - np.log(Sult)) / (np.log(N_max) - np.log(N_min))
            ) + (np.log(N_min))

            # Lifetime in number of cycles
            Nf = int(np.exp(i) / 4)  # divide by 4 to account for the coarseness of the approach

        # round Nf to nearest 100
        Nf = round(Nf, -2)

        # find exponential fit for the experimental data
        p = np.polyfit(np.log(data[:, 0]), np.log(data[:, 1]), 1)

        # plot the S-N curves
        plt.figure()
        plt.plot(np.exp(x), np.exp(S), "b-")
        plt.plot(data[:, 0], data[:, 1], "o", color="lightgreen")
        plt.plot(data[:, 0], np.exp(np.polyval(p, np.log(data[:, 0]))), "r-")
        if Smax is not None:
            # Plot intersection point
            plt.plot(np.exp(i), Smax, "ro")
            plt.plot([np.exp(i), np.exp(i)], [0, Smax], "r--")
            plt.plot([0, np.exp(i)], [Smax, Smax], "r--")
            # Plot the coodinate of the intersection point next to the point,
            # with a slight offset using numerical values
            plt.text(
                np.exp(i) * 1.2,
                Smax * 1.01,
                "({:.2e}, {:.2e})".format(np.exp(i), Smax),
                fontsize=12,
                verticalalignment="bottom",
            )
            # put the lifetime in a textbox on the top right corner
            plt.text(
                0.95,
                0.95,
                f"Nf: {Nf} flights",
                transform=plt.gca().transAxes,
                fontsize=14,
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(facecolor="white", alpha=0.5),
            )
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Number of cycles (Nf)")
        plt.ylabel("Stress")
        plt.legend(["Basquin law approximation", "Experimental data", "Exponential fit"])
        plt.savefig(
            os.path.join(os.path.dirname(__file__), "..", "..", "Figures", "Fatigue.pdf"),
            bbox_inches="tight",
            dpi=200,
        )
        plt.show()

        # # Plot the S-N curve
        # plt.figure()
        # plt.plot(np.exp(x), np.exp(S), "b-")
        # # plot the experimental data curve
        # plt.plot(data[:, 0], data[:, 1], "g-")
        # plt.xscale("log")
        # plt.yscale("log")
        # plt.xlabel("Number of cycles")
        # plt.ylabel("Stress")
        # plt.title("S-N curve")
        # plt.legend(["Basquin law approximation", "Experimental data"])
        # plt.show()


if __name__ == "__main__":
    stress_ult = aircraft_data["Materials"]["Aluminium"]["sigma_ult"] / 1e6
    nodes = 501
    load_factor = 1.5
    AoA = 0
    altitude = 3000

    wing_structure = WingStructure(aircraft_data, airfoil_shape, nodes)

    # t_skin, t_spar, no_stringers, A_stringer, spar_pos, h_frontspar, h_rearspar, chord_dist
    _, h_s1s2 = wing_structure.h_s1s2()
    # l_box_up, l_box_down = wing_structure.d_s1s2()

    h_s1s2 = h_s1s2[nodes // 2 :]
    # l_box_up, l_box_down = l_box_up[nodes//2:], l_box_down[nodes//2:]
    # ypts_halfspan =  wing_structure.ypts[nodes//2:]

    hmax = wing_structure.hmax_dist[nodes // 2 :]
    h_frontspar = h_s1s2[:, 0]  # height of the spar at 15% of the chord, as a function of y
    h_rearspar = h_s1s2[:, 1]  # height of the spar at 50% of the chord, as a function of y
    area = get_A(
        wing_structure.t_skin,
        wing_structure.t_spar_dist[nodes // 2 :],
        wing_structure.stringer_dist[nodes // 2 :],
        wing_structure.stringer_area,
        wing_structure.spar_pos,
        h_frontspar,
        h_rearspar,
        wing_structure.chord_distribution[nodes // 2 :],
    )

    Cl_DATA, Cdi_DATA, Cm_DATA = interpolate_Cl_Cd_Cm(Cl_data_wing, Cdi_data_wing, Cm_data_wing, wing_structure.ypts)

    material = aircraft_data["Geometry"]["wingbox_material"]

    L_cruise, D_cruise, M_cruise = get_force_distributions(
        AoA,
        altitude,
        aircraft_data["Performance"]["Vc_m/s"],
        Cl_DATA,
        Cdi_DATA,
        Cm_DATA,
        wing_structure.chord_distribution,
        ac_data=aircraft_data,
    )

    Vx, Vy, Vz, Mx, My, Mz = InternalLoads(
        L_cruise, D_cruise, M_cruise, wing_structure, ac_data=aircraft_data, load_factor=load_factor
    )

    axial_stresses = get_max_axial_stress(Mx, Vy, wing_structure.Ixx()[nodes // 2 :], hmax, area)
    stress_max = np.max(np.abs(axial_stresses)) / 1e6
    print(stress_max, stress_ult)

    fatigue_life(Sult=stress_ult, alpha=0.30, Smax=stress_max, verification=False, Experimental_SN=False)
    fatigue_life(Sult=stress_ult, alpha=0.30, Smax=stress_max, verification=False, Experimental_SN=True)
    fatigue_life(Sult=stress_ult, alpha=0.30, Smax=stress_max, verification=True)
