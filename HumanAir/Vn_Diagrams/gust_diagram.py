import numpy as np
import os
import matplotlib.pyplot as plt
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.isa import isa
from HumanAir.Vn_Diagrams.loading_diagram import calculate_manoeuvre_velocities
from HumanAir.aircraft_data import aircraft_data

COMMUTER = False

ft_to_m = 0.3048
FL_to_m = 0.3048 * 100.0
m_to_ft = 1 / 0.3048
m_to_FL = 1 / (0.3048 * 100.0)
lbs_to_kg = 0.45359237
kg_to_lbs = 1 / 0.45359237

G = 9.80665  # [m/s^2]

FONTSIZE = 40


def calculate_gust_diagram_loads(aircraft_data, Vc_ms, Vd_ms, V_S1, commuter_ac=False, h=0, temp_offset=0):
    """
    Calculate the applicable loads and velocities for the gust diagram

    Parameters
    ----------
    aircraft_data : dict
        Dictionary containing the aircraft data.
    Vc_ms, Vd_ms, V_S1 : float
        The cruise speed, dive speed, and stall speed in the landing configuration [m/s].
    commuter_ac : bool, default False
        Boolean indicating whether the aircraft is a commuter aircraft.
    h : float, default 0
        The altitude at which to calculate the ISA conditions [m].
        Note that this only affects the resulting velocities, not the loads.

    Returns
    -------
    n_max, n_min : float
        The maximum and minimum load factors, and the load factors at the respective velocities.
    V_B : float
        The velocity at which the gust diagram intersects the stall speed [m/s].
        Only relevant for commuter aircraft.
    n_cruise_pve, n_cruise_nve: float
        The load factors at the cruise speed.
    n_dive_pve, n_dive_nve: float
        The load factors at the dive speed.
    n_B_pve, n_B_nve: float
        The load factors at V_B.
    """

    _, _, rho = isa(h, temp_offset)
    _, _, rho_0 = isa(0)

    WS_Nm2 = aircraft_data["Performance"]["W/S_N/m2"]

    MGC_m = aircraft_data["Aero"]["MAC_wing"]
    clalpha = aircraft_data["Aero"]["CLalpha"]

    CLmax_clean = aircraft_data["Aero"]["CLmax_clean"]

    mu_g = 2 * WS_Nm2 / (rho * MGC_m * clalpha * G)
    k_g = 0.88 * mu_g / (5.3 + mu_g)
    # print(k_g)

    Ude_cruise_fps = 50
    Ude_cruise_ms = Ude_cruise_fps * ft_to_m

    Ude_dive_fps = 25
    Ude_dive_ms = Ude_dive_fps * ft_to_m

    Ude_B_fps = 66  # Only for commuter airplanes
    Ude_B_ms = Ude_B_fps * ft_to_m

    n_cruise_pve = 1 + (k_g * rho_0 * Ude_cruise_ms * Vc_ms * clalpha) / (2 * WS_Nm2)
    n_cruise_nve = 1 - (k_g * rho_0 * Ude_cruise_ms * Vc_ms * clalpha) / (2 * WS_Nm2)

    n_dive_pve = 1 + (k_g * rho_0 * Ude_dive_ms * Vd_ms * clalpha) / (2 * WS_Nm2)
    n_dive_nve = 1 - (k_g * rho_0 * Ude_dive_ms * Vd_ms * clalpha) / (2 * WS_Nm2)

    V_B_intersect = (
        k_g * rho_0 * Ude_B_ms * clalpha
        + np.sqrt((k_g * rho_0 * Ude_B_ms * clalpha) ** 2 + 8 * WS_Nm2 * rho * CLmax_clean)
    ) / (2 * rho * CLmax_clean)
    V_B_stall = V_S1 * np.sqrt(n_cruise_pve)

    V_B = min(V_B_intersect, V_B_stall)

    n_B_pve = 1 + (k_g * rho_0 * Ude_B_ms * V_B * clalpha) / (2 * WS_Nm2)  # Commuter only
    n_B_nve = 1 - (k_g * rho_0 * Ude_B_ms * V_B * clalpha) / (2 * WS_Nm2)

    n_max = max(n_B_pve, n_cruise_pve, n_dive_pve) if commuter_ac else max(n_cruise_pve, n_dive_pve)
    n_min = min(n_B_nve, n_cruise_nve, n_dive_nve) if commuter_ac else min(n_cruise_nve, n_dive_nve)

    return n_max, n_min, V_B, n_cruise_pve, n_cruise_nve, n_dive_pve, n_dive_nve, n_B_pve, n_B_nve


def plot_gust_diagram(
    Vc_ms,
    Vd_ms,
    n_cruise_pve,
    n_cruise_nve,
    n_dive_pve,
    n_dive_nve,
    V_B=None,
    n_B_pve=None,
    n_B_nve=None,
    commuter_ac=False,
    ac_name="aircraft",
    save_plot=True,
    show_plot=True,
):  # pragma: no cover
    """
    Plot the gust diagram

    Parameters
    ----------
    Vc_ms, Vd_ms : float
        The cruise and dive speed [m/s].
    n_cruise_pve, n_cruise_nve, n_dive_pve : float
        The load factors at the respective velocities.
    V_B : float, default None
        The velocity at which the gust diagram intersects the stall speed [m/s].
        Required for commuter aircraft, irrelevant for other aircraft.
    n_B_pve, n_B_nve : float, default None
        The load factors at V_B.
        Required for commuter aircraft, irrelevant for other aircraft.
    commuter_ac : bool, default False
        Boolean indicating whether the aircraft is a commuter aircraft.
    ac_name : str, default "aircraft"
        The name of the aircraft, for saving the plot.
    save_plot : bool, default True
        Boolean indicating whether to save the plot.
    show_plot : bool, default True
        Boolean indicating whether to show the plot.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The figure and axes objects of the plot.
    """

    plt.rcParams.update({"font.size": FONTSIZE})
    fig, ax = plt.subplots()
    fig.tight_layout()
    ax.axhline(linewidth=2, color="k")
    ax.axvline(linewidth=2, color="k")
    style = "r-"
    LIGHTCOLOUR = "dimgrey"

    if commuter_ac:
        ax.plot([0, V_B], [1, n_B_pve], style, linewidth=2, zorder=20)
        ax.plot([0, Vc_ms], [1, n_cruise_pve], color=LIGHTCOLOUR, linestyle="--", linewidth=2, zorder=20)
        ax.plot([V_B, Vc_ms], [n_B_pve, n_cruise_pve], style, linewidth=2, zorder=20)

        ax.plot([0, V_B], [1, n_B_nve], style)
        ax.plot([0, Vc_ms], [1, n_cruise_nve], color=LIGHTCOLOUR, linestyle="--", linewidth=2, zorder=20)
        ax.plot([V_B, Vc_ms], [n_B_nve, n_cruise_nve], style, linewidth=2, zorder=20)
    else:
        ax.plot([0, Vc_ms], [1, n_cruise_pve], style, linewidth=2, zorder=20)
        ax.plot([0, Vc_ms], [1, n_cruise_nve], style, linewidth=2, zorder=20)

    ax.plot([0, Vd_ms], [1, n_dive_pve], color=LIGHTCOLOUR, linestyle="--", linewidth=2, zorder=20)
    ax.plot([0, Vd_ms], [1, n_dive_nve], color=LIGHTCOLOUR, linestyle="--", linewidth=2, zorder=20)

    ax.plot([Vc_ms, Vd_ms], [n_cruise_pve, n_dive_pve], style, linewidth=2, zorder=20)
    ax.plot([Vc_ms, Vd_ms], [n_cruise_nve, n_dive_nve], style, linewidth=2, zorder=20)

    ax.plot([Vd_ms, Vd_ms], [n_dive_nve, n_dive_pve], style, linewidth=2, zorder=20)

    markersize = 15
    if commuter_ac:
        ax.plot(V_B, n_B_pve, "ro", ms=markersize, zorder=20)
        # ax.annotate(
        #     "B'",
        #     (V_B, n_B_pve),
        #     textcoords="offset points",
        #     xytext=(-2, 8),
        #     ha="center",
        #     fontweight="bold",
        #     fontsize=15,
        # )
        ax.plot(V_B, n_B_nve, "ro", ms=markersize, zorder=20)
        # ax.annotate(
        #     "G'",
        #     (V_B, n_B_nve),
        #     textcoords="offset points",
        #     xytext=(-2, -20),
        #     ha="center",
        #     fontweight="bold",
        #     fontsize=15,
        # )

    ax.plot(Vc_ms, n_cruise_pve, "ro", ms=markersize, zorder=20)
    # ax.annotate(
    #     "C'",
    #     (Vc_ms, n_cruise_pve),
    #     textcoords="offset points",
    #     xytext=(2, 8),
    #     ha="center",
    #     fontweight="bold",
    #     fontsize=15,
    # )
    ax.plot(Vc_ms, n_cruise_nve, "ro", ms=markersize, zorder=20)
    # ax.annotate(
    #     "F'",
    #     (Vc_ms, n_cruise_nve),
    #     textcoords="offset points",
    #     xytext=(2, -20),
    #     ha="center",
    #     fontweight="bold",
    #     fontsize=15,
    # )

    ax.plot(Vd_ms, n_dive_pve, "ro", ms=markersize, zorder=20)
    # ax.annotate(
    #     "D'",
    #     (Vd_ms, n_dive_pve),
    #     textcoords="offset points",
    #     xytext=(4, 8),
    #     ha="center",
    #     fontweight="bold",
    #     fontsize=15,
    # )
    ax.plot(Vd_ms, n_dive_nve, "ro", ms=markersize, zorder=20)
    # ax.annotate(
    #     "E'",
    #     (Vd_ms, n_dive_nve),
    #     textcoords="offset points",
    #     xytext=(4, -20),
    #     ha="center",
    #     fontweight="bold",
    #     fontsize=15,
    # )

    ax.plot(0, 1, "ro", ms=markersize, zorder=20)

    if commuter_ac:
        ax.plot([V_B, V_B], [n_B_nve, n_B_pve], linestyle="--", linewidth=2, color=LIGHTCOLOUR, zorder=-1)
        # ax.plot(Vc_ms, 0, marker="o", color="grey")
        ax.annotate(
            "$V_B$",
            (V_B, 0),
            textcoords="offset points",
            xytext=(30, 10),
            ha="center",
            fontsize=FONTSIZE,
            color=LIGHTCOLOUR,
        )

    ax.plot([Vc_ms, Vc_ms], [n_cruise_nve, n_cruise_pve], linestyle="--", linewidth=2, color=LIGHTCOLOUR, zorder=-1)
    # ax.plot(Vc_ms, 0, marker="o", color="grey")
    ax.annotate(
        "$V_C$",
        (Vc_ms, 0),
        textcoords="offset points",
        xytext=(30, 10),
        ha="center",
        fontsize=FONTSIZE,
        color=LIGHTCOLOUR,
    )

    ax.annotate(
        "$V_D$",
        (Vd_ms, 0),
        textcoords="offset points",
        xytext=(-30, 10),
        ha="center",
        fontsize=FONTSIZE,
        color=LIGHTCOLOUR,
    )

    ax.plot([0, Vd_ms], [1, 1], linestyle="--", linewidth=2, color=LIGHTCOLOUR, zorder=-1)

    ax.set_xlabel("Velocity (EAS) [m/s]")
    ax.set_ylabel("Load Factor [-]")

    ax.set_ylim(-2.2, 4.2)
    ax.set_yticks(np.arange(-2.0, 4.2, 1))

    ax.set_xticks(np.arange(0, 90.1, 20))
    ax.set_xticks(np.arange(0, 90.1, 10), minor=True)

    ax.tick_params("both", length=10, width=1, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    ax.grid()

    fig.set_size_inches(16, 9)
    fig.tight_layout()
    if save_plot:
        fig.savefig(
            os.path.join(os.path.dirname(__file__), "..", "..", "Figures", f"vn-gust-{ac_name}.pdf"),
            bbox_inches="tight",
            dpi=200,
        )
    if show_plot:
        plt.show()

    return fig, ax


if __name__ == "__main__":  # pragma: no cover
    Vc_ms, Vd_ms, _, V_S1, _, _ = calculate_manoeuvre_velocities(aircraft_data)

    (
        n_max,
        n_min,
        V_B,
        n_cruise_pve,
        n_cruise_nve,
        n_dive_pve,
        n_dive_nve,
        n_B_pve,
        n_B_nve,
    ) = calculate_gust_diagram_loads(aircraft_data, Vc_ms, Vd_ms, V_S1, COMMUTER, h=3000, temp_offset=18)

    plot_gust_diagram(
        Vc_ms,
        Vd_ms,
        n_cruise_pve,
        n_cruise_nve,
        n_dive_pve,
        n_dive_nve,
        V_B,
        n_B_pve,
        n_B_nve,
        COMMUTER,
        aircraft_data["name"],
        save_plot=False,
        show_plot=True,
    )

    print(f"n_max: {n_max:.2f}")
    print(f"n_min: {n_min:.2f}")
