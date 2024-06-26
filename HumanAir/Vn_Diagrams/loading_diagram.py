import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import warnings

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.isa import isa
from HumanAir.aircraft_data import aircraft_data

# PLOT_FLAPPED = False

FONT_SIZE = 40


def nmax_manoeuvre(MTOW_lbs):
    return min(2.1 + 24000 / (MTOW_lbs + 10000), 3.8)


def nmin_manoeuvre(MTOW_lbs):
    return -0.4 * nmax_manoeuvre(MTOW_lbs)


ft_to_m = 0.3048
FL_to_m = 0.3048 * 100.0
m_to_ft = 1 / 0.3048
m_to_FL = 1 / (0.3048 * 100.0)
lbs_to_kg = 0.45359237
kg_to_lbs = 1 / 0.45359237

G = 9.80665  # [m/s^2]


def calc_nmax_nmin_manoeuvre(MTOW_N):
    """
    Calculate the maximum and minimum load factors for the loading diagram.

    Parameters
    ----------
    MTOW_N : float
        The aircraft maximum take-off weight in Newton.

    Returns
    -------
    nmax, nmin: float
        The maximum and minimum load factors.
    """

    # MTOW_N = aircraft_data["MTOW_N"]
    MTOW_kg = MTOW_N / G
    MTOW_lbs = MTOW_kg * kg_to_lbs

    nmax = nmax_manoeuvre(MTOW_lbs)
    nmin = nmin_manoeuvre(MTOW_lbs)

    return nmax, nmin


def calculate_manoeuvre_velocities(
    aircraft_data,
    MTOW_N=None,
    WS_Nm2=None,
    CLmax_clean=None,
    CLmax_land=None,
    Vc_ms=None,
    rho=None,
    h=None,
    temp_offset=None,
):
    """
    Calculate the applicable velocities for the loading (manoeuvre) diagram.

    Parameters
    ----------
    aircraft_data : dict
        The aircraft data dictionary.
    MTOW_N : float, optional
        The aircraft maximum take-off weight in Newton. Taken from `aircraft_data` if not provided.
    WS_Nm2 : float, optional
        The wing loading (N/m^2). Taken from `aircraft_data` if not provided.
    CLmax_clean : float, optional
        The maximum lift coefficient in clean configuration. Taken from `aircraft_data` if not provided.
    CLmax_land : float, optional
        The maximum lift coefficient in landing configuration. Taken from `aircraft_data` if not provided.
    Vc_ms : float, optional
        The cruising speed in m/s. Taken from `aircraft_data` if not provided.
    rho : float, optional
        The air density. If not provided, it will be calculated using `h` or by default at sea-level.
    h : float, optional
        The altitude in meters. If not provided, it will be set to sea level.

    Returns
    -------
    Vc_ms, Vd_ms, V_A, V_S1, V_HH, V_S0
        The cruise, dive, manoeuvring, stall clean, point H (on diagram), and stall flapped speeds.
    """

    MTOW_N = MTOW_N if MTOW_N is not None else aircraft_data["Weights"]["MTOW_N"]

    nmax, nmin = calc_nmax_nmin_manoeuvre(MTOW_N)

    if rho is None:
        if temp_offset is not None:
            _, _, rho = isa(h, temp_offset) if h is not None else isa(0, temp_offset)
        else:
            _, _, rho = isa(h) if h is not None else isa(0)

    WS_Nm2 = WS_Nm2 if WS_Nm2 is not None else aircraft_data["Performance"]["W/S_N/m2"]
    WS_lbft2 = WS_Nm2 / G * kg_to_lbs / (m_to_ft**2)
    # print(f"W/S: {WS_Nm2} N/m^2 ({WS_lbft2:.3f} lb/ft^2)")
    # print(WS_lbft2)
    # print(33 * np.sqrt(WS_lbft2))

    Vc_ms = Vc_ms if Vc_ms is not None else aircraft_data["Performance"]["Vc_m/s"]

    # CS-23.335(b)(2)(i)
    k_vc_vd = 1.4
    if WS_lbft2 > 20:
        k_vc_vd = np.interp(WS_lbft2, [20, 100], [1.4, 1.35]) if WS_lbft2 < 100 else 1.35
    Vd_ms = k_vc_vd * Vc_ms
    # print(f"V_D: {Vd_ms:.2f} m/s")

    CLmax_clean = CLmax_clean if CLmax_clean is not None else aircraft_data["Aero"]["CLmax_clean"]
    CLmax_land = CLmax_land if CLmax_land is not None else aircraft_data["Aero"]["CLmax_Land"]

    V_A = np.sqrt(2 * nmax * WS_Nm2 / (rho * CLmax_clean))
    V_S1 = np.sqrt(2 * WS_Nm2 / (rho * CLmax_clean))
    V_HH = np.sqrt(
        2 * (-nmin) * WS_Nm2 / (rho * CLmax_clean)
    )  # NOTE: This is NOT!!! the same as V_H (which is the max level speed at sea level)

    if V_A > Vc_ms:
        V_A = Vc_ms
        warnings.warn("V_A is higher than V_C. Setting V_A = V_C.")

    V_S0 = np.sqrt(2 * WS_Nm2 / (rho * CLmax_land))

    # NOTE: Part below is for CS-25 aircraft that need to consider loading with flaps
    # nmax_flap = 2  # TODO: Check
    # V_I = np.sqrt(2 * 2 * WS_Nm2 / (rho * CLmax_land))
    # n_OI = q * CLmax_land / WS_Nm2
    # n_OI = np.concatenate((n_OI[np.where(n_OI < nmax_flap)], [nmax_flap]))  # needs to be changed
    # V_OI = np.concatenate((V_ms[: np.shape(n_OI)[0] - 1], [V_I]))

    # V_J = np.sqrt(2 * nmax_flap * WS_Nm2 / (rho * CLmax_clean))

    return Vc_ms, Vd_ms, V_A, V_S1, V_HH, V_S0


def plot_manoeuvre_diagram(
    Vc_ms,
    Vd_ms,
    V_A,
    V_S1,
    V_HH,
    V_S0,
    nmax,
    nmin,
    h,
    temp_offset,
    WS_Nm2,
    CLmax_clean,
    ac_name="aircraft",
    save_plot=False,
    show_plot=True,
):  # pragma: no cover
    """
    Plot the loading diagram.

    Parameters
    ----------
    Vc_ms, Vd_ms, V_A, V_S1, V_HH, V_S0 : float
        The applicable velocities.
    nmax, nmin : float
        The maximum and minimum load factors.
    h : float
        The altitude the aircraft flies at.
    WS_Nm2 : float
        The wing loading (N/m^2).
    CLmax_clean : float
        The maximum lift coefficient in clean configuration.
    ac_name : str, default "aircraft"
        The aircraft name used in the file name for the saved plot.
    save_plot : bool, default False
        Whether to save the plot or not.
    show_plot : bool, default True
        Whether to show the plot or not.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The figure and axes objects.
    """

    rho = isa(h, temp_offset)[2]

    V_ms = np.linspace(0, np.ceil(Vd_ms), 1000)
    q = 0.5 * rho * V_ms**2

    n_OA = q * CLmax_clean / WS_Nm2
    n_OA = np.concatenate((n_OA[np.where(V_ms < V_A)], [nmax]))
    V_OA = np.concatenate((V_ms[: np.shape(n_OA)[0] - 1], [V_A]))

    # NOTE: Part below is for CS-25 aircraft that need to consider loading with flaps
    # nmax_flap = 2  # TODO: Check
    # V_I = np.sqrt(2 * 2 * WS_Nm2 / (rho * CLmax_land))
    # n_OI = q * CLmax_land / WS_Nm2
    # n_OI = np.concatenate((n_OI[np.where(n_OI < nmax_flap)], [nmax_flap]))  # needs to be changed
    # V_OI = np.concatenate((V_ms[: np.shape(n_OI)[0] - 1], [V_I]))

    # V_J = np.sqrt(2 * nmax_flap * WS_Nm2 / (rho * CLmax_clean))
    # n_IJ = [nmax_flap, nmax_flap]
    # V_IJ = [V_I, V_J]

    n_AD = [nmax, nmax]
    V_AD = [V_A, Vd_ms]

    n_OH = -n_OA[np.where(n_OA < -nmin)]
    n_OH = np.concatenate((n_OH, [nmin]))
    V_OH = np.concatenate((V_ms[: np.shape(n_OH)[0] - 1], [V_HH]))

    n_HF = [nmin, nmin]
    V_HF = [V_HH, Vc_ms]

    n_FE = [nmin, 0]
    V_FE = [Vc_ms, Vd_ms]

    n_ED = [0, nmax]
    V_ED = [Vd_ms, Vd_ms]

    plt.rcParams.update({"font.size": FONT_SIZE})
    fig, ax = plt.subplots()
    fig.tight_layout()
    ax.axhline(linewidth=2, color="k")
    ax.axvline(linewidth=2, color="k")

    style = "r-"
    ax.plot(V_OA, n_OA, style, linewidth=2, zorder=20)
    ax.plot(V_AD, n_AD, style, linewidth=2, zorder=20)
    ax.plot(V_OH, n_OH, style, linewidth=2, zorder=20)
    ax.plot(V_HF, n_HF, style, linewidth=2, zorder=20)
    ax.plot(V_FE, n_FE, style, linewidth=2, zorder=20)
    ax.plot(V_ED, n_ED, style, linewidth=2, zorder=20)

    markersize = 15
    # if PLOT_FLAPPED:
    #     ax.plot(V_OI, n_OI, style, linewidth=2, zorder=20)
    #     ax.plot(V_IJ, n_IJ, style, linewidth=2, zorder=20)

    ax.plot(V_A, nmax, "ro", ms=markersize)
    # ax.annotate(
    #     "A",
    #     (V_A, nmax),
    #     textcoords="offset points",
    #     xytext=(-2, 8), ha="center",
    #     fontweight="bold",
    #     fontsize=FONT_SIZE
    # )

    ax.plot(Vd_ms, nmax, "ro", ms=markersize)
    # ax.annotate(
    #     "D",
    #     (Vd_ms, nmax),
    #     textcoords="offset points",
    #     xytext=(2, 8),
    #     ha="center",
    #     fontweight="bold",
    #     fontsize=FONT_SIZE,
    # )

    ax.plot(Vd_ms, 0, "ro", ms=markersize)
    # ax.annotate(
    #     "E", (Vd_ms, 0), textcoords="offset points", xytext=(8, 8), ha="center", fontweight="bold", fontsize=FONT_SIZE
    # )

    ax.plot(Vc_ms, nmin, "ro", ms=markersize)
    # ax.annotate(
    #     "F",
    #     (Vc_ms, nmin),
    #     textcoords="offset points",
    #     xytext=(2, -20),
    #     ha="center",
    #     fontweight="bold",
    #     fontsize=FONT_SIZE,
    # )

    ax.plot(V_HH, nmin, "ro", ms=markersize)
    # ax.annotate(
    #     "H",
    #     (V_H, nmin),
    #     textcoords="offset points",
    #     xytext=(-2, -20),
    #     ha="center",
    #     fontweight="bold",
    #     fontsize=FONT_SIZE,
    # )

    ax.plot(0, 0, "ro", ms=markersize)

    # if PLOT_FLAPPED:
    #     ax.plot(V_I, nmax_flap, "ro", ms=markersize)
    #     # ax.annotate(
    #     #     "I",
    #     #     (V_I, nmax_flap),
    #     #     textcoords="offset points",
    #     #     xytext=(-2, 8),
    #     #     ha="center",
    #     #     fontweight="bold",
    #     #     fontsize=15,
    #     # )

    #     ax.plot(V_J, nmax_flap, "ro", ms=markersize)
    #     # ax.annotate(
    #     #     "J",
    #     #     (V_J, nmax_flap),
    #     #     textcoords="offset points",
    #     #     xytext=(-2, 8),
    #     #     ha="center",
    #     #     fontweight="bold",
    #     #     fontsize=15,
    #     # )

    LIGHTCOLOUR = "dimgrey"

    ax.plot([Vc_ms, Vc_ms], [nmin, 0], linewidth=2, linestyle="--", color=LIGHTCOLOUR, zorder=-1)
    # ax.plot(Vc_ms, 0, marker="o", color="grey")
    ax.annotate(
        "$V_C$",
        (Vc_ms, 0),
        textcoords="offset points",
        xytext=(0, 11),
        ha="center",
        fontsize=FONT_SIZE,
        color=LIGHTCOLOUR,
    )

    ax.annotate(
        "$V_D$",
        (Vd_ms, 0),
        textcoords="offset points",
        xytext=(-38, 11),
        ha="center",
        fontsize=FONT_SIZE,
        color=LIGHTCOLOUR,
    )

    ax.plot([V_A, V_A], [0, nmax], linewidth=2, linestyle="--", color=LIGHTCOLOUR, zorder=-1)
    # ax.plot(Vc_ms, 0, marker="o", color="grey")
    ax.annotate(
        "$V_A$",
        (V_A, 0),
        textcoords="offset points",
        xytext=(0, -38),
        ha="center",
        fontsize=FONT_SIZE,
        color=LIGHTCOLOUR,
    )

    ax.plot([V_S1, V_S1], [-1, 1], linewidth=2, linestyle="--", color=LIGHTCOLOUR, zorder=-1)
    # ax.plot(V_S1, 0, marker="o", color="grey")
    ax.annotate(
        "$V_{S1}$",
        (V_S1, 0),
        textcoords="offset points",
        xytext=(40, 11),
        ha="center",
        fontsize=FONT_SIZE,
        color=LIGHTCOLOUR,
    )

    ax.annotate(
        f"$n_{{max}} = {nmax:.2f}$",
        (V_A + (Vd_ms - V_A) / 2, nmax),
        textcoords="offset points",
        xytext=(0, -38),
        ha="center",
        fontsize=FONT_SIZE,
    )
    ax.annotate(
        f"$n_{{min}} = {nmin:.2f}$",
        (V_HH + (Vc_ms - V_HH) / 2, nmin),
        textcoords="offset points",
        xytext=(0, 10),
        ha="center",
        fontsize=FONT_SIZE,
    )

    ax.plot([0, Vd_ms], [1, 1], linewidth=2, linestyle="--", color=LIGHTCOLOUR, zorder=5)

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
            os.path.join(os.path.dirname(__file__), "..", "..", "Figures", f"vn-loading-{ac_name}.pdf"),
            bbox_inches="tight",
            dpi=200,
        )
    if show_plot:
        plt.show()

    return fig, ax


if __name__ == "__main__":  # pragma: no cover
    h = 0  # [m]
    temp_offset = 0  # [deg C]

    nmax, nmin = calc_nmax_nmin_manoeuvre(aircraft_data["Weights"]["MTOW_N"])
    Vc_ms, Vd_ms, V_A, V_S1, V_HH, V_S0 = calculate_manoeuvre_velocities(aircraft_data, h=h, temp_offset=temp_offset)

    plot_manoeuvre_diagram(
        Vc_ms,
        Vd_ms,
        V_A,
        V_S1,
        V_HH,
        V_S0,
        nmax,
        nmin,
        h,
        temp_offset,
        aircraft_data["Performance"]["W/S_N/m2"],
        aircraft_data["Aero"]["CLmax_clean"],
        ac_name=aircraft_data["name"],
        save_plot=False,
        show_plot=True,
    )

    print("V_S1:", V_S1)
    print("V_S0:", V_S0)
    print("V_D:", Vd_ms)
    print("V_A:", V_A)
    print("nmax:", nmax)
    print("nmin:", nmin)
