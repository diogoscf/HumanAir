# import numpy as np
import json
import os
import matplotlib.pyplot as plt

import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

# from isa import isa

FILE = "design.json"
COMMUTER = True

ft_to_m = 0.3048
FL_to_m = 0.3048 * 100.0
m_to_ft = 1 / 0.3048
m_to_FL = 1 / (0.3048 * 100.0)
lbs_to_kg = 0.45359237
kg_to_lbs = 1 / 0.45359237
nm_to_km = 1.852
km_to_nm = 1 / 1.852

G = 9.80665  # [m/s^2]

# def power_required_W(W_N, MTOW_N, Vc_ms, WS_Nm2, CD0, AR, e, h = 0):
#     P_zerolift = 0.5 * isa(h)[2] * CD0 * MTOW_N * (Vc_ms**3) / WS_Nm2
#     P_induced = 2 * WS_Nm2 * (W_N**2) / (isa(h)[2] * (np.pi * AR * e) * Vc_ms * MTOW_N)
#     return P_zerolift + P_induced


def payload_range_points(filename):
    with open(os.path.join(os.path.dirname(__file__), "..", "Configurations", filename), "r", encoding="utf-8") as f:
        aircraft_data = json.load(f)

    contingency = aircraft_data["Contingency_C2W"]
    MTOW_N = aircraft_data["CL2Weight"]["MTOW_N"] / contingency
    OEW_N = aircraft_data["CL2Weight"]["OEW"] / contingency
    # Wpl_des_kg = aircraft_data["Wpl_des_kg"]
    # Wpl_des_N = Wpl_des_kg * G
    Wpl_des_N = aircraft_data["CL2Weight"]["Wpl_w/o_pilot"] / contingency
    Wpl_max_kg = aircraft_data["Weights"]["Wpl_max_kg"]
    Wpl_max_N = Wpl_max_kg * G
    MFW_N = aircraft_data["CL2Weight"]["Wfuel_N"] * 1.2 / contingency

    design_range_nm = aircraft_data["Performance"]["range_nm"]
    design_range_km = design_range_nm * nm_to_km

    Wf_des_N = MTOW_N - OEW_N - Wpl_des_N
    print(f"Wf_des_N: {Wf_des_N}")

    if Wf_des_N < 0:
        raise ValueError("OEW+Wpl exceeds MTOW.")
    if Wf_des_N > MFW_N:
        raise ValueError("Fuel weight for max payload exceeds MFW.")

    E_bat_available_Wh = (
        aircraft_data["Power_prop"]["E_bat_Wh"] * aircraft_data["Power_prop"]["eta_bat"]
    )  # * aircraft_data["DoD_bat"] * aircraft_data["η_electricmotor"] # uncommment for new version
    WP_NW = aircraft_data["Performance"]["W/P_N/W"]
    P_req_W = MTOW_N / WP_NW

    Vc_ms = aircraft_data["Performance"]["Vc_m/s"]
    Vc_kmh = Vc_ms * 3.6

    zero_fuel_endurance_h = E_bat_available_Wh / P_req_W
    zero_fuel_range_km = zero_fuel_endurance_h * Vc_kmh

    fuel_weight_maxpayload_N = MTOW_N - Wpl_max_N - OEW_N
    print(fuel_weight_maxpayload_N)
    E_fuel_maxpayload_Wh = (
        aircraft_data["Power_prop"]["E_fuel_Wh/kg"]
        * (fuel_weight_maxpayload_N / G)
        * aircraft_data["Power_prop"]["eta_generator"]
    )  # * aircraft_data["η_powertrain"] # uncommment for new version
    maxpayload_endurance_h = (E_bat_available_Wh + E_fuel_maxpayload_Wh) / P_req_W
    print(maxpayload_endurance_h)
    maxpayload_maxrange_km = maxpayload_endurance_h * Vc_kmh
    # print(fuel_weight_maxpayload_N/G, MTOW_N/G, Wpl_max_N/G, OEW_N/G)

    Wpl_maxfuel_N = MTOW_N - OEW_N - MFW_N
    E_fuel_available_Wh = (
        aircraft_data["Power_prop"]["E_fuel_Wh/kg"] * (MFW_N / G) * aircraft_data["Power_prop"]["eta_generator"]
    )  # * aircraft_data["η_powertrain"] # uncommment for new version
    max_fuel_endurance_h = (E_bat_available_Wh + E_fuel_available_Wh) / P_req_W
    max_fuel_range_km = max_fuel_endurance_h * Vc_kmh

    # CD0, AR, e = aircraft_data["CD0"], aircraft_data["AR"], aircraft_data["e"]
    # WS_Nm2 = aircraft_data["W/S_N/m2"]
    ferry_weight_N = OEW_N + MFW_N
    # ferry_power_req_W = power_required_W(ferry_weight_N, MTOW_N, Vc_ms, WS_Nm2, CD0, AR, e)
    ferry_power_req_W = ferry_weight_N / WP_NW
    ferry_endurance_h = (E_bat_available_Wh + E_fuel_available_Wh) / ferry_power_req_W
    ferry_range_km = ferry_endurance_h * Vc_kmh

    # print(f"Theoretical p_req: {P_req_W}")
    # p_req_calc = power_required_W(MTOW_N, MTOW_N, Vc_ms, WS_Nm2, CD0, AR, e, 0)
    # print(f"Calculated  p_req: {p_req_calc}")
    # print(f"Adjusted    p_req: {p_req_calc*((isa(1800)[2]/isa(0)[2])**(3/4))}")

    return (
        zero_fuel_range_km,
        maxpayload_maxrange_km,
        design_range_km,
        max_fuel_range_km,
        ferry_range_km,
        Wpl_max_N,
        Wpl_des_N,
        Wpl_maxfuel_N,
        aircraft_data["pretty_name"],
    )


if __name__ == "__main__":  # pragma: no cover
    plt.rcParams.update({"font.size": 30})
    fig, ax = plt.subplots()
    fig.tight_layout()
    ax.axhline(linewidth=2, color="k")
    ax.axvline(linewidth=2, color="k")
    style = "b-"

    for i, filename in enumerate(["design.json"]):
        # for i, filename in enumerate(["conventional.json"]):
        color = "brg"[i]
        style = color + "-"
        (
            zero_fuel_range_km,
            maxpayload_maxrange_km,
            design_range_km,
            max_fuel_range_km,
            ferry_range_km,
            Wpl_max_N,
            Wpl_des_N,
            Wpl_maxfuel_N,
            pretty_name,
        ) = payload_range_points(filename)

        ax.plot([0, zero_fuel_range_km * km_to_nm], [Wpl_max_N / G, Wpl_max_N / G], f"{color}--")
        ax.plot(
            [zero_fuel_range_km * km_to_nm, maxpayload_maxrange_km * km_to_nm], [Wpl_max_N / G, Wpl_max_N / G], style
        )
        ax.plot(
            [maxpayload_maxrange_km * km_to_nm, max_fuel_range_km * km_to_nm], [Wpl_max_N / G, Wpl_maxfuel_N / G], style
        )
        ax.plot([max_fuel_range_km * km_to_nm, ferry_range_km * km_to_nm], [Wpl_maxfuel_N / G, 0], style)

        ax.plot(zero_fuel_range_km * km_to_nm, Wpl_max_N / G, f"{color}o", ms=20)
        ax.plot(maxpayload_maxrange_km * km_to_nm, Wpl_max_N / G, f"{color}o", ms=20)
        print(maxpayload_maxrange_km * km_to_nm, Wpl_max_N / G)
        # ax.plot(design_range_km, Wpl_des_N, f"{color}o", ms=20)
        print(design_range_km * km_to_nm, Wpl_des_N / G)
        ax.plot(max_fuel_range_km * km_to_nm, Wpl_maxfuel_N / G, f"{color}o", ms=20)
        print(max_fuel_range_km * km_to_nm, Wpl_maxfuel_N / G)
        ax.plot(ferry_range_km * km_to_nm, 0, f"{color}o", ms=20, label=pretty_name)

    ax.plot(design_range_km * km_to_nm, Wpl_des_N / G, "X", mfc="purple", mec="purple", ms=20, label="Design Point")
    ax.plot([0, design_range_km * km_to_nm], [Wpl_des_N / G, Wpl_des_N / G], "--", color="purple")
    ax.plot([design_range_km * km_to_nm, design_range_km * km_to_nm], [0, Wpl_des_N / G], "--", color="purple")

    # ax.plot([0, zero_fuel_range_km], [Wpl_des_N, Wpl_des_N], "b--")
    # ax.plot([zero_fuel_range_km, design_range_km], [Wpl_des_N, Wpl_des_N], style)
    # ax.plot([design_range_km, max_fuel_range_km], [Wpl_des_N, Wpl_maxfuel_N], style)
    # ax.plot([max_fuel_range_km, ferry_range_km], [Wpl_maxfuel_N, 0], style)

    # ax.plot(zero_fuel_range_km, Wpl_des_N, "bo", ms=10)
    # ax.plot(design_range_km, Wpl_des_N, "bo", ms=10)
    # ax.plot(max_fuel_range_km, Wpl_maxfuel_N, "bo", ms=10)
    # ax.plot(ferry_range_km, 0, "bo", ms=10)

    # ax.set_xlabel("Range [km]")
    # ax.set_ylabel("Payload [N]")
    ax.set_xlabel("Range [NM]")
    ax.set_ylabel("Payload [kg]")

    ax.legend(loc="lower left", bbox_to_anchor=(0.01, 0.02))

    ax.grid()

    ax.set_xlim(-20, 1000)

    fig.set_size_inches(20, 7)
    # fig.savefig(
    #     os.path.join(os.path.dirname(__file__), "..", "..", "Figures", f"payload-range-{aircraft_data['name']}.pdf"),
    #     bbox_inches="tight",
    #     dpi=200,
    # )
    fig.savefig(
        os.path.join(os.path.dirname(__file__), "..", "..", "Figures", "payload-range.pdf"),
        bbox_inches="tight",
        dpi=200,
    )
    plt.show()
