import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.Class_I_Weight.Class_I_Weight_Estimation import WeightEstm
import numpy as np
import math

aircraft_data = {
    "name": "final_design",
    "pretty_name": "Final Design",
    "Contingency": 1.2,
    "Weights": {
        "MTOW_N": 25753,
        "W_L_N": 24500,
        "OEW_N": 17651,
        "MFW_N": 1390,
        "Wfuel_N": 1922,
        "Wbat_N": 11772,
        "Ww_N": 5000,
        "Wpl_des_kg": 630,
        "Wpl_max_kg": 720,
        "W_Pilot_N": 882.9,
    },
    "Iterations Class I": {"MTOW_kg": 2650.285, "A": 0.3632, "B": -81.896, "Aw": 0.09, "Bw": 14.018, "Wpl_des_kg": 630},
    "Power_prop": {
        "E_bat_Wh": 209445.55986030292,
        "E_fuel_Wh/kg": 11972,
        "E_bat_Wh/kg": 350,
        "P_ptr_kW/kg": 1.5,
        "eta_bat": 0.97,
        "eta_generator": 0.43,
        "DoD_bat": 0.8,
        "eta_electricmotor": 0.925,
        "eta_powertrain": 0.9216,
        "eta_p": 0.8500000000000001,
        "bat": 0.163,
        "P_req_cruise_W": 197683.3976973128,
        "K_n": 0.24,
        "int_fueltanks_fraction": 0,
        "N_e": 1,
        "N_t": 2,
        "P_req_TO_W": 247104.24712164098,
    },
    "Geometry": {
        "fus_height_m": 1.8,
        "fus_length_m": 10.09,
        "fus_width_m": 1.6,
        "MGC_m": 1.93,
        "XLEMAC_m": 3.7122108719556963,
        "tail_length_m": 4.5,
        "l_f_nonosecone": 9.5,
        "fuselage_max_perimeter": 6.8,
    },
    "Aero": {
        "CLmax_clean": 1.6,
        "CLmax_land": 2.5,
        "CLmax_TO": 2.6000000000000005,
        "CLalpha": 6.24,
        "CD0": 0.027999999999999997,
        "AR": 17.5,
        "AR_HS": 5,
        "AR_v": 2,
        "e": 0.82,
        "Taper_Wing": 0.4,
        "Taper_HS": 0.4,
        "QuarterChordSweep_Wing_deg": 0,
        "QuarterChordSweep_HS_deg": 30,
        "QuarterChordSweep_v_deg": 25,
        "HalfChordSweep_Wing_deg": 3.297,
        "HalfChordSweep_v_deg": 28.23,
        "deda": 0,
        "VhV": 1,
        "S_Wing": 41.67,
        "S_h": 5.1,
        "S_v": 3.168,
        "tc_m_Wing": 0.1367,
        "tc_m_HP": 0.12,
        "b_Wing": 17.61,
        "b_h": 5.051,
        "b_v": 2.518,
        "t_root_max_Wing": 0.46,
        "t_root_max_h": 0.3,
        "t_root_max_v": 0.3,
        "CLmax_Land": 2.0,
    },
    "Performance": {
        "W/S_N/m2": 836.4001000000001,
        "W/P_N/W": 0.1345174630075089,
        "Vc_m/s": 60.0,
        "Vh_m/s": 50,
        "M_D": 0.2,
        "range_nm": 600,
        "endurance": 5.2,
        "climb_rate": 4.5,
        "CO2": 46.89,
        "CO2_emissions_kg/kg": 3.16,
        "Altitude_Cruise_m": 3000,
        "Altitude_TO_m": 1800,
        "Altitude_Land_m:": 1800,
        "n_ult": 3.67,
        "P_cruise/P_TO": 0.8,
    },
    "Landing_gear": {
        "lm_m": 0.4137773297250364,
        "ln_m": 4.154111419961311,
        "ymin_m": 1.1454554584955832,
        "Hs_m": 0.2671098386386166,
        "Dwm_m": 0.551815,
        "Dwn_m": 0.329565,
        "Wtm_m": 0.1651,
        "Wtn_m": 0.127,
        "l_s_m": 1,
        "l_s_n": 1,
        "Retractable": True,
        "PMW_N": 11710.095983139854,
        "PNW_N": 2332.808033720294,
        "Xnw_m": 0.2,
        "Xmw_m": 4.7679236042696385,
    },
    "Stability": {
        "Cg_Aft": 0.45308741224228294,
        "Cg_Front": 0.27548662213803016,
        "Stability_Margin": 0.05,
        "XLEMAC_m": 3.7122108719556963,
        "X_cg_HS": 0.95,
        "C_L_h": -0.5,
        "C_L_AH": 1.72,
        "QCW_to_QCh": 6.9,
    },
    "General": {
        "Paint": True,
        "N_pax": 7,
        "N_row": 3,
        "PoweredFlightControls": True,
        "DuplicatedFlightControls": True,
        "APU": True,
    },
}


def test_Class_I_Weight_init():
    # Initialize the Class_I weight estimation object:
    class_I = WeightEstm(aircraft_data)

    # Check initialization:
    assert class_I.dict == {
        "name": "final_design",
        "pretty_name": "Final Design",
        "Contingency": 1.2,
        "Weights": {
            "MTOW_N": 25753,
            "W_L_N": 24500,
            "OEW_N": 17651,
            "MFW_N": 1390,
            "Wfuel_N": 1922,
            "Wbat_N": 11772,
            "Ww_N": 5000,
            "Wpl_des_kg": 630,
            "Wpl_max_kg": 720,
            "W_Pilot_N": 882.9,
        },
        "Iterations Class I": {
            "MTOW_kg": 2650.285,
            "A": 0.3632,
            "B": -81.896,
            "Aw": 0.09,
            "Bw": 14.018,
            "Wpl_des_kg": 630,
        },
        "Power_prop": {
            "E_bat_Wh": 209445.55986030292,
            "E_fuel_Wh/kg": 11972,
            "E_bat_Wh/kg": 350,
            "P_ptr_kW/kg": 1.5,
            "eta_bat": 0.97,
            "eta_generator": 0.43,
            "DoD_bat": 0.8,
            "eta_electricmotor": 0.925,
            "eta_powertrain": 0.9216,
            "eta_p": 0.8500000000000001,
            "bat": 0.163,
            "P_req_cruise_W": 197683.3976973128,
            "K_n": 0.24,
            "int_fueltanks_fraction": 0,
            "N_e": 1,
            "N_t": 2,
            "P_req_TO_W": 247104.24712164098,
        },
        "Geometry": {
            "fus_height_m": 1.8,
            "fus_length_m": 10.09,
            "fus_width_m": 1.6,
            "MGC_m": 1.93,
            "XLEMAC_m": 3.7122108719556963,
            "tail_length_m": 4.5,
            "l_f_nonosecone": 9.5,
            "fuselage_max_perimeter": 6.8,
        },
        "Aero": {
            "CLmax_clean": 1.6,
            "CLmax_land": 2.5,
            "CLmax_TO": 2.6000000000000005,
            "CLalpha": 6.24,
            "CD0": 0.027999999999999997,
            "AR": 17.5,
            "AR_HS": 5,
            "AR_v": 2,
            "e": 0.82,
            "Taper_Wing": 0.4,
            "Taper_HS": 0.4,
            "QuarterChordSweep_Wing_deg": 0,
            "QuarterChordSweep_HS_deg": 30,
            "QuarterChordSweep_v_deg": 25,
            "HalfChordSweep_Wing_deg": 3.297,
            "HalfChordSweep_v_deg": 28.23,
            "deda": 0,
            "VhV": 1,
            "S_Wing": 41.67,
            "S_h": 5.1,
            "S_v": 3.168,
            "tc_m_Wing": 0.1367,
            "tc_m_HP": 0.12,
            "b_Wing": 17.61,
            "b_h": 5.051,
            "b_v": 2.518,
            "t_root_max_Wing": 0.46,
            "t_root_max_h": 0.3,
            "t_root_max_v": 0.3,
            "CLmax_Land": 2.0,
        },
        "Performance": {
            "W/S_N/m2": 836.4001000000001,
            "W/P_N/W": 0.1345174630075089,
            "Vc_m/s": 60.0,
            "Vh_m/s": 50,
            "M_D": 0.2,
            "range_nm": 600,
            "endurance": 5.2,
            "climb_rate": 4.5,
            "CO2": 46.89,
            "CO2_emissions_kg/kg": 3.16,
            "Altitude_Cruise_m": 3000,
            "Altitude_TO_m": 1800,
            "Altitude_Land_m:": 1800,
            "n_ult": 3.67,
            "P_cruise/P_TO": 0.8,
        },
        "Landing_gear": {
            "lm_m": 0.4137773297250364,
            "ln_m": 4.154111419961311,
            "ymin_m": 1.1454554584955832,
            "Hs_m": 0.2671098386386166,
            "Dwm_m": 0.551815,
            "Dwn_m": 0.329565,
            "Wtm_m": 0.1651,
            "Wtn_m": 0.127,
            "l_s_m": 1,
            "l_s_n": 1,
            "Retractable": True,
            "PMW_N": 11710.095983139854,
            "PNW_N": 2332.808033720294,
            "Xnw_m": 0.2,
            "Xmw_m": 4.7679236042696385,
        },
        "Stability": {
            "Cg_Aft": 0.45308741224228294,
            "Cg_Front": 0.27548662213803016,
            "Stability_Margin": 0.05,
            "XLEMAC_m": 3.7122108719556963,
            "X_cg_HS": 0.95,
            "C_L_h": -0.5,
            "C_L_AH": 1.72,
            "QCW_to_QCh": 6.9,
        },
        "General": {
            "Paint": True,
            "N_pax": 7,
            "N_row": 3,
            "PoweredFlightControls": True,
            "DuplicatedFlightControls": True,
            "APU": True,
        },
    }


def test_OEW_prime():
    class_I = WeightEstm(aircraft_data)
    A = 0.3632
    MTOW_kg = 2650.285
    B = -81.896
    assert math.isclose(class_I.OEW_prime(), A * MTOW_kg + B, rel_tol=1e-3)


def test_PowertrainWeight():
    class_I = WeightEstm(aircraft_data)
    MTOW = 2650.285 * 9.81
    WP = 0.1345174630075089
    p_bat = 0.163
    eta_powertrain = 0.9216  # = eta_EM * eta
    p_ptr = 1.5 * 1000  # kW /kg

    assert math.isclose(
        class_I.PowertrainWeight(p_bat), (1 - p_bat) * (MTOW / WP) / (eta_powertrain * p_ptr), rel_tol=1e-3
    )


def test_BatteryWeight():
    class_I = WeightEstm(aircraft_data)
    p_bat = 0.163  # percentage of total aircraft energy provided by battery
    WP = 0.1345174630075089
    MTOW = 2650.285 * 9.81
    E = 5.2
    e_bat = 350
    eta_bat = 0.97
    eta_electricmotor = 0.925
    DOD_bat = 0.8

    assert math.isclose(
        class_I.BatteryWeight(p_bat),
        MTOW / WP * E * p_bat / (e_bat * eta_bat * eta_electricmotor * DOD_bat),
        rel_tol=1e-3,
    )


def test_FuelWeight():
    class_I = WeightEstm(aircraft_data)

    assert math.isclose(
        class_I.FuelWeight(class_I.dict["Power_prop"]["bat"]),
        1.15
        * (1 - class_I.dict["Power_prop"]["bat"])
        * 9.81
        * class_I.dict["Iterations Class I"]["MTOW_kg"]
        / class_I.dict["Performance"]["W/P_N/W"]
        * class_I.dict["Performance"]["endurance"]
        / (class_I.dict["Power_prop"]["eta_generator"] * class_I.dict["Power_prop"]["E_fuel_Wh/kg"]),
        rel_tol=1e-3,
    )


# 1.15 since MAF requires 15% reserve fuel.


def test_WingWeight():
    class_I = WeightEstm(aircraft_data)

    assert math.isclose(
        class_I.WingWeight(),
        class_I.dict["Iterations Class I"]["Aw"] * class_I.dict["Iterations Class I"]["MTOW_kg"]
        + class_I.dict["Iterations Class I"]["Bw"],
        rel_tol=1e-3,
    )


def test_Iterations():
    class_I = WeightEstm(aircraft_data)

    MTOW_new = 0
    MTOW_old = class_I.dict["Iterations Class I"]["MTOW_kg"]
    ok = False

    while (
        np.abs(
            (MTOW_new - class_I.dict["Iterations Class I"]["MTOW_kg"]) / class_I.dict["Iterations Class I"]["MTOW_kg"]
        )
        > 0.02
    ):
        if ok:
            class_I.dict["Iterations Class I"]["MTOW_kg"] = MTOW_new

        OEW_prime = class_I.OEW_prime()
        PowertrainWeight = class_I.PowertrainWeight(class_I.dict["Power_prop"]["bat"])
        BatteryWeight = class_I.BatteryWeight(class_I.dict["Power_prop"]["bat"])
        FuelWeight = class_I.FuelWeight(class_I.dict["Power_prop"]["bat"])
        WingWeight = class_I.WingWeight()

        MTOW_new = (
            OEW_prime
            + PowertrainWeight
            + BatteryWeight
            + FuelWeight
            + WingWeight
            + class_I.dict["Iterations Class I"]["Wpl_des_kg"]
        )

        if MTOW_new > 8000:
            break

        ok = True

    if MTOW_new < 4000:
        class_I.dict["Iterations Class I"]["MTOW_kg"] = MTOW_old
        assert np.allclose(
            class_I.Iterations(class_I.dict["Power_prop"]["bat"]),
            (
                MTOW_new,
                class_I.dict["Contingency"] * MTOW_new,
                class_I.dict["Contingency"] * OEW_prime,
                class_I.dict["Contingency"] * PowertrainWeight,
                class_I.dict["Contingency"] * BatteryWeight,
                class_I.dict["Contingency"] * FuelWeight,
                class_I.dict["Contingency"] * WingWeight,
                class_I.dict["Contingency"] * class_I.dict["Iterations Class I"]["Wpl_des_kg"],
            ),
            rtol=1e-3,
        )
        assert math.isclose(class_I.dict["Iterations Class I"]["MTOW_kg"], MTOW_old, rel_tol=1e-3)

    else:
        class_I.dict["Iterations Class I"]["MTOW_kg"] = MTOW_old
        assert np.allclose(
            class_I.Iterations(class_I.dict["Power_prop"]["bat"]),
            (
                0,
                class_I.dict["Contingency"] * MTOW_new,
                class_I.dict["Contingency"] * OEW_prime,
                class_I.dict["Contingency"] * PowertrainWeight,
                class_I.dict["Contingency"] * BatteryWeight,
                class_I.dict["Contingency"] * FuelWeight,
                class_I.dict["Contingency"] * WingWeight,
                class_I.dict["Contingency"] * class_I.dict["Iterations Class I"]["Wpl_des_kg"],
            ),
            rtol=1e-3,
        )
        assert math.isclose(class_I.dict["Iterations Class I"]["MTOW_kg"], MTOW_old, rel_tol=1e-3)
    # The reason for the if statement, is that if the design weighs less than 4000,
    # you actually care about the effect of the contingency, and how much the aircaft could
    # have weighed if the contingency was not applied.
    # If the weight is too large, it does not matter that much.


def test_PolynomialRegression():
    class_I = WeightEstm(aircraft_data)
    bat = np.arange(0, 0.18, 0.001)

    lst_P = []
    lst_bat = []

    for pbat in bat:
        row = class_I.Iterations(pbat)  # this is a tuple, not a row
        if row[0] != 0:
            lst_P.append(9.81 * row[1] / class_I.dict["Performance"]["W/P_N/W"])
            lst_bat.append(pbat)

    lst_bat = np.array(lst_bat)
    lst_P = np.array(lst_P)

    # Filter out non-positive values in lst_P
    valid_indices = lst_P > 0
    lst_bat = lst_bat[valid_indices]
    lst_P = lst_P[valid_indices]

    if len(lst_P) == 0:
        print([class_I.PolynomialRegression(bat)][0])
        assert np.allclose([class_I.PolynomialRegression(bat)][0][0], np.array([20, 20]), rtol=1e-3)
        assert np.allclose([class_I.PolynomialRegression(bat)][0][1], np.array([20, 20]), rtol=1e-3)

    coeff_exp = np.polyfit(lst_bat, np.log(lst_P), 1)  # assuming log of exponential is linear
    coeff_pol = np.polyfit(lst_bat, lst_P, 2)  # assuming second degree poly fit.

    # y_pol = coeff_pol[0] * lst_bat**2 + coeff_pol[1] * lst_bat + coeff_pol[2]
    # y_exp = np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * lst_bat)

    assert np.allclose([class_I.PolynomialRegression(bat)][0][0], coeff_exp, rtol=1e-3)
    assert np.allclose([class_I.PolynomialRegression(bat)][0][1], coeff_pol, rtol=1e-3)


if __name__ == "__main__":  # pragma: no cover
    test_Class_I_Weight_init()
    test_OEW_prime()
    test_BatteryWeight()
    test_FuelWeight()
    test_WingWeight()
    test_Iterations()
    test_PolynomialRegression()
    test_PowertrainWeight()
