import numpy as np
import math
import os
import sys

#
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from HumanAir.Class_II_Weight.Class_II_Weight import Class_II_Weight
from HumanAir.unit_conversions import m_to_ft, N_to_lbs, m_squared_to_ft_squared, m_s_to_kt, W_to_hp, lbs_to_N


aircraft_data = {
    "Contingency": 1.2,
    "Contingency_C2W": 1.12,
    "Weights": {
        "MTOW_N": 1000,
        "W_L_N": 200,
        "OEW_N": 100,
        "W_Pilot_N": 20,
        "Wfuel_N": 300,
        "Wbat_N": 0.5,
    },
    "Aero": {
        "S_Wing": 15,
        "AR": 6,
        "AR_HS": 4,
        "AR_v": 6.9,
        "QuarterChordSweep_Wing_deg": 4.2,
        "HalfChordSweep_Wing_deg": 50,
        "QuarterChordSweep_v_deg": 50,
        "Taper_Wing": 0.5,
        "tc_m_Wing": 0.1,
        "S_h": 13,
        "S_v": 11,
        "b_Wing": 89,
        "b_h": 90,
        "b_v": 91,
        "t_root_max_Wing": 121,
        "t_root_max_h": 122,
        "t_root_max_v": 123,
        "Taper_HS": 0.1,
        "QuarterChordSweep_HS_deg": 0,
        "deda": 0.1,
        "VhV": 1,
        "tc_m_HP": 0.1,
    },
    "Performance": {
        "n_ult": 6,
        "n_ult_l": 5.7,
        "Vh_m/s": 5,
        "Vc_m/s": 6,
        "N_pax": 4,
        "M_D": 0.8,
        "W/P_N/W": 0.1,
        "endurance": 4,
        "W/S_N/m2": 1,
        "Altitude_Cruise_m": 1,
        "Temp_offset_TO_Land_cruise": 0,
    },
    "Stability": {
        "QCW_to_QCh": 7,
        "C_L_h": 0.5,
        "C_L_AH": 0.5,
        "X_cg_HS": 0.1,
        "Cg_Aft": 2,
        "Cg_Front": 1,
        "Stability_Margin": 0.1,
    },
    "Geometry": {
        "l_f_nonosecone": 8,
        "fuselage_max_perimeter": 9,
        "fus_length_m": 10,
        "fus_width_m": 11,
        "fus_height_m": 12,
        "N_row": 2,
        "XLEMAC_m": 1,
    },
    "Power_prop": {
        "K_n": 0.24,
        "int": 0.5,
        "N_e": 5.5,
        "N_t": 6.6,
        "int_fueltanks_fraction": 0.5,
        "P_req_TO_W": 4200,
        "E_bat_Wh/kg": 200,
        "eta_bat": 0.8,
        "DoD_bat": 0.7,
        "eta_electricmotor": 0.9,
        "E_fuel_Wh/kg": 1000,
        "eta_generator": 1,
    },
    "Landing_gear": {"l_s_m": 0.7, "l_s_n": 0.7, "Retractable": "Yes", "Hs_m": 0.7},
    "General": {
        "Paint": True,
        "PoweredFlightControls": True,
        "DuplicatedFlightControls": True,
        "N_pax": 3,
        "N_row": 3,
        "APU": False,
        "HydrPne": True,
        "NacWght": True,
    },
    "CL2Weight": {},
    "Iterations Class I": {"Wpl_des_kg": 1},
}

"Structure Tests"


def test_class_II_weight_init():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    # Check if the attributes are correctly initialized
    assert weight_class.W_TO == N_to_lbs(1000) / 1.2  # Contingency
    assert weight_class.W_L == N_to_lbs(200) / 1.2
    assert weight_class.W_F == N_to_lbs(300) / 1.2
    assert weight_class.W_E == N_to_lbs(100) / 1.2  # OEW - Wpilot

    assert weight_class.S_Wing == m_squared_to_ft_squared(15)
    assert weight_class.S_h == m_squared_to_ft_squared(13)
    assert weight_class.S_v == m_squared_to_ft_squared(11)

    assert weight_class.n_ult == 6
    assert weight_class.n_ult_l == 5.7

    assert weight_class.AR_Wing == 6
    assert weight_class.AR_h == 4
    assert weight_class.AR_v == 6.9

    assert weight_class.QuarterChordSweep_Wing == np.deg2rad(4.2)
    assert weight_class.HalfChordSweep_Wing == np.deg2rad(50)
    assert weight_class.QuarterChordSweep_v == np.deg2rad(50)

    assert weight_class.Taper_Wing == 0.5
    assert weight_class.tc_m_Wing == m_to_ft(0.1)

    assert weight_class.b_Wing == m_to_ft(89)
    assert weight_class.b_h == m_to_ft(90)
    assert weight_class.b_v == m_to_ft(91)

    assert weight_class.t_root_max_Wing == m_to_ft(121)
    assert weight_class.t_root_max_h == m_to_ft(122)
    assert weight_class.t_root_max_v == m_to_ft(123)

    assert weight_class.V_H == m_s_to_kt(5)
    assert weight_class.V_c == m_s_to_kt(6)
    assert weight_class.M_D == 0.8
    assert weight_class.QCW_to_QCh == m_to_ft(7)
    assert weight_class.l_f_nonosecone == m_to_ft(8)
    assert weight_class.paint
    assert weight_class.p_max == m_to_ft(9)
    assert weight_class.N_pax == 3
    assert weight_class.N_row == 3
    assert weight_class.l_f == m_to_ft(10)
    assert weight_class.w_f == m_to_ft(11)
    assert weight_class.h_f == m_to_ft(12)

    assert weight_class.P_TO == W_to_hp(4200)
    assert weight_class.K_n == 0.24
    assert weight_class.K_p == 1.1
    assert weight_class.K_pg == 1.16
    assert weight_class.K_fsp == 6.65

    assert weight_class.int == 0.5

    assert weight_class.l_s_m == m_to_ft(0.7)
    assert weight_class.l_s_n == m_to_ft(0.7)
    assert weight_class.retractable == "Yes"

    assert weight_class.N_e == 5.5
    assert weight_class.N_t == 6.6

    assert weight_class.PoweredFlightControls
    assert weight_class.DuplicatedFlightControls


def test_wing_weight():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    # answer1 = 0.002933 * weight_class.S_Wing**1.018 * weight_class.AR_Wing**2.473 * weight_class.n_ult**0.611

    # answer2 = (
    #     96.948
    #     * (
    #         (weight_class.W_TO * weight_class.n_ult / (10**5)) ** 0.65
    #         * (weight_class.AR_Wing / np.cos(weight_class.QuarterChordSweep_Wing)) ** 0.57
    #         * (weight_class.S_Wing / 100) ** 0.61
    #         * ((1 + weight_class.Taper_Wing) / 2 * weight_class.tc_m_Wing) ** 0.36
    #         * (1 + weight_class.V_H / 500) ** 0.5
    #     )
    #     ** 0.993
    # )

    # answer3 = (
    #     0.00125
    #     * weight_class.W_TO
    #     * (weight_class.b_Wing / np.cos(weight_class.HalfChordSweep_Wing)) ** 0.75
    #     * (1 + (6.3 * np.cos(weight_class.HalfChordSweep_Wing) / weight_class.b_Wing) ** 0.5)
    #     * (weight_class.n_ult) ** 0.55
    #     * (
    #         weight_class.b_Wing
    #         * weight_class.S_Wing
    #         / (weight_class.W_TO * weight_class.t_root_max_Wing * np.cos(weight_class.HalfChordSweep_Wing))
    #     )
    #     ** 0.30
    # )

    assert math.isclose(weight_class.WingWeight()["Average"], 760, rel_tol=1e-3)


def test_empennage_weight():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = (3.184 * weight_class.W_TO**0.887 * weight_class.S_h**0.101 * weight_class.AR_h**0.138) / (
        174.04 * weight_class.t_root_max_h**0.223
    ) + (1.68 * weight_class.W_TO**0.567 * weight_class.S_v**0.1249 * weight_class.AR_v**0.482) / (
        639.95 * weight_class.t_root_max_v**0.747 * np.cos(weight_class.QuarterChordSweep_v) ** 0.882
    )
    result2 = (
        127
        * (
            (weight_class.W_TO * weight_class.n_ult / (10**5)) ** 0.87
            * (weight_class.S_h / 100) ** 1.2
            * 0.289
            * (weight_class.QCW_to_QCh / 10) ** 0.483
            * (weight_class.b_h / weight_class.t_root_max_h) ** 0.5
        )
        ** 0.458
        + 98.5
        * (
            (weight_class.W_TO * weight_class.n_ult / (10**5)) ** 0.87
            * (weight_class.S_v / 100) ** 1.2
            * 0.289
            * (weight_class.b_v / weight_class.t_root_max_v) ** 0.5
        )
        ** 0.458
    )
    result3 = 0.04 * (weight_class.n_ult * (weight_class.S_v + weight_class.S_h) ** 2) ** 0.75
    assert math.isclose(weight_class.EmpennageWeight()["Average"], (result1 + result2 + result3) / 3, rel_tol=1e-3)


def test_fuselage_weight():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = (
        14.86
        * weight_class.W_TO**0.144
        * (weight_class.l_f_nonosecone / weight_class.p_max) ** 0.778
        * (weight_class.l_f_nonosecone) ** 0.383
        * weight_class.N_pax**0.455
    )
    result2 = (
        200
        * (
            (weight_class.W_TO * weight_class.n_ult / (10**5)) ** 0.286
            * (weight_class.l_f / 10) ** 0.857
            * ((weight_class.w_f + weight_class.h_f) / 10)
            * (weight_class.V_c / 100) ** 0.338
        )
        ** 1.1
    )

    assert math.isclose(weight_class.FuselageWeight()["Average"], (result1 + result2) / 2, rel_tol=1e-3)


def test_nacelle_weight():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = weight_class.K_n * weight_class.P_TO
    result2 = 2.5 * weight_class.P_TO**0.5

    assert math.isclose(weight_class.NacelleWeight()["Average"], (result1 + result2) / 2, rel_tol=1e-3)


def test_landing_gear_weight(aircraft_data=aircraft_data):
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = (
        0.013 * weight_class.W_TO
        + 0.362 * weight_class.W_L**0.417 * weight_class.n_ult_l**0.95 * weight_class.l_s_m**0.183
        + 6.2
        + 0.0013 * weight_class.W_TO
        + 0.007157 * weight_class.W_L**0.749 * weight_class.n_ult_l * weight_class.l_s_n**0.788
    )
    if weight_class.retractable:
        result1 += 0.014 * weight_class.W_TO

    # result2 = 0.054 * weight_class.l_s_m**0.501 * (weight_class.W_L * weight_class.n_ult_l) ** 0.684

    assert math.isclose(weight_class.LandingGearWeight()["Average"], (result1 + result1) / 2, rel_tol=1e-3)


def test_structure_weight_total():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    assert math.isclose(
        weight_class.StructureWeight_Total(),
        weight_class.WingWeight()["Average"]
        + weight_class.EmpennageWeight()["Average"]
        + weight_class.FuselageWeight()["Average"]
        + weight_class.NacelleWeight()["Average"]
        + weight_class.LandingGearWeight()["Average"],
        rel_tol=1e-2,
    )


"PowerPlant Tests"


def test_fuel_system_weight():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = 0.4 * weight_class.W_F / weight_class.K_fsp
    result2 = (
        2.49
        * (
            (weight_class.W_F / weight_class.K_fsp) ** 0.6
            * (1 / (1 + weight_class.int)) ** 0.3
            * weight_class.N_t**0.2
            * weight_class.N_e**0.13
        )
        ** 1.21
    )
    result3 = (
        2 * (weight_class.W_F / 5.87) ** 0.667
    )  # TODO: I have no idea if we have single piston or multi piston in the end so please verify this(pg.91)

    assert math.isclose(weight_class.FuelSystemWeight()["Average"], (result1 + result2 + result3) / 3, rel_tol=1e-3)


def test_powerplant_weight_total():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = 2.575 * (weight_class.K_p * weight_class.P_TO) ** 0.922 * weight_class.N_e
    result2 = weight_class.K_pg * (weight_class.K_p * weight_class.P_TO + 0.24 * weight_class.P_TO)

    assert math.isclose(weight_class.PowerplantWeight_Total()["USAF"]["WeWaiWpropWp"], result1, rel_tol=1e-3)
    assert math.isclose(weight_class.PowerplantWeight_Total()["Torenbeek"]["Total"], result2, rel_tol=1e-3)


"Fixed Equipment Tests"


def test_flight_control_system():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = 0.0168 * (weight_class.W_TO - weight_class.W_bat)

    if not weight_class.PoweredFlightControls:
        result2 = 1.066 * (weight_class.W_TO - weight_class.W_bat) ** 0.626

        if not weight_class.DuplicatedFlightControls:
            result3 = 0.33 * (weight_class.W_TO - weight_class.W_bat) ** (2 / 3)

    else:
        result2 = 1.08 * (weight_class.W_TO - weight_class.W_bat) ** 0.7
        result3 = False

    assert math.isclose(weight_class.FlightControlSystem()["Cessna"], result1, rel_tol=1e-3)
    assert math.isclose(weight_class.FlightControlSystem()["USAF"], result2, rel_tol=1e-3)
    assert math.isclose(weight_class.FlightControlSystem()["Torenbeek"], result3, rel_tol=1e-3)


def test_hydraulics_pneumatics():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    assert math.isclose(
        weight_class.HydraulicsPneumatics()["Average"],
        (0.006 + 0.0120) / 2 * (weight_class.W_TO - weight_class.W_bat),
        rel_tol=1e-3,
    )


def test_instruments_avionics_electronics():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    assert math.isclose(
        weight_class.InstrumentsAvionicsElectronics()["Average"],
        40 + 0.008 * (weight_class.W_TO - weight_class.W_bat),
        rel_tol=1e-3,
    )


def test_electrical_system_weight():
    # initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = 0.0268 * weight_class.W_TO
    result2 = (
        426
        * (
            (weight_class.InstrumentsAvionicsElectronics()["Average"] + weight_class.FuelSystemWeight()["Average"])
            / 1000
        )
        ** 0.51
    )

    assert math.isclose(weight_class.ElectricalSystemWeight()["Average"], (result1 + result2) / 2, rel_tol=1e-3)


def test_aircon_pressurization_anti_decing_weight():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = (
        0.265
        * weight_class.W_TO**0.52
        * weight_class.N_pax**0.68
        * weight_class.InstrumentsAvionicsElectronics()["Average"] ** 0.17
        * weight_class.M_D**0.08
    )

    if weight_class.N_e <= 1:
        result2 = 2.5 * weight_class.N_pax

    else:
        result2 = 0.018 * weight_class.W_E

    assert math.isclose(
        weight_class.AirconPressurizationAntiDeicingWeight()["Average"], (result1 + result2) / 2, rel_tol=1e-2
    )


def test_oxygen_system():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    assert math.isclose(weight_class.OxygenSystem()["Average"], 20 + 0.5 * weight_class.N_pax, rel_tol=1e-3)


def test_APU_weight():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    if weight_class.APU:
        result1 = 0.0085 * weight_class.W_TO
    else:
        result1 = 0

    assert math.isclose(weight_class.APU_Weight()["Average"], result1, rel_tol=1e-3)


def test_furnishing():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    result1 = 0.412 * weight_class.N_pax**1.145 * (weight_class.W_TO - weight_class.W_bat) ** 0.489
    result2 = 5 + 13 * weight_class.N_pax + 25 * weight_class.N_row

    assert math.isclose(weight_class.Furnishings()["Average"], (result1 + result2) / 2, rel_tol=1e-3)


def test_auxiliary_gear():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    assert math.isclose(weight_class.AuxiliaryGear()["Average"], 0.01 * weight_class.W_E, rel_tol=1e-2)


def test_paint():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    if weight_class.paint:
        result1 = 0.0045 * (weight_class.W_TO - weight_class.W_bat)
    else:
        result1 = 0

    assert math.isclose(weight_class.Paint()["Average"], result1, rel_tol=1e-3)


def test_fixed_equipment_weight_total():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    assert math.isclose(
        weight_class.FixedEquipmentWeight_Total(),
        (
            weight_class.FlightControlSystem()["Average"]
            + weight_class.HydraulicsPneumatics()["Average"]
            + weight_class.InstrumentsAvionicsElectronics()["Average"]
            + weight_class.ElectricalSystemWeight()["Average"]
            + weight_class.AirconPressurizationAntiDeicingWeight()["Average"]
            + weight_class.OxygenSystem()["Average"]
            + weight_class.APU_Weight()["Average"]
            + weight_class.Furnishings()["Average"]
            + weight_class.AuxiliaryGear()["Average"]
            + weight_class.Paint()["Average"]
        ),
        rel_tol=1e-3,
    )


def test_new_empty_weight():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    assert math.isclose(
        lbs_to_N(weight_class.NewEmptyWeight(0)),
        lbs_to_N(
            weight_class.PowerplantWeight_Total()["Average"]
            + weight_class.StructureWeight_Total()
            + weight_class.FixedEquipmentWeight_Total()
        ),
        rel_tol=1,
    )


def test_NewOEW():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)
    bat_percent = 0.1

    assert math.isclose(
        weight_class.NewOEW(bat_percent),
        weight_class.NewEmptyWeight(bat_percent) + aircraft_data["Weights"]["W_Pilot_N"],
        rel_tol=1,
    )


def test_NewFuelWeight():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)
    bat_percent = 1
    assert math.isclose(0, weight_class.NewFuelWeight(bat_percent), rel_tol=1e-4)


def test_PolynomialRegression():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    bat_percent_list = np.arange(0, 1, 0.01)

    coefficent_list_exp = weight_class.PolynomialRegression(bat=bat_percent_list)[0]

    assert math.isclose(coefficent_list_exp[0], 10.64, rel_tol=1e-1)
    assert math.isclose(coefficent_list_exp[1], 12.9, rel_tol=1e-1)


def test_Iterations_C2W():
    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    bat_percent = 0

    assert math.isclose(weight_class.Iterarions_C2W(bat_percent=bat_percent)[0], 0, rel_tol=1e-3)
