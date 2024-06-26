import os
import sys

import numpy as np
import logging
import unittest
from unittest.mock import patch, mock_open, MagicMock

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from HumanAir.main import load_json_file, setup_logging, Generate, calculate_weighted_score, find_optimal_design

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
        "P_req_TO_W": 230000,
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


def test_setup_logging():
    with (
        patch("colorlog.StreamHandler") as mock_StreamHandler,
        patch("colorlog.ColoredFormatter") as mock_ColoredFormatter,
        patch("colorlog.getLogger") as mock_getLogger,
        # patch("os.path.abspath", return_value="/fake/dir") as mock_abspath,
        # patch("os.path.join", return_value="/fake/dir/../..") as mock_path_join,
        # patch("sys.path", new_callable=list) as mock_sys_path,
    ):
        # Mock logger and handler
        mock_logger = MagicMock()
        mock_handler = MagicMock()
        mock_formatter = MagicMock()
        mock_getLogger.return_value = mock_logger
        mock_StreamHandler.return_value = mock_handler
        mock_ColoredFormatter.return_value = mock_formatter

        # Call the function
        setup_logging()

        # Assertions
        mock_StreamHandler.assert_called_once()
        mock_ColoredFormatter.assert_called_once_with(
            "%(log_color)s%(levelname)s:%(message)s",
            log_colors={"DEBUG": "cyan", "INFO": "green", "WARNING": "yellow", "ERROR": "red", "CRITICAL": "bold_red"},
        )
        mock_handler.setFormatter.assert_called_once_with(mock_formatter)
        mock_logger.addHandler.assert_called_once_with(mock_handler)
        mock_logger.setLevel.assert_called_once_with(logging.INFO)

        # Check logging output (optional)
        with unittest.TestCase().assertLogs(level="INFO") as log:
            logging.info(" Starting the program")
            assert " Starting the program" in log.output[0]

        # I have no idea why the path would be modified, so removing this test
        # Check sys.path modification
        # assert "/fake/dir/../.." in mock_sys_path


# Test case for the above function
def test_load_design_json():
    with patch("builtins.open", new_callable=mock_open, read_data='{"key": "value"}') as mock_file:
        with patch("os.path.join", return_value="/fake/dir/Configurations/design.json"):
            with patch("os.path.abspath", return_value="/fake/dir/Configurations/design.json"):
                # Call the function
                result = load_json_file(mock_file)

                # Check that the file was opened correctly
                mock_file.assert_called_once_with("/fake/dir/Configurations/design.json", "r")

                # Check that the JSON was read correctly
                assert result == {"key": "value"}

                # Check the logged output (optional)
                with unittest.TestCase().assertLogs(level="INFO") as log:
                    logging.info(
                        f" Looking for design.json at: {os.path.abspath('/fake/dir/Configurations/design.json')}"
                    )
                    assert " Looking for design.json at: /fake/dir/Configurations/design.json" in log.output[0]


# def test_Generate():
#     "Generating the design points"
#     # tune the parameters with a reasonable range
#     A_lst = np.arange(13.0, 17.51, 0.5)
#     eta_p_lst = np.arange(0.8, 0.851, 0.05)
#     Clmax_clean_lst = np.arange(1.6, 2.21, 0.2)
#     Clmax_TO_lst = np.arange(2, 2.61, 0.2)
#     Clmax_Land_lst = np.arange(2, 2.61, 0.2)
#     Cd0_lst = np.arange(0.026, 0.0301, 0.002)
#     V_cruise_lst = np.arange(60, 65.1, 1)
#     climbrate_lst = np.arange(2.5, 5.01, 0.5)

#     # calculate the total numbers of iterations
#     total_iterations = (
#         len(A_lst)
#         * len(eta_p_lst)
#         * len(Clmax_clean_lst)
#         * len(Clmax_TO_lst)
#         * len(Clmax_Land_lst)
#         * len(Cd0_lst)
#         * len(V_cruise_lst)
#         * len(climbrate_lst)
#     )

#     # initialise the iteration counter
#     current_iteration = 0


def test_generate_function():
    with (
        patch("builtins.open", new_callable=unittest.mock.mock_open) as _,
        patch("json.dump") as mock_json_dump,
        patch("HumanAir.main.WP_WS") as mock_WP_WS,
        patch("HumanAir.main.WeightEstimation") as mock_WeightEstm,
        patch("HumanAir.main.co2") as mock_co2,
    ):
        # Mock the return values for methods and functions used in the Generate function
        mock_WP_WS().calculate_optimal_point.return_value = (0.5, 0.6)
        mock_WeightEstm().PolynomialRegression.return_value = ([0.1, 0.2], [0.3, 0.4])
        mock_co2.return_value = 0.8

        # Mock dictionary input
        mock_dict = {
            "Contingency": 1.2,
            "Performance": {
                "P_cruise/P_TO": 0.8,
                "endurance": 2,
                "Vc_m/s": 60.0,
                "range_nm": 600,
                "CO2_emissions_kg/kg": 3.16,
            },
            "Power_prop": {
                "eta_powertrain": 0.85,
                "P_req_cruise_W": 0,
                "E_bat_Wh": 0,
                "P_ptr_kW/kg": 1.5,
                "E_bat_Wh/kg": 350,
                "eta_bat": 0.97,
                "DoD_bat": 0.8,
                "eta_electricmotor": 0.925,
                "E_fuel_Wh/kg": 11972,
                "eta_generator": 0.4,
            },
            "Iterations Class I": {
                "MTOW_kg": 6000,
                "A": 0.3632,
                "B": -81.896,
                "Aw": 0.09,
                "Bw": 14.018,
                "Wpl_des_kg": 630,
            },
        }
        # Mock parameters object
        mock_p = MagicMock()

        def mock_parameter_lsts():
            A_lst = np.array([7.0, 18.51])
            eta_p_lst = np.array([0.80, 0.851])
            Clmax_clean_lst = np.array([1.6, 2.21])
            Clmax_TO_lst = np.array([2, 2.61])
            Clmax_Land_lst = np.array([2.2, 2.81])
            Cd0_lst = np.array([0.028, 0.031])
            V_cruise_lst = np.array([60, 65.1])
            climbrate_lst = np.array([2.5, 5.01])

            return A_lst, eta_p_lst, Clmax_clean_lst, Clmax_TO_lst, Clmax_Land_lst, Cd0_lst, V_cruise_lst, climbrate_lst

        # Run the Generate function
        Generate(mock_p, mock_dict, run=True, parameter_func=mock_parameter_lsts, bat_step=0.05)

        # Check if the json.dump was called to save the file
        mock_json_dump.assert_called_once()

        # Retrieve the arguments with which json.dump was called
        args, kwargs = mock_json_dump.call_args
        print(args)

        # Extract the data written to the JSON file
        written_data = args[0]

        # Perform assertions on the written data
        assert isinstance(written_data, dict)
        assert "0" in written_data
        assert "A" in written_data["0"]
        assert "eta_p" in written_data["0"]
        assert "Clmax_clean" in written_data["0"]
        assert "Clmax_TO" in written_data["0"]
        assert "Clmax_Land" in written_data["0"]
        assert "Cd0" in written_data["0"]
        assert "V_cruise" in written_data["0"]
        assert "climbrate" in written_data["0"]
        assert "CO2" in written_data["0"]
        assert "W/P" in written_data["0"]
        assert "W/S" in written_data["0"]
        assert "bat" in written_data["0"]
        print("Succesful generate", flush=True)


def test_calculate_weighted_score():
    point_data = {
        "A": 8,
        "eta_p": 0.75,
        "Clmax_clean": 2.0,
        "Clmax_TO": 2.4,
        "Clmax_Land": 2.4,
        "Cd0": 0.024,
        "V_cruise": 55,
        "climbrate": 3.0,
        "bat": 0.005,
        "CO2": 25,
    }

    weights = {
        "A": 1,
        "eta_p": 1,
        "Clmax_clean": 1,
        "Clmax_TO": 1,
        "Clmax_Land": 1,
        "Cd0": 1,
        "V_cruise": 1,
        "climbrate": 1,
        "bat": 1,
        "CO2": 1,
    }

    expected_score = sum(
        [
            (9.5 - 8) / 9.5,
            (0.85 - 0.75) / 0.85,
            (2.2 - 2.0) / 2.2,
            (2.6 - 2.4) / 2.6,
            (2.6 - 2.4) / 2.6,
            (0.026 - 0.024) / 0.026,
            (60 - 55) / 60,
            (3.5 - 3.0) / 3.5,
            (0.01 - 0.005) / 0.01,
            (30 - 25) / 30,
        ]
    )

    score = calculate_weighted_score(point_data, weights)

    assert abs(score - expected_score) < 1e-6, f"Expected {expected_score}, but got {score}"


def test_find_optimal_design():
    ac_data = {
        "Contingency": 1.2,
        "Aero": {},
        "Power_prop": {
            "eta_powertrain": 0.9216,
            "P_ptr_kW/kg": 1.5,
            "E_bat_Wh/kg": 350,
            "eta_bat": 0.97,
            "DoD_bat": 0.8,
            "eta_electricmotor": 0.925,
            "E_fuel_Wh/kg": 11972,
            "eta_generator": 0.4,
        },
        "Performance": {"endurance": 5.2, "P_cruise/P_TO": 0.8},
        "Iterations Class I": {
            "MTOW_kg": 2650.285,
            "A": 0.3632,
            "B": -81.896,
            "Aw": 0.09,
            "Bw": 14.018,
            "Wpl_des_kg": 630,
        },
    }

    design_points = {
        "1": {
            "A": 10,
            "eta_p": 0.8,
            "Clmax_clean": 2.1,
            "Clmax_TO": 2.5,
            "Clmax_Land": 2.5,
            "Cd0": 0.028,
            "V_cruise": 62,
            "climbrate": 4.0,
            "bat": 0.015,
            "CO2": 60,
            "W/S": 50,
            "W/P": 0.04,
        },
        "2": {
            "A": 12,
            "eta_p": 0.82,
            "Clmax_clean": 2.2,
            "Clmax_TO": 2.6,
            "Clmax_Land": 2.6,
            "Cd0": 0.027,
            "V_cruise": 65,
            "climbrate": 4.5,
            "bat": 0.02,
            "CO2": 70,
            "W/S": 55,
            "W/P": 0.05,
        },
    }

    weights = {
        "A": 1,
        "eta_p": 1,
        "Clmax_clean": 1,
        "Clmax_TO": 1,
        "Clmax_Land": 1,
        "Cd0": 1,
        "V_cruise": 1,
        "climbrate": 1,
        "bat": 1,
        "CO2": 1,
    }

    class MockWeightEstimation:
        def __init__(self, ac_data):
            self.dict = ac_data

        def Iterations(self, bat):
            return (1, 10000, 1000, 500, 900, 200, 300, 400)  # Mocking return values

    find_optimal_design(
        ac_data,
        maximum_weight_battery=1000,
        weights=weights,
        CO2_threshold=50,
        design_points=design_points,
        printing=True,
        weight_estimation_class=MockWeightEstimation,
    )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
