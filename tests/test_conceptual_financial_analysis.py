import math
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.FinancialAnalysis.conceptual_financial_analysis import hourly_operating_cost


def test_hourly_operating_cost():
    inp = {
        "mission_file": "maf_mission_graph.csv",
        "standard_aircraft_data": {
            "name": "c206",
            "pretty_name": "Cessna 206",
            "Vc_m/s": 73,
            "fuel_burn_L/h": 62.5,
            "fuel_density_kg/L": 0.717,
            "CO2_emissions_kg/kg": 3.05,
            "co2_fuel_%": 89,
            "fuel_per_hour": 144,
            "maintenance_per_hour": 129,
            "overhaul_per_hour": 55,
        },
        "aircraft_data": {
            "name": "final_design",
            "pretty_name": "Final Design",
            "MTOW_N": 15269,
            "OEW_N": 8816,
            "Wpl_des_kg": 540,
            "Wpl_max_kg": 630,
            "CLmax_clean": 1.6,
            "CLmax_land": 2.5,
            "W/S_N/m2": 618,
            "W/P_N/W": 0.118,
            "Vc_m/s": 60,
            "MGC_m": 1.93,
            "CLalpha": 6.24,
            "η_bat": 0.85,
            "DoD_bat": 0.8,
            "CD0": 0.02,
            "AR": 9.38,
            "e": 0.82,
            "Performance": {"Vc_m/s": 60, "range_nm": 600, "CO2_emissions_kg/kg": 3.16},
            "Weights": {"Wfuel_N": 1390},
            "Power_prop": {
                "E_bat_Wh": 186451,
                "P_req_cruise_W": 217309,
                "E_fuel_Wh/kg": 11972,
                "eta_generator": 0.45,
                "eta_electricmotor": 0.925,
                "eta_powertrain": 0.9216,
            },
        },
    }

    assert math.isclose(
        182.7,
        hourly_operating_cost(
            inp["mission_file"],
            inp["standard_aircraft_data"],
            inp["aircraft_data"],
            inp["aircraft_data"]["Weights"]["Wfuel_N"],
        ),
        abs_tol=0.05,
    )
