# import numpy as np
import sys
import os
import math

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data, c206_data
from HumanAir.CO2_Calculator.conceptual_co2 import _calculate_new_co2_flight_lengths, calculate_mission_freqs
from HumanAir.unit_conversions import nm_to_m

vol_jet_a1_price = 0.8  # US$/L
jet_a1_dens = 0.8025  # kg/m3


def hourly_operating_cost(mission_file, standard_aircraft_data=c206_data, ac_data=aircraft_data, fuel_weight=None):
    """
    Calculate the hourly operating cost of a new aircraft design based on the mission range. Overhaul taken
    from MAF values, so not yet dependent on the mission profile legs. Should be improved later.
    Maintenance _is_ dependent on flight time already though. Fuel cost based on inaccurate estimate of
    fuel burn (assuming cruise speed for all 600 nm), improved fuel burn estimate shall be used.

    Parameters
    ----------
    mission_file : str
        The name of the mission file to read.
    standard_ac_data : dict
        The standard aircraft data to use for the calculation. Will use C206 data by default.
    ac_data : dict
        The aircraft data to use for the calculation. Will use main new design data by default.

    Returns
    -------
    total_hourly_cost : float
        The total new aircraft operating cost per hour of flight in (US$)
    """
    overhaul_cost = (
        standard_aircraft_data["overhaul_per_hour"]
        * standard_aircraft_data["Vc_m/s"]
        / ac_data["Performance"]["Vc_m/s"]
    )

    mission_freqs = calculate_mission_freqs(mission_file)
    _, maintenance_cost = _calculate_new_co2_flight_lengths(
        mission_freqs,
        ac_data=ac_data,
        maintenance_standard_co2=None,
        V_standard_kts=None,
        standard_ac_data=standard_aircraft_data,
    )  # NOTE: This maintenance cost estimate is not really correct now

    endurance = nm_to_m(ac_data["Performance"]["range_nm"]) / ac_data["Performance"]["Vc_m/s"] / 3600

    fuel_burn = ac_data["Weights"]["Wfuel_N"] / 9.80665 / endurance
    fuel_cost = fuel_burn * vol_jet_a1_price / jet_a1_dens

    return overhaul_cost + maintenance_cost + fuel_cost


if __name__ == "__main__":  # pragma: no cover
    print(hourly_operating_cost("maf_mission_graph.csv", c206_data, ac_data=aircraft_data))

    # V&V
    input = {
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
            "Weights": {"MFW_N": 1390},
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
        hourly_operating_cost(input["mission_file"], input["standard_aircraft_data"], input["aircraft_data"]),
        abs_tol=0.05,
    )
