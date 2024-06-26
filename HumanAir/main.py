import sys
import os
import json
import numpy as np
import logging
import colorlog

# FINAL VERSION
"Dear Programmer Please do not remove this line, it is very important for the correct function of the main program"

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(script_dir, ".."))
sys.path.append(project_root)

from HumanAir.LoadingDiagram.Parameters import Parameters_ConvNoCanard as p
from HumanAir.Class_I_Weight.Class_I_Weight_Estimation import WeightEstm as WeightEstimation
from HumanAir.LoadingDiagram.Main import WP_WS
from HumanAir.Weights_and_CG.weight_fractions import iterate_cg_lg
from HumanAir.AerodynamicDesign.Aerodynamics_Main import aerodynamic_design
from HumanAir.AerodynamicDesign.FlapsDesign import flaps_design
from HumanAir.AerodynamicDesign.FullStability import TailIteration
from HumanAir.AerodynamicDesign.AileronSizing import AileronSizing, StickArm, AileronDerivatives
from HumanAir.AerodynamicDesign.WeathercockStability import VerticalTailSizing
from HumanAir.FinancialAnalysis.conceptual_financial_analysis import hourly_operating_cost
from HumanAir.Class_II_Weight.Class_II_Weight import RunClassII
from HumanAir.Class_II_Weight.Class_II_Weight import Class_II_Weight as ClassIIWeight
from HumanAir.Vn_Diagrams.design_values import calculate_load_design_values
from HumanAir.Vn_Diagrams.loading_diagram import calc_nmax_nmin_manoeuvre
from HumanAir.CO2_Calculator.co2v2 import calculate_co2_reduction_flightdist as co2
from HumanAir.CO2_Calculator.co2v2 import improvement_co2

# from HumanAir.StructuralAnalysis.LoadDistributions import load_distribution_diagram
from HumanAir.FuselageSizing.FuselageSizing import FuselageSizing


def setup_logging():
    handler = colorlog.StreamHandler()
    handler.setFormatter(
        colorlog.ColoredFormatter(
            "%(log_color)s%(levelname)s:%(message)s",
            log_colors={"DEBUG": "cyan", "INFO": "green", "WARNING": "yellow", "ERROR": "red", "CRITICAL": "bold_red"},
        )
    )
    logger = colorlog.getLogger()
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def load_json_file(file_name):
    "Getting the design.json file"
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the absolute path to the design.json file
    design_json_path = os.path.join(script_dir, "..", "HumanAir", "Configurations", file_name)

    # Print the absolute path for debugging
    logging.info(f" Looking for {file_name} at: {os.path.abspath(design_json_path)}")

    # Attempt to open the file
    with open(design_json_path, "r") as f:
        data = json.load(f)

    logging.info(f" Opening {file_name} successful")

    return data


"Generating the design points"


def create_parameters_lists():
    A_lst = np.arange(7.0, 18.51, 0.5)
    eta_p_lst = np.arange(0.80, 0.851, 0.05)
    Clmax_clean_lst = np.arange(1.6, 2.21, 0.2)
    Clmax_TO_lst = np.arange(2, 2.61, 0.2)
    Clmax_Land_lst = np.arange(2.2, 2.81, 0.2)
    Cd0_lst = np.arange(0.028, 0.031, 0.002)
    V_cruise_lst = np.arange(60, 65.1, 1)
    climbrate_lst = np.arange(2.5, 5.01, 0.5)

    return A_lst, eta_p_lst, Clmax_clean_lst, Clmax_TO_lst, Clmax_Land_lst, Cd0_lst, V_cruise_lst, climbrate_lst


def Generate(p, ac_data, run=False, parameter_func=create_parameters_lists, bat_step=0.001):
    # tune the parameters with a reasonable range
    # A_lst = np.arange(7.0, 18.51, 0.5)
    # eta_p_lst = np.arange(0.80, 0.851, 0.05)
    # Clmax_clean_lst = np.arange(1.6, 2.21, 0.2)
    # Clmax_TO_lst = np.arange(2, 2.61, 0.2)
    # Clmax_Land_lst = np.arange(2.2, 2.81, 0.2)
    # Cd0_lst = np.arange(0.028, 0.031, 0.002)
    # V_cruise_lst = np.arange(60, 65.1, 1)
    # climbrate_lst = np.arange(2.5, 5.01, 0.5)

    (
        A_lst,
        eta_p_lst,
        Clmax_clean_lst,
        Clmax_TO_lst,
        Clmax_Land_lst,
        Cd0_lst,
        V_cruise_lst,
        climbrate_lst,
    ) = parameter_func()

    # calculate the total numbers of iterations
    total_iterations = (
        len(A_lst)
        * len(eta_p_lst)
        * len(Clmax_clean_lst)
        * len(Clmax_TO_lst)
        * len(Clmax_Land_lst)
        * len(Cd0_lst)
        * len(V_cruise_lst)
        * len(climbrate_lst)
    )

    # initialise the iteration counter
    current_iteration = 0

    # run condition to not run the loop by mistake
    if run:
        idx = -1
        dict_iterations = {}
        for A in A_lst:
            for eta_p in eta_p_lst:
                for Clmax_clean in Clmax_clean_lst:
                    for Clmax_TO in Clmax_TO_lst:
                        Clmax_Land_lst = np.arange(Clmax_TO + 0.2, 2.81, 0.2)
                        for Clmax_Land in Clmax_Land_lst:
                            for Cd0 in Cd0_lst:
                                for V_cruise in V_cruise_lst:
                                    # update the endurance based on v cruise
                                    ac_data["Performance"]["endurance"] = 1111200 / V_cruise / 3600

                                    for climbrate in climbrate_lst:
                                        current_iteration += 1

                                        # print the iteration number every 500 steps
                                        if current_iteration % 500 == 0:
                                            logging.info(
                                                " Iteration: " + str(current_iteration) + "/" + str(total_iterations)
                                            )
                                            # print()

                                        # update the parameters
                                        p.A = A
                                        p.eta_p = eta_p
                                        p.Clmax_clean = Clmax_clean
                                        p.Clmax_TO = Clmax_TO
                                        p.Clmax_Land = Clmax_Land
                                        p.Cdo = Cd0
                                        p.V_cruise = V_cruise
                                        p.climbrate = climbrate

                                        # initialise the loading diagram to get W/P and W/S
                                        (
                                            ac_data["Performance"]["W/P_N/W"],
                                            ac_data["Performance"]["W/S_N/m2"],
                                        ) = WP_WS().calculate_optimal_point()

                                        # initialise the class I weight estimation
                                        WeightEstm = WeightEstimation(ac_data)

                                        # initialise the bat percentage
                                        bat = np.arange(0, 0.18, bat_step)
                                        coeff_exp, coeff_pol = WeightEstm.PolynomialRegression(bat)

                                        # run the regression to find power required cruise
                                        co2_ratio_max = 0

                                        # set the condition to find the first point with co2 ratio > 50
                                        ok = 0
                                        for step in range(len(bat)):
                                            # calculate the power required cruise
                                            ac_data["Power_prop"]["P_req_cruise_W"] = (
                                                ac_data["Performance"]["P_cruise/P_TO"]
                                                * np.exp(coeff_exp[1])
                                                * np.exp(coeff_exp[0] * bat[step])
                                            )
                                            ac_data["Power_prop"]["E_bat_Wh"] = (
                                                ac_data["Power_prop"]["P_req_cruise_W"]
                                                * ac_data["Performance"]["endurance"]
                                                / ac_data["Performance"]["P_cruise/P_TO"]
                                                * bat[step]
                                            )

                                            # calculate the co2 ratio for the specific combination of parameters
                                            co2_ratio = co2(ac_data=ac_data)

                                            print("----------")
                                            print(co2_ratio)
                                            print(ac_data["Power_prop"]["E_bat_Wh"])
                                            print(bat[step])
                                            print(ac_data["Power_prop"]["P_req_cruise_W"] / 0.8)
                                            print(
                                                np.sqrt(
                                                    ac_data["Power_prop"]["P_req_cruise_W"]
                                                    / ac_data["Performance"]["P_cruise/P_TO"]
                                                    * ac_data["Performance"]["W/P_N/W"]
                                                    / ac_data["Performance"]["W/S_N/m2"]
                                                    * A
                                                )
                                            )

                                            if co2_ratio * 100 > 15 and ok == 0:
                                                idx += 1
                                                ok = 1

                                            if (
                                                co2_ratio * 100 > co2_ratio_max
                                                and co2_ratio * 100 > 15
                                                and ac_data["Power_prop"]["E_bat_Wh"] < 200000
                                                and bat[step] > 0.1
                                                and ac_data["Power_prop"]["P_req_cruise_W"] / 0.8 < 330000
                                                and np.sqrt(
                                                    ac_data["Power_prop"]["P_req_cruise_W"]
                                                    / ac_data["Performance"]["P_cruise/P_TO"]
                                                    * ac_data["Performance"]["W/P_N/W"]
                                                    / ac_data["Performance"]["W/S_N/m2"]
                                                    * A
                                                )
                                                < 21.0
                                            ):  # the <250000 condition is for the battery to be able to be charged
                                                CO2 = co2_ratio

                                                # save the combination of parameters
                                                dict_iterations[str(idx)] = {}
                                                dict_iterations[str(idx)]["A"] = A
                                                dict_iterations[str(idx)]["eta_p"] = eta_p
                                                dict_iterations[str(idx)]["Clmax_clean"] = Clmax_clean
                                                dict_iterations[str(idx)]["Clmax_TO"] = Clmax_TO
                                                dict_iterations[str(idx)]["Clmax_Land"] = Clmax_Land
                                                dict_iterations[str(idx)]["Cd0"] = Cd0
                                                dict_iterations[str(idx)]["V_cruise"] = V_cruise
                                                dict_iterations[str(idx)]["climbrate"] = climbrate
                                                dict_iterations[str(idx)]["CO2"] = np.round(CO2 * 100, 2)
                                                dict_iterations[str(idx)]["W/P"] = ac_data["Performance"]["W/P_N/W"]
                                                dict_iterations[str(idx)]["W/S"] = ac_data["Performance"]["W/S_N/m2"]
                                                dict_iterations[str(idx)]["bat"] = bat[step]

                                                co2_ratio_max = co2_ratio

        # save the json file with all possible design options
        data_iterations_json_path = os.path.join(script_dir, "..", "HumanAir", "Configurations", "data_iterations.json")
        with open(data_iterations_json_path, "w+") as f:
            json.dump(dict_iterations, f, indent=4)


"Calculating the weighted score to find the best design point"


def calculate_weighted_score(point_data, weights):
    # setting the worst optimal data
    worst_optimal_data = {
        "A": 9.5,
        "eta_p": 0.85,
        "Clmax_clean": 2.2,
        "Clmax_TO": 2.6,
        "Clmax_Land": 2.6,
        "Cd0": 0.026,
        "V_cruise": 60,
        "climbrate": 3.5,
        "bat": 0.01,
        "CO2": 30,
    }

    # initialing the score
    score = 0

    # calculate the score
    for key, weight in weights.items():
        score += (
            (worst_optimal_data[key] - point_data.get(key, 0)) / worst_optimal_data[key] * weight
        )  # get a normalised value for the contribution of each key
    return score


"Finding the optimal design point with the highest score and the lowest battery weight"


def find_optimal_design(
    ac_data,
    maximum_weight_battery=1000,
    weights=None,
    CO2_threshold=50,
    design_points=None,
    printing=False,
    step=0,
    weight_estimation_class=WeightEstimation,
):
    optimum_design_points = {key: value for key, value in design_points.items() if value["CO2"] > CO2_threshold}
    scores = [(key, calculate_weighted_score(value, weights)) for key, value in optimum_design_points.items()]

    sorted_design_points = sorted(scores, key=lambda x: x[1], reverse=True)

    value = np.inf
    minimum = np.inf

    while value > maximum_weight_battery or design_points[sorted_design_points[step][0]]["Cd0"] < 0.027:
        step += 1

        optimum_design_option = design_points[sorted_design_points[step][0]]

        # Updating the design.json file with the optimum design option
        ac_data["Aero"]["AR"] = optimum_design_option["A"]
        ac_data["Power_prop"]["eta_p"] = optimum_design_option["eta_p"]
        ac_data["Aero"]["CLmax_clean"] = optimum_design_option["Clmax_clean"]
        ac_data["Aero"]["CLmax_TO"] = optimum_design_option["Clmax_TO"]
        ac_data["Aero"]["CLmax_Land"] = optimum_design_option["Clmax_Land"]
        ac_data["Aero"]["CD0"] = optimum_design_option["Cd0"]
        ac_data["Performance"]["Vc_m/s"] = optimum_design_option["V_cruise"]
        ac_data["Performance"]["climb_rate"] = optimum_design_option["climbrate"]
        ac_data["Power_prop"]["bat"] = optimum_design_option["bat"]
        ac_data["Performance"]["W/S_N/m2"] = optimum_design_option["W/S"]
        ac_data["Performance"]["W/P_N/W"] = optimum_design_option["W/P"]
        ac_data["Performance"]["CO2"] = optimum_design_option["CO2"]

        weight_estimation = WeightEstimation(ac_data)  # Correct instantiation of the WeightEstimation class
        iteration_results = weight_estimation.Iterations(ac_data["Power_prop"]["bat"])
        value = iteration_results[4]
        MTOW = iteration_results[1]
        ac_data["Power_prop"]["P_req_TO_W"] = 9.81 * MTOW / ac_data["Performance"]["W/P_N/W"]
        ac_data["Power_prop"]["P_req_cruise_W"] = 9.81 * 0.8 * MTOW / ac_data["Performance"]["W/P_N/W"]
        ac_data["Power_prop"]["E_bat_Wh"] = (
            ac_data["Power_prop"]["P_req_cruise_W"]
            * ac_data["Performance"]["endurance"]
            / ac_data["Performance"]["P_cruise/P_TO"]
            * ac_data["Power_prop"]["bat"]
        )

        if value < minimum:
            minimum = value

    # Printing option
    if printing:
        print(sorted_design_points[step])
        print(design_points[str(sorted_design_points[step][0])])
        print(weight_estimation.Iterations(ac_data["Power_prop"]["bat"]))


"Main Function to run the program"
if __name__ == "__main__":
    # set up the conditions to run the program
    run_generate = False
    run_classI = False
    run_classII = True
    run_fuselage_sizing = True

    # initialise the logging
    setup_logging()

    # Integration in progress v2
    logging.info(" Starting the program")

    # print which part of the main function will run
    print("\nWelcome to the HumanAir script. Hope you will have a good time!\n")
    print("The following parts of the program will run: \n")
    print("Generate design points: ", run_generate)
    print("Class I Weight Estimation: ", run_classI)
    print("Class II Weight Estimation: ", run_classII)
    print("Fuselage Sizing: ", run_fuselage_sizing)
    print("\n")

    # initialise the design.json file
    ac_data = load_json_file("design.json")

    M_D, V_H, n_ult_cruise, n_ult_land = calculate_load_design_values(ac_data)
    ac_data["Performance"]["n_ult"] = max(n_ult_land, n_ult_cruise)
    ac_data["Performance"]["M_D"] = M_D
    ac_data["Performance"]["Vh_m/s"] = V_H
    ac_data["Performance"]["n_ult_l"] = n_ult_land

    # run the generation code for the data iteration points
    if run_generate:
        logging.info(" Starting generating the new possible design points. This may take a while.")
        Generate(p, ac_data, run_generate)

    # run class I estimation
    if run_classI:
        logging.info(" Getting the data from the design point options")

        design_points = load_json_file("data_iterations.json")

        "Setting up the weights to get the optimal design point"
        # Weights for each key
        # increasing score: positive weight
        # decreasing score: negative weight
        # no influence: weight = 0

        importance_weights = {
            "A": +0.0,
            "eta_p": +0.1,
            "Clmax_clean": +0.05,
            "Clmax_TO": +0,
            "Clmax_Land": +0,
            "Cd0": -0.25,
            "V_cruise": -0.05,
            "climbrate": -0.1,
            "bat": -0.2,
            "CO2": -0.25,
        }

        maximum_weight_battery = 1000
        CO2_threshold = 20
        printing = False

        # set up that the optimal stability range is not yet set
        find_optimal_stability = False

        # initialise from which step to start to search from the design_iterations.json
        step = 0

        logging.info(" Starting the search for optimal stability range")

        while not find_optimal_stability:
            find_optimal_design(
                ac_data,
                maximum_weight_battery=maximum_weight_battery,
                weights=importance_weights,
                CO2_threshold=CO2_threshold,
                design_points=design_points,
                printing=printing,
                step=step,
            )
            logging.info(
                f"Finding the optimal design point with a maximum battery weight of {maximum_weight_battery}[kg]"
                + f"with a CO2 threshold of {CO2_threshold}[%] successful"
            )
            logging.info(" Calculating the weight components")

            weight_estimation = WeightEstimation(ac_data)
            component_weights = weight_estimation.Iterations(ac_data["Power_prop"]["bat"])

            print(
                f"Component weights:"
                f" OEW {round(component_weights[2], 2)}[kg],"
                f" Powertrain {round(component_weights[3], 2)}[kg],"
                f" Battery {round(component_weights[4], 2)}[kg],"
                f" Fuel {round(component_weights[5], 2)}[kg],"
                f" Wing {round(component_weights[6], 2)}[kg],"
                f" Wpl_des {round(component_weights[7], 2)}[kg]"
            )

            print("Total weight:", round(component_weights[1], 2), "[kg] including contingency")
            print("Contingency:", (round((ac_data["Contingency"] - 1) * 100, 0)), "%")

            # save the component weights to the dictionary in newtons
            ac_data["Weights"]["MTOW_N"] = 9.81 * round(component_weights[1], 2)
            ac_data["Weights"]["OEW_N"] = 9.81 * (
                round(component_weights[2], 2)
                + round(component_weights[3], 2)
                + round(component_weights[4], 2)
                + round(component_weights[6], 2)
            )
            ac_data["Weights"]["Wptr_N"] = 9.81 * round(component_weights[3], 2)
            ac_data["Weights"]["Wbat_N"] = 9.81 * round(component_weights[4], 2)
            ac_data["Weights"]["Wfuel_N"] = 9.81 * round(component_weights[5], 2)
            ac_data["Weights"]["Ww_N"] = 9.81 * round(component_weights[6], 2)
            ac_data["Weights"]["W_L_N"] = 9.81 * (round(component_weights[1], 2) - round(component_weights[5], 2))

            logging.info(" Calculating the weight components successful")

            # set up the condition to set up the range where the cg of the wing is with report of the mac
            xcg_location_percentage = np.arange(-0.1, 0.51, 0.1)
            logging.info(
                " Starting the search for the optimal stability range in terms of where to position the cg of the wing"
            )

            # iterate over the percentage to find the optimal stability range
            for pct in xcg_location_percentage:
                logging.info(" Calculating the Xcg excursion")

                wcg, CGlist, xlemac = iterate_cg_lg(ac_datafile=ac_data, PERCENTAGE=pct)

                # dont remove this line as it complies with nicholas's mood
                ac_data["Geometry"]["XLEMAC_m"] = xlemac
                mac_wing = aerodynamic_design(
                    ac_data,
                    checkwingplanform=False,
                    checkflowparameters=False,
                    checkstability=False,
                    checkhsplanform=False,
                )[0]

                ac_data["Stability"]["Cg_Aft"] = (round(max(CGlist), 2) - ac_data["Geometry"]["XLEMAC_m"]) / mac_wing
                ac_data["Stability"]["Cg_Front"] = (round(min(CGlist), 2) - ac_data["Geometry"]["XLEMAC_m"]) / mac_wing

                logging.info(" Prepare to check the stability")

                # dont remove this line as it complies with nicholas's mood
                aerodynamic_design(
                    ac_data,
                    checkwingplanform=False,
                    checkflowparameters=False,
                    checkstability=True,
                    checkhsplanform=False,
                )

                print(
                    "Is stability satisfied at a X_LEMAC "
                    + str(round(ac_data["Geometry"]["XLEMAC_m"], 2))
                    + " [m]"
                    + "|"
                    + "[Y/N]: "
                )
                answer = input()

                if answer.lower() == "y":
                    break

            if answer.lower() == "y":
                # print the range of the cg
                print(f"Xcg Range is between: {round(min(CGlist), 2)} and {round(max(CGlist), 2)} [m]")

                logging.info(" Calculating the Xcg excursion successful")
                logging.info(" Calculating the MAC")

                # dont remove this line as it complies with nicholas's mood
                mac_wing = aerodynamic_design(
                    ac_data,
                    checkwingplanform=False,
                    checkflowparameters=False,
                    checkstability=False,
                    checkhsplanform=False,
                )[0]
                print(f"MAC: {round(mac_wing, 2)} [m]")

                logging.info(" Calculating the MAC successful")
                logging.info(" Calculating the aerodynamic design")

                find_optimal_stability = True

                # get the aerodynamic specifications and save them to the dictionary
                (
                    mac_wing,
                    mac_HS,
                    c_root_wing,
                    c_tip_wing,
                    c_root_HS,
                    c_tip_HS,
                    S_Wing,
                    S_h,
                    b_Wing,
                    b_h,
                ) = aerodynamic_design(
                    ac_data,
                    checkwingplanform=True,
                    checkflowparameters=False,
                    checkstability=True,
                    checkhsplanform=True,
                )

                # updating the aerodynamic data in the dictionary
                ac_data["Aero"]["S_Wing"] = S_Wing
                ac_data["Aero"]["S_h"] = S_h
                ac_data["Aero"]["MAC_wing"] = mac_wing
                ac_data["Aero"]["MAC_HS"] = mac_HS
                ac_data["Aero"]["c_root_wing"] = c_root_wing
                ac_data["Aero"]["c_tip_wing"] = c_tip_wing
                ac_data["Aero"]["c_root_HS"] = c_root_HS
                ac_data["Aero"]["c_tip_HS"] = c_tip_HS
                ac_data["Aero"]["b_Wing"] = b_Wing
                ac_data["Aero"]["b_h"] = b_h

                logging.info(" Calculating the aerodynamic design successful")

                # calculate the loading distribution diagrams
                ac_data["Performance"]["n_max"], ac_data["Performance"]["n_min"] = calc_nmax_nmin_manoeuvre(
                    ac_data["Weights"]["MTOW_N"]
                )

                # logging.info(" Calculating the loading distribution diagram")

                # load_distribution_diagram(ac_data=ac_data)

                # logging.info(" Calculating the loading distribution diagram successful")
                logging.info(" Calculating the hourly price")

                # calculating the hourly cost
                cost = hourly_operating_cost(
                    "maf_mission_graph.csv", ac_data=ac_data, fuel_weight=ac_data["Weights"]["Wfuel_N"]
                )

                print(f"Cost: {round(cost, 2)} [US$]")

                print(ac_data["Stability"]["Cg_Aft"])

                logging.info(" Calculating the hourly price successful")
                logging.info(" Saving the modified design.json file")

                design_json_path = os.path.join(script_dir, "..", "HumanAir", "Configurations", "design.json")
                logging.info(" Design.json saved at: " + design_json_path)

                # save the updated dictionary
                with open(design_json_path, "w") as f:
                    json.dump(ac_data, f, indent=4)

                # calculating if we get the copium batteries how much co2 reduction would increase
                ac_data["Power_prop"]["E_bat_Wh"] = 685 / 350 * ac_data["Power_prop"]["E_bat_Wh"]

                print(
                    "Reduction with future expected battery technology: "
                    + str(round(co2(ac_data=ac_data) * 100, 2))
                    + "[%]"
                )
                ac_data["Power_prop"]["E_bat_Wh"] = 350 / 685 * ac_data["Power_prop"]["E_bat_Wh"]

    # run the class 2 estimation
    if run_classII:
        # initialise the class I weight estimation
        WeightEstm = ClassIIWeight(ac_data=ac_data)

        # initialise the bat percentage
        bat = np.arange(0, 0.3, 0.001)
        coeff_exp, coeff_pol = WeightEstm.PolynomialRegression(bat)

        # set the condition to find the first point with co2 ratio > 50
        ok = 0
        co2_ratio_max = 0
        pbat = ac_data["Power_prop"]["bat"]

        # update the power required and the battery energy
        old_P_cruise = ac_data["Power_prop"]["P_req_cruise_W"]
        old_E_bat = ac_data["Power_prop"]["E_bat_Wh"]

        # loop to find the best battery percentage after the class 2 weight estimation
        for step in range(len(bat)):
            # calculate the power required cruise
            ac_data["Power_prop"]["P_req_cruise_W"] = (
                ac_data["Performance"]["P_cruise/P_TO"] * np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * bat[step])
            )
            ac_data["Power_prop"]["E_bat_Wh"] = (
                ac_data["Power_prop"]["P_req_cruise_W"]
                * ac_data["Performance"]["endurance"]
                / ac_data["Performance"]["P_cruise/P_TO"]
                * bat[step]
            )

            # calculate the co2 ratio for the specific combination of parameters
            co2_ratio = co2(ac_data=ac_data)

            # set condition to find best battery percentage with highest co2 reduction and specific energy below
            if co2_ratio * 100 > co2_ratio_max and ac_data["Power_prop"]["E_bat_Wh"] < 189000:
                co2_ratio_max = co2_ratio
                pbat = bat[step]
                old_P_cruise = ac_data["Power_prop"]["P_req_cruise_W"]
                old_E_bat = ac_data["Power_prop"]["E_bat_Wh"]
                old_P_TO = ac_data["Power_prop"]["P_req_cruise_W"] / ac_data["Performance"]["P_cruise/P_TO"]

        # update the power required and the battery energy with the one found in the co2 loop
        ac_data["Power_prop"]["bat"] = pbat
        ac_data["Power_prop"]["P_req_cruise_W"] = old_P_cruise
        ac_data["Power_prop"]["E_bat_Wh"] = old_E_bat
        ac_data["Power_prop"]["P_req_TO_W"] = old_P_TO

        # calculate the class 2 weights components and print them
        logging.info(" Calculate Class II Weight Groups")
        class_2_dictionary = RunClassII(ac_data=ac_data, check=True, pbat=ac_data["Power_prop"]["bat"])

        # print the new MTOW and the battery weight with a 12% contingency
        print(f"Battery Weight: {round(class_2_dictionary['CL2Weight']['Wbat_N'] / 9.81, 2)} [kg]")
        print(f"Total Weight: {round(class_2_dictionary['CL2Weight']['MTOW_N'] / 9.81, 2)} [kg]")
        print(f"Contingency: {round((class_2_dictionary['Contingency_C2W'] - 1) * 100, 0)} %")

        logging.info(" Class II Weight Groups calculated successfully")

        # update the power required
        class_2_dictionary["Power_prop"]["P_req_TO_W"] = (
            class_2_dictionary["CL2Weight"]["MTOW_N"] / class_2_dictionary["Performance"]["W/P_N/W"]
        )
        class_2_dictionary["Power_prop"]["P_req_cruise_W"] = (
            class_2_dictionary["Performance"]["P_cruise/P_TO"]
            * class_2_dictionary["CL2Weight"]["MTOW_N"]
            / class_2_dictionary["Performance"]["W/P_N/W"]
        )
        class_2_dictionary["Power_prop"]["E_bat_Wh"] = (
            class_2_dictionary["CL2Weight"]["Wbat_N"] / 9.81 * class_2_dictionary["Power_prop"]["E_bat_Wh/kg"]
        )

        # print the co2 reduction with the current battery technology
        first_level = round(co2(ac_data=class_2_dictionary), 2)
        print("Reduction with current battery technology: " + str(round(first_level * 100, 2)) + "[%]")
        class_2_dictionary["Power_prop"]["E_bat_Wh"] = 685 / 350 * class_2_dictionary["Power_prop"]["E_bat_Wh"]

        # print the co2 reduction with the future expected battery technology
        second_level = round(co2(ac_data=class_2_dictionary), 2)
        print("Reduction with future expected battery technology: " + str(round(second_level * 100, 2)) + "[%]")
        logging.info(" Class II Weight Groups calculated successfully")

        # save the correct energy for current battery technology
        class_2_dictionary["Power_prop"]["E_bat_Wh"] = 350 / 685 * class_2_dictionary["Power_prop"]["E_bat_Wh"]

        logging.info(" Yearly prediction of CO2 reduction vs new Battery technology introduction year")
        improvement_co2(first_level=first_level, second_level=second_level, check_over_time=True)

        # calculate the loading distribution diagrams
        (
            class_2_dictionary["Performance"]["n_max"],
            class_2_dictionary["Performance"]["n_min"],
        ) = calc_nmax_nmin_manoeuvre(class_2_dictionary["Weights"]["MTOW_N"])

        # logging.info(" Calculating the loading distribution diagram")

        # # plot the load distribution diagrams
        # load_distribution_diagram(ac_data=class_2_dictionary)

        # logging.info(" Calculating the loading distribution diagram successful")
        logging.info(" Sizing the flaps")

        # sizing the flaps
        flaps_design(ac_data=class_2_dictionary)

        logging.info(" Calculating flap position and design succesfully")
        logging.info(" Sizing the ailerons")

        # sizing the ailerons and the stick arm forces
        AileronSizing(acd=class_2_dictionary)
        AileronDerivatives(acd=class_2_dictionary)
        StickArm(acd=class_2_dictionary, alpha=0.0, delta=14.0, h=3000.0, V=60.0)

        logging.info(" Calculating aileron position and design succesfully")
        logging.info(" Sizing the vertical tail")

        # size the vertical tail
        VerticalTailSizing(acd=class_2_dictionary)
        logging.info(" Sizing the vertical tail successful")
        logging.info(" Sizing the horizontal tail")

        # size the horizontal tail
        TailIteration(ac_datafile=class_2_dictionary)
        logging.info(" Sizing the horizontal tail successful")

        logging.info(" Calculating the cost")
        cost = hourly_operating_cost(
            "maf_mission_graph.csv", ac_data=class_2_dictionary, fuel_weight=class_2_dictionary["CL2Weight"]["Wfuel_N"]
        )
        print(f"Cost: {round(cost, 2)} [US$]")
        logging.info(" Calculating the cost successful")

        logging.info(" Sizing flaps")
        # sizing the flaps
        class_2_dictionary = flaps_design(ac_data=class_2_dictionary)

        logging.info(" Calculating flap position and design succesfully")

        # save the updated dictionary
        design_json_path = os.path.join(script_dir, "Configurations", "design.json")
        logging.info(" Design.json saved at: " + design_json_path)

        with open(design_json_path, "w") as f:
            json.dump(class_2_dictionary, f, indent=4)

        logging.info(" Program finished successfully")

    # run the fuselage sizing script
    if run_fuselage_sizing:
        # initialise the gear clearance
        s_gear = 0.2

        # load the design.json file
        fuselage_sizing_dict = load_json_file("design.json")

        xcg_bat = fuselage_sizing_dict["Stability"]["Xcg_battery_m"] / fuselage_sizing_dict["Geometry"]["fus_length_m"]

        logging.info(" Calculating the fuselage dimension")

        # get the fuselage sizing class
        fuselage_size = FuselageSizing(ac_data=fuselage_sizing_dict, bat_xcg=xcg_bat)

        # fuselage_size = FuselageSizing(ac_data=fuselage_sizing_dict, bat_xcg=xcg_bat)

        # show all of the dimensions of the fuselage
        print("Top Width", round(fuselage_size.top_width(), 2), "[m]")
        print("Bottom Width", round(fuselage_size.bottom_width(s_gear=s_gear), 2), "[m]")
        print("Fuselage Height", round(fuselage_size.height(), 2), "[m]")
        print(
            "Fuselage Length Without Engine",
            round(fuselage_size.length_fus() - fuselage_size.l_engine - fuselage_size.l_enbu, 2),
            "[m]",
        )
        print("Fuselage Length", round(fuselage_size.length_fus(), 2), "[m]")
        print("Maximum Perimeter", round(fuselage_size.maximum_perimeter(s_gear=s_gear), 2), "[m]")
        print("Fuselage Wetted Area", round(fuselage_size.fuselage_wetted(s_gear=s_gear), 2), "[m^2]")
        print("Main Strut Length", round(fuselage_size.length_main_strut(s_gear=s_gear), 2), "[m]")
        print("Nose Strut Length", round(fuselage_size.h_nose_strut, 2), "[m]")

        # plot the side and front view
        print(fuselage_size.below_position(s_gear=0.1))
        fuselage_size.plot_side_drawing(s_gear=0.2, ac_data=fuselage_sizing_dict)
        fuselage_size.plot_front_view(s_gear=0.2)

        # save the updated dictionary
        design_json_path = os.path.join(script_dir, "Configurations", "design.json")
        logging.info(" Design.json saved at: " + design_json_path)

        with open(design_json_path, "w") as f:
            json.dump(fuselage_sizing_dict, f, indent=4)

        logging.info(" Program finished successfully")

        logging.info(" Calculating the fuselage dimension successful")

print("Finally, I am free!")
