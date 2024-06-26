import sys
import os
import json
import numpy as np
from tqdm import tqdm

# Integration in progress
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
sys.path.append(project_root)

from HumanAir.LoadingDiagram.Parameters import Parameters_ConvNoCanard as p
from HumanAir.Class_I_Weight.Class_I_Weight_Estimation import WeightEstm as WeightEstimation
from HumanAir.LoadingDiagram.Main import WP_WS
from HumanAir.CO2_Calculator.conceptual_co2 import calculate_co2_reduction_average_flight as co2  # NOTE: Probably wrong

# Hyperparameters
learning_rate = 0.01
max_iterations = 1000
tolerance = 1.1e-5

# Load the json file once
with open("../Configurations/design.json", "r") as f:
    base_dict = json.load(f)

# Initial guesses for the parameters
params = {
    "A": 8.5,
    "eta_p": 0.85,
    "Clmax_clean": 1.8,
    "Clmax_TO": 2.4,
    "Clmax_Land": 2.4,
    "Cd0": 0.025,
    "V_cruise": 63,
    "climbrate": 4.0,
}

# define the lower bound of the values
MIN_VALUES = {
    "A": 6.0,
    "eta_p": 0.8,
    "Clmax_clean": 1.5,
    "Clmax_TO": 2.0,
    "Clmax_Land": 2.0,
    "Cd0": 0.025,
    "V_cruise": 50,
    "climbrate": 3.2,
}


# get the cost function and save the best value depending on the battery percentage
def objective_function(params):
    p.A = params["A"]
    p.eta_p = params["eta_p"]
    p.Clmax_clean = params["Clmax_clean"]
    p.Clmax_TO = params["Clmax_TO"]
    p.Clmax_Land = params["Clmax_Land"]
    p.Cdo = params["Cd0"]
    p.V_cruise = params["V_cruise"]
    p.climbrate = params["climbrate"]

    working_dict = base_dict.copy()
    working_dict["endurance"] = 1111200 / p.V_cruise / 3600
    working_dict["W/P"], working_dict["W/S"] = WP_WS().calculate_optimal_point()
    weight_estimation = WeightEstimation(working_dict)

    bat = np.arange(0, 0.18, 0.001)

    coeff_exp, coeff_pol = weight_estimation.PolynomialRegression(bat)

    best_co2_ratio = 0
    best_battery_percentage = 0

    for step in range(len(bat)):
        working_dict["P_req_cruise_W"] = (
            working_dict["P_cruise/P_TO"] * np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * bat[step])
        )
        working_dict["E_bat_Wh"] = (
            working_dict["P_req_cruise_W"] * working_dict["endurance"] / working_dict["P_cruise/P_TO"] * bat[step]
        )

        co2_ratio = co2(ac_data=working_dict)

        if co2_ratio > best_co2_ratio and working_dict["E_bat_Wh"] < 250000:
            best_co2_ratio = co2_ratio
            best_battery_percentage = bat[step]

    return best_co2_ratio, best_battery_percentage


# compute the gradient of the cost function for a specific parameter
def compute_gradient_for_parameter(param_name, param_value):
    params_up = params.copy()
    params_down = params.copy()

    # Ensure that parameters don't go below certain values
    params_up[param_name] = max(params_up[param_name], MIN_VALUES[param_name])
    params_down[param_name] = max(params_down[param_name], MIN_VALUES[param_name])

    # updating the upper and lower bound paramaters by a given tolerance
    params_up[param_name] += 0.001 * params_up[param_name]
    params_down[param_name] -= 0.001 * params_down[param_name]

    # find the cost function for the upper and lower bound parameters
    f_up, bat_up = objective_function(params_up)
    f_down, bat_down = objective_function(params_down)

    # calculate the gradient
    grad = (f_up - f_down) / (2 * tolerance)

    # saving the parameters for the best case
    if f_up * 100 or f_down * 100 > 50:
        saving_dict = params.copy()
        saving_dict["CO2"] = f_up
        saving_dict["Battery"] = bat_up
        # if f_up > f_down:
        #     maxi = f_up
        # else:
        #     maxi = f_down

        with open("optimized_params_possible.json", "w") as file:
            json.dump(saving_dict, file, indent=4)

    # Choose the battery percentage based on which objective value is better
    if f_up > f_down:
        battery_percentage = bat_up
    else:
        battery_percentage = bat_down

    return grad, battery_percentage


# Perform gradient descent for each parameter
for iteration in tqdm(range(max_iterations), desc="Gradient Descent"):
    for param_name in params:
        param_value = params[param_name]

        # Compute gradient for the current parameter
        grad, battery_percentage = compute_gradient_for_parameter(param_name, param_value)

        # Update parameter value
        params[param_name] -= learning_rate * grad

    # Check convergence
    converged = True
    for param_name in params:
        grad, _ = compute_gradient_for_parameter(param_name, params[param_name])
        if abs(grad) >= 0.01:
            converged = False
            break

    if converged:
        print(f"Converged after {iteration} iterations")
        break

# Save the optimized parameters
with open("optimized_params.json", "w") as file:
    json.dump(params, file, indent=4)

print("Optimization complete.")
