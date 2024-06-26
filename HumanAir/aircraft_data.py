import json
import os
import pandas as pd
import re

# from types import SimpleNamespace as Namespace

with open(os.path.join(os.path.dirname(__file__), "Configurations", "design.json"), "r", encoding="utf-8") as f:
    aircraft_data = json.load(f)

with open(os.path.join(os.path.dirname(__file__), "Configurations", "c206.json"), "r", encoding="utf-8") as f:
    c206_data = json.load(f)


def save_ac_data_to_json(ac_data=aircraft_data, filename="design.json"):  # pragma: no cover
    with open(os.path.join(os.path.dirname(__file__), "Configurations", filename), "w", encoding="utf-8") as f:
        json.dump(ac_data, f, indent=4)


def import_aero_data(file_name):  # pragma: no cover
    data = {}
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the absolute path to the design.json file
    file_path = os.path.join(script_dir, file_name)

    with open(file_path, "r") as f:
        lines = f.readlines()
        angle_file = lines[0]
        angles = re.findall(r"VLM1 -\s*-?\d+\.?\d*", angle_file)
        angles = [float(re.search(r"VLM1 -\s*(-?\d+\.?\d*)", angle).group(1)) for angle in angles]

        for line in lines[1:]:
            values = line.split()
            for i in range(0, len(values), 2):
                angle_index = i // 2
                angle = angles[angle_index]
                y_positions = float(values[i])
                lift_distribution = float(values[i + 1])

                if angle not in data:
                    data[angle] = {"y_span": [], "coefficient": []}

                data[angle]["y_span"].append(y_positions)
                data[angle]["coefficient"].append(lift_distribution)

    return data


airfoil_filepath = os.path.join(os.path.dirname(__file__), "AerodynamicData", "wing_airfoil.txt")
airfoil_shape = pd.read_csv(airfoil_filepath, sep="\\s+", header=None, names=["x", "y"], skiprows=1)

cl_wing_filepath = os.path.join(os.path.dirname(__file__), "AerodynamicData", "26-Wing_LocalLiftCoeff.txt")
cm_wing_filepath = os.path.join(os.path.dirname(__file__), "AerodynamicData", "26-Wing_TotalCm.txt")
cdi_wing_filepath = os.path.join(os.path.dirname(__file__), "AerodynamicData", "26-Wing_InducedDrag.txt")

Cl_data_wing = import_aero_data(cl_wing_filepath)
Cm_data_wing = import_aero_data(cm_wing_filepath)
Cdi_data_wing = import_aero_data(cdi_wing_filepath)
