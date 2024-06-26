import numpy as np

# import csv
import os

# import json
import sys

# import time
# import pandas as pd
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

# from aircraft_data import aircraft_data, c206_data
from HumanAir.CO2_Calculator.conceptual_co2 import calculate_mission_freqs

# https://en.wikipedia.org/wiki/Probability_distribution_fitting


if __name__ == "__main__":
    mission_freqs = calculate_mission_freqs("maf_mission_graph.csv")
    print(mission_freqs)
    # plt.stairs(mission_freqs[:,2]*100, np.concatenate((np.array([0]), mission_freqs[:,0].flatten())))
    plt.scatter(np.log(mission_freqs[:, 0]), mission_freqs[:, 2] * 100)

    plt.show()
