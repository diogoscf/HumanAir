import os
import sys

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
sys.path.append(project_root)


import HumanAir.FlightPerformance.aircraft as aircraft

acf = aircraft.Aircraft()

"""========== Aircraft Design Parameters =========="""

name = "Our own current design"
A = acf.AR
e = acf.e_clean
eta_p = acf.max_eff_prop

Clmax_clean = acf.CLmax_clean
Clmax_TO = acf.CLmax_TO
Clmax_Land = acf.CLmax_land


Cdo = acf.CD0_clean  # Cessna :0.028 B2: 0.0065, 0.0165

"""========== Mission Parameters =========="""
TOP = 60  # function of airstrip surface characteristics and length
h_TO = 750
h_Cruise = 3000  # lowest is most efficient for prop, we have min req, so this parameter is = min allowed (no change)
h_Land = 750
s_land = 424  # 500 corrected for grass landing
f = 1  # W_to/W_land - changed to 1 since if aircraft flies on batteries W_to = W_land - used to be 0.95
temp_offset = 18

V_stall = 27  # Cessna 206 with STOL kit
V_climb = 1.2 * V_stall  # not used
V_cruise = 60  # Cessna 206

climbrate = 5
climbgradient = 0.083
Cl_SafetyFactor = 1.2
# nmax=4.5 not used

CruisePower = 0.8
CruiseWeight = 1
