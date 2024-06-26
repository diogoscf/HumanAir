"""========== Aircraft Design Parameters =========="""

name = "Conventional Aircraft, No Canard"
shortname = "conventional"
A = 11.5  # Cessna 206: 9.38, B2: 5.75
e = 0.82
eta_p = 0.85

Clmax_clean = 1.8
Clmax_TO = 2.6
Clmax_Land = 2.8


Cdo = 0.028  # Cessna :0.028 B2: 0.0065, 0.0165

"""========== Mission Parameters =========="""
TOP = 60  # ftion of airstrip surface characteristics and length
h_TO = 750
h_Cruise = 3000  # lowest is most efficient for prop, we have min req, so this parameter is = min allowed (no change)
h_Land = 750
s_land = 424  # 500 corrected for grass landing
f = 1  # W_to/W_land - changed to 1 since if aircraft flies on batteries W_to = W_land - used to be 0.95
temp_offset = 18

V_stall = 25  # Cessna 206 with STOL kit
V_climb = 1.2 * V_stall
V_cruise = 60  # Cessna 206

climbrate = 4.5
climbgradient = 0.083
Cl_SafetyFactor = 1.2
# nmax=4.5 not used

CruisePower = 0.8
CruiseWeight = 1
