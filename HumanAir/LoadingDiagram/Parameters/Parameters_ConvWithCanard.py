# import numpy as np

"""========== Aircraft Design Parameters =========="""
name = "Conventional Aircraft, With Canard"
shortname = "canard"
A = 12  # Cessna 206: 9.38, B2: 5.75
e = 0.82
TOP = 66.91
eta_p = 0.85

Clmax_clean = 1.9
Clmax_TO = 1.9
Clmax_Land = 2.5


Cdo = 0.02  # Cessna :0.028 B2: 0.0065, 0.0165 https://www.icas.org/ICAS_ARCHIVE/ICAS1982/ICAS-82-1.2.2.pdf

"""========== Mission Parameters =========="""
h_TO = 1800
h_Cruise = 3000
h_Land = 1800

V_stall = 31.38  # Cessna 206
V_climb = 1.2 * V_stall
V_cruise = 60  # Cessna 206

s_land = 500
f = 0.997  # W_to/W_land
climbrate = 5
Cl_SafetyFactor = 1.2
nmax = 4.5

CruisePower = 0.8
CruiseWeight = 0.9
