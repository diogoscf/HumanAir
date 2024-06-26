# import numpy as np
# import json
import os
import sys

# import matplotlib.pyplot as plt
# sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.AerodynamicDesign.ISA import ISA
from HumanAir.AerodynamicDesign.PlanformDesign import Planform
from HumanAir.AerodynamicDesign.FlowParameters import Flow
from HumanAir.AerodynamicDesign.LongitudinalStability import LongitudinalStability


from HumanAir.aircraft_data import aircraft_data


def aerodynamic_design(
    ac_data=aircraft_data,
    checkwingplanform=False,
    checkflowparameters=False,
    checkstability=False,
    checkhsplanform=False,
):
    MTOW = ac_data["Weights"]["MTOW_N"]
    WS = ac_data["Performance"]["W/S_N/m2"]
    AR_Wing = ac_data["Aero"]["AR"]
    AR_HS = ac_data["Aero"]["AR_HS"]
    Taper_Wing = ac_data["Aero"]["Taper_Wing"]
    Taper_HS = ac_data["Aero"]["Taper_HS"]
    QuarterChordSweep_Wing = ac_data["Aero"]["QuarterChordSweep_Wing_deg"]
    QuarterChordSweep_HS = ac_data["Aero"]["QuarterChordSweep_HS_deg"]
    CruiseHeight = ac_data["Performance"]["Altitude_Cruise_m"]
    Temp_offset = ac_data["Performance"]["Temp_offset_TO_Land_cruise"]
    TemperatureGradient = -0.0065
    CruiseVelocity = ac_data["Performance"]["Vc_m/s"]
    CLh = ac_data["Stability"]["C_L_h"]
    CLah = ac_data["Stability"]["C_L_AH"]
    Xcgh = ac_data["Stability"]["X_cg_HS"]
    XLEMAC = ac_data["Geometry"]["XLEMAC_m"]
    CgAft = ac_data["Stability"]["Cg_Aft"]
    CgFwd = ac_data["Stability"]["Cg_Front"]
    SM = ac_data["Stability"]["Stability_Margin"]
    deda = ac_data["Aero"]["deda"]
    VhV = ac_data["Aero"]["VhV"]
    FuselageLength = ac_data["Geometry"]["fus_length_m"]
    tc_wing = ac_data["Aero"]["tc_m_Wing"]
    tc_HP = ac_data["Aero"]["tc_m_HP"]

    WingPlanform = Planform(AR_Wing, Taper_Wing, QuarterChordSweep_Wing, tc_wing, MTOW=MTOW, WS=WS)

    c_root_wing = WingPlanform.RootChord()
    c_tip_wing = WingPlanform.TipChord()
    S_Wing = WingPlanform.WingSurfaceArea()
    b_Wing = WingPlanform.WingSpan()

    ISACruise = ISA(CruiseHeight, Temp_offset, TemperatureGradient)

    if checkwingplanform:
        WingPlanform.PlotWingPlanform()

    Flowvariables = Flow(CruiseVelocity, ISACruise, WingPlanform)

    if checkflowparameters:
        print("MACH = ", Flowvariables.Mach())
        print("Reynolds = ", Flowvariables.Reynolds())
        print("Beta =", Flowvariables.Beta())

    """========== Stability Analysis =========="""
    # replace with 'FX_63-137.json' and for pycharm 'NACA0012.json'
    # replace with 'c:\\Users\\nicho\\Documents\\GitHub\\ae3200-dse-g15\\HumanAir\\Configurations\\FX_63-137.json'
    # and same for the other airfoil for vscode (change nicho to your username)
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the absolute path to the design.json file
    airfoil_wing_json_path = os.path.join(script_dir, "Airfoils", "FX_63-137.json")
    airfoil_elevator_json_path = os.path.join(script_dir, "Airfoils", "NACA0012.json")

    # Stab = LongitudinalStability(
    #     CLh,
    #     CLah,
    #     Xcgh,
    #     XLEMAC,
    #     CgAft,
    #     CgFwd,
    #     SM,
    #     deda,
    #     VhV,
    #     FuselageLength,
    #     WingPlanform,
    #     "../AerodynamicDesign/Airfoils/FX_63-137.json",
    #     "../AerodynamicDesign/Airfoils/NACA0012.json",
    # )

    Stab = LongitudinalStability(
        CLh,
        CLah,
        Xcgh,
        XLEMAC,
        CgAft,
        CgFwd,
        SM,
        deda,
        VhV,
        FuselageLength,
        WingPlanform,
        airfoil_wing_json_path,
        airfoil_elevator_json_path,
    )

    if checkstability:
        Stab.Plotting()

    HSPlanform = Planform(AR_HS, Taper_HS, QuarterChordSweep_HS, tc_HP, S=WingPlanform.WingSurfaceArea() * Stab.ShS())

    c_root_HS = HSPlanform.RootChord()
    c_tip_HS = HSPlanform.TipChord()

    if checkhsplanform:
        HSPlanform.PlotWingPlanform()

    mac_wing = WingPlanform.MAC()
    mac_HS = HSPlanform.MAC()
    S_h = HSPlanform.WingSurfaceArea()
    b_h = HSPlanform.WingSpan()

    return mac_wing, mac_HS, c_root_wing, c_tip_wing, c_root_HS, c_tip_HS, S_Wing, S_h, b_Wing, b_h


if __name__ == "__main__":  # pragma: no cover
    # replace with 'design.json' for pycharm
    # replace with 'c:\\Users\\nicho\\Documents\\GitHub\\ae3200-dse-g15\\HumanAir\\Configurations\\design.json'
    # for vscode (change nicho to your username)
    # with open("design.json",'r') as f:
    #     dict = json.load(f)
    aerodynamic_design(
        aircraft_data, checkwingplanform=True, checkflowparameters=True, checkstability=True, checkhsplanform=True
    )
