import pandas as pd
import numpy as np
import sys
from math import tan, sqrt, floor, ceil
import time
import os
import matplotlib.pyplot as plt
from sympy import symbols, solve

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, ".."))
sys.path.append(project_root)

from HumanAir.aircraft_data import aircraft_data
from HumanAir.isa import isa
from HumanAir.FuselageSizing.FuselageSizing import FuselageSizing


def component_mass(ac_datafile=aircraft_data):
    # Convert weights to kg and with contingency, put them in a list (wcg)
    wcg = np.zeros((3, 9))

    MTOW_cont = ac_datafile["CL2Weight"]["MTOW_N"] / 9.81
    wcg[1, -1] = (
        ac_datafile["CL2Weight"]["OEW"] / 9.81 - ac_datafile["CL2Weight"]["W_pilot"] / 9.81
    )  # Empty mass (without pilot)
    wcg[1, 0] = ac_datafile["CL2Weight"]["Wing Weight"] / 9.81  # Wing mass
    wcg[1, 1] = ac_datafile["CL2Weight"]["Landing Gear Weight"] / 9.81 * 0.772371311  # Main gear mass
    wcg[1, 2] = ac_datafile["CL2Weight"]["Total Powerplant Weight"] / 9.81  # Powertrain mass
    wcg[1, 3] = ac_datafile["CL2Weight"]["Landing Gear Weight"] / 9.81 * 0.227628689  # Nose wheel mass
    wcg[1, 4] = ac_datafile["CL2Weight"]["Fuselage Weight"] / 9.81  # Fuselage mass
    wcg[1, 5] = ac_datafile["CL2Weight"]["Empennage Weight"] / 9.81  # Empennage mass
    wcg[1, 6] = ac_datafile["CL2Weight"]["Total Fixed Equipment Weight"] / 9.81  # Fixed equipment mass
    wcg[1, 7] = ac_datafile["CL2Weight"]["Wbat_N"] / 9.81  # Battery mass

    # Set up weight fractions based on MTOW
    for weight in range(len(wcg[0, :])):
        wcg[0, weight] = wcg[1, weight] / MTOW_cont

    # Check whether component fractions add up to EW fraction (margin kept for discrepancy in Class 2 weight)
    fracsum = np.sum(wcg[0, 0:-1])
    if abs(fracsum / wcg[0, -1] - 1) > 0.03:
        print("WARNING: WEIGHT FRACTIONS DIFFER MORE THAN 3%")
        print("Expected OEW/MTOW:", wcg[0, -1])
        print("Summed OEW/MTOW:", fracsum)
        con = input("Continue? (y/n): ")
        if con.lower() == "n":
            sys.exit("Weight fractions do not add up")

    # Return fractions and masses of each component: Wing, MLG, pwtr, NLG, fus, emp, FE, bat, EW
    return wcg


def find_lg(nose_loading, aftcg, wcg, ac_datafile=aircraft_data, max_load_nose=None, fwdcg=None):
    # Import compatible tyre database, in ascending order of maximum load
    tyre_file = os.path.join(os.path.dirname(__file__), "tiredata.csv")
    tyres = pd.read_csv(tyre_file, index_col=0).to_numpy()

    # Determine weight (in kg) on each tyre
    Pmw = (1 - nose_loading) * ac_datafile["CL2Weight"]["MTOW_N"] / (2 * 9.81)
    Pnw = nose_loading * ac_datafile["CL2Weight"]["MTOW_N"] / 9.81
    if max_load_nose:
        Pnw = max_load_nose  # Accounts for additional load from front CG, done after iterations

    # Determine weight on entire main gear (necessary for stability calculations)
    Pmg = (1 - nose_loading) * ac_datafile["CL2Weight"]["MTOW_N"] / (9.81)

    # Choose smallest available main and nose tyres from the list: Ascend list until load is satisfied
    for tyre in range(len(tyres[:, 0])):
        Wt_m = tyres[tyre, 0]
        Dw_m = tyres[tyre, 1]
        if tyres[tyre, 2] >= Pmw:
            break

    for tyre in range(len(tyres[:, 0])):
        Wt_n = tyres[tyre, 0]
        Dw_n = tyres[tyre, 1]
        if tyres[tyre, 2] >= Pnw:
            break

    if Pmw > tyres[-1, 2]:
        print("WARNING: NO TYRE AVAILABLE")
        con = input("Continue? (y/n): ")
        if con == "n":
            sys.exit("Too heavy for landing gear")

    # Find empty weight CG height w.r.t fuselage bottom: Weighted average
    # TODO: Import these from structures or dict
    h_f = ac_datafile["Geometry"]["fus_height_m"]
    Hcg_wing = 1.0 * h_f
    Hcg_gear = 0 * h_f
    Hcg_pwtr = 0.5 * h_f
    Hcg_fus = 0.5 * h_f
    Hcg_emp = 1.1 * h_f
    # Hcg_FE = 0.5 * h_f
    Hcg_bat = 0.1 * h_f
    Hcg_OEW = (
        Hcg_wing * wcg[1, 0]
        + Hcg_gear * ac_datafile["CL2Weight"]["Landing Gear Weight"]
        + Hcg_pwtr * wcg[1, 2]
        + Hcg_fus * wcg[1, 4]
        + Hcg_emp * wcg[1, 5]
        + Hcg_bat * wcg[1, 6]
        + Hcg_bat * wcg[1, 7]
    ) / np.sum(wcg[1, 0:8])

    # Find necessary strut height for tail strike-avoidance
    TB = ac_datafile["Landing_gear"]["Tipback_angle_deg"]
    SP = ac_datafile["Landing_gear"]["Scrape_angle_deg"]

    H_s = 1.5 * Dw_m  # initial guess
    l_m = tan(np.radians(TB)) * (Hcg_OEW + H_s + 0.5 * Dw_m)  # Distance from most aft CG
    H_strike = (
        ac_datafile["Geometry"]["fus_length_m"] - ac_datafile["Geometry"]["tail_length_m"] - (aftcg + l_m)
    ) * np.tan(np.radians(SP)) - 0.5 * Dw_m

    # Iterate over strut height so both tip-back and scrape angle requirements are met
    i = 1.0
    while i > 0.0001:
        H_s = H_strike
        l_m = tan(np.radians(TB)) * (Hcg_OEW + H_s + 0.5 * Dw_m)
        H_strike = (
            ac_datafile["Geometry"]["fus_length_m"] - ac_datafile["Geometry"]["tail_length_m"] - (aftcg + l_m)
        ) * np.tan(np.radians(SP)) - 0.5 * Dw_m
        i = abs(H_s / H_strike - 1)
        # Set minimum strut height as 60% of wheel diameter for clearance
        if H_strike < 0.6 * Dw_m:
            H_strike = 0.6 * Dw_m
            i = 0

    H_s = H_strike
    # Determine minimum wheel base required for tip-over protection
    Pnw = nose_loading * ac_datafile["CL2Weight"]["MTOW_N"] / 9.81
    l_n = l_m * Pmg / Pnw
    ymin = (l_m + l_n) / (sqrt(l_n**2 * tan(np.radians(55)) ** 2 / (Hcg_OEW + H_s + 0.5 * Dw_m) ** 2 - 1))

    # Determine minimum wheel base for most forward
    if fwdcg:
        l_mfwd = wcg[2, 1] - fwdcg
        l_nfwd = fwdcg - wcg[2, 3]
        ymin_fwd = (l_mfwd + l_nfwd) / (
            sqrt(l_nfwd**2 * tan(np.radians(55)) ** 2 / (Hcg_OEW + H_s + 0.5 * Dw_m) ** 2 - 1)
        )

        # Select limiting wheelbase
        ymin = max([ymin, ymin_fwd])

    # Write geometric landing gear values to dict
    ac_datafile["Landing_gear"]["lm_m"] = l_m
    ac_datafile["Landing_gear"]["ln_m"] = l_n
    ac_datafile["Landing_gear"]["ymin_m"] = ymin
    ac_datafile["Landing_gear"]["Hs_m"] = H_s
    ac_datafile["Landing_gear"]["Dwm_m"] = Dw_m
    ac_datafile["Landing_gear"]["Dwn_m"] = Dw_n
    ac_datafile["Landing_gear"]["Wtm_m"] = Wt_m
    ac_datafile["Landing_gear"]["Wtn_m"] = Wt_n
    ac_datafile["Landing_gear"]["PMW_N"] = Pmw * 9.81
    ac_datafile["Landing_gear"]["PNW_N"] = Pnw * 9.81

    return l_m, l_n, Pmg, Pnw, H_s, Pmw, ymin, Dw_m, Wt_m, Dw_n, Wt_n


def potato_diagrams(ftb_list, btf_list, Xcg_OEW, xlemac, name, ac_datafile=aircraft_data, plot=False):
    # Initiate list of masses and Xcg with aircraft empty mass (without pilot)
    masslist_1 = [ac_datafile["CL2Weight"]["OEW"] / 9.81 - ac_datafile["CL2Weight"]["W_pilot"] / 9.81]
    CGlist_1 = [Xcg_OEW]

    # Determine CG excursion for the front-to-back loading of imported loading case: serial weighted average
    for i in range(ftb_list.shape[1]):
        CGlist_1.append(
            (masslist_1[i] * CGlist_1[i] + ftb_list[0, i] * ftb_list[1, i]) / (masslist_1[i] + ftb_list[1, i])
        )
        masslist_1.append(masslist_1[i] + ftb_list[1, i])

    # Initiate list of masses and Xcg with aircraft empty mass (without pilot)
    masslist_2 = [ac_datafile["CL2Weight"]["OEW"] / 9.81 - ac_datafile["CL2Weight"]["W_pilot"] / 9.81]
    CGlist_2 = [Xcg_OEW]

    # Determine CG excursion for the back-to-front loading of imported loading case: serial weighted average
    for i in range(btf_list.shape[1]):
        CGlist_2.append(
            (masslist_2[i] * CGlist_2[i] + btf_list[0, i] * btf_list[1, i]) / (masslist_2[i] + btf_list[1, i])
        )
        masslist_2.append(masslist_2[i] + btf_list[1, i])

    # Find most forward and aft CG from all points, apply 2% margin for movement, and convert to %MAC
    CGfwd = min(CGlist_1 + CGlist_2) - 0.02 * ac_datafile["Aero"]["MAC_wing"]
    CGaft = max(CGlist_1 + CGlist_2) + 0.02 * ac_datafile["Aero"]["MAC_wing"]
    CGfwdMAC = (CGfwd - xlemac) / ac_datafile["Aero"]["MAC_wing"]
    CGaftMAC = (CGaft - xlemac) / ac_datafile["Aero"]["MAC_wing"]

    # Plot the potato diagram for the selected loading case: the name is imported
    if plot:
        plt.plot(
            ((np.array(CGlist_1) - xlemac) / ac_datafile["Aero"]["MAC_wing"]),
            masslist_1,
            marker="o",
            color="red",
            markersize=4,
        )
        plt.plot(
            ((np.array(CGlist_2) - xlemac) / ac_datafile["Aero"]["MAC_wing"]),
            masslist_2,
            marker="o",
            color="blue",
            markersize=4,
        )
        plt.plot([CGfwdMAC, CGfwdMAC], [0, max(masslist_1)], color="dimgrey", linewidth=1, linestyle="dashdot")
        plt.plot([CGaftMAC, CGaftMAC], [0, max(masslist_1)], color="dimgrey", linewidth=1, linestyle="dashdot")
        plt.ylabel(r"Mass [kg]")
        plt.xlabel(r"$X_{cg}$ [MAC]")
        plt.title(name)
        plt.xlim(floor(CGfwdMAC * 20) / 20, ceil(CGaftMAC * 20) / 20)
        plt.ylim(masslist_1[0], 1.1 * (max(masslist_1 + masslist_2) - masslist_1[0]) + masslist_1[0])
        plt.tight_layout()
        plt.show()

    return [CGfwd, CGaft], (masslist_1 + masslist_2), (CGlist_1 + CGlist_2)


def cg_excursion(Xcg_OEW, xlemac, ac_datafile=aircraft_data, plot=False):
    # TODO: Update Xcg with fuselage sizing values and fuel tank position
    # Define passenger Xcg positions from the nose, 2 passengers per row
    Xcg_row1 = 4.02
    Xcg_row2 = 5.11
    Xcg_row3 = 6.20

    # Define pilot/cargo/fuel Xcg from the nose
    Xcg_pilot = 2.922
    Xcg_luggage = 7
    Xcg_cargomax = 6.5  # Maximum aft position of full cargo load -> dangerous goods
    Xcg_fuel = xlemac + 0.4 * ac_datafile["Aero"]["MAC_wing"]  # ADSEE tells me fuel is at 40% MAC

    # Import mass values for the various types of payload
    m_pax = ac_datafile["CL2Weight"]["Passenger Mass"]  # person without luggage
    m_lug = ac_datafile["CL2Weight"]["Luggage Mass"]  # luggage portion of the passenger mass budget
    m_pil = ac_datafile["CL2Weight"]["W_pilot"] / 9.81  # no luggage, all 90kg on seat
    m_pl = ac_datafile["CL2Weight"]["Wpl_w/o_pilot"] / 9.81  # Payload mass without pilot
    m_f = ac_datafile["CL2Weight"]["Wfuel_N"] / 9.81  # Fuel mass

    # Set up front-to-back and back-to-front passenger loading mass and CG sequence
    Pax_ftb = np.array(
        [[Xcg_row1, Xcg_row1, Xcg_row2, Xcg_row2, Xcg_row3, Xcg_row3], [m_pax, m_pax, m_pax, m_pax, m_pax, m_pax]]
    )
    Pax_btf = np.array(
        [[Xcg_row3, Xcg_row3, Xcg_row2, Xcg_row2, Xcg_row1, Xcg_row1], [m_pax, m_pax, m_pax, m_pax, m_pax, m_pax]]
    )

    # Set up mass and CG sequence for the remaining payload items
    pil = np.array([[Xcg_pilot], [m_pil]])
    lug = np.vstack((Xcg_luggage * np.ones(6), m_lug * np.ones(6)))
    pl = np.array([[Xcg_cargomax], [m_pl]])
    fl = np.array([[Xcg_fuel], [m_f]])

    # Define all possible (limiting) loading situations, pilot always goes with passengers
    fcp_1 = np.hstack((fl, lug, pil, Pax_ftb))
    fcp_2 = np.hstack((fl, lug, Pax_btf, pil))

    fpc_1 = np.hstack((fl, pil, Pax_ftb, lug))
    fpc_2 = np.hstack((fl, Pax_btf, pil, lug))

    cpf_1 = np.hstack((lug, pil, Pax_ftb, fl))
    cpf_2 = np.hstack((lug, Pax_btf, pil, fl))

    cfp_1 = np.hstack((lug, fl, pil, Pax_ftb))
    cfp_2 = np.hstack((lug, fl, Pax_btf, pil))

    pcf_1 = np.hstack((pil, Pax_ftb, lug, fl))
    pcf_2 = np.hstack((Pax_btf, pil, lug, fl))

    pfc_1 = np.hstack((pil, Pax_ftb, fl, lug))
    pfc_2 = np.hstack((Pax_btf, pil, fl, lug))

    pmpl_1 = np.hstack((pil, pl))
    pmpl_2 = np.hstack((pl, pil))

    # Get CG excursion for each loading situation from potato diagrams
    CGlist_fcp, masslist1, fullCG1 = potato_diagrams(
        fcp_1, fcp_2, Xcg_OEW, xlemac, "Fuel/Luggage/Passengers", ac_datafile, plot
    )
    CGlist_fpc, masslist2, fullCG2 = potato_diagrams(
        fpc_1, fpc_2, Xcg_OEW, xlemac, "Fuel/Passengers/Luggage", ac_datafile, plot
    )
    CGlist_cpf, masslist3, fullCG3 = potato_diagrams(
        cpf_1, cpf_2, Xcg_OEW, xlemac, "Luggage/Passengers/Fuel", ac_datafile, plot
    )
    CGlist_cfp, masslist4, fullCG4 = potato_diagrams(
        cfp_1, cfp_2, Xcg_OEW, xlemac, "Luggage/Fuel/Passengers", ac_datafile, plot
    )
    CGlist_pcf, masslist5, fullCG5 = potato_diagrams(
        pcf_1, pcf_2, Xcg_OEW, xlemac, "Passengers/Luggage/Fuel", ac_datafile, plot
    )
    CGlist_pfc, masslist6, fullCG6 = potato_diagrams(
        pfc_1, pfc_2, Xcg_OEW, xlemac, "Passengers/Fuel/Luggage", ac_datafile, plot
    )
    CGlist_pmpl, masslist7, fullCG7 = potato_diagrams(
        pmpl_1, pmpl_2, Xcg_OEW, xlemac, "Full Cargo Load", ac_datafile, plot
    )

    # Combine lists for easy exporting: CG is only forward and aft, fullCG is all encountered CG positions
    Combined_CG = CGlist_fcp + CGlist_fpc + CGlist_cpf + CGlist_cfp + CGlist_pcf + CGlist_pfc + CGlist_pmpl
    Combined_mass = masslist1 + masslist2 + masslist3 + masslist4 + masslist5 + masslist6 + masslist7
    Combined_FullCG = fullCG1 + fullCG2 + fullCG3 + fullCG4 + fullCG5 + fullCG6 + fullCG7

    return Combined_CG, Combined_mass, Combined_FullCG


def TailAero_copy(l_H, acd=aircraft_data):  # pragma: no cover
    # (copy from FullStability.py)
    # Get and convert values to imp*rial
    taper_h = acd["Aero"]["Taper_HS"]
    Cr_h = 3 / 2 * acd["Aero"]["MAC_HS"] * ((1 + taper_h) / (1 + taper_h * taper_h**2))  # Root chord of stabiliser
    Cpw = Cr_h * (1 - (1 - taper_h) * acd["Power_prop"]["Dp_m"] / acd["Aero"]["b_h"])  # Chord at prop wash edge
    Shslip = (Cr_h + Cpw) * acd["Power_prop"]["Dp_m"] / 2 * 10.7639104  # Tail area affected by prop wash
    Sh = acd["Aero"]["S_h"] * 10.7639104
    U1 = acd["Performance"]["Vc_m/s"] * 3.2808399
    Pav = 0.98 * 400 * acd["Power_prop"]["eta_p"]
    Dp = acd["Power_prop"]["Dp_m"] / 0.3048
    qbar = (
        0.5
        * isa(acd["Performance"]["Altitude_Cruise_m"], acd["Performance"]["Temp_offset_TO_Land_cruise"])[2]
        * acd["Performance"]["Vc_m/s"] ** 2
    ) * 0.020885

    # Magic Roskam equation, only god knows how this works
    eta_H = 1 + (Shslip / Sh) * (2200 * Pav / (qbar * U1 * np.pi * Dp**2))
    VhVcorr = acd["Aero"]["VhV"] * eta_H

    # Downwash Gradient supposedly from Slingerland whoever that may be
    r = float(l_H * 2 / acd["Aero"]["b_Wing"])
    mtv = float(0 * 2 / acd["Aero"]["b_Wing"])
    deda = (
        (
            (r / (r**2 + mtv**2)) * 0.4876 / np.sqrt(r**2 + 0.6319 + mtv**2)
            + (1 + (r**2 / (r**2 + 0.7915 + 5.0734 * mtv**2)) ** 0.3113) * (1 - np.sqrt(mtv**2 / (1 + mtv**2)))
        )
        * acd["Aero"]["CLalpha"]
        / (np.pi * acd["Aero"]["AR"])
    )

    return VhVcorr, deda


def Rotation_check(Weightlist, wcg, xlemac, FullCGlist, acd=aircraft_data):
    # Get air density at 750m take-off altitude and ISA offset
    rho = isa(acd["Performance"]["Altitude_TO_m"], acd["Performance"]["Temp_offset_TO_Land_cruise"])[2]

    # Get required values for magic formula that is 100% verified correct, trust me
    Vto = acd["Performance"]["V_takeoff_m/s"]  # Take-off rotation velocity
    S = acd["Aero"]["S_Wing"]  # Wing surface area
    Sh = acd["Aero"]["S_h"]  # Tail surface area
    l_w = wcg[2, 1] - (xlemac + 1 / 4 * acd["Aero"]["MAC_wing"])  # Wing C.P. distance to aft gear (point of rotation)
    lH_lg = acd["Stability"]["QCW_to_QCh"] - l_w  # Stabiliser distance to aft gear
    H_pwtrT = (
        0.5 * acd["Geometry"]["fus_height_m"] + acd["Landing_gear"]["Hs_m"] + 0.5 * acd["Landing_gear"]["Dwm_m"]
    )  # Height of thrust application
    Tto = 6003

    # Find tail effectiveness
    VhVcorr, deda = TailAero_copy(acd["Stability"]["QCW_to_QCh"])

    # Find moment coefficient of the wing + flaps deployed for take-off
    Cmacw = acd["Aero"]["Cm_0_wing"] * (
        acd["Aero"]["AR"]
        * np.cos(acd["Aero"]["QuarterChordSweep_Wing_deg"]) ** 2
        / (acd["Aero"]["AR"] + 2 * np.cos(acd["Aero"]["QuarterChordSweep_Wing_deg"]))
    )
    Cmacflap = acd["Flaps"]["mu2_takeoff"] * (
        -acd["Flaps"]["mu1_takeoff"] * acd["Flaps"]["deltaCLmax_takeoff"] * acd["Flaps"]["cprime_c_takeoff"]
        - (
            acd["Flaps"]["CL_AoA0_takeoff"]
            + acd["Flaps"]["deltaCLmax_takeoff"] * (1 - acd["Flaps"]["Swf"] / acd["Aero"]["S_Wing"])
        )
        * 1
        / 8
        * acd["Flaps"]["cprime_c_takeoff"]
        * (acd["Flaps"]["cprime_c_takeoff"] - 1)
    )
    Cmw_TO = Cmacw + Cmacflap

    # Tail lift rate
    SweepHS_05 = np.tan(np.radians(acd["Aero"]["QuarterChordSweep_HS_deg"])) - 4 / acd["Aero"]["AR_HS"] * (
        0.25 * (1 - acd["Aero"]["Taper_HS"]) / (1 + acd["Aero"]["Taper_HS"])
    )
    ClaH = (
        2
        * np.pi
        * acd["Aero"]["AR_HS"]
        / (2 + sqrt(4 + (acd["Aero"]["AR_HS"] / 0.95) ** 2 * (1 + tan(SweepHS_05) ** 2)))
    )

    # Magic formula that is 100% correct: outputs maximum allowable nose loading for weight conditions
    # Just pray it works, derive it yourself if you don't believe me
    tail_lift = []
    # Weightlist in [kg]
    for i in range(len(Weightlist)):
        tail_lift.append(
            1
            / lH_lg
            * (
                acd["Flaps"]["CL_AoA0_takeoff"] * 1 / 2 * rho * Vto**2 * S * l_w
                - Tto * H_pwtrT
                - Weightlist[i] * 9.81 * (wcg[2, 1] - FullCGlist[i])
                + Cmw_TO * 1 / 2 * rho * Vto**2 * S * acd["Aero"]["MAC_wing"]
            )
        )

    # Find Elevator effectiveness
    CLH_corr = 2 * np.array(tail_lift) / (rho * Sh * (VhVcorr * Vto) ** 2)
    e0 = deda * 6.5  # Downwash at zero angle of attack [deg]
    a_H = acd["Stability"]["i_H_deg"] - e0
    tau_e = np.abs(
        (np.min(CLH_corr) + ClaH * np.radians(a_H)) / (ClaH * np.radians(acd["Stability"]["Elevator_deflection_deg"]))
    )

    # Find sizing through fitted polynomial
    x = symbols("x")
    CeC = solve(-6.624 * x**4 + 12.07 * x**3 - 8.292 * x**2 + 3.295 * x + 0.004942 - tau_e, x)
    CeC = min([n for n in CeC if n.is_real])
    Elev = True
    if CeC > 0.4:
        Elev = False
        # print("WARNING: Your aircraft did not take off and crashed into a school")

    # Update elevator sizing
    acd["Stability"]["Hinge_chord_elevator"] = ceil((1 - CeC) * 20) / 20
    acd["Stability"]["Hinge_chord_elevator"]

    return Elev


def iterate_cg_lg(ac_datafile=aircraft_data, PERCENTAGE=0.2, bat_xcg=0.5, plot=False):
    # Set desired distance of nosewheel from nose [m]
    nose_distance = ac_datafile["Landing_gear"]["Nose_gear_distance"]
    nose_loading = 0.08  # Initial guess and minimum required

    # Get components fractions and weights
    wcg = component_mass(ac_datafile)

    # Get preliminary wing CG location with ol X LEMAC
    CGw_MAC = 0.4 * ac_datafile["Aero"]["MAC_wing"]
    wcg[2, 0] = ac_datafile["Geometry"]["XLEMAC_m"] + CGw_MAC  # Wing distance from tip

    # Find powertrain CG from the nose: Weighted average
    Weng = ac_datafile["CL2Weight"]["Engine Weight"]
    Wprop = ac_datafile["CL2Weight"]["Propeller Weight"]
    Wmotor = ac_datafile["CL2Weight"]["Electric Motor Weight"]

    wcg[2, 2] = (
        Weng * ac_datafile["Stability"]["Xcg_engine_m"]
        + Wprop * ac_datafile["Stability"]["Xcg_prop_m"]
        + Wmotor * ac_datafile["Stability"]["Xcg_motor_m"]
    ) / (Weng + Wprop + Wmotor)

    # Guesstimate of other component CG locations
    wcg[2, 4] = 0.4 * ac_datafile["Geometry"]["fus_length_m"]  # Fuselage distance from tip
    wcg[2, 5] = 0.9 * ac_datafile["Geometry"]["fus_length_m"]  # Empennage distance from tip
    wcg[2, 6] = 0.4 * ac_datafile["Geometry"]["fus_length_m"]  # Fixed equipment distance from tip
    wcg[2, 7] = bat_xcg * ac_datafile["Geometry"]["fus_length_m"]  # Battery distance from tip

    # Weighted average calculation for the empty weight CG, now still excludes landing gear
    Xcg_OEW = (
        (ac_datafile["Geometry"]["XLEMAC_m"] + CGw_MAC) * wcg[1, 0]
        + np.average(wcg[2, [2, 4, 5, 6, 7]], weights=wcg[1, [2, 4, 5, 6, 7]]) * np.sum(wcg[1, [2, 4, 5, 6, 7]])
    ) / (wcg[1, 0] + np.sum(wcg[1, [2, 4, 5, 6, 7]]))
    wcg[2, -1] = Xcg_OEW

    # Get preliminary CG excursion from the nose
    xlemac = ac_datafile["Geometry"]["XLEMAC_m"]
    aftcg = max(cg_excursion(Xcg_OEW, xlemac)[0])

    # First guess for landing gear placement before going into iteration: CG locations placed at wheel location
    l_m, l_n, Pmg, Pnw, H_s = find_lg(nose_loading, aftcg, wcg, ac_datafile)[0:5]  # Nose loading of 8% initially
    wcg[2, 1] = aftcg + l_m
    wcg[2, 3] = aftcg - l_n
    # Iterate on CG and X LEMAC positions
    iter = 1.0

    while iter > 0.0001:  # convergence criterion
        # Save old X LEMAC for convergence calculation
        xlemacold = xlemac

        # Determine new empty weight CG including landing gear and updated X LEMAC
        Xcg_OEW = np.average(wcg[2, 0:8], weights=wcg[1, 0:8])

        # Find most aft CG for current configuration and X LEMAC
        CGlist = cg_excursion(Xcg_OEW, xlemac)[0]
        aftcg = max(CGlist)

        # Revise nosewheel loading and placement in case wheel is too far forward; in front of the nose
        if wcg[2, 3] < nose_distance:
            wcg[2, 3] = nose_distance
            l_n = aftcg - wcg[2, 3]
            nose_loading = 1 / (l_n / l_m + 1)

        # Place nosewheel at set distance from nose in case it's behind; Makes it easier to fit in fuselage
        l_m, l_n, Pmg, Pnw, H_s = find_lg(nose_loading, aftcg, wcg, ac_datafile)[0:5]
        wcg[2, 1] = aftcg + l_m
        wcg[2, 3] = aftcg - l_n
        if wcg[2, 3] > nose_distance:
            wcg[2, 3] = nose_distance
            l_n = aftcg - wcg[2, 3]
            l_m = l_n * Pnw / Pmg
            wcg[2, 1] = aftcg + l_m

        # Update X LEMAC with current landing gear location
        wcg[2, -1] = Xcg_OEW
        cgwg = np.average(wcg[2, 0:2] - xlemac, weights=wcg[1, 0:2])  # wing group cg location
        xlemac = np.average(wcg[2, 2:8], weights=wcg[1, 2:8]) + ac_datafile["Aero"]["MAC_wing"] * (
            (cgwg / ac_datafile["Aero"]["MAC_wing"]) * np.sum(wcg[1, 0:2]) / np.sum(wcg[1, 2:8])
            - PERCENTAGE * (1 + np.sum(wcg[1, 0:2]) / np.sum(wcg[1, 2:8]))
        )

        # Update wing CG and convergence parameter
        wcg[2, 0] = CGw_MAC + xlemac
        iter = abs(xlemacold / xlemac - 1)

    # Check nose wheel loading for forward CG location (not larger than rotation maximum)
    # Import full data on CG excursion for converged configuration
    CGlist, masslist, FullCGlist = cg_excursion(Xcg_OEW, xlemac, ac_datafile, plot=plot)
    nose_combined = np.vstack((np.array(FullCGlist), np.array(masslist)))

    # Check whether wheel pressure is below the maximum pressure for all conditions
    Pressure_list = []
    for i in range(nose_combined.shape[1]):
        l_mfwd = wcg[2, 1] - nose_combined[0, i]
        l_nfwd = nose_combined[0, i] - wcg[2, 3]

        # Find and compare pressure fraction: This formula works, trust
        nose_load_fwd = 1 / (l_nfwd / l_mfwd + 1)
        Pressure_list.append(nose_load_fwd * nose_combined[1, i])

    # Update nose gear dimensions for additional pressure
    Wheel_pressure = max(Pressure_list)
    find_lg(nose_loading, max(CGlist), wcg, ac_datafile, max_load_nose=Wheel_pressure, fwdcg=min(FullCGlist))[0:5]

    # Check Rotation possibility
    Rotate = np.array(Rotation_check(masslist, wcg, xlemac, FullCGlist))
    # Update dictionary values
    ac_datafile["Geometry"]["XLEMAC_m"] = xlemac
    ac_datafile["Landing_gear"]["Xmw_m"] = wcg[2, 1]
    ac_datafile["Landing_gear"]["Xnw_m"] = wcg[2, 3]

    return wcg, CGlist, xlemac, Rotate


# iterate though lemac such that the x_lemac is larger than 3.2m and that the landing gear can fit with the batteries
def optimised_xlemac_landing_gears(ac_data=aircraft_data, percentage=0.2, bat_xcg_init=0.2, lemac_limit=3.2):
    # initialise the sizing parameter
    sizing = False
    bat_xcg = bat_xcg_init

    # loop to get the optimised lemac and landing gear data with the batteries
    # in the right position and the lemac larger than 3.2m
    while not sizing:
        _, CGlist, xlemac, rotate = iterate_cg_lg(ac_datafile=ac_data, PERCENTAGE=percentage, bat_xcg=bat_xcg)

        # check if the nose gear can be placed
        if not rotate:
            # print("Nose gear collapsed, increasing battery xcg")
            bat_xcg += 0.01

        # check if the xlemac is acceptable
        elif xlemac > lemac_limit and rotate:
            # TODO: call the function from fuselage sizing and check if the boxes overlap and returns the lemac

            # get the position of the batteries and the landing gears
            fuselage_sizing = FuselageSizing(ac_data=ac_data, bat_xcg=bat_xcg)
            bellow_position = fuselage_sizing.below_position(s_gear=0.1)  # s_gear is the clearance
            # print(bellow_position)
            # check if the sizings are not overlapping
            if (
                bellow_position["nose landing gear"][1]
                < bellow_position["battery"][0]
                # and bellow_position["battery"][1] < bellow_position["main landing gear"][0]
            ):
                # print("DA")
                # if the sizing is correct, break the loop and return the optimised xlemac
                # and update the xcg of the batteries
                sizing = True

                # update the lemac and the xcg of the batteries
                ac_data["Geometry"]["XLEMAC_m"] = xlemac
                ac_data["Stability"]["Xcg_battery_m"] = bat_xcg * ac_data["Geometry"]["fus_length_m"]

                # update the cg limits
                ac_data["Stability"]["Cg_Aft"] = (max(CGlist) - xlemac) / ac_data["Aero"]["MAC_wing"]
                ac_data["Stability"]["Cg_Front"] = (min(CGlist) - xlemac) / ac_data["Aero"]["MAC_wing"]

            # check if the battery fits after the main landing gear
            # elif bellow_position["main landing gear"][1] + 0.05 < bellow_position["battery"][0]:
            #     #print("NU")
            #     # if the sizing is correct, break the loop and return the optimised xlemac
            #     # and update the xcg of the batteries
            #     sizing = True
            #
            #     # update the lemac and the xcg of the batteries
            #     ac_data["Geometry"]["XLEMAC_m"] = xlemac
            #     ac_data["Stability"]["Xcg_battery_m"] = bat_xcg * ac_data["Geometry"]["fus_length_m"]
            #
            #     # update the cg limits
            #     ac_data["Stability"]["Cg_Aft"] = (max(CGlist) - xlemac) / ac_data["Aero"]["MAC_wing"]
            #     ac_data["Stability"]["Cg_Front"] = (min(CGlist) - xlemac) / ac_data["Aero"]["MAC_wing"]

            # if the sizing is not correct, increase the battery xcg
            else:
                bat_xcg += 0.01

        # if the xlemac is not acceptable, increase the battery xcg
        elif xlemac < lemac_limit:
            bat_xcg += 0.01

        # if the xlemac is acceptable, but no rotation, increase the battery xcg
        else:
            bat_xcg += 0.01

        # if it gets over an unrealistic value, break the loop
        if bat_xcg > 0.8:
            break

    return bat_xcg


def calculate_lh(ac_data=aircraft_data):
    # lh is defined as the distance from quarter chord location of the wing to the quarter chord
    # location of the horizontal tail
    QCW_mac = ac_data["Geometry"]["XLEMAC_m"] + 0.25 * ac_data["Aero"]["MAC_wing"]
    hinge_chord_percentage = ac_data["Stability"]["Hinge_chord_elevator"]

    # get the horizontal stabiliser data from the aircraft data
    AR_h = ac_data["Aero"]["AR_HS"]
    taper_h = ac_data["Aero"]["Taper_HS"]
    # c_root_h = ac_data["Aero"]["c_root_HS"]
    b_h = ac_data["Aero"]["b_h"]

    # calculate the leading edge angle of the horizontal stabiliser and the x lemac
    tan_LE_sweep = tan(0) - 4 / AR_h * ((-hinge_chord_percentage) * (1 - taper_h) / (1 + taper_h))

    # calculate where the mac of the horizontal stabiliser wrt the leading edge
    y_mac_h = b_h / 6 * (1 + 2 * taper_h) / (1 + taper_h)
    x_mac_h = y_mac_h * tan_LE_sweep

    # get the quarter chord location of the horizontal stabiliser
    QCH_mac = (
        x_mac_h + 0.25 * ac_data["Aero"]["MAC_HS"] + ac_data["Geometry"]["fus_length_m"] - ac_data["Aero"]["c_root_HS"]
    )

    # update the aircraft data with the new lh

    ac_data["Stability"]["QCW_to_QCh"] = QCH_mac - QCW_mac
    return QCH_mac - QCW_mac


if __name__ == "__main__":  # pragma: no cover
    init = time.process_time()
    # print(iterate_cg_lg(aircraft_data, PERCENTAGE=0.2, bat_xcg=0.5, plot=False))
    # calculate_lh(ac_data=aircraft_data)
    bat_xcg = optimised_xlemac_landing_gears(ac_data=aircraft_data, percentage=0.1, bat_xcg_init=0.1, lemac_limit=3.2)
    total = time.process_time() - init
    print(total)
