import sys
import os

import numpy as np
from matplotlib import pyplot as plt

"Unit test in progress"

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from HumanAir.AerodynamicDesign.Aerodynamics_Main import aerodynamic_design
from HumanAir.aircraft_data import aircraft_data
from HumanAir.unit_conversions import m_to_ft, N_to_lbs, m_squared_to_ft_squared, m_s_to_kt, W_to_hp, lbs_to_N, ft_to_m


class Class_II_Weight:
    def __init__(self, ac_data):
        self.dict = ac_data
        self.W_TO = N_to_lbs(ac_data["Weights"]["MTOW_N"]) / ac_data["Contingency"]
        self.W_L = N_to_lbs(ac_data["Weights"]["W_L_N"]) / ac_data["Contingency"]
        self.W_F = N_to_lbs(ac_data["Weights"]["Wfuel_N"]) / ac_data["Contingency"]
        self.W_E = N_to_lbs(ac_data["Weights"]["OEW_N"]) / ac_data["Contingency"]
        self.W_bat = N_to_lbs(ac_data["Weights"]["Wbat_N"]) / ac_data["Contingency"]

        self.S_Wing = m_squared_to_ft_squared(ac_data["Aero"]["S_Wing"])
        self.S_h = m_squared_to_ft_squared(ac_data["Aero"]["S_h"])
        self.S_v = m_squared_to_ft_squared(ac_data["Aero"]["S_v"])

        self.n_ult = ac_data["Performance"]["n_ult"]
        self.n_ult_l = ac_data["Performance"]["n_ult_l"]

        self.AR_Wing = ac_data["Aero"]["AR"]
        self.AR_h = ac_data["Aero"]["AR_HS"]
        self.AR_v = ac_data["Aero"]["AR_v"]

        self.QuarterChordSweep_Wing = np.deg2rad(ac_data["Aero"]["QuarterChordSweep_Wing_deg"])
        self.HalfChordSweep_Wing = np.deg2rad(ac_data["Aero"]["HalfChordSweep_Wing_deg"])
        self.QuarterChordSweep_v = np.deg2rad(ac_data["Aero"]["QuarterChordSweep_v_deg"])

        self.Taper_Wing = ac_data["Aero"]["Taper_Wing"]
        self.tc_m_Wing = m_to_ft(ac_data["Aero"]["tc_m_Wing"])

        self.b_Wing = m_to_ft(ac_data["Aero"]["b_Wing"])
        self.b_h = m_to_ft(ac_data["Aero"]["b_h"])
        self.b_v = m_to_ft(ac_data["Aero"]["b_v"])

        self.t_root_max_Wing = m_to_ft(ac_data["Aero"]["t_root_max_Wing"])
        self.t_root_max_h = m_to_ft(ac_data["Aero"]["t_root_max_h"])
        self.t_root_max_v = m_to_ft(ac_data["Aero"]["t_root_max_v"])

        self.V_H = m_s_to_kt(ac_data["Performance"]["Vh_m/s"])
        self.V_c = m_s_to_kt(ac_data["Performance"]["Vc_m/s"])
        self.M_D = ac_data["Performance"]["M_D"]
        self.QCW_to_QCh = m_to_ft(ac_data["Stability"]["QCW_to_QCh"])
        self.l_f_nonosecone = m_to_ft(ac_data["Geometry"]["l_f_nonosecone"])
        self.paint = ac_data["General"]["Paint"]
        self.p_max = m_to_ft(ac_data["Geometry"]["fuselage_max_perimeter"])
        self.N_pax = ac_data["General"]["N_pax"]
        self.N_row = ac_data["General"]["N_row"]
        self.l_f = m_to_ft(ac_data["Geometry"]["fus_length_m"])
        self.w_f = m_to_ft(ac_data["Geometry"]["fus_width_m"])
        self.h_f = m_to_ft(ac_data["Geometry"]["fus_height_m"])

        self.P_TO = W_to_hp(ac_data["Power_prop"]["P_req_TO_W"])
        self.K_n = ac_data["Power_prop"]["K_n"]  # 0.37 for radial, 0.24 for horizontally opposed
        self.K_p = 1.1  # TODO: check if supercharging is used (pg.84)
        self.K_pg = 1.16
        self.K_fsp = 6.65  # TODO: check if this is the correct fuel in lbs/gal (pg.91)

        self.int = ac_data["Power_prop"]["int_fueltanks_fraction"]

        ac_data["Landing_gear"]["l_s_m"] = ac_data["Landing_gear"]["Hs_m"]
        ac_data["Landing_gear"]["l_s_n"] = ac_data["Landing_gear"]["Hs_m"]
        self.l_s_m = m_to_ft(ac_data["Landing_gear"]["l_s_m"])
        self.l_s_n = m_to_ft(ac_data["Landing_gear"]["l_s_n"])
        self.retractable = ac_data["Landing_gear"]["Retractable"]

        self.N_e = ac_data["Power_prop"]["N_e"]
        self.N_t = ac_data["Power_prop"]["N_t"]

        self.PoweredFlightControls = ac_data["General"]["PoweredFlightControls"]
        self.DuplicatedFlightControls = ac_data["General"]["DuplicatedFlightControls"]

        self.APU = ac_data["General"]["APU"]
        self.HydrPne = ac_data["General"]["HydrPne"]
        self.NacWght = ac_data["General"]["NacWght"]

        """========== Structure Weight =========="""

    def WingWeight(self):
        results = {}

        results["Cessna"] = 0.002933 * self.S_Wing**1.018 * self.AR_Wing**2.473 * self.n_ult**0.611
        results["USAF"] = (
            96.948
            * (
                (self.W_TO * self.n_ult * 10 ** (-5)) ** 0.65
                * (self.AR_Wing / np.cos(self.QuarterChordSweep_Wing)) ** 0.57
                * (self.S_Wing / 100) ** 0.61
                * ((1 + self.Taper_Wing) / 2 * self.tc_m_Wing) ** 0.36
                * (1 + self.V_H / 500) ** 0.5
            )
            ** 0.993
        )
        results["Torenbeek"] = (
            0.00125
            * self.W_TO
            * (self.b_Wing / np.cos(self.HalfChordSweep_Wing)) ** 0.75
            * (1 + (6.3 * np.cos(self.HalfChordSweep_Wing) / self.b_Wing) ** 0.5)
            * self.n_ult**0.55
            * (self.b_Wing * self.S_Wing / (self.t_root_max_Wing * self.W_TO * np.cos(self.HalfChordSweep_Wing)))
            ** 0.30
        )

        # results["Average"] = np.average([results["Cessna"], results["USAF"], results["Torenbeek"]])
        results["Average"] = 760
        return results

    def EmpennageWeight(self):
        results = {}
        """ Cessna """
        results["Cessna"] = {}
        results["Cessna"]["W_h"] = (3.184 * self.W_TO**0.887 * self.S_h**0.101 * self.AR_h**0.138) / (
            174.04 * self.t_root_max_h**0.223
        )
        results["Cessna"]["W_v"] = (1.68 * self.W_TO**0.567 * self.S_v**0.1249 * self.AR_v**0.482) / (
            639.95 * self.t_root_max_v**0.747 * (np.cos(self.QuarterChordSweep_v)) ** 0.882
        )
        results["Cessna"]["W_c"] = 0
        results["Cessna"]["Total"] = results["Cessna"]["W_h"] + results["Cessna"]["W_v"] + results["Cessna"]["W_c"]

        """ USAF """
        results["USAF"] = {}
        results["USAF"]["W_h"] = (
            127
            * (
                (self.W_TO * self.n_ult * 10 ** (-5)) ** 0.87
                * (self.S_h / 100) ** 1.2
                * 0.289
                * (self.QCW_to_QCh / 10) ** 0.483
                * (self.b_h / self.t_root_max_h) ** 0.5
            )
            ** 0.458
        )
        results["USAF"]["W_v"] = (
            98.5
            * (
                (self.W_TO * self.n_ult * 10 ** (-5)) ** 0.87
                * (self.S_v / 100) ** 1.2
                * 0.289
                * (self.b_v / self.t_root_max_v) ** 0.5
            )
            ** 0.458
        )
        results["USAF"]["W_c"] = 0
        results["USAF"]["Total"] = results["USAF"]["W_h"] + results["USAF"]["W_v"] + results["USAF"]["W_c"]

        """ Torenbeek """
        results["Torenbeek"] = 0.04 * (self.n_ult * (self.S_v + self.S_h) ** 2) ** 0.75

        results["Average"] = np.average([results["Cessna"]["Total"], results["USAF"]["Total"], results["Torenbeek"]])

        return results

    def FuselageWeight(self):
        results = {}

        results["Cessna"] = (
            14.86
            * self.W_TO**0.144
            * (self.l_f_nonosecone / self.p_max) ** 0.778
            * self.l_f_nonosecone**0.383
            * self.N_pax**0.455
        )
        results["USAF"] = (
            200
            * (
                (self.W_TO * self.n_ult * 10 ** (-5)) ** 0.286
                * (self.l_f / 10) ** 0.857
                * ((self.w_f + self.h_f) / 10)
                * (self.V_c / 100) ** 0.338
            )
            ** 1.1
        )

        results["Average"] = np.average([results["Cessna"], results["USAF"]])

        return results

    def NacelleWeight(self):
        results = {}

        results["Cessna"] = self.P_TO * self.K_n
        results["Torenbeek"] = 2.5 * self.P_TO**0.5

        if self.NacWght:
            results["Average"] = np.average([results["Cessna"], results["Torenbeek"]])
        else:
            results["Average"] = 0
        return results

    def LandingGearWeight(self):
        results = {}

        results["Cessna"] = (
            0.013 * self.W_TO
            + 0.362 * self.W_L**0.417 * self.n_ult_l**0.950 * self.l_s_m**0.183
            + 6.2
            + 0.0013 * self.W_TO
            + 0.007157 * self.W_L**0.749 * self.n_ult_l * self.l_s_n**0.788
        )
        if self.retractable:
            results["Cessna"] += 0.014 * self.W_TO

        results["USAF"] = 0.054 * self.l_s_m**0.501 * (self.W_L * self.n_ult_l) ** 0.684

        #  Can easily change below between Cessna twice or Cessna and USAF
        results["Average"] = np.average([results["Cessna"], results["Cessna"]])

        self.l_s_n = ft_to_m(self.l_s_n)
        self.l_s_m = ft_to_m(self.l_s_m)
        return results

    def StructureWeight_Total(self):
        self.dict["CL2Weight"]["Total Structures Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.WingWeight()["Average"]
            + self.EmpennageWeight()["Average"]
            + self.FuselageWeight()["Average"]
            + self.NacelleWeight()["Average"]
            + self.LandingGearWeight()["Average"]
        )
        self.dict["CL2Weight"]["Wing Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(self.WingWeight()["Average"])
        self.dict["CL2Weight"]["Empennage Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.EmpennageWeight()["Average"]
        )
        self.dict["CL2Weight"]["Fuselage Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.FuselageWeight()["Average"]
        )
        self.dict["CL2Weight"]["Nacelle Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.NacelleWeight()["Average"]
        )
        self.dict["CL2Weight"]["Landing Gear Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.LandingGearWeight()["Average"]
        )
        return (
            self.WingWeight()["Average"]
            + self.EmpennageWeight()["Average"]
            + self.FuselageWeight()["Average"]
            + self.NacelleWeight()["Average"]
            + self.LandingGearWeight()["Average"]
        )

    """========== Powerplant Weight =========="""

    def FuelSystemWeight(self):
        results = {}

        results["Cessna"] = 0.40 * self.W_F / self.K_fsp
        results["USAF"] = (
            2.49
            * ((self.W_F / self.K_fsp) ** 0.6 * (1 / (1 + self.int)) ** 0.3 * self.N_t**0.20 * self.N_e**0.13) ** 1.21
        )
        results["Torenbeek"] = 2 * (self.W_F / 5.87) ** 0.667

        results["Average"] = np.average([results["Cessna"], results["USAF"], results["Torenbeek"]])
        return results

    def PowerplantWeight_Total(self):
        self.dict["CL2Weight"]["Total Powerplant Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(1102)
        results = {}

        results["USAF"] = {}
        results["USAF"]["WeWaiWpropWp"] = 2.575 * (self.K_p * self.P_TO) ** 0.922 * self.N_e
        results["USAF"]["Total"] = (
            results["USAF"]["WeWaiWpropWp"] + self.FuelSystemWeight()["USAF"] - self.NacelleWeight()["Average"]
        )

        results["Torenbeek"] = {}
        results["Torenbeek"]["Total"] = self.K_pg * (self.K_p * self.P_TO + 0.24 * self.P_TO)
        results["Average"] = np.average([results["USAF"]["Total"], results["Torenbeek"]["Total"]])
        results["Average"] = 1102

        return results

    """========== Fixed Equipment Weight =========="""

    def FlightControlSystem(self):
        results = {}

        results["Cessna"] = 0.0168 * (self.W_TO - self.W_bat)
        if not self.PoweredFlightControls:
            results["USAF"] = 1.066 * (self.W_TO - self.W_bat) ** 0.626
            if not self.DuplicatedFlightControls:
                results["Torenbeek"] = 0.33 * (self.W_TO - self.W_bat) ** (2 / 3)
        else:
            results["USAF"] = 1.08 * (self.W_TO - self.W_bat) ** 0.7
            results["Torenbeek"] = False

        finalresults = []
        for i in results.values():
            if i is not False:
                finalresults.append(i)
        results["Average"] = np.average(finalresults)
        return results

    def HydraulicsPneumatics(self):
        constant = (0.006 + 0.0120) / 2  # TODO: check which value to take more accurately(pg.101)
        if self.HydrPne:
            return {"Average": constant * (self.W_TO - self.W_bat)}
        else:
            return {"Average": 0.0}

    def InstrumentsAvionicsElectronics(self):
        return {"Average": 40 + 0.008 * (self.W_TO - self.W_bat)}

    def ElectricalSystemWeight(self):
        results = {}
        results["Cessna"] = 0.0268 * self.W_TO
        results["USAF"] = (
            426
            * ((self.FuelSystemWeight()["Average"] + self.InstrumentsAvionicsElectronics()["Average"]) / 1000) ** 0.51
        )

        results["Average"] = np.average([results["Cessna"], results["USAF"]])
        return results

    def AirconPressurizationAntiDeicingWeight(self):
        results = {}

        results["USAF"] = (
            0.265
            * (self.W_TO - self.W_bat) ** 0.52
            * self.N_pax**0.68
            * self.InstrumentsAvionicsElectronics()["Average"] ** 0.17
            * self.M_D**0.08
        )
        if self.N_e > 1:
            results["Torenbeek"] = 0.018 * (self.W_E - self.W_bat)
        else:
            results["Torenbeek"] = 2.5 * self.N_pax

        results["Average"] = np.average([results["USAF"], results["Torenbeek"]])
        return results

    def OxygenSystem(self):
        return {"Average": 20 + 0.5 * self.N_pax}

    def APU_Weight(self):
        results = {}
        if self.APU:
            results["Average"] = 0.0085 * (self.W_TO - self.W_bat)
        else:
            results["Average"] = 0
        return results

    def Furnishings(self):
        results = {}
        results["Cessna"] = 0.412 * self.N_pax**1.145 * (self.W_TO - self.W_bat) ** 0.489
        results["Torenbeek"] = (
            5 + 13 * self.N_pax + 25 * self.N_row
        )  # TODO: check if there is multiengine or not (pg.108)

        results["Average"] = np.average([results["Cessna"], results["Torenbeek"]])
        return results

    def AuxiliaryGear(self):
        return {"Average": 0.01 * (self.W_E - self.W_bat)}

    def Paint(self):
        if self.paint:
            results = {"Average": 0.0045 * (self.W_TO - self.W_bat)}
        else:
            results = {"Average": 0}
        return results

    def NewBatteryWeight(self, bat):
        return (
            lbs_to_N(self.W_TO)
            / self.dict["Performance"]["W/P_N/W"]
            * self.dict["Performance"]["endurance"]
            * bat
            / self.dict["Power_prop"]["E_bat_Wh/kg"]
            / self.dict["Power_prop"]["eta_bat"]
            / self.dict["Power_prop"]["DoD_bat"]
            / self.dict["Power_prop"]["eta_electricmotor"]
        )

    def NewFuelWeight(self, bat):
        return (
            1.15
            * lbs_to_N(self.W_TO)
            / self.dict["Performance"]["W/P_N/W"]
            * (1 - bat)
            * self.dict["Performance"]["endurance"]
            / self.dict["Power_prop"]["E_fuel_Wh/kg"]
            / self.dict["Power_prop"]["eta_generator"]
        )

    def FixedEquipmentWeight_Total(self):
        self.dict["CL2Weight"]["Total Fixed Equipment Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.FlightControlSystem()["Average"]
            + self.HydraulicsPneumatics()["Average"]
            + self.InstrumentsAvionicsElectronics()["Average"]
            + self.ElectricalSystemWeight()["Average"]
            + self.AirconPressurizationAntiDeicingWeight()["Average"]
            + self.OxygenSystem()["Average"]
            + self.APU_Weight()["Average"]
            + self.Furnishings()["Average"]
            + self.AuxiliaryGear()["Average"]
            + self.Paint()["Average"]
        )
        self.dict["CL2Weight"]["Flight Control Systems Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.FlightControlSystem()["Average"]
        )
        self.dict["CL2Weight"]["Hydraulics and/or Pneumatics Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.HydraulicsPneumatics()["Average"]
        )
        self.dict["CL2Weight"]["Instruments, Avionics and Electronics Weight"] = self.dict[
            "Contingency_C2W"
        ] * lbs_to_N(self.InstrumentsAvionicsElectronics()["Average"])
        self.dict["CL2Weight"]["Electrical System Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.ElectricalSystemWeight()["Average"]
        )
        self.dict["CL2Weight"]["Airconditioning, Pressurization and Anti or Deicing Weight"] = self.dict[
            "Contingency_C2W"
        ] * lbs_to_N(self.AirconPressurizationAntiDeicingWeight()["Average"])
        self.dict["CL2Weight"]["Oxygen System Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.OxygenSystem()["Average"]
        )
        self.dict["CL2Weight"]["APU Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(self.APU_Weight()["Average"])
        self.dict["CL2Weight"]["Furnishings Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.Furnishings()["Average"]
        )
        self.dict["CL2Weight"]["Auxiliary Gear Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(
            self.AuxiliaryGear()["Average"]
        )
        self.dict["CL2Weight"]["Paint Weight"] = self.dict["Contingency_C2W"] * lbs_to_N(self.Paint()["Average"])

        return (
            self.FlightControlSystem()["Average"]
            + self.HydraulicsPneumatics()["Average"]
            + self.InstrumentsAvionicsElectronics()["Average"]
            + self.ElectricalSystemWeight()["Average"]
            + self.AirconPressurizationAntiDeicingWeight()["Average"]
            + self.OxygenSystem()["Average"]
            + self.APU_Weight()["Average"]
            + self.Furnishings()["Average"]
            + self.AuxiliaryGear()["Average"]
            + self.Paint()["Average"]
        )

    def NewEmptyWeight(self, bat_percent):
        return lbs_to_N(
            (
                self.PowerplantWeight_Total()["Average"]
                + self.StructureWeight_Total()
                + self.FixedEquipmentWeight_Total()
            )
        ) + 9.81 * self.NewBatteryWeight(bat_percent)

    def NewOEW(self, bat_percent):
        return self.NewEmptyWeight(bat_percent) + self.dict["Weights"]["W_Pilot_N"]

    # iteration function between class I and class II
    def Iterarions_C2W(self, bat_percent):
        # set up the old MTOW and the new one for the iteration loop
        MTOW_new = 0

        # initialise the class 2 weight subdictionary
        # self.dict["CL2Weight"] = {}
        self.dict["CL2Weight"]["MTOW_N"] = (
            self.dict["Weights"]["MTOW_N"] / self.dict["Contingency"]
        )  # the first value shall be the one returned from the class 1 weight estimation

        ok = False

        # iterate until the difference between the old and new MTOW is less than 2%
        while np.abs((MTOW_new - self.dict["CL2Weight"]["MTOW_N"]) / self.dict["CL2Weight"]["MTOW_N"]) > 0.02:
            # ok condition so that it doesnt update for the first step as it is needed to be saved later
            if ok:
                self.dict["CL2Weight"]["MTOW_N"] = MTOW_new

            # calculate the OEW, Battery Weight, and Fuel Weight
            OEW = self.NewOEW(bat_percent)
            BatteryWeight = self.NewBatteryWeight(bat_percent)
            FuelWeight = self.NewFuelWeight(bat_percent)

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
            ) = aerodynamic_design(ac_data=self.dict)
            self.dict["Aero"]["S_Wing"] = S_Wing
            self.dict["Aero"]["S_h"] = S_h
            self.dict["Aero"]["MAC_wing"] = mac_wing
            self.dict["Aero"]["MAC_HS"] = mac_HS
            self.dict["Aero"]["c_root_wing"] = c_root_wing
            self.dict["Aero"]["c_tip_wing"] = c_tip_wing
            self.dict["Aero"]["c_root_HS"] = c_root_HS
            self.dict["Aero"]["c_tip_HS"] = c_tip_HS
            self.dict["Aero"]["b_Wing"] = b_Wing
            self.dict["Aero"]["b_h"] = b_h

            # trasnform the new MTOW to N
            MTOW_new = (
                OEW
                + 9.81 * FuelWeight
                + 9.81 * self.dict["Iterations Class I"]["Wpl_des_kg"]
                - self.dict["Weights"]["W_Pilot_N"]
            )
            self.W_TO = N_to_lbs(MTOW_new)
            self.W_E = N_to_lbs(OEW)
            self.W_F = N_to_lbs(9.81 * FuelWeight)
            self.W_bat = N_to_lbs(9.81 * BatteryWeight)
            self.W_L = N_to_lbs(MTOW_new - 9.81 * FuelWeight)

            if MTOW_new > 80000:
                break

            ok = True

        # if MTOW is less than 60000, return the values for each group weight
        if MTOW_new < 80000:
            # self.W_E = self.W_E - N_to_lbs(9.81 * BatteryWeight)
            self.dict["Iterations Class I"]["MTOW_kg"] = 2650
            self.dict["CL2Weight"]["MTOW_N"] = self.dict["Contingency_C2W"] * MTOW_new
            self.dict["CL2Weight"]["OEW"] = self.dict["Contingency_C2W"] * OEW
            self.dict["CL2Weight"]["Wbat_N"] = self.dict["Contingency_C2W"] * 9.81 * BatteryWeight
            self.dict["CL2Weight"]["Wfuel_N"] = self.dict["Contingency_C2W"] * 9.81 * FuelWeight
            self.dict["CL2Weight"]["Wpl_w/o_pilot"] = self.dict["Contingency_C2W"] * (
                9.81 * self.dict["Iterations Class I"]["Wpl_des_kg"] - self.dict["Weights"]["W_Pilot_N"]
            )
            self.dict["CL2Weight"]["W_pilot"] = self.dict["Contingency_C2W"] * self.dict["Weights"]["W_Pilot_N"]
            # print MTOW w/o cont, MTOW w cont, OEW w cont, Bat weight w cont, Fuel weight w cont,
            # Payload w contingency, Structures w contingency, Fuel system w contingency, Powerplant w contingency,
            # Fixed equipment w contingency
            return (MTOW_new, self.dict["CL2Weight"]["MTOW_N"])
        else:
            # if MTOW is greater than 60000, the program shall return invalid values
            self.dict["Iterations Class I"]["MTOW_kg"] = 2650
            # print MTOW w/o cont, MTOW w cont, OEW w cont, Bat weight w cont,
            # Fuel weight w cont, Payload w contingency, Structures w contingency, Fuel system w contingency,
            # Powerplant w contingency, Fixed equipment w contingency
            return (0, self.dict["CL2Weight"]["MTOW_N"])

    def PolynomialRegression(self, bat):
        lst_P = []
        lst_bat = []

        for pbat in bat:
            row = self.Iterarions_C2W(pbat)
            if row[0] != 0:
                lst_P.append(row[1] / self.dict["Performance"]["W/P_N/W"])
                lst_bat.append(pbat)

        lst_bat = np.array(lst_bat)
        lst_P = np.array(lst_P)

        # Filter out non-positive values in lst_P
        valid_indices = lst_P > 0
        lst_bat = lst_bat[valid_indices]
        lst_P = lst_P[valid_indices]

        if len(lst_P) == 0:
            return np.array([20, 20]), np.array([20, 20])

        coeff_exp = np.polyfit(lst_bat, np.log(lst_P), 1)
        coeff_pol = np.polyfit(lst_bat, lst_P, 2)

        # y_pol = coeff_pol[0] * lst_bat**2 + coeff_pol[1] * lst_bat + coeff_pol[2]
        # y_exp = np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * lst_bat)

        return coeff_exp, coeff_pol

    def SensitivityAnalysis(self, components, percentage_changes):  # pragma: no cover
        """
        Perform sensitivity analysis on the specified component weights.

        Parameters:
        components (list of str): List of component weights to vary (e.g., "Wing", "Empennage", "Fuselage", etc.).
        percentage_changes (list of float): List of percentage changes to apply to the component weights.

        Returns:
        dict: Dictionary containing component names, percentage changes, and corresponding MTOW values.
        """
        original_weights = {
            "Wing": self.WingWeight()["Average"],
            "Empennage": self.EmpennageWeight()["Average"],
            "Fuselage": self.FuselageWeight()["Average"],
            "LandingGear": self.LandingGearWeight()["Average"],
            "Total Structures": self.StructureWeight_Total(),
            "Powerplant": self.PowerplantWeight_Total()["Average"],
            "FlightControl": self.FlightControlSystem()["Average"],
            "Hydraulics": self.HydraulicsPneumatics()["Average"],
            "Instruments": self.InstrumentsAvionicsElectronics()["Average"],
            "Electrical": self.ElectricalSystemWeight()["Average"],
            "Aircon": self.AirconPressurizationAntiDeicingWeight()["Average"],
            "Oxygen": self.OxygenSystem()["Average"],
            "Furnishings": self.Furnishings()["Average"],
            "Auxiliary Gear": self.AuxiliaryGear()["Average"],
            "Paint": self.Paint()["Average"],
            "Fixed Equipment": self.FixedEquipmentWeight_Total(),
            "Fuel": self.NewFuelWeight(bat=0.144),
            "Battery": self.NewBatteryWeight(bat=0.144),
        }

        results = {component: {} for component in components}

        for component in components:
            original_weight = original_weights[component]
            for pct_change in percentage_changes:
                # Apply percentage change
                new_weight = original_weight * (1 + pct_change / 100)

                # Update the corresponding method to return the new weight
                if component == "Wing":
                    self.WingWeight = lambda: {"Average": new_weight}
                elif component == "Empennage":
                    self.EmpennageWeight = lambda: {"Average": new_weight}
                elif component == "Fuselage":
                    self.FuselageWeight = lambda: {"Average": new_weight}
                elif component == "LandingGear":
                    self.LandingGearWeight = lambda: {"Average": new_weight}
                elif component == "Total Structures":
                    self.StructureWeight_Total = lambda: new_weight
                elif component == "Powerplant":
                    self.PowerplantWeight_Total = lambda: {"Average": new_weight}
                elif component == "FlightControl":
                    self.FlightControlSystem = lambda: {"Average": new_weight}
                elif component == "Hydraulics":
                    self.HydraulicsPneumatics = lambda: {"Average": new_weight}
                elif component == "Instruments":
                    self.InstrumentsAvionicsElectronics = lambda: {"Average": new_weight}
                elif component == "Electrical":
                    self.ElectricalSystemWeight = lambda: {"Average": new_weight}
                elif component == "Aircon":
                    self.AirconPressurizationAntiDeicingWeight = lambda: {"Average": new_weight}
                elif component == "Oxygen":
                    self.OxygenSystem = lambda: {"Average": new_weight}
                elif component == "Furnishings":
                    self.Furnishings = lambda: {"Average": new_weight}
                elif component == "Auxiliary Gear":
                    self.AuxiliaryGear = lambda: {"Average": new_weight}
                elif component == "Paint":
                    self.Paint = lambda: {"Average": new_weight}
                elif component == "Total Fixed Equipment":
                    self.FixedEquipmentWeight_Total = lambda: new_weight
                elif component == "Fuel":
                    self.NewFuelWeight = lambda bat_percent: new_weight
                elif component == "Battery":
                    self.NewBatteryWeight = lambda bat_percent: new_weight
                # Add other components as necessary

                # Recalculate MTOW
                MTOW_new = self.Iterarions_C2W(bat_percent=0.144)[0]

                # Store result
                results[component][pct_change] = MTOW_new

        # Restore original methods
        self.WingWeight = lambda: {"Average": original_weights["Wing"]}
        self.EmpennageWeight = lambda: {"Average": original_weights["Empennage"]}
        self.FuselageWeight = lambda: {"Average": original_weights["Fuselage"]}
        self.LandingGearWeight = lambda: {"Average": original_weights["LandingGear"]}
        self.TotalStructuresWeight = lambda: original_weights["Total Structures"]
        self.FuelSystemWeight = lambda: original_weights["Fuel System"]
        self.PowerplantWeight = lambda: original_weights["Powerplant"]
        self.FlightControlSystemWeight = lambda: original_weights["FlightControl System"]
        self.HydraulicsWeight = lambda: original_weights["Hydraulics"]
        self.InstrumentsWeight = lambda: original_weights["Instruments"]
        self.ElectricalSystemWeight = lambda: original_weights["Electrical System"]
        self.AirconWeight = lambda: original_weights["Aircon"]
        self.OxygenSystemWeight = lambda: original_weights["Oxygen System"]
        self.FurnishingsWeight = lambda: original_weights["Furnishings"]
        self.AuxiliaryGearWeight = lambda: original_weights["Auxiliary Gear"]
        self.PaintWeight = lambda: original_weights["Paint"]
        self.FixedEquipmentWeight = lambda: original_weights["Total Fixed Equipment"]
        self.NewFuelWeight = lambda: original_weights["Fuel"]
        self.NewBatteryWeight = lambda: original_weights["Battery"]
        # Add other components as necessary

        return results

    def plot_sensitivity_analysis(self, sensitivity_results):  # pragma: no cover
        """
        Plot the sensitivity analysis results.

        Parameters:
        sensitivity_results (dict): Results from the SensitivityAnalysis method.
        original_mtow (float): The original MTOW value.
        """
        components = list(sensitivity_results.keys())
        percentage_changes = [pct for pct in sensitivity_results[components[0]].keys() if pct != 0]  # Exclude 0% change

        fig, axes = plt.subplots(3, 2, figsize=(22, 15), sharey=True)
        # fig.tight_layout(pad=5.0)

        for i, pct_change in enumerate(percentage_changes):
            row = i // 2
            col = i % 2

            relative_errors = [
                (sensitivity_results[component][pct_change] - sensitivity_results[component][0])
                / sensitivity_results[component][0]
                * 100
                for component in components
            ]

            # Sort components and relative errors based on the relative errors
            sorted_indices = np.argsort(relative_errors)[::-1]
            sorted_components = [components[i] for i in sorted_indices]
            sorted_relative_errors = [relative_errors[i] for i in sorted_indices]

            axes[row, col].bar(sorted_components, sorted_relative_errors, color="skyblue")
            axes[row, col].set_title(f"Percentage Change: {pct_change}%", fontsize=20)  # Adjust font size as needed
            # axes[row, col].set_xlabel("Component", fontsize=20)  # Adjust font size as needed
            axes[row, col].set_ylabel("RelError MTOW (%)", fontsize=20)  # Adjust font size as needed
            axes[row, col].set_ylim([-20, 30])  # Adjust the limits as needed for better visualization
            axes[row, col].grid(True, axis="y", linestyle="--", alpha=0.7)
            axes[row, col].tick_params(axis="x", rotation=90)
            axes[row, col].tick_params(axis="both", labelsize=20)  # Adjust tick label font size

        plt.tight_layout()
        plt.savefig("Sensitivity_Analysis.pdf")
        plt.show()


def RunClassII(ac_data=aircraft_data, check=None, pbat=0.0):  # pragma: no cover
    # initialise the class II weight object
    p = Class_II_Weight(ac_data)
    p.Iterarions_C2W(pbat)

    dict = p.dict

    if check:
        print("========== Structures Weight ==========")
        print("\nWing Weight = ", round(dict["CL2Weight"]["Wing Weight"] / 9.81, 2), " [kg]")
        print("Empennage Weight =", round(dict["CL2Weight"]["Empennage Weight"] / 9.81, 2), " [kg]")
        print("Fuselage Weight = ", round(dict["CL2Weight"]["Fuselage Weight"] / 9.81, 2), " [kg]")
        print("Nacelle Weight = ", round(dict["CL2Weight"]["Nacelle Weight"] / 9.81, 2), " [kg]")
        print("Landing Gear Weight = ", round(dict["CL2Weight"]["Landing Gear Weight"] / 9.81, 2), " [kg]")
        print("\nTotal Structures Weight = ", round(dict["CL2Weight"]["Total Structures Weight"] / 9.81, 2), " [kg]")
        print("\n\n ========== Powerplant Weight ==========")
        # print('\n Fuel System Weight = ', lbs_to_N(p.FuelSystemWeight()["Average"]), " [N]")
        print("\n Total Powerplant Weight = ", round(dict["CL2Weight"]["Total Powerplant Weight"] / 9.81, 2), " [kg]")
        print("\n\n ========== Fixed Equipment Weight ==========")
        print(
            "\n Flight Control Systems Weight = ",
            round(dict["CL2Weight"]["Flight Control Systems Weight"] / 9.81, 2),
            " [kg]",
        )
        print(
            "Hydraulics and/or Pneumatics Weight = ",
            round(dict["CL2Weight"]["Hydraulics and/or Pneumatics Weight"] / 9.81, 2),
            " [kg]",
        )
        print(
            "Instruments, Avionics and Electronics Weight = ",
            round(dict["CL2Weight"]["Instruments, Avionics and Electronics Weight"] / 9.81, 2),
            " [kg]",
        )
        print("Electrical System Weight = ", round(dict["CL2Weight"]["Electrical System Weight"] / 9.81, 2), " [kg]")
        print(
            "Airconditioning, Pressurization and Anti or Deicing Weight = ",
            round(dict["CL2Weight"]["Airconditioning, Pressurization and Anti or Deicing Weight"] / 9.81, 2),
            " [kg]",
        )
        print("Oxygen System Weight = ", round(dict["CL2Weight"]["Oxygen System Weight"] / 9.81, 2), " [kg]")
        print("APU Weight = ", round(dict["CL2Weight"]["APU Weight"] / 9.81, 2), " [kg]")
        print("Furnishings Weight = ", round(dict["CL2Weight"]["Furnishings Weight"] / 9.81, 2), " [kg]")
        print("Auxiliary Gear Weight = ", round(dict["CL2Weight"]["Auxiliary Gear Weight"] / 9.81, 2), " [kg]")
        print("Paint Weight = ", round(dict["CL2Weight"]["Paint Weight"] / 9.81, 2), " [kg]")
        print(
            "\nTotal Fixed Equipment Weight = ",
            round(dict["CL2Weight"]["Total Fixed Equipment Weight"] / 9.81, 2),
            " [kg]",
        )
        print("\n\n ========== Total Operating Weights ==========")
        # print("\nTotal Empty Weight = ", round(p.NewEmptyWeight(pbat) / 9.81, 2), " [kg]")
        print("Total Operating Empty Weight = ", round(dict["CL2Weight"]["OEW"] / 9.81, 2), " [kg]")

    return dict


if __name__ == "__main__":  # pragma: no cover
    # Populate with actual aircraft data
    class_ii_weight = Class_II_Weight(ac_data=aircraft_data)
    percentage_changes = [-10, -5, 0, 5, 10, 15, 20]  # Example percentage changes
    components = [
        "Wing",
        "Empennage",
        "Fuselage",
        "LandingGear",
        "Total Structures",
        "Powerplant",
        "FlightControl",
        "Hydraulics",
        "Instruments",
        "Electrical",
        "Aircon",
        "Oxygen",
        "Furnishings",
        "Auxiliary Gear",
        "Paint",
        "Fixed Equipment",
        "Fuel",
        "Battery",
    ]
    sensitivity_results = class_ii_weight.SensitivityAnalysis(components, percentage_changes)

    # print(sensitivity_results)
    # for pct_change, mtow in sensitivity_results.items():
    #     print(mtow)
    #     print(f"Percentage Change: {pct_change}%, MTOW: {mtow} N")

    class_ii_weight.plot_sensitivity_analysis(sensitivity_results)

    RunClassII(ac_data=aircraft_data, check=True, pbat=0)
