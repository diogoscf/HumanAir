import numpy as np
from scipy.optimize import root_scalar, minimize_scalar  # type: ignore[import-untyped]

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir import aircraft_data

from HumanAir.FlightPerformance.helper import density
from HumanAir.FlightPerformance import thrust_power, takeoff_landing


class Aircraft:
    """
    Aircraft object for evaluating performance.

    Since the electric engine is more powerful than the fuel engine it assumes
    the fuel engine is used.

    IMPORTANT: the thrust function stores intermediate values for faster
    calculations. If aircraft parameters change it should be ensured that
    these intermediate calculations are ereased.
    """

    def __init__(self, FILE="design.json"):
        dat = aircraft_data.aircraft_data

        #
        # performance
        #

        self.V_cruise = dat["Performance"]["Vc_m/s"]
        self.h_cruise = dat["Performance"]["Altitude_Cruise_m"]
        self.h_cruise_max = dat["Performance"]["Altitude_Max_Cruise_m"]
        self.h_TO = dat["Performance"]["Altitude_TO_m"]
        self.h_land = dat["Performance"]["Altitude_Land_m"]
        self.dT_default = dat["Performance"]["Temp_offset_TO_Land_cruise"]

        #
        # aerodynamic data
        #

        self.CD0_clean = dat["Aero"]["CD0"]
        self.AR = dat["Aero"]["AR"]
        self.e_clean = dat["Aero"]["e"]
        self.S_clean = dat["Aero"]["S_Wing"]
        self.b = dat["Aero"]["b_Wing"]
        self.CL_Dmin = 0.0  # 0.146#0.17 voor matchen point met Cdmin=Cd0 # 0.16 voor L/D van 19 #0.866

        self.CLmax_clean = dat["Aero"]["CLmax_clean"]
        self.CLmax_TO = dat["Aero"]["CLmax_TO"]
        self.CLmax_land = dat["Aero"]["CLmax_Land"]

        self.CL_ground_TO = dat["Flaps"]["CL_AoA0_landing"]
        self.CL_ground_land = dat["Flaps"]["CL_AoA0_takeoff"]

        #
        # flaps
        #

        # if any of these statements are not true flap drag calculation has to
        # be revised using roskam pt 6
        if (
            dat["Aero"]["QuarterChordSweep_Wing_deg"] != 0
            or not np.isclose(dat["Flaps"]["cf_c"], 0.25)
            or dat["Flaps"]["deflection_takeoff"] != 45
            or dat["Flaps"]["deflection_landing"] != 45
            or dat["Flaps"]["flap_start"] > 0.3
            or dat["Flaps"]["flap_start"] < 0.1
            or dat["Flaps"]["flap_end"] > 0.75
            or dat["Flaps"]["flap_end"] < 0.65
            or not np.isclose(self.AR, 11)
        ):
            raise Exception("Flap drag calculations need to be revised")

        self.S_wf = dat["Flaps"]["Swf"]
        self.Sprime_S_ld = dat["Flaps"]["Sprime_S_landing"]
        self.Sprime_S_TO = dat["Flaps"]["Sprime_S_takeoff"]
        self.flap_defl_TO = dat["Flaps"]["deflection_takeoff"]
        self.flap_defl_ld = dat["Flaps"]["deflection_landing"]

        #
        # weights
        #

        self.W_OE = dat["CL2Weight"]["OEW"]
        self.W_MF = dat["CL2Weight"]["Wfuel_N"]  # max fuel weight
        self.W_pl_no_pilot = dat["CL2Weight"]["Wpl_w/o_pilot"]  # class ii
        # self.W_pl_des       = dat["Weights"]["Wpl_des_kg"] * 9.80655 # class i
        # self.W_pl_max       = dat["Weights"]["Wpl_max_kg"] * 9.80655 # class i
        self.W_MTO = dat["CL2Weight"]["MTOW_N"]

        #
        # propulsion
        #

        self.max_cont_power_sealevel = dat["Power_prop"]["Fuel_engine_P_max_cont_W"]
        self.takeoff_power_sealevel = dat["Power_prop"]["Fuel_engine_P_TO_W"]
        self.eff_powertrain = dat["Power_prop"]["eta_powertrain"]
        self.max_eff_prop = dat["Power_prop"]["eta_p"]
        self.fuel_cons_TO = dat["Power_prop"]["Newton_fuel_per_W_per_S_TO"]
        self.fuel_cons_flight = dat["Power_prop"]["Newton_fuel_per_W_per_S_flight"]

        self.electric_takeoff_power = dat["Power_prop"]["Electric_engine_P_TO_W"]
        self.electric_max_cont_power = dat["Power_prop"]["Electric_engine_P_max_cont_W"]
        self.eff_electric_motor = dat["Power_prop"]["eta_electricmotor"]
        self.eff_battery = dat["Power_prop"]["eta_bat"]
        self.max_bat_cap = dat["Power_prop"]["E_bat_Wh"]
        self.max_DoD_bat = dat["Power_prop"]["DoD_bat"]

        self.number_of_engines = 1
        self.propeller_diameter = dat["Power_prop"]["Dp_m"]
        self.spinner_diameter = 0.4

        #
        # misc
        #

        self.wing_h_above_ground = (
            dat["Geometry"]["fus_height_m"] + dat["Landing_gear"]["Hs_m"] + dat["Landing_gear"]["Dwm_m"] / 2
        )

        self.retractable_gear = dat["Landing_gear"]["Retractable"]

        x_cg = dat["Stability"]["Cg_Front"] * dat["Aero"]["MAC_wing"] + dat["Geometry"]["XLEMAC_m"]
        x_main, x_nose = dat["Landing_gear"]["Xmw_m"], dat["Landing_gear"]["Xnw_m"]
        self.weight_on_MLG = (x_cg - x_nose) / (
            x_main - x_nose
        )  # fraction of weight on MLG with front most possible CG position

        #
        # calculated data - normal equations cannot be used since CL_CDmin =/= 0
        #

        # CDmin to get CD0 at CL=0
        self.CDmin = self.CD0_clean - self.CL_Dmin**2 / (np.pi * self.AR * self.e_clean)
        print(f"CDmin: {self.CDmin:.4f}")

        # find CL/CD max for clean config
        def CL_CD(CL):
            return -(CL / self.CD(CL))  # minus so that minimize_scalar() will maximise

        res = minimize_scalar(CL_CD, bounds=[0, self.CLmax_clean])
        self.LDmax = -res.fun
        self.CL_LDmax = res.x
        self.CD_LDmax = self.CD(res.x)

        # find CL^3/CD^2 max for clean config
        def CL3_CD2(CL):
            return -(CL**3 / self.CD(CL) ** 2)  # minus so that minimize_scalar() will maximise

        res = minimize_scalar(CL3_CD2, bounds=[0, self.CLmax_clean])
        self.L3D2max = -res.fun
        self.CL_L3D2max = res.x
        self.CD_L3D2max = self.CD(res.x)

        # self.LDmax    = 0.5 * np.sqrt(np.pi*self.AR*self.e_clean/self.CD0_clean) # ruijgrok p107
        # self.CD_LDmax = 2 * self.CD0_clean # ruijgrok p106
        # self.CL_LDmax = self.CD_LDmax * self.LDmax

        # self.L3D2max    = 3 * np.sqrt(3) / 16 * np.pi * self.AR * self.e_clean \
        #     * np.sqrt(np.pi * self.AR * self.e_clean / self.CD0_clean) # for P_r_min for max climbrate; ruijgrok p107
        # self.CD_L3D2max = 4 * self.CD0_clean # ruijgrok p107
        # self.CL_L3D2max = (self.L3D2max * self.CD_L3D2max**2)**(1/3)

    def S(self, flaps="up"):
        if flaps == "up":
            return self.S_clean
        elif flaps == "TO":
            return self.S_clean * self.Sprime_S_TO
        elif flaps == "land":
            return self.S_clean * self.Sprime_S_ld

    def CD(
        self, CL, gear="up", flaps="up", ground_effect_h=None
    ):  # this parabolic approximation is about accurate until CL=1 (ruijgrok p226)
        """Parabolic estimation of C_D for given C_L for given gear/flaps conditions (default=clean)"""

        CDmin = self.CDmin
        e = self.e_clean

        # ground effect reduces induced drag, using nicolai section 10.2
        # (credit to Wieselberger for the equation), according to gudmundson
        # a good approximation for 0.033 < h_b < 0.25
        sigma = 0
        if ground_effect_h is not None:
            h_b = (self.wing_h_above_ground + ground_effect_h) / self.b
            sigma = (1 - 1.32 * h_b) / (1.05 + 7.4 * h_b)
            sigma = max(sigma, 0)  # for around h_b > .75 sigma will become negative

        if gear == "down" and self.retractable_gear:
            CDmin += 0.011  # p264 of Nicolai: Fundamentals of Aircraft and Airship Design: Volume 1

        if flaps == "TO" or flaps == "land":
            # roskam pt 6 p82 and further
            delta_CL = (
                self.CLmax_TO - self.CLmax_clean
            )  # assume constant delta CL (not according to Roskam but close enough)
            delta_CDpro = 0.065 * self.S_wf / self.S_clean  # induced drag
            delta_CDind = 0.2**2 * delta_CL**2  # profile drag
            delta_CDint = 0.4 * delta_CDpro  # inteference drag

            return (
                CDmin
                + delta_CDint
                + delta_CDpro
                + (delta_CDind + (CL - self.CL_Dmin - delta_CL) ** 2 / (np.pi * self.AR * e)) * (1 - sigma)
            )
        else:
            return CDmin + (CL - self.CL_Dmin) ** 2 / (np.pi * self.AR * e) * (1 - sigma)

    def P_shaft(self, h, dT, use_takeoff_power=False, electric=False):
        return thrust_power.P_shaft(self, h, dT, use_takeoff_power=use_takeoff_power)

    def P_a(self, h, dT, use_takeoff_power=False, V=None):
        return thrust_power.P_a(self, h, dT, use_takeoff_power=use_takeoff_power, V=V)

    def T(self, V_ms, h, dT, use_takeoff_power=False, electric=False):
        return thrust_power.T(self, V_ms, h, dT, use_takeoff_power=use_takeoff_power)

    def prop_eff(self, V, h, dT, use_takeoff_power=False):
        return thrust_power.prop_eff(self, V, h, dT, use_takeoff_power=use_takeoff_power)

    def fuel_rate(self, P_shaft, is_takeoff):
        return thrust_power.fuel_rate(self, P_shaft, is_takeoff)

    def bat_cap_rate(self, P_shaft):
        return thrust_power.bat_cap_rate(self, P_shaft)

    def P_r(self, W, V, h, dT):
        rho = density(h, dT)
        CL = W / (0.5 * rho * V**2 * self.S())
        return V * 0.5 * rho * V**2 * self.S() * self.CD(CL)

    def V_Dmin(self, W, h, dT):
        """Calculates velocity corresponding with minimum drag in clean configuration"""
        a = 1 / np.sqrt(self.CD0_clean * np.pi * self.AR * self.e_clean)
        b = np.sqrt(W / self.S() * 2 / density(h, dT) * a)  # ruijgrok p224, same as V at (L/D)max
        c = np.sqrt(W / (0.5 * density(h, dT) * self.S * self.CL_LDmax))
        if not np.isclose(b, c):
            print(f"V_D min error: {b:.2f} - {c:.2f}")
        return b

    def D_min(self, W):
        """Minimum drag (in Newtons) for given weight in clean configuration. Corresponds to V_Dmin"""
        return W / self.LDmax  # ruijgrok p220

    def P_r_min(self, W, h, dT):
        """Minimum possible value for  power required for given conditions in clean configuration,
        note that this assumes a small flight path angle."""
        return W * np.sqrt(W / self.S() * 2 / density(h, dT) * 1 / self.L3D2max)  # ruijgrok p221

    def V_Prmin(self, W, h, dT):
        return np.sqrt(W / (0.5 * density(h, dT) * self.S() * self.CL_L3D2max))  # assumes small flight path angle

    def RC_max(self, W, h, dT):
        """
        Calculates maximum possible climb rate for clean configuration in given
        condition. Assumes constant V.

        Uses small angle approximation for L = W.

        Note that for high altitudes the value is not accurate due to bad
        high-CL drag estimation.

        Parameters
        ----------
        W : float
            Aircraft gross weight [N].
        h : float
            Geopotential altitude [m].
        dT : TYPE
            ISA temperature offset [deg C].

        Returns
        -------
        float
            Maximum climb rate or 0 if CL>1 [m/s].
        float
            Speed at which max climb rate is reached or 0 if CL>1 [m/s]

        """
        rho = density(h, dT)

        min_speed = self.stall_speed(W, h, dT)
        max_speed = self.V_max(W, h, dT)
        # V_opt = 0
        # CL_opt = 0
        # RC_max = 0

        # # rc_list = []
        # def RC(V): # calculates climbrate *-1
        #     CL = W / (0.5 * rho * V**2 * self.S()) # small angle assumption, L=W cos(gamma) -> L=W
        #     CD = self.CD(CL)
        #     P_r = V * (0.5 * rho * V**2 * self.S() * CD)
        #     P_a = self.P_a(h, dT, V=V)
        #     return (P_a - P_r) / W, CL

        # # this iterates over V speeds in reverse order, as soon as CL becomes
        # # greater than 1 it stops and returns zero. This way this function
        # # return RC=0 V_opt=0
        # for V in np.arange(min_speed, max_speed+1, 0.5)[::-1]:
        #     rc, cl = RC(V)
        #     # rc_list.append(rc)
        #     if rc > RC_max:
        #         V_opt = V
        #         CL_opt = cl
        #         RC_max = rc
        # print(f"RC CL = {CL_opt:.3f}")

        # code for plotting RC-V graph
        # plt.figure(figsize=(10,7))
        # plt.plot(np.arange(min_speed, max_speed+1, 1), rc_list)
        # plt.xlabel("Velocity [m/s]")
        # plt.ylabel("Rate of climb")
        # plt.show()

        def RC(V):  # calculates climbrate *-1
            CL = W / (0.5 * rho * V**2 * self.S())  # small angle assumption, L=W cos(gamma) -> L=W
            CD = self.CD(CL)
            P_r = V * (0.5 * rho * V**2 * self.S() * CD)
            P_a = self.P_a(h, dT, V=V)
            return -(P_a - P_r) / W  # minus so that minimize() will maximize

        bounds = [min_speed, max_speed]
        if bounds[0] > bounds[1]:
            return 0

        res = minimize_scalar(RC, bounds=bounds, method="bounded")  # minus to convert back to positive
        RC_max = -res.fun
        V_opt = res.x

        #
        # alternate method, does not take into account varying prop efficiency,
        # result is almost the same
        #
        # rc2 = (self.P_a(h, dT) - self.P_r_min(W, h, dT)) / W
        # print(f"RC {res.x:.2f} - {rc_max:.2f} - {rc2:.2f}")

        return np.array([RC_max, V_opt])

    def climb_slope_max(self, W, h, dT, gear="up", flaps="up"):
        """
        Calculates climb slope for given conditions, assuming constant V. The
        drag is multiplied by 1.5 to account for high drag at high CL, however
        this is not accurate at all.

        Note that CS23 requires the speed at which this is measured to be
        >= 1.2 * the stall speed in this config.

        It uses the small angle approximation L = W.

        Takes into account lower power available at low speeds, and since the
        max climb slope is achieved at high CL (=low speeds) this method gives
        signifiantly lower results than when constant P_a is assumed.

        Parameters
        ----------
        W : float
            Aircraft gross weight.
        h : float
            Altitude.
        dT : float
            ISA temperature offset.
        gear : string, optional
            "up" or "down". The default is "up".
        flaps : float, optional
            "up", "TO", or "landing. The default is "up".

        Returns
        -------
        float
            Maximum climb slope ranging 0-100%.

        """
        rho = density(h, dT)
        S = self.S(flaps=flaps)

        min_speed = 1.2 * self.stall_speed(W, h, dT, flaps=flaps)  # according to CS23
        # min_speed = np.sqrt(W / (0.5 * rho * S * 1)) # CLmax = 1, beyond that no accurate results
        max_speed = self.V_max(W, h, dT)

        V_opt = 0
        CL_opt = 0
        slope_max = 0

        def slope(V):  # calculates slope *-1
            CL = W / (0.5 * rho * V**2 * S)  # small angle assumption, L=W cos(gamma) -> L=W
            CD = 1.5 * self.CD(
                CL, gear=gear, flaps=flaps
            )  # inaccurate with current CD calculations, currently multiplied by 1.5 as adjustment
            D = 0.5 * rho * V**2 * S * CD
            T = self.T(V, h, dT)
            return (T - D) / W, CL

        for V in np.arange(min_speed, max_speed + 1, 1):
            sl, cl = slope(V)
            if sl > slope_max:
                V_opt = V
                CL_opt = cl
                slope_max = sl * 100

        print(f"Slope CL = {CL_opt:.2f}; at speed {V_opt:.2f}; min speed {min_speed:.2f}")
        # def slope(V): # calculates slope *-1
        #     CL = min(1, W / (0.5 * rho * V**2 * self.S)) # small angle assumption, L=W cos(gamma) -> L=W
        #     CD = self.CD(CL, gear=gear, flaps=flaps) # inaccurate with current CD calculations
        #     D = 0.5 * rho * V**2 * self.S * CD
        #     T = self.T(V, h, dT)
        #     return -(T - D) / W # minus so that minimize() will maximize

        # min_speed = 1.2 * self.stall_speed(W, h, dT, flaps=flaps) # according to CS23
        # max_speed = self.V_max(W, h, dT)

        # if min_speed > max_speed:
        #     return 0 # this happends when h > ceiling, just return 0

        # minus to convert back to positive
        # res = minimize_scalar(slope, bounds=[min_speed, max_speed], method="bounded")
        # slope_max = np.tan(np.arcsin(-res.fun)) * 100 # result is sin(gamma), we need tan(gamma)

        # if res.x < min_speed:
        #     raise Exception(f"Speed for climb slope = {res.x:.2f}, which less than minimum {min_speed:.2f}.")

        #
        # alternate method for testing (does not take into account varying prop eff and thus returns too high max slope)
        #

        # if flaps == "up":
        #     CL = self.CLmax_clean / self.CL_climb_safetyfactor
        # elif flaps == "TO":
        #     CL = self.CLmax_TO / self.CL_climb_safetyfactor
        # elif flaps == "land":
        #     CL = self.CLmax_land / self.CL_climb_safetyfactor
        # slope_max2 = (
        #     100 * self.P_a(h, dT) / W * 1 / np.sqrt(W / self.S * 2 / density(h, dT) * 1 / CL)
        #     - self.CD(CL, gear=gear, flaps=flaps) / CL
        # )
        # print(f"SL {res.x:.2f} - {slope_max:.2f} - {slope_max2:.2f}")

        return np.array([slope_max, V_opt])

    def stall_speed(self, W, h, dT, flaps="up"):
        """
        Calcualtes stall speed for given conditions in given config (default
        is clean).

        Parameters
        ----------
        W : float
            Aircraft gross weight [N].
        h : float
            Geopotential altitude [m].
        dT : float
            ISA temperature offset [deg C].
        flaps : string, optional
            "up", "TO", or "landing". The default is "up".

        Returns
        -------
        float
            Stall speed [m/s].

        """
        if flaps == "up":
            CL = self.CLmax_clean
        elif flaps == "TO":
            CL = self.CLmax_TO
        elif flaps == "land":
            CL = self.CLmax_land

        return np.sqrt(W / (0.5 * density(h, dT) * self.S(flaps=flaps) * CL))

    def V_max(self, W, h, dT, use_max_prop_eff=False):
        """
        Calculates the maximum true airspeed in clean configuration.

        It does so by finding the point where P_a() equals the drag force times
        velocity.

        When use_max_prop_eff is set to true it will not take into account
        low-speed losses. This is to prevent an infinite loop when used by
        the thrust function.

        Note that it is inaccurate for high altitudes because of bad high-CL
        drag prediction. Also the root-scalar method is unpredictable for high
        altitudes.


        Parameters
        ----------
        W : float
            Aircraft weight [N].
        h : float
            Geopotential altitude [m].
        dT : float
            ISA temperature offset [deg C].
        use_max_prop_efficiency : boolean
            Whether to assume max prop efficiency

        Returns
        -------
        float
            The maximum true airspeed in m/s.
        """
        rho = density(h, dT)

        def diff(V):  # returns difference between power available and power required
            if type(V) is not np.float64:
                V = V[0]
            CL = W / (0.5 * rho * self.S() * V**2)
            CD = self.CD(CL)
            D = 0.5 * rho * self.S() * CD * V**2
            if use_max_prop_eff:
                return self.P_a(h, dT) - D * V
            else:
                return self.P_a(h, dT, V=V) - D * V

        return root_scalar(diff, x0=150).root

    def power_setting(self, V, W, h, dT):
        """
        Returns power setting (0-1 and >1 if impossible speed) to maintain
        steady horizontal flight at given speed in clean configuration.

        Returns 0 if C_L > 1 because with current drag estimation CD is
        inaccurate for C_L>1.

        Parameters
        ----------
        V : float
            True airpseed [m/s].
        W : float
            Aircraft gross weight [N].
        h : float
            Altitude [m].
        dT : float
            ISA temperature offset [deg C].

        Returns
        -------
        float
            Power setting [0-1 and >1 if impossible].

        """
        CL = W / (0.5 * density(h, dT) * V**2 * self.S())
        D = 0.5 * density(h, dT) * V**2 * self.S() * self.CD(CL)

        if CL > 1:
            print(h)
            return 0

        return V * D / self.P_a(h, dT, V=V)

    #
    # Take-off and landing
    #

    def takeoff_ground_run(self, W, h, dT, slope, surface, electric=False, calc_time=False):
        return takeoff_landing.takeoff_ground_run(
            self, W, h, dT, slope, surface, electric=electric, calc_time=calc_time
        )

    def landing_ground_distance(self, W, h, dT, slope, surface, reversible_pitch=False, calc_time=False):
        return takeoff_landing.landing_ground_distance(
            self, W, h, dT, slope, surface, reversible_pitch=reversible_pitch, calc_time=calc_time
        )
