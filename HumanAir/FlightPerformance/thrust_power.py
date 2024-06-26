import numpy as np
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import HumanAir.unit_conversions as conv
from HumanAir.FlightPerformance.helper import density


def P_shaft(acf, h, dT, use_takeoff_power=False, electric=False):
    """
    Returns the maximum available shaft power (therefore excluding propeller
    efficiency) for either continuous operation or short operation.

    It is assumed that shaft power decreases with air density but is constant
    with velocity. This is a good assumption according to Ruijgrok p132-p133
    and Gudmundson. The decrease with altitude is computed using the Gagg and
    Ferrar model (sourced from Gudmundson).

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    h : float
        The geopotential altitude [m].
    dT : float
        ISA temperature offset [deg C].
    use_takeoff_power : boolean, optional
        Whether to use takeoff power. The default is False, which will result
        in the maximum continuous power to be returned.

    Returns
    -------
    float
        The maximum available shaft power [W].

    """

    if electric:
        if use_takeoff_power:
            return acf.electric_takeoff_power * acf.eff_powertrain
        else:
            return acf.electric_max_cont_power * acf.eff_powertrain
    else:
        # temperature correction is included here (temperature only influences
        # P_shaft by changing the density - see gudmundson example 7-3)
        alt_correction = 1.132 * (density(h, dT) / 1.225) - 0.132  # Gagg and Ferrar model, Gudmundsen eq 7-16

        if use_takeoff_power:
            return acf.takeoff_power_sealevel * acf.eff_powertrain * alt_correction
        else:
            return acf.max_cont_power_sealevel * acf.eff_powertrain * alt_correction


def P_a(acf, h, dT, use_takeoff_power=False, V=None):
    """
    Max continuous power available for given conditions. If no airspeed is
    given then the maximum propeller efficiency is used instead of the speed
    dependent one.

    Constant propeller efficiency at max value is a good assumption for high
    velocities according to Ruijgrok p176.

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    h : float
        The geopotential altitude [m].
    dT : float
        The ISA temperature offset.
    use_takeoff_power : boolean, optional
        Wether to use takeoff power instead of max continuous power. The
        default is False.
    V : float, optional
        Airspeed. The default is None, which results in max prop efficiency
        being used.

    Returns
    -------
    float
        Power available [W].

    """

    if V is None:
        return acf.max_eff_prop * P_shaft(acf, h, dT, use_takeoff_power=use_takeoff_power)
    else:
        return acf.prop_eff(V, h, dT) * P_shaft(acf, h, dT, use_takeoff_power=use_takeoff_power)


"""
Thrust calculations for a certain h and dT and use_takeoff_power are stored
here so that the matrix does not need to be recalculated for every call to the
thrust function.

It is an array of dicts, each dict containing:

    "h": float, altitude

    "dT": float, ISA temperature offset

    "use_takeoff_power": boolean, whether max constant or max TO power was used

    "electric": boolean, whether the electric engine was used

    V_C_ms : float, maximum velocity this calculation can be used for in m/s,
        beyond this point T = P_a/V is used.

    "abcd_vector": 1D array with four floats, this is multiplied (dot product)
        with the vector [V**3, V**2, V, 1] (with V in kt) to get the thrust in
        lbf.
"""
thrust_calculations = []  # type: ignore[var-annotated]


def T(acf, V_ms, h, dT, use_takeoff_power=False, electric=False):
    """
    Approximates the thrust for given conditions. It uses the cubic spline
    interpolation as derived in Gudmundsun section 14.4.2 (method #3).

    https://www.sciencedirect.com/science/article/pii/B9780123973085000143#s0235

    Note that it stores intermediate calculations for faster repeated thrust
    calculations with the same h, dT and use_takeoff power. To force a 'reset'
    (for example if aircraft parameters change) thrust_calculations must be set
    to an empty array.

    The main purpose of this function is for calculating thrust at low
    airspeeds; for cruise speed and beyond one can simply use T=P/V*eff_prop.

    Basically it interpolates between the static thrust at zero airspeed,
    P/V*Max_prop_efficiency at cruise and P/V*Max_prop_efficiency at maximum
    airspeed. This is done since the efficiency for constant-speed propellers
    can be assumed to be constant for high air speeds. The limitation is that
    we do not know when the maximum efficiency is reached: that specific speed
    may be lower than the cruise speed, giving an underestimated thrust for
    V < V_cruise.

    This function will throw an exception for V>V_max as beyond the maximum
    speed the interpolation is no longer valid (see Gudmundson fig 14-42).

    For calculating V_max it uses the OEW.

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    V : float
        True airspeed [m/s].
    h : float
        Geopotential altitude [m].
    dT : float
        ISA temperature offset.
    use_takeoff_power : boolean, optional
        Whether to use takeoff power instead of maximum continuous power. The
        default is False.

    Returns
    -------
    thrust_N : float
        Thrust in N.

    """

    # retrieve previously calculated abcd_vector if it exists
    abcd_vector = None
    for calc in thrust_calculations:
        if (
            np.isclose(calc["h"], h)
            and np.isclose(calc["dT"], dT)
            and calc["use_takeoff_power"] == use_takeoff_power
            and calc["electric"] == electric
        ):
            abcd_vector = calc["abcd_vector"]
            V_C_ms = calc["V_C_ms"]

    # if no previously calculated abcd_vector exists, calculate it
    if abcd_vector is None:

        rho = density(h, dT)

        # max prop efficiency
        eff_prop_max = acf.max_eff_prop

        # get shaft power
        P_sh_W = P_shaft(acf, h, dT, use_takeoff_power, electric=electric)
        P_sh_hp = conv.W_to_hp(P_sh_W)

        # calculate static thrust
        A_prop = np.pi * (acf.propeller_diameter / 2) ** 2
        A_spinner = np.pi * (acf.spinner_diameter / 2) ** 2
        T_static_N = (
            0.85 * P_sh_W ** (2 / 3) * (2 * rho * A_prop) ** (1 / 3) * (1 - A_spinner / A_prop)
        )  # Gudmundson eq 14-64
        T_static_lb = conv.N_to_lbs(T_static_N)

        # cruise speed
        V_C_ms = acf.V_cruise
        V_C_kt = conv.m_s_to_kt(V_C_ms)

        # max speed
        V_H_ms = acf.V_max(acf.W_OE, h, dT, use_max_prop_eff=True)
        V_H_kt = conv.m_s_to_kt(V_H_ms)

        # cruise thrust
        T_cruise_N = P_sh_W * eff_prop_max / V_C_ms
        T_cruise_lb = conv.N_to_lbs(T_cruise_N)

        # max speed thrust
        T_Vmax_N = P_sh_W * eff_prop_max / V_H_ms
        T_Vmax_lb = conv.N_to_lbs(T_Vmax_N)

        # print(f"For h={h:.2f} dT={dT:.2f} use_takeoff_power={use_takeoff_power}:\n \
        #         Static thrust: {T_static_N:.2f}\n \
        #         Cruise thrust: {T_cruise_N:.2f}\n \
        #         V_max thrust:  {T_Vmax_N:.2f}")

        #
        # now performing the actual calculations according to Gudmundson eq 14-41 and eq 14-42
        #

        # matrix
        MTX = [
            [0, 0, 0, 1],
            [V_C_kt**3, V_C_kt**2, V_C_kt, 1],
            [3 * V_C_kt**2, 2 * V_C_kt, 1, 0],
            [V_H_kt**3, V_H_kt**2, V_H_kt**2, 1],
        ]

        # vecotr
        VCT = [T_static_lb, T_cruise_lb, -eff_prop_max * 325.8 * P_sh_hp / V_C_kt**2, T_Vmax_lb]

        # first multiplication
        abcd_vector = np.dot(np.linalg.inv(MTX), VCT)

        # store calculated abcd_vector
        thrust_calculations.append(
            {
                "h": h,
                "dT": dT,
                "use_takeoff_power": use_takeoff_power,
                "electric": electric,
                "V_C_ms": V_C_ms,
                "abcd_vector": abcd_vector,
            }
        )

    # check if V_ms is out of bounds
    if V_ms > V_C_ms:
        # this method is not accurate for V > V_cruise, it was found that if
        # V_max is too far from V_cruise the interpolation gets off, so beyond
        # V_cruise P_a/V is returned. This works since one of the constraints
        # of the interpolation is that dT/dV=0 at V_cruise
        return P_a(acf, h, dT, use_takeoff_power=use_takeoff_power) / V_ms

    # velocity vector
    V_kt = conv.m_s_to_kt(V_ms)
    VCT_2 = [V_kt**3, V_kt**2, V_kt, 1]

    # calculate thrust
    thrust_lbs = np.dot(abcd_vector, VCT_2)
    thrust_N = conv.lbs_to_N(thrust_lbs)

    return thrust_N


def prop_eff(acf, V, h, dT, use_takeoff_power=False):
    """
    Calculates prop efficiency using eff_p = T*V/P_shaft.

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    V : float
        True airspeed [m/s].
    h : float
        Geopoential altitude [m].
    dT : float
        ISA temperature offset [deg C].
    use_takeoff_power : boolean, optional
        Whether to use takoeff or max continuous power. The default is False.

    Returns
    -------
    float
        Prop efficiency [0-1].

    """
    eff = (
        T(acf, V, h, dT, use_takeoff_power=use_takeoff_power)
        * V
        / P_shaft(acf, h, dT, use_takeoff_power=use_takeoff_power)
    )
    if eff > 1 and V > 0:  # warn user, V>0 condition to ignore warnings for nonsensical values
        print(
            f"Propeller efficiency > 1. V={V} m/s, h={h} m, dT={dT} degC, takeoffpower={use_takeoff_power}, eff={eff}"
        )
    return eff


def fuel_rate(acf, P_shaft, is_takeoff):
    """
    Calculates the fuel consumed for given shaft power required in N/s.

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    P_shaft : float
        Shaft power required [W].
    takeoff : boolean
        Whether this is during take-off

    Returns
    -------
    float
        Fuel consumption rate [N/s].

    """
    if P_shaft < 0:
        raise Exception("P_shaft supplied is negative")

    # from engine specs:
    # T/O fuel consumption 0.17 kg/hp/hr -> 0,00000062101634094798 N/W/s
    # averaged non-TO fuel consumption 0.16 kg/hp -> 0,00000058448596795104 N/W/s
    if is_takeoff:
        return P_shaft / acf.eff_powertrain * acf.fuel_cons_TO
    else:
        return P_shaft / acf.eff_powertrain * acf.fuel_cons_flight


def bat_cap_rate(acf, P_shaft):
    """
    Calculates the rate at which battery capacity is depleted in Wh/s.

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    P_shaft : float
        Shaft power required [W].

    Returns
    -------
    float
        Rate at which battery capacity is consumed [Wh/s].

    """
    if P_shaft < 0:
        raise Exception("P_shaft supplied is negative")

    return P_shaft / acf.eff_powertrain / acf.eff_electric_motor / acf.eff_battery / 3600  # conv to Wh/s
