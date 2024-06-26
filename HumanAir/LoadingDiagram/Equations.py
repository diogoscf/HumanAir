import numpy as np

"""========== General Equations =========="""


def Density(h, temp_offset, rho0=1.225, constlambda=-0.0065, g0=9.80665, T0=288.15, R=287.05):
    Temp = T0 + constlambda * h  # temperature correction is for temperature offset
    Correction = Temp / (Temp + temp_offset)
    return Correction * rho0 * (1 + (constlambda * h / T0)) ** (-(g0 / (R * constlambda) + 1))


# def Power(Pto, h, temp_offset, rho0=1.225, T0=288.15, constlambda=-0.0065):
#    return Pto*(Density(h, temp_offset)/rho0)**(3/4)


def Cd(cdo, A, e, Cl):
    return cdo + Cl**2 / (np.pi * A * e)


"""========== Flight Operations =========="""


def Stallspeedx(h, temp_offset, Vs, Clmax, T0=288.15, constlambda=-0.0065):
    # Temp = T0 + constlambda * h  # temperature correction is for temperature offset
    return 0.5 * Density(h, temp_offset) * Vs**2 * Clmax


def Takeoff(TOP, WS, h, temp_offset, ClmaxTO, rho0=1.225, T0=288.15, constlambda=-0.0065):
    constsigma = Density(h, temp_offset) / rho0
    Clto = ClmaxTO  # / 1.21  # roskam pt1 p95
    return TOP * Clto * constsigma / WS


def Landingx(Clmax_land, h, temp_offset, sland, f, T0=288.15, constlambda=-0.0065):
    return (Clmax_land * Density(h, temp_offset) * sland / 0.305199384478051392) / (
        2 * f
    )  # 0.5915=including airborne phase, 0.305199384478051392=only ground run


def Cruise(
    etap,
    h,
    temp_offset,
    Cd0,
    Vcruise,
    WS,
    A,
    e,
    cruisepowersetting,
    cruiseweightwrtMTOW,
    rho0=1.225,
    T0=288.15,
    constlambda=-0.0065,
):
    a = (cruisepowersetting / cruiseweightwrtMTOW) * etap * (Density(h, temp_offset) / rho0) ** 0.75
    b = (
        (Cd0 * 0.5 * Density(h, temp_offset) * Vcruise**3) / (cruiseweightwrtMTOW * WS)
        + (cruiseweightwrtMTOW * WS) / (np.pi * A * e * 0.5 * Density(h, temp_offset) * Vcruise)
    ) ** (-1)
    return a * b


def Climbrate(etap, A, e, h, temp_offset, Cdo, climbrate, WS, T0=288.15, constlambda=-0.0065):
    a = climbrate + (np.sqrt(WS) * np.sqrt(2 / (Density(h, temp_offset)))) / (1.345 * (A * e) ** 0.75 / (Cdo**0.25))
    return etap / a


def Climbgradient(etap, WS, climbgradient, V, A, e, Clmax, Clsafetyfactor, Cdo, h, temp_offset):
    return etap / (
        np.sqrt(WS)
        * (climbgradient + Cd(Cdo, A, e, Clmax / Clsafetyfactor) / (Clmax / Clsafetyfactor))
        * np.sqrt(2 / (Density(h, temp_offset) * (Clmax / Clsafetyfactor)))
    )


# def Manouvering(Cdo, h, V, WS, nmax, A, e, etap, T0=288.15, constlambda=-0.0065):
#    return ((Cdo*0.5*Density(h, temp_offset)*V**3/WS+WS*nmax**2/(np.pi*A*e*0.5*Density(h, temp_offset)*V))/etap)**(-1)
