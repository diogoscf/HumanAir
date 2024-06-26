nautical_mile = 1852  # meter
foot = 0.3048  # meter
G = 9.80665  # m/s^2


def ft_to_m(a):
    return a * foot


def m_to_ft(a):
    return a / foot


def nm_to_m(a):
    return a * nautical_mile


def m_to_nm(a):
    return a / nautical_mile


def ft_min_to_m_s(a):
    return a * foot / 60


def m_s_to_ft_min(a):
    return a / foot * 60


def kt_to_m_s(a):
    return a * nautical_mile / 60 / 60


def m_s_to_kt(a):
    return a / nautical_mile * 60 * 60


def N_to_lbs(a):
    return 0.224808943 * a


def lbs_to_N(a):
    return a * 4.4482216153


def m_squared_to_ft_squared(a):
    return 10.7639104 * a


def W_to_hp(a):
    return 0.00134102209 * a
