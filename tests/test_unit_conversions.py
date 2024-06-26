import os
import sys
import math

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from HumanAir.unit_conversions import (
    ft_to_m,
    m_to_ft,
    nm_to_m,
    m_to_nm,
    ft_min_to_m_s,
    m_s_to_ft_min,
    kt_to_m_s,
    m_s_to_kt,
    N_to_lbs,
    lbs_to_N,
    m_squared_to_ft_squared,
    W_to_hp,
)


def test_m_to_ft():
    assert math.isclose(m_to_ft(10), 32.81, rel_tol=1e-2)


def test_ft_to_m():
    assert math.isclose(ft_to_m(10), 3.05, rel_tol=1e-2)


def test_nm_to_m():
    assert math.isclose(nm_to_m(10), 18520, rel_tol=1e-2)


def test_m_to_nm():
    assert math.isclose(m_to_nm(18520), 10, rel_tol=1e-2)


def test_ft_min_to_m_s():
    assert math.isclose(ft_min_to_m_s(100), 0.508, rel_tol=1e-2)


def test_m_s_to_ft_min():
    assert math.isclose(m_s_to_ft_min(100), 19685.03937, rel_tol=1e-2)


def test_kt_to_m_s():
    assert math.isclose(kt_to_m_s(10), 5.144, rel_tol=1e-2)


def test_m_s_to_kt():
    assert math.isclose(m_s_to_kt(10), 19.4384, rel_tol=1e-2)


def test_N_to_lbs():
    assert math.isclose(N_to_lbs(10), 2.248, rel_tol=1e-2)


def test_lbs_to_N():
    assert math.isclose(lbs_to_N(10), 44.482, rel_tol=1e-2)


def test_m_squared_to_ft_squared():
    assert math.isclose(m_squared_to_ft_squared(100), 1076.39, rel_tol=1e-2)


def test_w_to_hp():
    assert math.isclose(W_to_hp(100), 0.134, rel_tol=1e-2)
