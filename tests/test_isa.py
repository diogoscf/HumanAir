import os
import sys
import math

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from HumanAir.isa import (
    round_sig,
    t1,
    p1,
    rho1,
    isa,
    FL_to_m,
)


def test_round_sig():
    assert round_sig(12345, 3) == 12300
    assert round_sig(0.0012345, 3) == 0.00123
    assert round_sig(12345, 1) == 10000
    assert round_sig(0.00012345, 2) == 0.00012


def test_t1():
    assert math.isclose(t1(288.15, -0.0065, 0, 11000), 216.65, rel_tol=1e-2)
    assert math.isclose(t1(216.65, 0.001, 20000, 32000), 228.65, rel_tol=1e-2)


def test_p1():
    assert math.isclose(p1(101325, 288.15, -0.0065, 0, 11000), 22632, rel_tol=1e-2)
    assert math.isclose(p1(22632, 216.65, 0.0, 11000, 20000), 5474.9, rel_tol=1e-2)


def test_rho1():
    assert math.isclose(rho1(101325, 288.15), 1.225, rel_tol=1e-2)
    assert math.isclose(rho1(22632, 216.65), 0.36391, rel_tol=1e-2)


def test_isa():
    t, p, rho = isa(0)
    assert math.isclose(t, 288.15, rel_tol=1e-2)
    assert math.isclose(p, 101325, rel_tol=1e-2)
    assert math.isclose(rho, 1.225, rel_tol=1e-2)

    t, p, rho = isa(11000)
    assert math.isclose(t, 216.65, rel_tol=1e-2)
    assert math.isclose(p, 22632, rel_tol=1e-2)
    assert math.isclose(rho, 0.36391, rel_tol=1e-2)


def test_FL_to_m():
    assert math.isclose(FL_to_m(100), 3048, rel_tol=1e-2)
    assert math.isclose(FL_to_m(200), 6096, rel_tol=1e-2)
