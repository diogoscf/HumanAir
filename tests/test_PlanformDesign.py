import sys
import os
from math import isclose
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.AerodynamicDesign.PlanformDesign import Planform


def test_WingSurfaceArea():
    AR = 8.0
    Taper = 0.3
    QuarterChordSweep = 25.0
    tc = 0.1
    MTOW = 5000.0
    WS = 100.0
    S = 1

    planform = Planform(AR, Taper, QuarterChordSweep, tc, MTOW=MTOW, WS=WS, S=S)

    assert isclose(planform.WingSurfaceArea(), S, rel_tol=1e-3)


def test_WingSpan():
    AR = 8.0
    Taper = 0.3
    QuarterChordSweep = 25.0
    tc = 0.1
    MTOW = 5000.0
    WS = 100.0
    S = 1

    planform = Planform(AR, Taper, QuarterChordSweep, tc, MTOW=MTOW, WS=WS, S=S)

    assert isclose(planform.WingSpan(), np.sqrt(1 * 8), rel_tol=1e-3)


def test_RootChord():
    AR = 8.0
    Taper = 0.3
    QuarterChordSweep = 25.0
    tc = 0.1
    MTOW = 5000.0
    WS = 100.0
    S = 1

    planform = Planform(AR, Taper, QuarterChordSweep, tc, MTOW=MTOW, WS=WS, S=S)

    assert isclose(planform.RootChord(), 2 * 1 / ((1 + 0.3) * np.sqrt(1 * 8)), rel_tol=1e-3)


def test_TipChord():
    AR = 8.0
    Taper = 0.3
    QuarterChordSweep = 25.0
    tc = 0.1
    MTOW = 5000.0
    WS = 100.0
    S = 1

    planform = Planform(AR, Taper, QuarterChordSweep, tc, MTOW=MTOW, WS=WS, S=S)

    assert isclose(planform.TipChord(), 0.3 * 2 * 1 / ((1 + 0.3) * np.sqrt(1 * 8)), rel_tol=1e-3)


def test_MAC():
    AR = 8.0
    Taper = 0.3
    QuarterChordSweep = 25.0
    tc = 0.1
    MTOW = 5000.0
    WS = 100.0
    S = 1

    planform = Planform(AR, Taper, QuarterChordSweep, tc, MTOW=MTOW, WS=WS, S=S)

    assert isclose(
        planform.MAC(), (2 * 2 * 1 / ((1 + 0.3) * np.sqrt(1 * 8)) / 3) * (1 + 0.3 + 0.3**2) / (1 + 0.3), rel_tol=1e-3
    )


def test_MAC_y():
    AR = 8.0
    Taper = 0.3
    QuarterChordSweep = 25.0
    tc = 0.1
    MTOW = 5000.0
    WS = 100.0
    S = 1

    planform = Planform(AR, Taper, QuarterChordSweep, tc, MTOW=MTOW, WS=WS, S=S)

    assert isclose(planform.MAC_y(), np.sqrt(1 * 8) / 6 * (1 + 2 * 0.3) / (1 + 0.3), rel_tol=1e-3)


def test_HalfChordSweep():
    AR = 8.0
    Taper = 0.3
    QuarterChordSweep = 25.0
    tc = 0.1
    MTOW = 5000.0
    WS = 100.0
    S = 1

    planform = Planform(AR, Taper, QuarterChordSweep, tc, MTOW=MTOW, WS=WS, S=S)

    assert isclose(planform.HalfChordSweep(), 21.75, rel_tol=1e-2)


def test_t_root_max():
    AR = 8.0
    Taper = 0.3
    QuarterChordSweep = 25.0
    tc = 0.1
    MTOW = 5000.0
    WS = 100.0
    S = 1

    planform = Planform(AR, Taper, QuarterChordSweep, tc, MTOW=MTOW, WS=WS, S=S)

    assert isclose(planform.t_root_max(), 0.1 * 2 * 1 / ((1 + 0.3) * np.sqrt(1 * 8)), rel_tol=1e-3)
