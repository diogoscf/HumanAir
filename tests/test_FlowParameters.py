import sys
import os
from math import isclose

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.AerodynamicDesign.ISA import ISA
from HumanAir.AerodynamicDesign.FlowParameters import Flow


class MockWing:
    @staticmethod
    def MAC():
        return 10.0


flow_init = Flow(100, ISA(0, 0, -0.0065), MockWing)


def test_mach():
    assert isclose(flow_init.Mach(), 100 / ISA(0, 0, -0.0065).SpeedOfSound(), rel_tol=1e-2)


def test_reynolds():
    assert isclose(
        flow_init.Reynolds(),
        100 * 10.0 * ISA(0, 0, -0.0065).Density() / ISA(0, 0, -0.0065).DynamicViscosity(),
        rel_tol=1e-2,
    )


def test_beta():
    assert isclose(flow_init.Beta(), (1 - (100 / ISA(0, 0, -0.0065).SpeedOfSound()) ** 2) ** 0.5, rel_tol=1e-2)
