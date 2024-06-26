import sys
import numpy as np
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.AerodynamicDesign.ISA import ISA


def test_temperature():
    # Test the Temperature method of ISA class
    isa_sea_level = ISA(0, 0, -0.0065)  # ISA model at sea level

    result_sea_level = isa_sea_level.Temperature()

    # Expected temperatures
    expected_sea_level = 288.15  # K

    assert np.isclose(result_sea_level, expected_sea_level, rtol=0.01)


def test_density():
    # Test the Density method of ISA class
    isa_sea_level = ISA(0, 0, -0.0065)  # ISA model at sea level

    result_sea_level = isa_sea_level.Density()

    # Expected densities (approximate values)
    expected_sea_level = 1.225  # kg/m^3 (standard sea level density)

    assert np.isclose(result_sea_level, expected_sea_level, rtol=0.01)


def test_dynamic_viscosity():
    # Test the DynamicViscosity method of ISA class
    isa_sea_level = ISA(0, 0, -0.0065)  # ISA model at sea level

    result_sea_level = isa_sea_level.DynamicViscosity()

    # Expected dynamic viscosities
    expected_sea_level = 1.8

    assert np.isclose(result_sea_level, expected_sea_level, rtol=1)


def test_speed_of_sound():
    # Test the SpeedOfSound method of ISA class
    isa_sea_level = ISA(0, 0, -0.0065)  # ISA model at sea level

    result_sea_level = isa_sea_level.SpeedOfSound()

    # Expected speed of sounds (approximate values)
    expected_sea_level = 340.3  # m/s (speed of sound at sea level)

    assert np.isclose(result_sea_level, expected_sea_level, rtol=0.1)
