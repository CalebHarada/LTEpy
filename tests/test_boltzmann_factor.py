import numpy as np
import pytest

import sys
sys.path.append('../LTEpy')

from LTEpy import lte, atom
from LTEpy.constants import EVOLT



def test_bfact():
    """ Check that the boltzman factors are calculated correctly for the
    hydrogen atom at 1000 K.
    
    """

    # hydrogen boltzmann factors at 1000 K
    hbf_1000K = np.array([3.475e68, 1.365e17, 4.127e7, 1.922e4, 551.6])
    # hydrogen atom
    hydrogen=atom.Hydrogen(levels=np.arange(1,6))

    temp = 1000 # Kelvin
    bfact = lte.Boltzmann_Factor(temp, hydrogen).bfact
    good = np.isclose(bfact, hbf_1000K, rtol=1e-3)
    bads = np.logical_not(good)
    if np.any(bads):
        err = f"Levels {hydrogen.levels[bads]=} with {bfact[bads]=} should have boltzmann factors {hbf_1000K[bads]}"
        raise ValueError(err)

def test_zero_Kelvin():
    """ Make sure the code raises a ValueError if you try to use 0 Kelvin temperature.
    
    """
    hydrogen=atom.Hydrogen(levels=np.arange(1,6))
    temp = 0 # Kelvin
    
    with pytest.raises(ValueError):
        bfact = lte.Boltzmann_Factor(temp, hydrogen).bfact

    
if __name__ == "__main__":
    test_zero_Kelvin()
    test_bfact()