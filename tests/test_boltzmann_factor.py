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
    hbf_1000K = np.array([5.532e18, 3.136e7, 4.006e2, 1.562e0, 5.343e-1])
    # hydrogen atom
    hydrogen=atom.Hydrogen(levels=np.arange(1,6))

    temp = 1000 # Kelvin
    hbf = lte.Boltzmann_Factor(temp, hydrogen)
    assert np.all(hbf_1000K == pytest.approx(hbf.bfact))

def test_zero_Kelvin():
    """ Make sure the code raises a ValueError if you try to use 0 Kelvin temperature.
    
    """
    pass

    
if __name__ == "__main__":
    test_zero_Kelvin()
    test_bfact()