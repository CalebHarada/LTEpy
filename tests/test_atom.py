import numpy as np
import pytest

import sys
sys.path.append('../')

from LTEpy import lte, atom

def test_atom_constructor():
    """ Test building an atom with the Atom() class for user input energy levels
    and degeneracies.
    
    """
    pass

def test_negative_energy():
    """ Make sure a ValueError is raised if you try to use positive energy levels.
    
    """
    pass

def test_positive_gdegen():
    """ Make sure a ValueError is raised if you try to use negative or 0 degeneracies.
    
    """
    pass

