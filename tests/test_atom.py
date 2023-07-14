import numpy as np
import pytest

import sys
sys.path.append('../LTEpy')

from LTEpy import atom
from LTEpy.constants import EVOLT

def test_atom_constructor():
    """ Test building an atom with the Atom() class for user input levels, energy,
    and degeneracies.
    
    """

    levels = np.array([1,2,3,4])
    gdegen = 2*levels**2
    energy = -13.6*EVOLT/levels**2

    atom1=atom.Atom(gdegen, energy)
    atom2=atom.Atom(gdegen, energy, levels)

    # make sure levels match
    if np.any(atom1.levels != atom2.levels):
        err = 'Atom levels not generated correctly'
        raise ValueError(err)

    # check energy
    assert np.all(atom1.energy == energy)
    assert np.all(atom2.energy == energy)

    # check degeneracy
    assert np.all(atom1.gdegen == gdegen)
    assert np.all(atom2.gdegen == gdegen)


def test_negative_energy():
    """ Make sure a ValueError is raised if you try to use positive energies.
    
    """
    levels = np.array([1,2,3,4])
    gdegen = 2*levels**2

    # try a positive energy
    energy = -13.6*EVOLT/levels**2
    energy[0] = 1
    with pytest.raises(ValueError):
        atom1=atom.Atom(gdegen, energy)


def test_positive_gdegen():
    """ Make sure a ValueError is raised if you try to use negative or 0 degeneracies.
    
    """
    levels = np.array([1,2,3,4])
    energy = -13.6/levels**2

    # try a 0 degeneracy
    gdegen = 2*levels**2
    gdegen[0] = 0
    with pytest.raises(ValueError):
        atom1=atom.Atom(gdegen, energy)

    # try a negative degeneracy
    gdegen = 2*levels**2
    gdegen[0] = -1
    with pytest.raises(ValueError):
        atom2=atom.Atom(gdegen, energy)


def test_hydrogen_atom():
    """ Test building an atom with the Atom() class for user input energy levels.
    
    """

    levels = np.array([1,2,3,4])
    gdegen = 2*levels**2
    energy = -13.6*EVOLT/levels**2

    hydrogen=atom.Hydrogen(levels)

    # make sure levels match
    assert np.all(hydrogen.levels == levels)

    # check energy levels
    assert np.all(hydrogen.energy == energy)

    # check degeneracy
    assert np.all(hydrogen.gdegen == gdegen)

    
if __name__ == "__main__":
    test_atom_constructor()
    test_negative_energy()
    test_positive_gdegen()
    test_hydrogen_atom()