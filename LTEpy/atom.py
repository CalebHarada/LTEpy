
import numpy as np

import atom
from constants import EVOLT, KBOLTZ,


class Atom():
    """ Class for storing energy levels and degeneracies for any atom.
    
    """

    def __init__(self, name, degen, energy, levels=None):
        """ 
        Parameters
        ----------
        name : str
            name of atom
        degen : dict
            degeneracies
        elevels : dict
            energy levels
        
        """
        self.name = name
        self.degen = degen
        self.energy = energy

        if levels is None:
            levels = np.arange(1,len(degen))

        self.levels = levels


    def add_level(self, degen, energy, level):
        pass 

class Hydrogen(Atom):
    """ Class for a hydrogen atom, including hydrogen-specific energy level functions.
    
    """

    _ZNUM = 1

    def energy_at_level(level):
        """ Calculate the energy level of a hydrogen atom.

        Parameters
        ----------
        level : arraylike
            Energy level, n
        """
        energy = -13.6*EVOLT/level**2
        return energy