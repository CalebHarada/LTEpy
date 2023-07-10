
import numpy as np

from LTEpy.constants import EVOLT, KBOLTZ


class Atom():
    """ Class for storing energy levels and degeneracies for any atom.
    
    """

    def __init__(self, name, gdegen, energy, levels=None):
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
        self.gdegen = gdegen
        self.energy = energy

        if levels is None:
            levels = np.arange(1,len(gdegen))

        self.levels = levels


    def add_level(self, gdegen, energy, level):
        raise NotImplementedError("'add_level' is not yet implemented")

class Hydrogen(Atom):
    """ Class for a hydrogen atom, including hydrogen-specific energy level functions.
    
    """

    _ZNUM = 1


    def energy_at_level(self, levels):
        """ Calculate the energy at each level of a hydrogen atom.

        Parameters
        ----------
        level : arraylike
            Energy level, n

        Returns
        -------
        energy : arraylike
            Energy of each energy level, in cgs units (ergs)

        TODO: Use more precise version of this eq. using Z and rydberg const
        E = -13.6 eV / n^2
        """
        energy = -13.6*EVOLT/levels**2
        return energy
    
    def gdegen_at_level(self, levels):
        """ Calculate the degeneracy at each level of a hydrogen atom.

        Parameters
        ----------
        level : arraylike
            Energy level, n

        Returns
        -------
        gdegen : arraylike
            Degeneracy of each energy level
        
        g = 2 n^2
        """
        gdegen = 2*levels**2
        self.gdegen=gdegen
        return gdegen

    
    def partition_function():
        """ Calculate the 
        
        
        """
        raise NotImplementedError("'partition_function() not yet implemented")