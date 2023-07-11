import numpy as np

from LTEpy.constants import EVOLT, KBOLTZ


class Atom():
    """ Class for storing energy levels and degeneracies for any atom.
    
    Attributes : 
    levels : NDarray of ints
        Energy levels
    energy : NDarray of floats
        Energy of each energy level, in ergs
    gdegen : NDarray of ints
        Degeneracy of each energy level
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
    
    def boltzmann_factor(self, levii, temp):
        """ Calculate the Boltzmann factor for a given energy level.
        gi
        Parameters
        ----------
        levii : int
            Energy level, ii. Must be included in levels
        temp : arraylike


        Returns
        -------
        boltzfact : scalar
            Boltzmann factor, proportional to the probability of being in state ii

        """
        ii = list(self.levels).index(levii)
        Eii = self.energy[ii]

        boltzfact = np.exp(-Eii/KBOLTZ/temp)
        return boltzfact
    
    def partition_function(self, temp):
        """ Calculate the partition function, the sum of all the Boltzmann factors,
        using all levels belonging to the atom.
        
        
        """
        sum = 0
        for lev in self.levels:
            sum += self.boltzmann_factor(lev, temp)
        return sum


class Hydrogen(Atom):
    """ Class for a hydrogen atom, including hydrogen-specific energy level functions.
    
    """

    def __init__(self, name='hydrogen', levels=np.arange(1,11)):
        """ 
        Parameters
        ----------
        name : str
            name of atom
        levels : NDarray of integers
            Levels at which to calculate energy and degeneracy, 
            default 1 to 10.
        
        """
        self.name = name
        self.levels = levels
        self.gdegen = self.gdegen_at_level(levels)
        self.energy = self.energy_at_level(levels)

    def energy_at_level(self, levels):
        """ Calculate the energy at each level of a hydrogen atom.

        Parameters
        ----------
        levels : arraylike
            Energy level(s), n

        Returns
        -------
        energy : arraylike
            Energy of each energy level, in cgs units (ergs)

        TODO: Use more precise/generic version of this eq. using Z and rydberg const
        E = -13.6 eV / n^2
        """
        energy = -13.6*EVOLT/levels**2
        return energy
    
    def gdegen_at_level(self, levels):
        """ Calculate the degeneracy at each level of a hydrogen atom.

        Parameters
        ----------
        level : arraylike
            Energy level(s), n

        Returns
        -------
        gdegen : arraylike
            Degeneracy of each energy level
        
        g = 2 n^2
        """
        gdegen = 2*levels**2
        return gdegen

