import numpy as np

from LTEpy.constants import EVOLT, KBOLTZ


class Atom():
    """ Class for storing energy levels and degeneracies for any atom.
    
    Attributes : 
        levels (NDarray of ints) energy levels
        energy (NDarray of floats): Energy of each energy level, in ergs
        gdegen (NDarray of ints): Degeneracy of each energy level
    """

    def __init__(self, gdegen, energy, levels=None):
        """ 
        Args:
            gdegen (array): degeneracies.
            energy (array): energy at each level.
            levels (array or None): number of each energy level corresponding to gdegen and energy.
        
        """

        if np.any(gdegen<=0):
            raise ValueError(f"All {gdegen=} must be positive.")
        self.gdegen = gdegen

        if np.any(energy>=0):
            raise ValueError(f"All {energy=} must be negative for bound states.")
        self.energy = energy

        if levels is None:
            levels = np.arange(1,len(gdegen)+1)

        self.levels = levels
    
    def boltzmann_factor(self, levii, temp):
        """ Calculate the Boltzmann factor for a given energy level.
        
        Args:
            levii (int): Energy level, ii. Must be included in levels.
            temp (float): Temperature in Kelvin.


        Returns:
            bfact (float): Boltzmann factor, proportional to the probability of being in state ii.
        
        NOTE: Not tested
        """
        ii = list(self.levels).index(levii)
        Eii = self.energy[ii]

        bfact = np.exp(-Eii/KBOLTZ/temp)
        return bfact
    
    def partition_function(self, temp):
        """ Calculate the partition function, the sum of all the Boltzmann factors,
        using all levels belonging to the atom.
        
        Args:
            temp (float): Temperature in Kelvin.

        NOTE: Not tested
        """
        sum = 0
        for lev in self.levels:
            sum += self.boltzmann_factor(lev, temp)
        return sum


class Hydrogen(Atom):
    """ Class for hydrogen atoms.
    
    """

    def __init__(self, levels=np.arange(1,11)):
        """ 
        Args:
            levels (NDarray of integers): Levels at which to calculate energy and degeneracy, default 1 to 10.
        
        """
        self.levels = levels
        self.gdegen = self.gdegen_at_level(levels)
        self.energy = self.energy_at_level(levels)

    def energy_at_level(self, levels):
        """ Calculate the energy at each level of a hydrogen atom.

        Args: 
            levels (arraylike of integers): Energy level(s).

        Returns
            energy (arraylike): Energy at (each) level, in ergs.
        
        ..math::
            E_n = (-13.6\mathrm{eV}) / n^2

        TODO: Use more precise/generalizable version of this eq. 
        """

        energy = -13.6*EVOLT/levels**2
        return energy
    
    def gdegen_at_level(self, levels):
        """ Calculate the degeneracy at each level of a hydrogen atom.

        Args:
            level (arraylike of integers): Energy level(s)

        Returns:
            gdegen (arraylike): Degeneracy of each energy level
        
        .. math::
            g(n) = 2 n^2
        """
        gdegen = 2*levels**2
        return gdegen

