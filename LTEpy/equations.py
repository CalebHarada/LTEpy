import abc
import numpy as np

import constants

class _LTE(abc.ABC):
    """
    
    """

    def __init__(self, temp):
        self.temp = temp



class Planck(_LTE):
    """ 
    
    """


class Maxwell_Boltzmann(_LTE):
    """ 
    
    """

class Boltzmann_Gibbs(_LTE): 
    """ Class for calculating Boltzmann Distribution (aka Gibbs Distribution) for two energy 
    levels of a given atom.
    
    """ 
    def __init__(self, temp, atom, lev1, lev2):
        """ Initialize

        """
        self.temp = temp
        self.atom = atom
        self.lev1 = lev1
        self.lev2 = lev2

    def pipi(self):
        """ Calculate probability ratio p_i/p_j between energy levels i and j
     

        p_i/p_j = exp[(E_j - E_i)/kT]
        """
        atom = self.atom
        ii = atom.levels.index(self.levi)
        jj = atom.levels.index(self.levj)

        deltaE = atom.energy[ii] - atom.energy[jj] 

        pipj = np.exp(deltaE/constants.)


        

    
    
        

class atom():
    """ Class for storing energy levels and degeneracies for a given atom.
    
    """

    def __init__(self, name, degen, energy):
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





