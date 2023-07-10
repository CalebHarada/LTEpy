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

    def __init__(self, temp):
        super().__init__(temp)

        self.temp = temp


class Maxwell_Boltzmann(_LTE):
    """ 
    
    """

class Boltzmann_Gibbs(_LTE): 
    """ 
    
    """ 
    def __init__(self, temp, atom1, atom2):
    
        

class atom():
    """
    
    """
    def __init__(self, name, degen, energy, levels=None):
        """ 
        Parameters
        ----------
        name : str

        degen : arraylike
            degeneracies
        elevels : arraylike
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





