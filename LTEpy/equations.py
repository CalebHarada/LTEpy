import abc
import numpy as np

from LTEpy.constants import KBOLTZ

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

class Boltzmann_Factor(_LTE): 
    """ Class for calculating Boltzmann Factor 
    between two energy levels of a given atom.
    
    TODO: Add subclass, Gibbs factor, that also allows for a mu term.
    """ 
    def __init__(self, temp, atom): # , levjj):
        """ Initialize

        Parameters
        ----------
        temp : scalar
            Temperature in K
        atom : equations.Atom Object
            Atom, contains levels, energy levels, and degeneracies.

        """
        self.temp = temp
        self.atom = atom
        self.bfact = None


    def boltzmann_factor(self):
        """ Calculate probability ratio of probabilities p_i/p_j between energy levels i and j.
     

        factor = exp[ - E_i)/kT]
        """
        if self.bfact is None:
            atom = self.atom
            bfact = np.zeros_like(atom.levels)
            for ii, lev in enumerate(atom.levels):
                energy = atom.energy[ii]
                bfact[ii] = (np.exp(-np.float64(energy)/KBOLTZ/self.temp))
        self.bfact = bfact

        return self.bfact


    # def ninj(self):
    #     """ Calculate the ratio of number densities n_i/n_j between energy levels i and j
    #     from the probability ratio and degeneracies.
        
        
    #     n_i/n_j = (g_1/g_2) * exp[(E_j - E_i)/kT] = (g_1/g_2) * (p_i/p_j)
    #     """

    #     atom = self.atom
    #     gigj = self.gii / self.gjj
    #     pipj = self._pipj()
    #     self.ninj = gigj * pipj

    #     return self.ninj



        

    
        







