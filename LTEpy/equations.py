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


        # These values are calculated as needed by the class when the corresponding methods are called
        # self._pipj = None          #: Probability ratio of ith to jth energy levels
        # self._deltaE = None        #: Difference between ith and jth energy levels)
        

    # @property
    # def _deltaE(self):
    #     """ Calculate the difference between energy levels.
        
    #     """
    #     if self._deltaE is None:
    #         deltaE = self.Ejj - self.Eii
    #         self._deltaE = deltaE

    #     return self._deltaE

    def _pipj(self):
        """ Calculate probability ratio of probabilities p_i/p_j between energy levels i and j.
     

        p_i/p_j = exp[(E_j - E_i)/kT]
        """
        pipj = np.exp(-(self.Eii - self.Ejj)
                        /KBOLTZ/self.temp)
        self._pipj = pipj

        return self._pipj


    def ninj(self):
        """ Calculate the ratio of number densities n_i/n_j between energy levels i and j
        from the probability ratio and degeneracies.
        
        
        n_i/n_j = (g_1/g_2) * exp[(E_j - E_i)/kT] = (g_1/g_2) * (p_i/p_j)
        """

        atom = self.atom
        gigj = self.gii / self.gjj
        pipj = self._pipj()
        self.ninj = gigj * pipj

        return self.ninj



        

    
        







