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

class Boltzmann_Gibbs(_LTE): 
    """ Class for calculating Boltzmann Distribution (aka Gibbs Distribution, if mu is included) 
    between two energy levels of a given atom.
    
    """ 
    def __init__(self, temp, atom, levii, levjj):
        """ Initialize

        Parameters
        ----------
        temp : scalar
            Temperature in K
        atom : equations.Atom Object
            Atom
        levii : integer
            First energy level, i. Must be included in atom.levels
        levjj : integer
            Second energy level, j. Must be included in atom.levels


        """
        self.temp = temp
        self.atom = atom
        self.levii = levii
        self.levjj = levjj

        # ---- Parameters from atom
        ii = list(atom.levels).index(levii)
        jj = list(atom.levels).index(levjj)
        self.ii = ii
        self.jj = jj

        self.gii = atom.gdegen[ii]
        self.gjj = atom.gdegen[jj]

        self.Eii = atom.energy[ii]
        self.Ejj = atom.energy[jj]

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



        

    
        







