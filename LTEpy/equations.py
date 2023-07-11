import abc
import numpy as np
import matplotlib.pyplot as plt

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

        BUG: This only works for very large temperatures because the bfact is too large otherwise.
        """
        if self.bfact is None:
            atom = self.atom
            bfact = np.zeros_like(atom.levels)
            for ii, lev in enumerate(atom.levels):
                energy = atom.energy[ii]
                bfact[ii] = np.float64(np.exp(-(energy)/KBOLTZ/self.temp))
        self.bfact = bfact

        return self.bfact
    
    def plot_bfact(self, levmin=None, levmax=None):
        xx = self.atom.levels
        yy = self.bfact

        fig, ax = plt.subplots()
        ax.plot(xx, yy)
        ax.set_xlabel('Energy Level, $n$')
        ax.set_ylabel('Boltzmann Factor, $\exp(-E_n/kT)$')
        ax.set_xscale('log')
        ax.set_yscale('log')



        

    
        







