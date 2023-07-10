import abc
import numpy as np

from constants import KBOLTZ, 

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
        ii = atom.levels.index(levii)
        jj = atom.levels.index(levjj)
        self.ii = ii
        self.jj = jj

        self.gii = atom.degen(ii)
        self.gjj = atom.degen(jj)

        self.Eii = atom.energy(ii)
        self.Ejj = atom.energy(jj)

        # ---- Derived and internal parameters
        self._pipj = None

    @property
    def _pipj(self):
        """ Calculate probability ratio of probabilities p_i/p_j between energy levels i and j.
     

        p_i/p_j = exp[(E_j - E_i)/kT]
        """
        if self._pipj is None:
            deltaE = self.Ejj - self.Eii
            pipj = np.exp(deltaE/KBOLTZ/self.temp)
            self._pipj = pipj

        return self._pipj

    def _ninj(self):
        """ Calculate the ratio of number densities n_i/n_j between energy levels i and j.
        
        """
        atom = self.atom




        

    
    
        

class Atom():
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





