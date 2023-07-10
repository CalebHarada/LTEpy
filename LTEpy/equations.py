import abc
import numpy as np
import matplotlib.pyplot as plt


from constants import *


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

    def compute_B_nu(self, nu):
        """
        
        """

        spec_intensity_nu = (2 * HPLANCK * nu**3 / SPLC**2) / (np.exp(HPLANCK * nu / (KBOLTZ * self.temp)) - 1)

        return spec_intensity_nu
    
    def compute_B_lambda(self, wl):
        """
        
        """

        spec_intensity_lamb = (2 * HPLANCK * SPLC**2 / wl**5) / (np.exp(HPLANCK * SPLC / (wl * KBOLTZ * self.temp)) - 1)

        return spec_intensity_lamb
    

    def plot_B_nu(self, nu_1, nu_2, N_nu=500, lw=1, **kwargs):
        """
        
        """

        nus = np.linspace(nu_1, nu_2, N_nu)

        fig, ax = plt.subplots()
        ax.plot(nus, self.compute_B_nu(nus), lw=lw, **kwargs)
        
        plt.show()


    def plot_B_lambda(self, wl_1, wl_2, N_wl=500, lw=1, **kwargs):
        """
        
        """

        wls = np.linspace(wl_1, wl_2, N_wl)

        fig, ax = plt.subplots()
        ax.plot(wls, self.compute_B_lambda(wls), lw=lw, **kwargs)
        
        plt.show()







class Maxwell_Boltzmann(_LTE):
    """ 
    
    """

class Boltzmann_Gibbs(_LTE): 
    """ 
    
    """ 
    #def __init__(self, temp, atom1, atom2):
    
        

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





