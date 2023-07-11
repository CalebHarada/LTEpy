import abc
import numpy as np
import matplotlib.pyplot as plt

from LTEpy.constants import *
from LTEpy.constants import KBOLTZ, EVOLT



class _LTE(abc.ABC):
    """
    
    """

    def __init__(self, temp):
        self.temp = temp



class Planck(_LTE):
    """ Class for calculating the Planck Function (either in wavelength or frequency) given an LTE object

    """

    def __init__(self, temp):
        """ Initialize

        Parameters
        ----------
        temp : scalar
            Temperature in K

        """

        super().__init__(temp)

        self.temp = temp

    def compute_B_nu(self, nu):
        """ Calculate the specific intensity B_nu (in cgs units) for a given frequency

        Parameters
        ----------
        nu : scalar
            frequency in Hz

        
        B_nu = (2*h * nu^3 / c^2) * 1/(exp[h*nu / k*T] - 1)
            
        """

        spec_intensity_nu = (2 * HPLANCK * nu**3 / SPLC**2) / (np.exp(HPLANCK * nu / (KBOLTZ * self.temp)) - 1)

        return spec_intensity_nu
    
    def compute_B_lambda(self, wl):
        """ Calculate the specific intensity B_lambda (in cgs units) for a given wavelength

        Parameters
        ----------
        wl : scalar
            wavelength in cm

        
        B_lambda = (2*h * c^2 / lambda^5) * 1/(exp[h*c / lambda*k*T] - 1)
            
        """

        spec_intensity_lamb = (2 * HPLANCK * SPLC**2 / wl**5) / (np.exp(HPLANCK * SPLC / (wl * KBOLTZ * self.temp)) - 1)

        return spec_intensity_lamb
    

    def plot_B_nu(self, nu_1, nu_2, N_nu=500, lw=1, **kwargs):
        """ Plot specific intensity B_nu between two frequencies
        
        Parameters
        ----------
        nu_1 : scalar
            first frequency in Hz
        nu_2 : scalar
            second frequency in Hz
        N_nu : scalar
            number of frequency points to plot
        lw : scalar
            line width for plotting
        
        """

        nus = np.linspace(nu_1, nu_2, N_nu)

        fig, ax = plt.subplots()
        ax.plot(nus, self.compute_B_nu(nus), lw=lw, **kwargs)
        
        plt.show()


    def plot_B_lambda(self, wl_1, wl_2, N_wl=500, lw=1, **kwargs):
        """ Plot specific intensity B_lambda between two wavelengths
        
        Parameters
        ----------
        wl_1 : scalar
            first wavelength in cm
        wl_2 : scalar
            second wavelength in cm
        N_wl : scalar
            number of wavelength points to plot
        lw : scalar
            line width for plotting
        
        """

        wls = np.linspace(wl_1, wl_2, N_wl)

        fig, ax = plt.subplots()
        ax.plot(wls, self.compute_B_lambda(wls), lw=lw, **kwargs)
        
        plt.show()







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
        label = '$T$=%.2eK' % self.temp

        fig, ax = plt.subplots()
        hh, = ax.plot(xx, yy, label=label)
        ax.set_xlabel('Energy Level, $n$')
        ax.set_ylabel('Boltzmann Factor, $\exp(-E_n/kT)$')
        # ax.set_xscale('log')
        ax.set_yscale('log')

        return fig, hh



        

    
        







