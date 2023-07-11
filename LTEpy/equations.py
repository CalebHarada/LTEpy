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



        

    
        







