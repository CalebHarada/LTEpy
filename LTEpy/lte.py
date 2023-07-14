import abc
import numpy as np
import matplotlib.pyplot as plt

from LTEpy.constants import *
import LTEpy.plot as plot


class _LTE(abc.ABC):
    """
    
    """

    def __init__(self, temp):
        self.temp = temp



class Planck(_LTE):
    """Planck's Law.

    Class for calculating Planck's Law for a blackbody
    in Local Thermodynamic Equilibrium (either in wavelength or frequency).

    Args:
        temp (float): Temperature in K.
    """

    def __init__(self, temp):

        super().__init__(temp)

        if temp > 0:
            self.temp = temp
        elif temp == 0:
            raise ValueError('Temperature must be greater than 0.')
        else:
            raise ValueError('Temperature cannot be negative.')
        

    def set_temp(self, temp):
        """Set temperature.

        Change the temperature of this Planck LTE object.

        Args:
            temp (float): Temperature in K.
        """

        if temp > 0:
            self.temp = temp
        elif temp == 0:
            raise ValueError('Temperature must be greater than 0.')
        else:
            raise ValueError('Temperature cannot be negative.')


    def compute_B_nu(self, nu):
        """Spectral radiance in frequency space.

        Calculate the spectral radiance (or specific intensity) :math:`B_\\nu`
        for a given frequency in cgs units.

        .. math::
            B_\\nu(T) = \\frac{2 h \\nu^3}{c^2} \\frac{1}{e^{h \\nu / k T} - 1}

        Args:
            nu (float or :obj:`np.array`): Frequency in Hz.

        Returns:
            float or :obj:`np.array`: Spectral radiance in erg/s/cm^2/sr/Hz.
        """

        spec_intensity_nu = (2 * HPLANCK * nu**3 / SPLC**2) / (np.exp(HPLANCK * nu / (KBOLTZ * self.temp)) - 1)

        return spec_intensity_nu
    
    
    def compute_B_lambda(self, wl):
        """Spectral radiance in wavelength space.

        Calculate the spectral radiance (or specific intensity) :math:`B_\\lambda`
        for a given wavelength.

        .. math::
            B_\\lambda(T) = \\frac{2 h c^2}{\\lambda^5} \\frac{1}{e^{h c / \\lambda k T} - 1}

        Args:
            wl (float or :obj:`np.array`): Wavelength in cm.

        Returns:
            float or :obj:`np.array`: Spectral radiance in erg/s/cm^2/sr/nm.
        """

        spec_intensity_lamb = (2 * HPLANCK * SPLC**2 / wl**5) / (np.exp(HPLANCK * SPLC / (wl * KBOLTZ * self.temp)) - 1)
        spec_intensity_lamb *= 1e-7

        return spec_intensity_lamb
    

    def compute_lambda_max(self):
        """Wien's Displacement Law.

        Calculate the wavelength of maximum spectral radiance :math:`\\lambda_{max}`
        from Wien's Law in cm.

        .. math::
            \\lambda_{max} = \\frac{b}{T}

        Returns:
            float: Peak wavelength in cm.
        """

        wl_max = B_WIEN / self.temp
        
        return wl_max

        
    def plot_B_nu(self, nu_1, nu_2, N_nu=500, lw=1, log_scale=True, ax=None, **ax_kwargs):
        """Plot spectral radiance in frequency space.

        Plot the spectral radiance :math:`B_\\nu` between two frequencies.

        Args:
            nu_1 (float): First frequency in Hz.
            nu_2 (float): Second frequency in Hz.
            N_nu (int, optional): Number of frequency points to plot. Defaults to 500.
            lw (int, optional): Plot line width. Defaults to 1.
            log_scale (bool, optional): Option to plot with a log scale. Defaults to True.
            ax (:obj:`matplotlib.pyplot.Axes`, optional): Matplotlib axis for plotting. Defaults to None.
            **ax_kwargs: Keyword arguments passed to :obj:`matplotlib.pyplot.Axes` object.

        Returns:
            :obj:`matplotlib.pyplot.Axes`: Matplotlib axis
        """

        nus = np.linspace(nu_1, nu_2, N_nu)  # define frequency array

        if not ax:
            _, ax = plt.subplots()

        ax.plot(nus, self.compute_B_nu(nus), lw=lw, label="{:} K".format(self.temp), **ax_kwargs)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("$B_\\nu$ (erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$ Hz$^{-1}$)")
        ax.legend(loc=1)
        
        if log_scale:
            ax.set_xscale("log")
            ax.set_yscale("log")
        
        return ax


    def plot_B_lambda(self, wl_1, wl_2, N_wl=500, lw=1, log_scale=True, ax=None, **ax_kwargs):
        """Plot spectral radiance in wavelength space.

        Plot the spectral radiance :math:`B_\\lambda` between two wavelengths.

        Args:
            wl_1 (float): First wavelength in nm. 
            wl_2 (float): Second wavelength in nm.
            N_wl (int, optional): Number of wavelength points to plot. Defaults to 500.
            lw (int, optional): Plot line width. Defaults to 1.
            log_scale (bool, optional): Option to plot with a log scale. Defaults to True.
            ax (:obj:`matplotlib.pyplot.Axes`, optional): Matplotlib axis for plotting. Defaults to None.
            **ax_kwargs: Keyword arguments passed to :obj:`matplotlib.pyplot.Axes` object.

        Returns:
            :obj:`matplotlib.pyplot.Axes`: Matplotlib axis
        """

        wls = np.linspace(wl_1, wl_2, N_wl)  # define wavelength array

        if log_scale:
            wls = np.geomspace(wl_1, wl_2, N_wl)

        if not ax:
            _, ax = plt.subplots()

        ax.plot(wls, self.compute_B_lambda(wls * 1e-7), lw=lw, label="{:} K".format(self.temp), **ax_kwargs)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("$B_\lambda$ (erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$ nm$^{-1}$)")
        ax.legend(loc=1)

        if log_scale:
            ax.set_xscale("log")
            ax.set_yscale("log")
        
        return ax
    
        
    def plot_lambda_max(self, ax, **vline_kwargs):
        """Plot peak wavelength.

        Plot a vertical line at the peak wavelength :math:`\\lambda_{max}` in nm according to Wien's Law.

        Args:
            ax (:obj:`matplotlib.pyplot.Axes`): Matplotlib axis for plotting.
            **vline_kwargs: Keyword arguments passed to :obj:`matplotlib.pyplot.axvline`
        """

        wl_max = self.compute_lambda_max()

        ax.axvline(wl_max * 1e7, **vline_kwargs)
        







class Maxwell_Boltzmann(_LTE):
    """Maxwell-Boltzmann Distribution.

    Class for calculating the Maxwell-Boltzmann Distribution for a system
    in Local Thermodynamic Equilibrium, given a particle mass.

    Args:
        temp (float): Temperature in K.
        mass (float): Particle mass in AMU.
    """

    def __init__(self, temp, mass):

        super().__init__(temp)
        self.temp = temp
        self.mass = mass


    def set_mass(self, mass):
        """Set mass.

        Change the mass of this Maxwell-Boltzmann LTE object.

        Args:
            mass (float): Mass in AMU.
        """

        self.mass = mass
    

    def set_temp(self, temp):
        """Set temperature.

        Change the temperature of this Maxwell-Boltzmann LTE object.

        Args:
            temp (float): Temperature in K.
        """

        self.temp = temp
    

    def compute_maxwell_boltzmann(self, speed):
        """Maxwell-Boltzmann speed distribution.

        Compute the probabilty density of a particle with a given mass having a given speed.

        .. math::
            f(v) = \\bigg(\\frac{m}{2 \pi k T}\\bigg)^{3/2} 4 \pi v^2 e^{-m v^2 / 2 k T}

        Args:
            speed (float or :obj:`np.array`): Speed in cm/s.

        Returns:
            float or :obj:`np.array`: Probability density s/cm
        """

        mass = self.mass * AMU

        f_v = 4 * np.pi * speed**2 * (mass / (2 * np.pi * KBOLTZ * self.temp))**1.5 * np.exp(-mass * speed**2 / (2 * KBOLTZ * self.temp))

        return f_v
    
    
    def plot_fv(self, speed_1, speed_2, N_speed=500, lw=1, ax=None, **ax_kwargs):
        """Plot the Maxwell-Boltzmann Distribution.

        Plot the distributino of particle speeds defined by the Maxwell-Boltzmann
        Distribution between two speeds.

        Args:
            speed_1 (float): First speed in cm/s.
            speed_2 (float): Second speed in cm/s.
            N_speed (int, optional): Number of speed points to plot. Defaults to 500.
            lw (int, optional): Plot line width. Defaults to 1.
            ax (:obj:`matplotlib.pyplot.Axes`, optional): Matplotlib axis for plotting. Defaults to None.
            **ax_kwargs: Keyword arguments passed to :obj:`matplotlib.pyplot.Axes` object.

        Returns:
            :obj:`matplotlib.pyplot.Axes`: Matplotlib axis
        """

        speeds = np.linspace(speed_1, speed_2, N_speed)  # define speed array

        if not ax:
            _, ax = plt.subplots()

        ax.plot(speeds, self.compute_maxwell_boltzmann(speeds), lw=lw, **ax_kwargs)
        ax.set_xlabel("Speed (cm s$^{-1}$)")
        ax.set_ylabel("Probability density (s cm$^{-1}$)")
        
        return ax



class Boltzmann_Factor(_LTE): 
    """ Class for calculating Boltzmann Factor 
    between two energy levels of a given atom.
    
    """ 
    def __init__(self, temp, atom): # , levjj):
        """ Initialize

        Args:
            temp (float): Temperature in K
            atom (:obj:'atom.Atom'): Atom, contains levels, energy levels, and degeneracies.

        """

        # check that temperature is positive
        if temp == 0: 
            err = f"{temp=} must be > 0 Kelvin."
            raise ValueError(err)
        
        # set attributes
        self.temp = temp
        self.atom = atom
        self._bfact = None

    @property
    def bfact(self):
        """ Boltzmann factor
        
        Calculate Boltzmann factor for each energy level.

        Returns:
            _bfact (:obj:`np.array`) : Boltzmann factor

     
        .. math::
            \\beta(n, T) = \exp \\bigg(\\frac{-E_n}{k_B T}\\bigg)

        """
        if self._bfact is None:
            atom = self.atom
            bfact = np.zeros_like(atom.levels, dtype=float)
            for ii, lev in enumerate(atom.levels):
                eng = atom.energy[ii]
                bfact[ii] = np.float64(np.exp(-(eng)/KBOLTZ/self.temp))
            self._bfact = bfact
        return self._bfact


    def draw_bfact(self, ax, levmin=None, levmax=None, color=None):
        """ Draw the Boltzmann factors
        
        Plot the Boltzmann factors for each energy level at a fixed temperature, 
        and return the line handle.
        
        Args:
            levmin (float or None): Lowest level to plot.
            levman (float or None): Highest level to plot.
            
        Returns:
            :obj:'matplotlib.Line2D.Line': Line handle
        """
        if levmin is not None:
            iimin = list(self.atom.levels).index[levmin]
        else:
            iimin = 0
        if levmax is not None:
            iimax = list(self.atom.levels).index[levmax]
        else:
            iimax = -1
        xx = self.atom.levels[iimin:iimax]
        yy = self.bfact[iimin:iimax]
        label = '$T$=%.2eK' % self.temp

        hh, = ax.plot(xx, yy, '-o', label=label, color=color)
        ax.set_xlabel(plot.LABEL_LEVEL)
        ax.set_ylabel(plot.LABEL_BFACT)

        return hh,

    def plot_bfact(self, levmin=None, levmax=None,):
        """ Plot the Boltzmann factor 
        
        Plot the Boltzmann factor for each energy level at a fixed temperature.

        Args:
            levmin (float or None): Lowest level to plot.
            levman (float or None): Highest level to plot.
            
        Returns:
            :obj:`matplotlib.pyplot.Figure`: Matplotlib figure
            :obj:'matplotlib.Line2D.Line': Line handle
        
        """

        fig, ax = plot.figax(xscale='linear')
        hh, = self.draw_bfact(ax, levmin, levmax)
        ax.set_xlabel(plot.LABEL_LEVEL)
        ax.set_ylabel(plot.LABEL_BFACT)

        return fig, hh



        

    
        







