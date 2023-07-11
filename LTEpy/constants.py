""" Numerical Constants

Constants in CGQ units. 

-----
* [cm] = centimeter
* [g] = gram
* [s] = second
* [erg] = cm^2 * g / s^2
* [Jy] = jansky = [erg/s/cm^2/Hz]
* [fr] franklin = statcoulomb = electro-static unit [esu]
* [K] Kelvin


(source: https://github.com/nanograv/holodeck/blob/main/holodeck/constants.py)
"""

import numpy as np
import astropy as ap
import astropy.constants  # noqa

# ---- Fundamental Constants
NWTG = ap.constants.G.cgs.value             #: Newton's Gravitational Constant [cm^3/g/s^2]
SPLC = ap.constants.c.cgs.value             #: Speed of light [cm/s]
MELC = ap.constants.m_e.cgs.value           #: Electron Mass [g]
MPRT = ap.constants.m_p.cgs.value           #: Proton Mass [g]
AMU = ap.constants.u.cgs.value              #: Atomic Mass Unit [g]
QELC = ap.constants.e.gauss.value           #: Fundamental unit of charge (electron charge) [fr]
KBOLTZ = ap.constants.k_B.cgs.value         #: Boltzmann constant [erg/K]
HPLANCK = ap.constants.h.cgs.value          #: Planck constant [erg/s]
SIGMA_SB = ap.constants.sigma_sb.cgs.value  #: Stefan-Boltzmann constant [erg/cm^2/s/K^4]
SIGMA_T = ap.constants.sigma_T.cgs.value    #: Thomson/Electron -Scattering cross-section [cm^2]
B_WIEN = ap.constants.b_wien.cgs.value      #: Wien's Displacement Constant [cm K]


# ---- Typical astronomy units
MSOL = ap.constants.M_sun.cgs.value                                #: Solar Mass [g]
LSOL = ap.constants.L_sun.cgs.value                                #: Solar Luminosity [erg/s]
RSOL = ap.constants.R_sun.cgs.value                                #: Solar Radius [cm]
PC = ap.constants.pc.cgs.value                                     #: Parsec [cm]
AU = ap.constants.au.cgs.value                                     #: Astronomical Unit [cm]
YR = ap.units.year.to(ap.units.s)                                  #: year [s]
EVOLT = ap.units.eV.to(ap.units.erg)                               #: Electronvolt in ergs
JY = ap.units.jansky.to(ap.units.g/ap.units.s**2)                  #: Jansky [erg/s/cm^2/Hz]

