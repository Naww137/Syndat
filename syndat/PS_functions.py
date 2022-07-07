#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:18:04 2022

@author: nwx
"""

import numpy as np
import matplotlib.pyplot as plt


# assuming boundary condition selected s.t. shift factor is eliminated for s wave but not others!


ac = 6.7e-15 # m or 6.7 femtometers
M = 63
m = 1
hbar = 6.582119569e-16 # eV-s
c = 2.99792458e8 # m/s
m_eV = 939.565420e6 # eV/c^2

# !!! energy in eV !!!


def k_wavenumber(E):
    """
    Calculate the wavenumber of the compound state nucleus.

    This function calculates the wavenumber of the compound state nucleus
    created by an incident neutron and a Cu-63 atom. Some nucelar parameters 
    are housed within this function, this could be updated for other nuclei.

    Parameters
    ----------
    E : float or numpy.ndarray
        Energy of the incident neutron.

    Returns
    -------
    float or numpy.ndarray
        Returns either the scalar k or vector of k at different energies.
    
    Examples
    --------
    >>> from sample_resparm import PS_functions
    >>> PS_functions.k_wavenumber(10.4)
    697379673687.8605
    >>> PS_functions.k_wavenumber(np.array([1,10,100,1000]))
    array([2.16248257e+11, 6.83837032e+11, 2.16248257e+12, 6.83837032e+12])
    """
    k = 1/hbar * M/(m+M) * np.sqrt(2*m_eV*E) * 1/c
    return k

def PS_recursive(E, ac, orbital_angular_momentum):
    """
    Calculates penetrability and shift functions using recursion.

    This function calculates the centifugal barrier penetrability as well as
    the shift factor for a neutron incident on Cu-63. The recursive implementation
    will allow for high l values. Some nucelar parameters are housed within 
    this function, this could be updated for other nuclei.

    Parameters
    ----------
    E : numpy.ndarray
        Energy of the incident neutron.
    ac : float 
        The channel radius in meters.
    orbital_angular_momentum : int
        Orbital angular momentum of the pair, describes waveform (l).

    Returns
    -------
    S_array : numpy.ndarray
        Array shift factors, each row is an l value, columns traverse energy vector given.
    P_array : numpy.ndarray
        Array penetrability, each row is an l value, columns traverse energy vector given.

    See Also
    --------
    PS_explicit : Calculates penetrability and shift functions using explicit definitions.
    
    Examples
    --------
    >>> from sample_resparm import PS_functions
    >>> PS_functions.PS_recursive(np.array([10.4]), 6.7e-15, 2)
    [array([[ 0.        ],
            [-0.99997817],
            [-1.99999272]]),
     array([[4.67244381e-03],
            [1.02005310e-07],
            [2.47442770e-13]])]
    """
    rho = k_wavenumber(E)*ac

    S_array = np.zeros([orbital_angular_momentum+1,len(E)])
    P_array = np.ones([orbital_angular_momentum+1,len(E)])
    P_array[0,:] *= rho
    
    for l in range(1,orbital_angular_momentum+1):
        S = (rho**2*(l-S_array[l-1]) / ((l-S_array[l-1])**2 + P_array[l-1]**2)) - l
        P = rho**2*P_array[l-1] / ((l-S_array[l-1])**2 + P_array[l-1]**2)
        
        S_array[l,:]=S; P_array[l,:]=P
        
    return S_array, P_array


def PS_explicit(E, orbital_angular_momentum):
    """
    Calculates penetrability and shift functions using explicit definitions.

    This function calculates the centifugal barrier penetrability as well as
    the shift factor for a neutron incident on Cu-63. The explicit implementation
    only allows for l-values up to 3. Some nucelar parameters are housed within 
    this function, this could be updated for other nuclei.

    Parameters
    ----------
    E : float or numpy.ndarray
        Energy of the incident neutron.
    ac : float 
        The channel radius in meters.
    orbital_angular_momentum : int
        Orbital angular momentum of the pair, describes waveform (l).

    Returns
    -------
    S_array : array-like
        shift factor(s) at given energy.
    P_array : array-like
        Penetrability at given energy.

    See Also
    --------
    PS_recursive : Calculates penetrability and shift functions using recursion.
    
    Examples
    --------
    >>> from sample_resparm import PS_functions
    >>> PS_functions.PS_explicit(np.array([10.4, 10.5]), 6.7e-15, 2)
    [array([-1.99999272, -1.99999265]), array([2.4744277e-13, 2.5343386e-13])]
    """
    rho = k_wavenumber(E)*ac
    if orbital_angular_momentum == 0:
        P = rho
        S = np.zeros(len(E))
    if orbital_angular_momentum == 1:
        P = rho**3/(1+rho**2)
        S = -1/(1+rho**2)
    if orbital_angular_momentum == 2:
        P = rho**5/(9+3*rho**2+rho**4)
        S = -(18+3*rho**2)/(9+3*rho**2+rho**4)
    if orbital_angular_momentum == 3:
        P = rho**7/(255+45*rho**2+6*rho**4+rho**6)
        S = -(675+90*rho**2+6*rho**4)/(225+45*rho**2+6*rho**4+rho**6)
        
    return S, P
        
        
        


        
    



