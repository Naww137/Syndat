#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 12:54:52 2022

@author: noahwalton
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats as stats


def gstat(J, I, i):
    """
    Calculates the spin statistical factor for a given $J_{\pi,\alpha}$.

    The spin statistical factor (g) is a weigting factor describing the probability of the different total angular momenta, given by:
    
    $ g_{J\alpha} = \frac{2J+1}{(2i+1)(2I+1)} $ 


    Parameters
    ----------
    J : Float or int
        Total angular momentum of the channel.
    I : float or int
        Spin and parity of the target particle.
    i : float or int
        Spin and parity of the incident particle.

    Returns
    -------
    _type_
        _description_

    """
    return (2*J+1)/((2*i+1)*(2*I+1))


def k_wavenumber(E, M, m):
    """
    Calculates the angular wavenumber of the compound state nucleus at a given incident energy, units .

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

    """
    # constants 
    hbar = 6.582119569e-16 # eV-s
    c = 2.99792458e8 # m/s
    mn_eV = 939.565420e6 # eV/c^2
    constant = (np.sqrt(2*mn_eV)/c/hbar)*(1e-14) # 0.002197 #sqrt(2Mn)/hbar 

    k = (M/(M+m))*constant*np.sqrt(E)
    return k
    

def PS_recursive(E, ac, M, m, orbital_angular_momentum):
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
            Scattering channel radius in meters.
        M : float or int
            Mass of the target nucleus.
        m : float or int
            Mass of the incident particle.
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
        rho = k_wavenumber(E, M, m)*ac

        S_array = np.zeros([orbital_angular_momentum+1,len(E)])
        P_array = np.ones([orbital_angular_momentum+1,len(E)])
        P_array[0,:] *= rho
        
        for l in range(1,orbital_angular_momentum+1):
            S = (rho**2*(l-S_array[l-1]) / ((l-S_array[l-1])**2 + P_array[l-1]**2)) - l
            P = rho**2*P_array[l-1] / ((l-S_array[l-1])**2 + P_array[l-1]**2)
            
            S_array[l,:]=S; P_array[l,:]=P
            
        return S_array, P_array


def FofE_explicit(E, ac, M, m, orbital_angular_momentum):
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
        Scattering channel radius in meters.
    M : float or int
            Mass of the target nucleus.
    m : float or int
        Mass of the incident particle.
    orbital_angular_momentum : int
        Orbital angular momentum of the pair, describes waveform (l).

    Returns
    -------
    S_array : array-like
        Shift factor(s) at given energy.
    P_array : array-like
        Penetrability at given energy.
    psi : array-like
        Potential scattering phase shift.

    See Also
    --------
    PS_recursive : Calculates penetrability and shift functions using recursion.
    
    Examples
    --------
    >>> from sample_resparm import PS_functions
    >>> PS_functions.PS_explicit(np.array([10.4, 10.5]), 6.7e-15, 2)
    [array([-1.99999272, -1.99999265]), array([2.4744277e-13, 2.5343386e-13])]
    """
    
    assert(orbital_angular_momentum == 0, "Phase shift function in syndat.scattering theory needs to be updated for higher-order waveforms")

    k = k_wavenumber(E, M, m)
    rho = k*ac
    phi = rho
    
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

    if orbital_angular_momentum > 3:
        raise ValueError("PS_explicit cannot handle orbital_angular_momenta > 3, use PS_recursive")
        
    return S, P, phi, k




def reduced_width_square_2_partial_width(E, ac, M, m, reduced_widths_square, orbital_angular_momentum):
    S, P, psi, k = FofE_explicit(np.array(E), ac, M, m, orbital_angular_momentum)
    partial_widths = 2*P*reduced_widths_square 
    return partial_widths


def SLBW(E, pair, resonance_ladder):
    """
    _summary_

    _extended_summary_

    Parameters
    ----------
    E : _type_
        _description_
    pair : _type_
        _description_
    resonance_ladder : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """

    # TODO: check types in resonance ladder

    if resonance_ladder.empty:
        xs_cap = 0
        xs_scat = 1/np.sqrt(E)
        xs_tot = xs_scat + xs_cap
        return xs_tot, xs_scat, xs_cap
    # TODO: if resonance ladder is empty, calculate nominal cross section - the scattering phase shift still exists, but what is the Jpi to calculate it? 1/np.sqrt(E)?
        
    xs_cap = 0; xs_scat = 0
    group_by_J = dict(tuple(resonance_ladder.groupby('J')))

    for Jpi in group_by_J:
        
        J_df = group_by_J[Jpi]
        # assert J > 0 
        Jpi = abs(Jpi)

        # orbital_angular_momentum = J_df.lwave.unique()
        orbital_angular_momentum = np.unique(J_df.lwave.values[0]) # just takes the first level's l-wave vector, all should be the same in each group
        assert len(orbital_angular_momentum) == 1, "Cannot handle different l-waves contributing to multichannel widths"
        
        # calculate functions of energy -> shift, penetrability, phase shift
        g = gstat(Jpi, pair.I, pair.i) #(2*J+1)/( (2*ii+1)*(2*I+1) );   # spin statistical factor g sub(j alpha)
        shift, penetration, phi, k = FofE_explicit(E, pair.ac, pair.M, pair.m, orbital_angular_momentum[0])

        # calculate capture
        sum1 = 0
        for index, row in J_df.iterrows():
            E_lambda = row.E
            Gg = row.Gg * 1e-3
            # gnx2 = sum([row[ign] for ign in range(2,len(row))]) * 1e-3  # Not sampling multiple, single-channel particle widths
            gnx2 = row.gnx2 * 1e-3 
            Gnx = 2*penetration*gnx2

            d = (E-E_lambda)**2 + ((Gg+Gnx)/2)**2 
            sum1 += (Gg*Gnx) / ( d )

        xs_cap += (np.pi*g/(k**2))*sum1


        # calculate scatter
        sum1 = 0
        sum2 = 0
        sum3 = 0
        for index, row in J_df.iterrows():
            E_lambda = row.E
            Gg = row.Gg * 1e-3
            # gnx2 = sum([row[ign] for ign in range(2,len(row))]) * 1e-3 # Not sampling multiple, single-channel particle widths
            gnx2 = row.gnx2 * 1e-3 
            Gnx = 2*penetration*gnx2

            Gtot = Gnx+Gg
            d = (E-E_lambda)**2 + ((Gg+Gnx)/2)**2 
            sum1 += Gnx*Gtot/d
            sum2 += Gnx*(E-E_lambda)/d
            sum3 += (Gnx*(E-E_lambda)/d)**2 + (Gnx*Gtot/d/2)**2

        xs_scat += (np.pi*g/(k**2))* ( (1-np.cos(2*phi))*(2-sum1) + 2*np.sin(2*phi)*sum2 + sum3 )

    # calculate total
    xs_tot = xs_cap+xs_scat

    return xs_tot, xs_scat, xs_cap






def SLBW_capture(g, k, E, resonance_ladder):
    """
    Calculates a multi-level capture cross section using SLBW formalism.

    _extended_summary_

    Parameters
    ----------
    g : float
        Spin statistical factor $g_{J,\alpha}$.
    k : float or array-like
        Angular wavenumber or array of angular wavenumber values corresponding to the energy vector.
    E : float or array-like
        KE of incident particle or array of KE's.
    resonance_ladder : DataFrame
        DF with columns for 

    Returns
    -------
    xs
        SLBW capture cross section.
    """

    if len(k) != len(E):
        raise ValueError("Vector of angular wavenumbers, k(E), does not match the length of vector E")

    xs = 0
    constant = (np.pi*g/(k**2))
    for index, row in resonance_ladder.iterrows():
        E_lambda = row.E
        gnx2 = sum([row[ign] for ign in range(2,len(row))]) * 1e-3
        Gnx = 2*P(E)*gnx2
        Gg = row.Gg * 1e-3
        d = (E-E_lambda)**2 + ((Gg+Gn)/2)**2 
        xs += (Gg*Gn) / ( d )
    xs = constant*xs
    return xs
    


