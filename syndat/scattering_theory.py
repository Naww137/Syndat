

import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats as stats


def g(J, I, i):
    """
    Calculates the spin statistical factor for a given $J_{\pi,\alpha}$.

    The spin statistical factor (g) is a weigting factor describing the probability of the different total angular momenta, give by:
    </br> $ g_{J\alpha} = \frac{2J+1}{(2i+1)(2I+1)} $ </br>

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
    Calculates the angular wavenumber of the compound state nucleus at a given incident energy.

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
    Constant = 0.002197 #sqrt(2Mn)/hbar
    hbar = 6.582119569e-16 # eV-s
    c = 2.99792458e8 # m/s
    m_eV = 939.565420e6 # eV/c^2

    # k = 1/hbar * M/(m+M) * np.sqrt(2*m_eV*E) * 1/c
    k = Constant*(M/(M+m))*np.sqrt(E)
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


def P_S_psi_explicit(E, ac, M, m, orbital_angular_momentum):
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

    rho = k_wavenumber(E, M, m)*ac
    psi = rho
    
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
        
    return S, P, psi




def reduced_width_square_2_partial_width(E, ac, M, m, reduced_widths_square, orbital_angular_momentum):
    S,P,psi = P_S_psi_explicit(np.array(E), ac, M, m, orbital_angular_momentum)
    partial_widths = 2*P*reduced_widths_square 
    return partial_widths


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
    resonance_ladder : _type_
        _description_

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
        Gn = sum([row[ign] for ign in range(2,len(row))]) * 1e-3
        Gg = row.Gg * 1e-3
        E_lambda = row.E
        d = (E-E_lambda)**2 + ((Gg+Gn)/2)**2 
        xs += (Gg*Gn) / ( d )
    xs = constant*xs
    return xs
    
