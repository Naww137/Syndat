#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:28:37 2022

@author: noahwalton
"""

import numpy as np
import matplotlib.pyplot as plt
import os


def sample_wigner_invCDF(N_samples):
    """
    Sample the wigner distribution.

    This function simply samples from the wigner distribution using inverse
    CDF sampling and is used by other functions for generating resonance level spacing.

    Parameters
    ----------
    N_samples : int
        Number of samples and/or length of sample vector.

    Returns
    -------
    numpy.ndarray or float
        Array of i.i.d. samples from wigner distribution.

    Notes
    -----
    See preliminary methods for sampling resonance level spacing from GOE.
    
    Examples
    --------
    >>> from sample_resparm import sample_levels
    >>> np.random.seed(7)
    >>> sample_levels.sample_wigner_invCDF(2,10)
    array([1.7214878 , 1.31941784])
    """
    samples = np.sqrt(-4/np.pi*np.log(np.random.default_rng().uniform(low=0.0,high=1.0,size=N_samples)))
    #samples = np.sqrt(-4*np.log(np.random.default_rng().uniform(low=0.0,high=1.0,size=N_samples)))      # remove the pi terms to match GOE
    if N_samples == 1:
        samples = np.ndarray.item(samples)
    return samples

def generate_GOE(N):
    A = np.random.default_rng().standard_normal(size=(N,N))/np.sqrt(2*N)
    X = A + np.transpose(A)
    return X

def wigner_PDF(x, avg_level_spacing):
    x = x/avg_level_spacing
    y = (np.pi/2) * x * np.exp(-np.pi*(x**2)/4)
    #y = (1/2) * x * np.exp(-(x**2)/4)   # remove the pi terms to match GOE
    y = y/avg_level_spacing
    return y

def sample_resonance_levels(E0, N_levels, avg_level_spacing, method):
    
    E0 = E0+avg_level_spacing*np.random.default_rng().uniform(low=0.0,high=1.0) # offset starting point so we are not just finding the distribution each time
    if method == 'invCDF':
        level_spacing = avg_level_spacing*sample_wigner_invCDF(N_levels)
        levels = [E0+level_spacing[0]]
        
        for ilevel in range(1,N_levels):
            levels.append(levels[ilevel-1]+level_spacing[ilevel])
            
    elif method == 'GOE':
        # Use resonance_generator from FUDGE:
        import sys
        sys.path.insert(0, '/Users/colefritsch/ENCORE/syndat/nuc_syndat/fudge/')
        from brownies.BNL.restools.resonance_generator import getFakeResonanceSet
        from xData import XYs1d

        EMax = (N_levels+5)*avg_level_spacing
        levelDensity = XYs1d.XYs1d(data=[[E0, 1.0/avg_level_spacing], [EMax, 1.0/avg_level_spacing]],
                                                axes=XYs1d.XYs1d.defaultAxes(
                                                    labelsUnits={
                                                            XYs1d.yAxisIndex: ('level_density', '1/eV'),
                                                            XYs1d.xAxisIndex: ('excitation_energy', 'eV')}))
        resonances = getFakeResonanceSet(E0=E0, aveD=avg_level_spacing, style='goe', L=0, J=0, levelDensity=levelDensity, DOFs=1, domainMin=E0,
                                         domainMax=EMax, widthKeys=())
        levels = np.array(resonances.data)[:N_levels,0]
        level_spacing = np.diff(levels)

    else:
        print('method for sampling resonance levels is not recognized')
        os.sys.exit()
            
    return levels, level_spacing



def sample_RRR_levels(E_range, avg_level_spacing, method='GOE'):
    """
    Sample the resonance energy levels.

    This function samples the wigner distribution using invCDF method in order 
    to get a ladder of resonance energy levels within a specified range. The energy range given
    is expanded by 5-6 times the average level spacing, a resonance ladder is sampled over that, 
    then the ladder is filtered to the energy range of interest.

    Parameters
    ----------
    E_range : array-like
        Array energies in RRR, only need min/max.
    avg_level_spacing : float
        Average level spacing value to scale wigner distribution.

    Returns
    -------
    levels : numpy.ndarray
        Array of resonance energy levels.
    spacings : numpy.ndarray
        Array of i.i.d. resonance level spacing samples.

    See Also
    --------
    sample_resonance_levels : Samples a specified number of resonance energy levels.

    Notes
    -----
    
    
    Examples
    --------
    >>> from sample_resparm import sample_levels
    >>> np.random.seed(7)
    >>> sample_levels.sample_RRR_levels([0.1,10], 2)
    ([array([6.31322239]),
      array([6.56223504]),
      array([8.65279185]),
      array([10.27692974])],
     [array([6.21322239]),
      array([0.24901265]),
      array([2.09055681]),
      array([1.62413789])])
    """
    
    if method == 'invCDF':
        # randomly offset starting point so we are not just finding the distribution fixed to this point with ML
        # is this necessary?
        # sample a ladder 5-6 average level spacings before and 5-6 average level spacings after window
        E0 = min(E_range)-avg_level_spacing*np.random.default_rng().uniform(low=5.0,high=6.0)     
        E_end = max(E_range)+avg_level_spacing*np.random.default_rng().uniform(low=5.0,high=6.0)   
        
        levels = []; spacings = []
        spacing = avg_level_spacing*sample_wigner_invCDF(1)
        spacings.append(spacing)
        level = E0+spacing
        
        while level < E_end:
            levels.append(level)
            spacing = avg_level_spacing*sample_wigner_invCDF(1)
            spacings.append(spacing)
            level = levels[-1] + spacing

        levels = list(filter(lambda l: l<max(E_range) and l>min(E_range), levels))

    elif method == 'GOE':
        # Use resonance_generator from FUDGE:
        import sys
        sys.path.insert(0, '/Users/colefritsch/ENCORE/syndat/nuc_syndat/fudge/')
        from brownies.BNL.restools.resonance_generator import getFakeResonanceSet
        from xData import XYs1d

        levelDensity = XYs1d.XYs1d(data=[[min(E_range), 1.0/avg_level_spacing], [max(E_range), 1.0/avg_level_spacing]],
                                                axes=XYs1d.XYs1d.defaultAxes(
                                                    labelsUnits={
                                                            XYs1d.yAxisIndex: ('level_density', '1/eV'),
                                                            XYs1d.xAxisIndex: ('excitation_energy', 'eV')}))
        resonances = getFakeResonanceSet(E0=min(E_range), aveD=avg_level_spacing, style='goe', L=0, J=0, levelDensity=levelDensity, DOFs=1, domainMin=min(E_range),
                                         domainMax=max(E_range), widthKeys=())
        levels = np.array(resonances.data)[:,0]
        spacings = np.diff(levels)
        print('Hello!')
            
    return levels, spacings


def compare_pdf_to_samples(level_spacing_vector, avg_level_spacing, method):
    
    fig = plt.figure(num=1,frameon=True); ax = fig.gca()
    
    x = np.linspace(0,max(level_spacing_vector),10000)
    plt.plot(x, wigner_PDF(x,avg_level_spacing), color='r', label='Wigner PDF', zorder=10)
        
    if method == 'GOE':
        print(); print('WARNING: ')
        print('GOE sampling does not match wigner pdf exactly'); print()
        plt.hist(level_spacing_vector, bins=75, density=True, ec='k', linewidth=0.75,color='cornflowerblue', zorder=2, label='GOE')

    elif method == 'invCDF':
        plt.hist(level_spacing_vector, bins=75, density=True, ec='k', linewidth=0.75,color='cornflowerblue', zorder=2, label='invCDF')
        
    else:
        print(); print('WARNING: ')
        print('no appropriate method selected for pdf comparison'); print()
    
    ax.set_facecolor('whitesmoke'); ax.grid(color='w', linestyle='-', linewidth=2, zorder=1, alpha=1)
    ax.set_xlabel('Level Spacing'); ax.set_ylabel('Normalized Frequency'); plt.title('Distribution of Level Spacing Samples')
    plt.legend()
    plt.show(); plt.close()
    
    return


