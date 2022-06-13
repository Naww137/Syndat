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
    samples = np.sqrt(-4/np.pi*np.log(np.random.default_rng().uniform(low=0.0,high=1.0,size=N_samples)))
    #samples = np.sqrt(-4*np.log(np.random.default_rng().uniform(low=0.0,high=1.0,size=N_samples)))      # remove the pi terms to match GOE
    return samples

def generate_GOE(N):
    A = np.random.default_rng().standard_normal(size=(N,N))/np.sqrt(2*N)
    X = A + np.transpose(A)
    return X

def wigner_PDF(x):
    y = (np.pi/2) * x * np.exp(-np.pi*(x**2)/4)
    #y = (1/2) * x * np.exp(-(x**2)/4)   # remove the pi terms to match GOE
    return y

def sample_resonance_levels(E0, N_levels, avg_level_spacing, method):
    
    if method == 'invCDF':
        level_spacing = avg_level_spacing*sample_wigner_invCDF(N_levels)
            
    elif method == 'GOE':
        level_spacing = []
        for ilevel in range(N_levels):
            X = generate_GOE(2)
            eigenvalues = np.linalg.eigvals(X)
            spacing = avg_level_spacing*abs(np.diff(eigenvalues))
            level_spacing.append(spacing.item())
    else:
        print('method for sampling resonance levels is not recognized')
        os.sys.exit()
            
    E0 = E0+avg_level_spacing*np.random.default_rng().uniform(low=0.0,high=1.0) # offset starting point so we are not just finding the distribution each time
    levels = [E0+level_spacing[0]]
    
    for ilevel in range(1,N_levels):
        levels.append(levels[ilevel-1]+level_spacing[ilevel])
            
    return levels, level_spacing


def compare_pdf_to_samples(level_spacing_vector, avg_level_spacing, method):
    
    fig = plt.figure(num=1,frameon=True); ax = fig.gca()
    
    x = np.linspace(0,max(level_spacing_vector),10000)
    plt.plot(x, wigner_PDF(x), color='r', label='Wigner PDF', zorder=10)
    
    if avg_level_spacing != 1:
        print(); print('WARNING: ')
        print('pdf has not been transformed for <D> other than 1 - will not match samples'); print()
        
    if method == 'GOE':
        print(); print('WARNING: ')
        print('GOE sampling has not been transformed by factor of pi - will not match pdf or invCDF sampling'); print()
        plt.hist(level_spacing_vector, bins=75, density=True, ec='k', linewidth=0.75,color='cornflowerblue', zorder=2, label='GOE')

    elif method == 'invCDF':
        plt.hist(level_spacing_vector, bins=75, density=True, ec='k', linewidth=0.75,color='pink', zorder=2, label='invCDF')
        
    else:
        print(); print('WARNING: ')
        print('no appropriate method selected for pdf comparison'); print()
    
    ax.set_facecolor('whitesmoke'); ax.grid(color='w', linestyle='-', linewidth=2, zorder=1, alpha=1)
    ax.set_xlabel('Level Spacing'); ax.set_ylabel('Normalized Frequency'); plt.title('Distribution of Level Spacing Samples')
    plt.legend()
    
    return




