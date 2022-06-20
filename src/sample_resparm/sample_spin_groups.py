#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 07:48:31 2022

@author: nwx
"""

import os
import numpy as np
from sample_resparm import sample_widths
from sample_resparm import sample_levels
import pandas as pd

#%%






def quant_vec_sum(a,b):
    """
    Calculates a quantum vector sum.

    This function performs a quantum vector sum, a.k.a. it maps the quantum 
    triangular relationship between two integers or half integers.

    Parameters
    ----------
    a : float or int
        a variable.
    b : float or int
        a variable.

    Returns
    -------
    numpy.ndarray
        Array of all possible quantum values.
    """
    vec = np.arange(abs(a-b), a+b+1, 1)
    return vec


def map_quantum_numbers(I,i,l_wave_max, print_out):
    """
    Maps the possible quantum numbers for pair.

    This function maps out the possible quantum spin numbers (Jpi) for a given
    particle pair up to some maximum considered incident waveform (l-wave).

    Parameters
    ----------
    I : float or int
        Intrinsic spin of the target particle.
    i : float or int
        Intrinsic spin of the incident paritcle.
    l_wave_max : int
        Maximum considered incident waveform (l-wave).
    print_out : bool
        User option to print out quantum spin (J) mapping to console.

    Returns
    -------
    Jn : array-like
        List containing possible J, # of contibuting channels, and contibuting 
        waveforms for negative parity. Formatted as (J,#chs,[l-wave, l-wave])
    Jp : array-like
        List containing possible J, # of contibuting channels, and contibuting 
        waveforms for positive parity. Formatted as (J,#chs,[l-wave, l-wave])
    Notes
    -----
    Currently this function is explicitly coded for target particles
    with negative parity and incident particles with positive parity.
    
    Examples
    --------
    >>> from sample_resparm import sample_spin_groups
    >>> sample_spin_groups.map_quantum_numbers(3/2,1/2,2, False)
    ([(1.0, 1, [0.0]), (2.0, 1, [0.0])],
     [(0.0, 1, [1.0]),
      (1.0, 2, [1.0, 1.0]),
      (2.0, 2, [1.0, 1.0]),
      (3.0, 1, [1.0])])
    """
    Jn = []; Jp = []
    S = quant_vec_sum(I,i)
    L = range(l_wave_max+1)
    
    possible_Jpi = {}
    J_negative = []; J_positive = []
    for i_l, l in enumerate(L):
        this_l = {}
        
        if l==0 or l==2:
            parity = -1
        elif l==1 or l==3:
            parity = 1
        
        for i_s, s in enumerate(S):
            js = quant_vec_sum(s,l)
            this_l[f's={s}'] = js
            for j in js:
                if parity == 1:
                    J_positive.append([l,s,j])
                if parity == -1:
                    J_negative.append([l,s,j])
            
        possible_Jpi[f'l={l}'] = this_l
        
    Jn_total = np.array(J_negative)
    Jn_unique = np.unique(Jn_total[:,2])

    for j in Jn_unique:
        entrance_channels = np.count_nonzero(Jn_total[:,2] == j)
        
        ls = []; ss = []
        for i, jtot in enumerate(Jn_total[:,2]):
            if jtot == j:
                ls.append(Jn_total[i,0])
                ss.append(Jn_total[i,1])
                
        Jn.append((j,entrance_channels,ls))
        
    # condition for if there are no positive parity (l=0 only)
    if L[-1] > 0:
        Jp_total = np.array(J_positive)
        Jp_unique = np.unique(Jp_total[:,2]) 
        
        for j in Jp_unique:
            entrance_channels = np.count_nonzero(Jp_total[:,2] == j)
            
            ls = []; ss = []
            for i, jtot in enumerate(Jp_total[:,2]):
                if jtot == j:
                    ls.append(Jp_total[i,0])
                    ss.append(Jp_total[i,1])
                
            Jp.append((j,entrance_channels,ls))
    
    
    if print_out:
        print()
        print('The following arrays describe all possible spin groups for a each parity.\n\
        The data is given as a tuple where the first value is the integer \n\
        or half integer total quantum spin J and the second value is the \n\
        number of entrance channels for that spin group. \n\
        * See the dictionary "possible_Jpi" for a nested packing structure.')
        print()
        print('Spin group data for negative parity\n(J-, #Chs, l-waves)')
        for each in Jn:
            print(each)
        print()
        print('Spin group data for positive parity\n(J+, #Chs, l-waves)')
        for each in Jp:
            print(each)

    return Jn, Jp




def sample_each_Jpi(I, i, l_wave_max,  
                    Erange, 
                    Davg, Gavg, 
                    Gavg_swave, 
                    print_out = True,
                    save_csv = True, 
                    sammy_run_folder = os.getcwd()):
    """
    Samples a full resonance parameter ladder for each possible spin group.

    This function samples resonance parameters (Energy and widths) for each 
    possible spin group (Jpi) of a given particle pair. The results can be 
    printed to the console and/or saved to a csv.

    Parameters
    ----------
    I : float or int
        Intrinsic spin of the target particle.
    i : float or int
        Intrinsic spin of the incident paritcle.
    l_wave_max : int
        Maximum considered incident waveform (l-wave).
    Erange : array-like
        Array of resolve resonance range energy, only requires min/max.
    Davg : array-like
        Nested list of average level spacing for each spin group number. First 
        list is for negative parity (J-) second is for positive parity (J+).
    Gavg : array-like
        Nested list of average widths for each spin group number. First 
        list is for negative parity (J-) second is for positive parity (J+).
    Gavg_swave : float
        Average width used to sample agregate capture widths. **Unsure of this value.
    print_out : bool
        User option to print out quantum spin (J) mapping to console.
    save_csv : bool
        User option to save resonance ladders to csv.
    sammy_run_folder : str
        Folder in which the csv(s) containing resparm ladders will be saved.

    Notes
    -----
    Unsure of the average capture width for Gg sampling.
    
    Returns
    -------
    Jn_df : DataFrame
        Pandas DataFrame conatining a resonance parameter ladder for each 
        quantum spin group with negative parity (all J-).
    Jp_df : DataFrame
        Pandas DataFrame conatining a resonance parameter ladder for each 
        quantum spin group with positive parity (all J+).
    """
    
    [Jn, Jp] = map_quantum_numbers(I,i,l_wave_max, print_out)
    
# =============================================================================
#     negative parity Js
# =============================================================================
    Jn_df = pd.DataFrame()
    for ij, j in enumerate(Jn):
        
        [levels, level_spacing] = sample_levels.sample_RRR_levels(Erange, Davg[0][ij])
        
        Gnx = []; gnx = []
        for ichannel, lwave in enumerate(j[2]):
            [red_nwidth_2, nwidth] = sample_widths.sample_RRR_widths(levels, Gavg[0][ij], 1, lwave)
            Gnx.append(nwidth); gnx.append(red_nwidth_2)
              
        # what is the orbital angular momentum of the gamma width?
        [gg, Gg] = sample_widths.sample_RRR_widths(levels, Gavg_swave, 100, 0)
        
        Gn = pd.DataFrame(Gnx)
        E_Gg = pd.DataFrame([levels, Gg], index=['E','Gg'])
        
        parm_df = pd.concat([E_Gg,Gn], axis=0)
        parm_dfv = parm_df.transpose()
        
        Jn_df[f'{j}'] = ""
        Jn_df = pd.concat([Jn_df,parm_dfv], axis=1)
        
        if save_csv:
            parm_dfv.to_csv(os.path.join(sammy_run_folder, f'J_neg_{j[0]}.csv'))
        
# =============================================================================
#         if print_out:
#             print(); print(j)
#             print(parm_dfv); print()
# =============================================================================
            
# =============================================================================
#       positive parity Js
# =============================================================================
    Jp_df = pd.DataFrame()
    for ij, j in enumerate(Jp):
        
        [levels, level_spacing] = sample_levels.sample_RRR_levels(Erange, Davg[1][ij])
        
        Gnx = []; gnx = []
        for ichannel, lwave in enumerate(j[2]):
            [red_nwidth_2, nwidth] = sample_widths.sample_RRR_widths(levels, Gavg[1][ij], 1, lwave)
            Gnx.append(nwidth); gnx.append(red_nwidth_2)
              
        # what is the orbital angular momentum of the gamma width?
        [gg, Gg] = sample_widths.sample_RRR_widths(levels, Gavg_swave, 100, 0)
        
        Gn = pd.DataFrame(Gnx)
        E_Gg = pd.DataFrame([levels, Gg], index=['E','Gg'])
        
        parm_df = pd.concat([E_Gg,Gn], axis=0)
        parm_dfv = parm_df.transpose()
        
        Jp_df[f'{j}'] = ""
        Jp_df = pd.concat([Jp_df,parm_dfv], axis=1)
        
        if save_csv:
            parm_dfv.to_csv(os.path.join(sammy_run_folder, f'J_pos_{j[0]}.csv'))
        
# =============================================================================
#         if print_out:
#             print(); print(j)
#             print(parm_dfv); print()
# =============================================================================
            
# =============================================================================
#    save both
# =============================================================================
    if save_csv:
        Jn_df.to_csv(os.path.join(sammy_run_folder,'Jn_all.csv'))
        Jp_df.to_csv(os.path.join(sammy_run_folder,'Jp_all.csv'))
        
    return Jn_df, Jp_df
  